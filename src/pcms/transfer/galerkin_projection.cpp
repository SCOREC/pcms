#include "pcms/transfer/galerkin_projection.hpp"
#include "pcms/field/detail/field_ops.hpp"
#include "pcms/field/field_ops.hpp"
#include "pcms/utility/arrays.hpp"
#include <Kokkos_Core.hpp>

namespace pcms::transfer
{

namespace
{

// samples -= g. A free function: nvcc does not accept an extended lambda
// inside a private member function (ApplyImpl).
void SubtractInPlace(Kokkos::View<Real**, DeviceMemorySpace> samples,
                     Kokkos::View<Real**, DeviceMemorySpace> g)
{
  const auto num_points = static_cast<LO>(samples.extent(0));
  const auto num_components = static_cast<int>(samples.extent(1));
  Kokkos::parallel_for(
    "projection_residual",
    Kokkos::RangePolicy<DefaultExecutionSpace>(0, num_points),
    KOKKOS_LAMBDA(const LO i) {
      for (int c = 0; c < num_components; ++c) {
        samples(i, c) -= g(i, c);
      }
    });
}

const Geometry& SourceGeometry(const PointEvaluatorFactory& source)
{
  const auto geometry = source.GetGeometry();
  if (geometry == nullptr) {
    throw pcms_error("GalerkinProjection: the source has no geometry");
  }
  return *geometry;
}

Quadrature MakeQuadrature(const FunctionSpace& target_space,
                          const ProjectionOptions& options)
{
  if (options.quadrature) {
    return *options.quadrature;
  }
  return CreateGaussQuadrature(target_space.GetGeometry(), 2);
}

} // namespace

GalerkinProjection::GalerkinProjection(
  std::shared_ptr<const PointEvaluatorFactory> source_ptr,
  std::shared_ptr<const FunctionSpace> target_space_ptr,
  ProjectionOptions options)
  : source_(std::move(source_ptr)),
    target_space_(std::move(target_space_ptr)),
    right_hand_side_(target_space_, MakeQuadrature(*target_space_, options)),
    mass_(target_space_, target_space_->CreateQuadrature()),
    functional_(detail::ScratchFunction(*target_space_)),
    coefficients_(detail::ScratchFunction(*target_space_))
{
  const Geometry& source_geometry = SourceGeometry(*source_);
  const FunctionSpace& target_space = *target_space_;
  // The source is bound to the points in the cells of its geometry that the
  // quadrature names (no search), else searched for.
  source_at_points_ =
    source_->CreatePointEvaluator<Real>(detail::BindQuadrature(
      source_geometry, right_hand_side_.GetQuadrature(), options.to_source));
  if (options.mass_matrix_type == MassMatrixType::Lumped) {
    lumped_ =
      std::make_unique<Function<Real>>(detail::ScratchFunction(target_space));
    mass_.AssembleLumped(*lumped_);
  } else {
    solver_ = std::make_unique<detail::LinearSolver>(mass_, options.solve);
  }
  if (options.control_variate) {
    control_variate_.emplace(
      source_, target_space_,
      source_geometry.Bind(target_space.GetLayout()->GetDOFHolderCoordinates(),
                           OutOfBoundsPolicy{}, options.to_source));
  }
  samples_ = Kokkos::View<Real**, DeviceMemorySpace>(
    "projection_samples", right_hand_side_.GetQuadrature().points->NumPoints(),
    target_space.GetLayout()->GetNumComponents());
}

void GalerkinProjection::ApplyImpl(const Field<Real>& source,
                                   Field<Real>& target) const
{
  const int num_components = target_space_->GetLayout()->GetNumComponents();
  if (source.GetLayout().GetNumComponents() != num_components ||
      source.GetData().GetValueType() != target.GetData().GetValueType()) {
    throw pcms_error("GalerkinProjection::Apply: the source and the target "
                     "must have the same value type and number of components");
  }
  if (&target.GetLayout() != target_space_->GetLayout().get()) {
    throw pcms_error("GalerkinProjection::Apply: target field layout mismatch");
  }
  const ValueBasis written =
    source_at_points_->OutputBasis(source.GetData().GetValueBasis());
  pcms::detail::CheckWrittenValueBasis("GalerkinProjection::Apply", written,
                                       target.GetData().GetValueBasis());
  auto samples = samples_;

  // 1. f at the points, minus the control variate g there; the control
  //    variate is the interpolation of the source onto the target space.
  source_at_points_->Evaluate(source, ValueView<Real, DeviceMemorySpace>(
                                        written, pcms::MakeRank2View(samples)));
  if (control_variate_) {
    if (!pcms::detail::SameValueBasis(
          control_variate_->values.GetData().GetValueBasis(),
          target.GetData().GetValueBasis())) {
      control_variate_->values = target_space_->CreateFunction<Real>(
        "", target.GetData().GetValueBasis());
    }
    control_variate_->interpolator.Apply(source, control_variate_->values);
    Kokkos::View<Real**, DeviceMemorySpace> g(
      "projection_control_variate", samples.extent(0), samples.extent(1));
    right_hand_side_.Evaluator().Evaluate(
      control_variate_->values,
      ValueView<Real, DeviceMemorySpace>(
        control_variate_->values.GetData().GetValueBasis(),
        pcms::MakeRank2View(g)));
    SubtractInPlace(samples, g);
  }

  // 2. b = E^T (W f), 3. M x = b (x = M^-1 b for the lumped mass)
  right_hand_side_.Assemble(
    ValueView<const Real, DeviceMemorySpace>(
      functional_.GetData().GetValueBasis(), MakeConstRank2View(samples_)),
    functional_);
  if (lumped_ != nullptr) {
    detail::Divide(functional_, *lumped_, target);
  } else {
    solver_->Solve(functional_, target);
  }
  // 4. x = g + M^-1 r
  if (control_variate_) {
    Axpy(Real(1), control_variate_->values, target);
  }
}

bool GalerkinProjection::SupportsTranspose() const noexcept
{
  return !control_variate_ && source_at_points_->SupportsTranspose();
}

void GalerkinProjection::ApplyTransposeImpl(const Field<Real>& source,
                                            Field<Real>& target) const
{
  if (!SupportsTranspose()) {
    throw pcms_error("GalerkinProjection::ApplyTranspose: needs a source "
                     "evaluator with the transpose and no control variate");
  }
  if (&source.GetLayout() != target_space_->GetLayout().get()) {
    throw pcms_error("GalerkinProjection::ApplyTranspose: the source must be "
                     "a field of the target space's layout");
  }
  const ValueBasis written =
    source_at_points_->OutputBasis(target.GetData().GetValueBasis());
  pcms::detail::CheckWrittenValueBasis("GalerkinProjection::ApplyTranspose",
                                       source.GetData().GetValueBasis(),
                                       written);
  // z = M^-1 b, E z at the points, W, then S^T into the functional
  if (lumped_ != nullptr) {
    detail::Divide(source, *lumped_, coefficients_);
  } else {
    solver_->Solve(source, coefficients_);
  }
  right_hand_side_.Evaluator().Evaluate(
    coefficients_,
    ValueView<Real, DeviceMemorySpace>(coefficients_.GetData().GetValueBasis(),
                                       MakeRank2View(samples_)));
  const auto weighted =
    right_hand_side_.Weight(ValueView<const Real, DeviceMemorySpace>(
      written, MakeConstRank2View(samples_)));
  pcms::detail::SetZero(target);
  source_at_points_->EvaluateTranspose(weighted, target);
}

} // namespace pcms::transfer
