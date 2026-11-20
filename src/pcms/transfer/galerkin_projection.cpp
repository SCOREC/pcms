#include "pcms/transfer/galerkin_projection.hpp"
#include "pcms/field/field_ops.hpp"
#include "pcms/utility/arrays.hpp"
#include <Kokkos_Core.hpp>

namespace pcms::transfer
{

namespace
{

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
    functional_(target_space_->CreateFunction<Real>())
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
      std::make_unique<Function<Real>>(target_space.CreateFunction<Real>());
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

void GalerkinProjection::Apply(const Field<Real>& source,
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
  const auto num_points = static_cast<LO>(samples_.extent(0));
  // Shallow copy so the device lambda captures the view, not `this`.
  auto samples = samples_;

  // 1. f at the points, minus the control variate g there; the control
  //    variate is the interpolation of the source onto the target space.
  source_at_points_->Evaluate(
    source, ValueView<Real, DeviceMemorySpace>(source.GetData().GetValueBasis(),
                                               pcms::MakeRank2View(samples)));
  if (control_variate_) {
    control_variate_->interpolator.Apply(source, control_variate_->values);
    Kokkos::View<Real**, DeviceMemorySpace> g(
      "projection_control_variate", samples.extent(0), samples.extent(1));
    right_hand_side_.Evaluator().Evaluate(
      control_variate_->values,
      ValueView<Real, DeviceMemorySpace>(
        control_variate_->values.GetData().GetValueBasis(),
        pcms::MakeRank2View(g)));
    Kokkos::parallel_for(
      "projection_residual",
      Kokkos::RangePolicy<DefaultExecutionSpace>(0, num_points),
      KOKKOS_LAMBDA(const LO i) {
        for (int c = 0; c < num_components; ++c) {
          samples(i, c) -= g(i, c);
        }
      });
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

} // namespace pcms::transfer
