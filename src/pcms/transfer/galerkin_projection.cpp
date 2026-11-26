#include "pcms/transfer/galerkin_projection.hpp"
#include "pcms/field/detail/field_ops.hpp"
#include "pcms/field/field_ops.hpp"

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
  std::shared_ptr<const PointEvaluatorFactory> source,
  std::shared_ptr<const FunctionSpace> target_space, ProjectionOptions options)
  : target_space_(std::move(target_space)),
    right_hand_side_(
      target_space_, source,
      detail::FormOptions{.quadrature = MakeQuadrature(*target_space_, options),
                          .to_trial = options.to_source}),
    mass_(target_space_),
    functional_(target_space_->CreateFunction<Real>()),
    coefficients_(target_space_->CreateFunction<Real>()),
    correction_(target_space_->CreateFunction<Real>())
{
  if (!mass_.SupportsApply()) {
    throw pcms_error("GalerkinProjection: the target space's evaluator does "
                     "not support the transpose");
  }
  const Geometry& source_geometry = SourceGeometry(*source);
  const FunctionSpace& target = *target_space_;
  if (options.mass_matrix_type == MassMatrixType::Lumped) {
    lumped_ = std::make_unique<Function<Real>>(target.CreateFunction<Real>());
    mass_.AssembleLumped(*lumped_);
  } else {
    solver_ = std::make_unique<detail::LinearSolver>(mass_, options.solve);
  }
  if (options.control_variate) {
    if (target.HasNodalDOFs()) {
      interpolation_ = std::make_unique<pcms::detail::Interpolator<Real>>(
        source, target_space_,
        source_geometry.Bind(target.GetLayout()->GetDOFHolderCoordinates(),
                             OutOfBoundsPolicy{}, options.to_source));
    } else {
      ProjectionOptions inner = options;
      inner.quadrature = target.GetQuadrature();
      inner.control_variate = false;
      interpolation_ = std::make_unique<GalerkinProjection>(
        source, target_space_, std::move(inner));
    }
    control_variate_.emplace(right_hand_side_.WithTrial(target_space_));
  }
}

Function<Real>& GalerkinProjection::InBasis(Function<Real>& scratch,
                                            const ValueBasis& basis) const
{
  if (scratch.GetData().GetValueBasis() != basis) {
    scratch = target_space_->CreateFunction<Real>("", basis.variance,
                                                  basis.system, basis.scaling);
  }
  return scratch;
}

void GalerkinProjection::ApplyImpl(const Field<Real>& source,
                                   Field<Real>& target) const
{
  if (source.GetLayout().GetNumComponents() !=
        target_space_->GetLayout()->GetNumComponents() ||
      source.GetData().GetValueType() != target.GetData().GetValueType()) {
    throw pcms_error("GalerkinProjection::Apply: the source and the target "
                     "must have the same value type and number of components");
  }
  if (&target.GetLayout() != target_space_->GetLayout().get()) {
    throw pcms_error("GalerkinProjection::Apply: target field layout mismatch");
  }
  const ValueBasis basis =
    TrialBasis(right_hand_side_, source.GetData().GetValueBasis());
  auto& functional = InBasis(functional_, basis);
  right_hand_side_.Apply(source, functional);
  if (control_variate_) {
    // b -= B_VV I u
    auto& interpolant = InBasis(coefficients_, basis);
    auto& correction = InBasis(correction_, basis);
    interpolation_->Apply(source, interpolant);
    control_variate_->Apply(interpolant, correction);
    Axpy(Real(-1), correction, functional);
  }
  if (lumped_ != nullptr) {
    detail::Divide(functional, *lumped_, target);
  } else {
    solver_->Solve(functional, target);
  }
  if (control_variate_) {
    Axpy(Real(1), coefficients_, target);
  }
}

bool GalerkinProjection::SupportsTranspose() const noexcept
{
  return right_hand_side_.SupportsTranspose() &&
         (interpolation_ == nullptr || interpolation_->SupportsTranspose());
}

void GalerkinProjection::ApplyTransposeImpl(const Field<Real>& source,
                                            Field<Real>& target) const
{
  if (!SupportsTranspose()) {
    throw pcms_error("GalerkinProjection::ApplyTranspose: needs a source "
                     "evaluator and a control variate with the transpose");
  }
  if (&source.GetLayout() != target_space_->GetLayout().get()) {
    throw pcms_error("GalerkinProjection::ApplyTranspose: the source must be "
                     "a field of the target space's layout");
  }
  const ValueBasis& basis = source.GetData().GetValueBasis();
  // z = M^-1 v
  auto& z = InBasis(coefficients_, basis);
  if (lumped_ != nullptr) {
    detail::Divide(source, *lumped_, z);
  } else {
    solver_->Solve(source, z);
  }
  if (control_variate_) {
    // I^T (v - B_VV z)
    auto& residual = InBasis(functional_, basis);
    auto& correction = InBasis(correction_, basis);
    control_variate_->Apply(z, correction);
    residual.SetDOFHolderDataUnchecked(source.GetDOFHolderData().GetValues());
    residual.Synchronize();
    Axpy(Real(-1), correction, residual);
    interpolation_->ApplyTranspose(residual, target);
  } else {
    pcms::detail::SetZero(target);
  }
  right_hand_side_.ApplyTranspose(z, target);
}

} // namespace pcms::transfer
