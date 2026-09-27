#include "pcms/transfer/omega_h_control_variate_projection.hpp"
#include "pcms/transfer/conservative_projection_solver.hpp"
#include "pcms/transfer/omega_h_mass_integrator.hpp"
#include "pcms/transfer/omega_h_mc_rhs_integrator.hpp"
#include "pcms/utility/arrays.h"
#include <Kokkos_Core.hpp>
#include <Omega_h_array.hpp>

namespace pcms
{

OmegaHControlVariateProjection::OmegaHControlVariateProjection(
  const PointEvaluatorFactory& source, const FunctionSpace& target_space,
  int samples_per_element, MonteCarloSampling sampling, uint64_t seed,
  MassMatrixType mass_matrix_type)
  : target_layout_(std::dynamic_pointer_cast<const OmegaHLagrangeLayout>(
      target_space.GetLayout())),
    rhs_integrator_(std::make_unique<OmegaHMonteCarloRHSIntegrator>(
      target_layout_, target_space.GetCoordinateSystem(), samples_per_element,
      sampling, seed)),
    interpolator_(source, target_space),
    control_variate_(target_space.CreateFunction<Real>())
{
  const auto sample_coords = rhs_integrator_->GetIntegrationPoints();
  source_at_samples_ = source.CreatePointEvaluator<Real>(
    EvaluationRequest::FromCoordinates(sample_coords));
  control_variate_at_samples_ = target_space.CreatePointEvaluator<Real>(
    EvaluationRequest::FromCoordinates(sample_coords));

  // Mass integrator is only needed to build the solver; PETSc reference-counts
  // the matrix so it remains alive inside the KSP after this scope ends.
  OmegaHMassIntegrator mass_integrator(
    target_layout_, target_space.GetCoordinateSystem(), mass_matrix_type);
  solver_ = std::make_unique<GalerkinProjectionSolver>(mass_integrator,
                                                       *rhs_integrator_);

  target_values_ = Kokkos::View<Real**, DeviceMemorySpace>(
    "cv_target_values", target_layout_->GetNumOwnedDofHolder(),
    target_layout_->GetNumComponents());
  const auto num_samples = sample_coords.GetValues().extent(0);
  f_samples_ =
    Kokkos::View<Real**, DeviceMemorySpace>("cv_f_samples", num_samples, 1);
  residual_ =
    Kokkos::View<Real**, DeviceMemorySpace>("cv_residual", num_samples, 1);
}

// Defined here so that GalerkinProjectionSolver (forward-declared in the
// header) is a complete type when unique_ptr's destructor is instantiated.
OmegaHControlVariateProjection::~OmegaHControlVariateProjection() = default;

void OmegaHControlVariateProjection::Apply(const Field<Real>& source,
                                           Field<Real>& target) const
{
  const int num_dof_holders = target_layout_->GetNumOwnedDofHolder();
  const int num_components = target_layout_->GetNumComponents();

  // 1. Control variate: interpolate the source field onto the target space.
  interpolator_.Apply(source, control_variate_);

  // 2. Sample the source field and the control variate at the fixed Monte
  //    Carlo sample points; the stochastic RHS integrates the residual.
  const std::size_t num_samples = f_samples_.extent(0);
  // Shallow copies so the device lambda captures the views, not `this`.
  auto f_samples = f_samples_;
  auto residual = residual_;
  source_at_samples_->Evaluate(
    source, ValueView<Real, DeviceMemorySpace>(source.GetData().GetValueBasis(),
                                               MakeRank2View(f_samples)));
  control_variate_at_samples_->Evaluate(
    control_variate_,
    ValueView<Real, DeviceMemorySpace>(
      control_variate_.GetData().GetValueBasis(), MakeRank2View(residual)));
  Kokkos::parallel_for(
    "cv_residual", Kokkos::RangePolicy<DefaultExecutionSpace>(0, num_samples),
    KOKKOS_LAMBDA(int i) {
      residual(i, 0) = f_samples(i, 0) - residual(i, 0);
    });

  // 3. Project the residual: M * delta = r.
  const auto delta = solver_->Solve(MakeConstRank2View(residual));

  // 4. x = g + delta. delta is indexed by active PETSc row, while the field is
  //    indexed by local DOF holder.
  const auto g_nodal = control_variate_.GetDOFHolderData().GetValues();
  const auto global_to_local = target_layout_->GetGlobalToLocalPermutation();

  if (static_cast<int>(g_nodal.extent(0)) != num_dof_holders ||
      static_cast<int>(g_nodal.extent(1)) != num_components) {
    throw pcms_error("OmegaHControlVariateProjection::Apply: control variate "
                     "extents do not match the target layout");
  }

  const auto values = MakeRank2View(target_values_);
  Kokkos::parallel_for(
    "cv_add_correction",
    Kokkos::RangePolicy<DefaultExecutionSpace>(0, num_dof_holders),
    KOKKOS_LAMBDA(int i) {
      for (int c = 0; c < num_components; ++c) {
        values(i, c) =
          g_nodal(i, c) + delta[global_to_local(i) * num_components + c];
      }
    });
  Kokkos::fence();
  target.SetDOFHolderDataUnchecked(MakeConstRank2View(target_values_));
}

} // namespace pcms
