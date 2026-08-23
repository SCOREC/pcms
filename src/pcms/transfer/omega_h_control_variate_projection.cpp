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
  const FunctionSpace& source_space, const FunctionSpace& target_space,
  int samples_per_element, MonteCarloSampling sampling, uint64_t seed,
  MassMatrixType mass_matrix_type)
  : target_layout_(std::dynamic_pointer_cast<const OmegaHLagrangeLayout>(
      target_space.GetLayout())),
    rhs_integrator_(std::make_unique<OmegaHMonteCarloRHSIntegrator>(
      target_layout_, target_space.GetCoordinateSystem(), samples_per_element,
      sampling, seed)),
    interpolator_(source_space, target_space),
    control_variate_(target_space.CreateFunction<Real>())
{
  const auto sample_coords = rhs_integrator_->GetIntegrationPoints();
  source_at_samples_ = source_space.CreatePointEvaluator<Real>(
    EvaluationRequest::FromCoordinates(sample_coords));
  control_variate_at_samples_ = target_space.CreatePointEvaluator<Real>(
    EvaluationRequest::FromCoordinates(sample_coords));

  // Mass integrator is only needed to build the solver; PETSc reference-counts
  // the matrix so it remains alive inside the KSP after this scope ends.
  OmegaHMassIntegrator mass_integrator(
    target_layout_, target_space.GetCoordinateSystem(), mass_matrix_type);
  solver_ = std::make_unique<GalerkinProjectionSolver>(mass_integrator,
                                                       *rhs_integrator_);
}

// Defined here so that GalerkinProjectionSolver (forward-declared in the
// header) is a complete type when unique_ptr's destructor is instantiated.
OmegaHControlVariateProjection::~OmegaHControlVariateProjection() = default;

void OmegaHControlVariateProjection::Apply(const Field<Real>& source,
                                           Field<Real>& target) const
{
  // 1. Control variate: interpolate the source field onto the target space.
  interpolator_.Apply(source, control_variate_);

  // 2. Sample the source field and the control variate at the fixed Monte
  //    Carlo sample points; the stochastic RHS integrates the residual.
  const std::size_t num_samples =
    rhs_integrator_->GetIntegrationPoints().GetValues().extent(0);
  Kokkos::View<Real**, DeviceMemorySpace> f_samples("cv_f_samples", num_samples,
                                                    1);
  Kokkos::View<Real**, DeviceMemorySpace> residual("cv_residual", num_samples,
                                                   1);
  source_at_samples_->Evaluate(source, MakeRank2View(f_samples));
  control_variate_at_samples_->Evaluate(control_variate_,
                                        MakeRank2View(residual));
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
  const int nverts = static_cast<int>(g_nodal.extent(0));

  Kokkos::View<Real**, DeviceMemorySpace> result("cv_result", nverts, 1);
  Kokkos::parallel_for(
    "cv_add_correction", Kokkos::RangePolicy<DefaultExecutionSpace>(0, nverts),
    KOKKOS_LAMBDA(int i) {
      result(i, 0) = g_nodal(i, 0) + delta[global_to_local(i)];
    });
  Kokkos::fence();

  target.SetDOFHolderDataUnchecked(MakeConstRank2View(result));
}

} // namespace pcms
