#ifndef PCMS_TRANSFER_OMEGA_H_CONTROL_VARIATE_PROJECTION_HPP
#define PCMS_TRANSFER_OMEGA_H_CONTROL_VARIATE_PROJECTION_HPP

#include "pcms/field/field.h"
#include "pcms/field/function_space.h"
#include "pcms/field/layout/omega_h_lagrange.h"
#include "pcms/field/point_evaluator.h"
#include "pcms/transfer/interpolator.h"
#include "pcms/transfer/mass_matrix_type.hpp"
#include "pcms/transfer/monte_carlo_sampling.hpp"
#include "pcms/transfer/transfer_operator.hpp"
#include <Kokkos_Core.hpp>
#include <cstdint>
#include <memory>

namespace pcms
{

class GalerkinProjectionSolver;
class OmegaHMonteCarloRHSIntegrator;

// Variance-reduced Monte Carlo Galerkin projection using a control variate.
//
// The control variate g is the interpolation of the source field onto the
// target space (Interpolator: evaluate the source field at the target nodal
// locations, interpolate with the target basis). Since g lives in the target
// space, its L2 projection is g itself, so only the residual f - g needs the
// stochastic RHS:
//   M x = M g + \int phi (f - g)  =>  x = g + M^{-1} r,
//   r_j ~= (|T| / N) * sum_i phi_j(x_i) (f(x_i) - g(x_i))
// Because g approximates f, the sampled residual has far smaller variance
// than sampling f directly (and is exactly zero when f is in the target
// space).
//
// Apply() per-call cost: one source-field evaluation at the target nodes
// (interpolation), one source-field and one target-field evaluation at the
// fixed sample points, RHS assembly, and a KSPSolve with the cached
// factorization.
class OmegaHControlVariateProjection : public TransferOperator<Real>
{
public:
  OmegaHControlVariateProjection(
    const FunctionSpace& source_space, const FunctionSpace& target_space,
    int samples_per_element, MonteCarloSampling sampling,
    uint64_t seed = 8675309,
    MassMatrixType mass_matrix_type = MassMatrixType::Consistent);

  // Defined in the .cpp where GalerkinProjectionSolver is a complete type.
  ~OmegaHControlVariateProjection() override;

  void Apply(const Field<Real>& source, Field<Real>& target) const override;

  void Apply(TransferKey, const Field<Real>& source,
             Rank2View<Real, DeviceMemorySpace> out) const override;

private:
  std::shared_ptr<const OmegaHLagrangeLayout> target_layout_;
  std::unique_ptr<OmegaHMonteCarloRHSIntegrator> rhs_integrator_;
  Interpolator<Real> interpolator_;
  // Scratch field on the target space holding the interpolated control
  // variate; overwritten on every Apply.
  mutable Field<Real> control_variate_;
  std::unique_ptr<PointEvaluator<Real>> source_at_samples_;
  std::unique_ptr<PointEvaluator<Real>> control_variate_at_samples_;
  std::unique_ptr<GalerkinProjectionSolver> solver_;

  mutable Kokkos::View<Real**, DeviceMemorySpace> target_values_;
  mutable Kokkos::View<Real**, DeviceMemorySpace> f_samples_;
  mutable Kokkos::View<Real**, DeviceMemorySpace> residual_;
};

} // namespace pcms

#endif // PCMS_TRANSFER_OMEGA_H_CONTROL_VARIATE_PROJECTION_HPP
