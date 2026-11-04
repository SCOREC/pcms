#ifndef PCMS_TRANSFER_OMEGA_H_CONSERVATIVE_PROJECTION_HPP
#define PCMS_TRANSFER_OMEGA_H_CONSERVATIVE_PROJECTION_HPP

#include "pcms/field/function_space.hpp"
#include "pcms/backends/omega_h/lagrange_layout.hpp"
#include "pcms/field/point_evaluator.hpp"
#include "pcms/transfer/mass_matrix_type.hpp"
#include "pcms/transfer/transfer_operator.hpp"
#include <Kokkos_Core.hpp>
#include <memory>

namespace pcms::transfer
{

class GalerkinProjectionSolver;
class OmegaHIntersectionRHSIntegrator;

// Conservative Galerkin projection between Omega_h order-1 Lagrange spaces.
//
// Construction (one-time cost):
//   - Computes mesh intersections and quadrature data
//   (OmegaHIntersectionRHSIntegrator).
//   - Localizes integration points in the source mesh (PointEvaluator).
//   - Assembles the target mass matrix and factors it via KSP
//     (GalerkinProjectionSolver). With MassMatrixType::Lumped the mass
//     matrix is row-sum lumped and the KSP degenerates to an exact diagonal
//     scaling; conservation is unaffected, but linear fields are no longer
//     reproduced exactly on order-1 targets.
//
// Apply() (per-call cost):
//   - Evaluates the source field at the fixed integration points.
//   - Assembles a fresh RHS vector.
//   - Calls KSPSolve reusing the cached factorization.
//   No mesh intersection, no matrix assembly, no refactorization per call.
class OmegaHConservativeProjection : public TransferOperator<Real>
{
public:
  OmegaHConservativeProjection(
    const FunctionSpace& source_space, const FunctionSpace& target_space,
    MassMatrixType mass_matrix_type = MassMatrixType::Consistent);

  // Defined in the .cpp where GalerkinProjectionSolver is a complete type.
  ~OmegaHConservativeProjection() override;

  void Apply(const Field<Real>& source, Field<Real>& target) const override;

private:
  std::shared_ptr<const pcms::backend::omega_h::OmegaHLagrangeLayout> source_layout_;
  std::shared_ptr<const pcms::backend::omega_h::OmegaHLagrangeLayout> target_layout_;
  std::unique_ptr<OmegaHIntersectionRHSIntegrator> rhs_integrator_;
  std::unique_ptr<PointEvaluator<Real>> evaluator_;
  std::unique_ptr<GalerkinProjectionSolver> solver_;
  mutable Kokkos::View<Real**, DeviceMemorySpace> target_values_;
};

} // namespace pcms::transfer

#endif // PCMS_TRANSFER_OMEGA_H_CONSERVATIVE_PROJECTION_HPP
