#ifndef PCMS_TRANSFER_OMEGA_H_MC_RHS_INTEGRATOR_HPP
#define PCMS_TRANSFER_OMEGA_H_MC_RHS_INTEGRATOR_HPP

#include "pcms/field/coordinate_view.hpp"
#include "pcms/field/function_space.h"
#include "pcms/field/layout/omega_h_lagrange.h"
#include "pcms/transfer/linear_form_integrator.hpp"
#include "pcms/transfer/monte_carlo_sampling.hpp"
#include <Kokkos_Core.hpp>
#include <cstdint>
#include <memory>

namespace pcms
{

// Monte Carlo RHS integrator for conservative L2 projection onto order-1
// Lagrange spaces on Omega_h simplex meshes (triangles in 2D, tetrahedra in
// 3D).
//
// Instead of intersection-based quadrature, each load-vector entry
//   b_j = \int phi_j f dx
// is estimated per target element T from a uniform sample x_1..x_N over T:
//   b_j ~= (|T| / N) * sum_i phi_j(x_i) f(x_i)
// Unit-square samples (pseudo-random) are mapped to each triangle
// with the area-preserving sqrt parametrization. Sample locations are fixed
// at construction and exposed through GetIntegrationPoints(); Assemble()
// consumes source-field values sampled there, exactly like
// OmegaHIntersectionRHSIntegrator, so this drops into
// GalerkinProjectionSolver::Solve.
class OmegaHMonteCarloRHSIntegrator : public LinearFormIntegrator
{
public:
  OmegaHMonteCarloRHSIntegrator(const FunctionSpace& target_space,
                                int samples_per_element,
                                MonteCarloSampling sampling,
                                uint64_t seed = 8675309);
  OmegaHMonteCarloRHSIntegrator(
    std::shared_ptr<const OmegaHLagrangeLayout> target_layout,
    std::shared_ptr<const CoordinateSystem> target_coordinate_system,
    int samples_per_element, MonteCarloSampling sampling,
    uint64_t seed = 8675309);
  ~OmegaHMonteCarloRHSIntegrator();

  OmegaHMonteCarloRHSIntegrator(const OmegaHMonteCarloRHSIntegrator&) = delete;
  OmegaHMonteCarloRHSIntegrator& operator=(
    const OmegaHMonteCarloRHSIntegrator&) = delete;

  Vec GetVector() const noexcept override;

  CoordinateView<DeviceMemorySpace> GetIntegrationPoints()
    const noexcept override;

  void Assemble(
    Rank2View<const Real, DeviceMemorySpace> sampled_values) override;

private:
  Vec vec_ = nullptr;
  Kokkos::View<Real**, DeviceMemorySpace>
    coords_; // [num_pts][dim] sample coordinates
  Kokkos::View<PetscInt*, DeviceMemorySpace> node_gids_; // COO indices
  Kokkos::View<Real*, DeviceMemorySpace> coeffs_;        // basis * |T| / N
  int nbary_ = 3; // barycentric DOFs per sample (Dim+1): 3 in 2D, 4 in 3D
};

} // namespace pcms

#endif // PCMS_TRANSFER_OMEGA_H_MC_RHS_INTEGRATOR_HPP
