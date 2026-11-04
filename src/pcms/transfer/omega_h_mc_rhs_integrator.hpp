#ifndef PCMS_TRANSFER_OMEGA_H_MC_RHS_INTEGRATOR_HPP
#define PCMS_TRANSFER_OMEGA_H_MC_RHS_INTEGRATOR_HPP

#include "pcms/utility/coordinate_view.hpp"
#include "pcms/backends/omega_h/geometry.hpp"
#include "pcms/field/function_space.hpp"
#include "pcms/backends/omega_h/lagrange_layout.hpp"
#include "pcms/transfer/linear_form_integrator.hpp"
#include <Kokkos_Core.hpp>
#include <cstdint>
#include <memory>

namespace pcms
{

namespace detail
{
struct OmegaHScalarSpace;
} // namespace detail

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
// at construction and exposed through GetIntegrationPoints(), with the target
// element each lies in (IntegrationPointOwners(), for binding the target space
// without a search);
// |T| comes from the target geometry's differential. Assemble() consumes
// source-field values sampled there, exactly like
// OmegaHIntersectionRHSIntegrator, so this drops into
// GalerkinProjectionSolver::Solve.
class OmegaHMonteCarloRHSIntegrator : public LinearFormIntegrator
{
public:
  /// @param target_space order-1 scalar Cartesian Lagrange space on an Omega_h
  ///        geometry
  OmegaHMonteCarloRHSIntegrator(const FunctionSpace& target_space,
                                int samples_per_element,
                                uint64_t seed = 8675309);
  ~OmegaHMonteCarloRHSIntegrator();

  OmegaHMonteCarloRHSIntegrator(const OmegaHMonteCarloRHSIntegrator&) = delete;
  OmegaHMonteCarloRHSIntegrator& operator=(
    const OmegaHMonteCarloRHSIntegrator&) = delete;

  Vec GetVector() const noexcept override;

  CoordinateView<DeviceMemorySpace> GetIntegrationPoints()
    const noexcept override;

  /// @return the target element each sample lies in
  [[nodiscard]] OwnerHints IntegrationPointOwners() const;

  void Assemble(
    Rank2View<const Real, DeviceMemorySpace> sampled_values) override;

private:
  OmegaHMonteCarloRHSIntegrator(detail::OmegaHScalarSpace target,
                                int samples_per_element, uint64_t seed);

  Vec vec_ = nullptr;
  std::shared_ptr<const OmegaHGeometry> geometry_;
  Kokkos::View<LO*, DeviceMemorySpace> elements_; // target element per sample
  Kokkos::View<Real**, DeviceMemorySpace>
    coords_; // [num_pts][dim] sample coordinates
  Kokkos::View<PetscInt*, DeviceMemorySpace> node_gids_; // COO indices
  Kokkos::View<Real*, DeviceMemorySpace> coeffs_;        // basis * |T| / N
  int nbary_ = 3; // barycentric DOFs per sample (Dim+1): 3 in 2D, 4 in 3D
};

} // namespace pcms

#endif // PCMS_TRANSFER_OMEGA_H_MC_RHS_INTEGRATOR_HPP
