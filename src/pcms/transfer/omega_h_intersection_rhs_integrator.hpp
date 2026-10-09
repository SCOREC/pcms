#ifndef PCMS_TRANSFER_OMEGA_H_INTERSECTION_RHS_INTEGRATOR_HPP
#define PCMS_TRANSFER_OMEGA_H_INTERSECTION_RHS_INTEGRATOR_HPP

#include "pcms/field/coordinate_view.hpp"
#include "pcms/field/function_space.h"
#include "pcms/field/layout/omega_h_lagrange.h"
#include "pcms/transfer/linear_form_integrator.hpp"
#include "pcms/transfer/mesh_intersection.hpp"
#include <Kokkos_Core.hpp>
#include <memory>

namespace pcms
{

// Conservative-projection RHS integrator for Lagrange spaces on Omega_h 2D
// simplex meshes. Source and target orders are handled independently: the
// source enters only through pointwise samples (any supported order), while the
// target order selects the basis and DOF layout (P0: one element-constant DOF;
// P1: three vertex DOFs). The quadrature order is source_order + target_order.
//
// Construction validates both spaces, computes mesh intersections, and runs a
// two-pass device kernel (count then fill) to build the quadrature data and
// allocate the owned RHS vector sized to the target DOF count.
//
// Assemble(sampled_values):
//   Zeros the owned vector, then for each integration point i and local target
//   DOF j:
//     contribution = stored_coeff[i][j] * sampled_values(i, 0)
//   Scatter-adds via VecSetValuesCOO with ADD_VALUES.
// GetVector() returns the assembled vector (non-owning handle).
class OmegaHIntersectionRHSIntegrator : public LinearFormIntegrator
{
public:
  OmegaHIntersectionRHSIntegrator(const FunctionSpace& source_space,
                                  const FunctionSpace& target_space);
  OmegaHIntersectionRHSIntegrator(
    std::shared_ptr<const OmegaHLagrangeLayout> source_layout,
    std::shared_ptr<const CoordinateSystem> source_coordinate_system,
    std::shared_ptr<const OmegaHLagrangeLayout> target_layout,
    std::shared_ptr<const CoordinateSystem> target_coordinate_system);
  ~OmegaHIntersectionRHSIntegrator();

  Vec GetVector() const noexcept override;

  CoordinateView<DeviceMemorySpace> GetIntegrationPoints()
    const noexcept override;

  void Assemble(
    Rank2View<const Real, DeviceMemorySpace> sampled_values) override;

private:
  Vec vec_ = nullptr;
  Kokkos::View<Real**, DeviceMemorySpace>
    coords_; // [num_pts][dim] integration point coordinates
  Kokkos::View<PetscInt*, DeviceMemorySpace>
    node_gids_;                                   // COO idx, num_pts*ndof
  Kokkos::View<Real*, DeviceMemorySpace> coeffs_; // weights, num_pts*ndof
  int ndof_per_elem_ = 0;                         // target DOFs/element
};

// Builds an OmegaHIntersectionRHSIntegrator for conservative L2 projection.
// Both spaces must use OmegaHLagrangeLayout, scalar, Cartesian, order-0 or
// order-1 (source and target orders may differ).
std::unique_ptr<LinearFormIntegrator> BuildOmegaHConservativeRHSIntegrator(
  const FunctionSpace& source_space, const FunctionSpace& target_space);

} // namespace pcms

#endif // PCMS_TRANSFER_OMEGA_H_INTERSECTION_RHS_INTEGRATOR_HPP
