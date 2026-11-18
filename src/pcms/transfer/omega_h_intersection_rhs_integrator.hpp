#ifndef PCMS_TRANSFER_OMEGA_H_INTERSECTION_RHS_INTEGRATOR_HPP
#define PCMS_TRANSFER_OMEGA_H_INTERSECTION_RHS_INTEGRATOR_HPP

#include "pcms/utility/coordinate_view.hpp"
#include "pcms/field/function_space.hpp"
#include "pcms/backends/omega_h/lagrange_layout.hpp"
#include "pcms/backends/omega_h/geometry.hpp"
#include "pcms/geometry/bound_source.hpp"
#include "pcms/transfer/linear_form_integrator.hpp"
#include "pcms/backends/omega_h/mesh_intersection.hpp"
#include <Kokkos_Core.hpp>
#include <memory>

namespace pcms::transfer
{

// Conservative-projection RHS integrator for Lagrange spaces on Omega_h 2D
// simplex meshes. Source and target orders are handled independently: the
// source enters only through pointwise samples (any supported order), while the
// target order selects the basis and DOF layout (P0: one element-constant DOF;
// P1: three vertex DOFs). The quadrature order is source_order + target_order.
//
// Construction validates both spaces (scalar Cartesian Lagrange on simplex
// meshes of the same dimension), computes mesh intersections, and runs a
// two-pass device kernel (count then fill) to build the quadrature data and
// allocate the owned RHS vector sized to the target DOF count. The source space
// must be on the Omega_h geometry of its layout's mesh: the intersection reuses
// that geometry's point search.
//
// The integration points are generated in source elements, so the source is
// bound to them without a search:
// source_geometry->Bind(GetIntegrationPoints(), IntegrationPointOwners()).
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
  ~OmegaHIntersectionRHSIntegrator();

  Vec GetVector() const noexcept override;

  CoordinateView<DeviceMemorySpace> GetIntegrationPoints()
    const noexcept override;

  /// @return the source element each integration point lies in
  [[nodiscard]] OwnerHints IntegrationPointOwners() const;

  void Assemble(
    Rank2View<const Real, DeviceMemorySpace> sampled_values) override;

private:
  Vec vec_ = nullptr;
  std::shared_ptr<const backend::omega_h::OmegaHGeometry> source_geometry_;
  Kokkos::View<Real**, DeviceMemorySpace>
    coords_; // [num_pts][dim] integration point coordinates
  // source element of each integration point
  Kokkos::View<LO*, DeviceMemorySpace> source_elements_;
  Kokkos::View<PetscInt*, DeviceMemorySpace>
    node_gids_;                                   // COO idx, num_pts*ndof
  Kokkos::View<Real*, DeviceMemorySpace> coeffs_; // weights, num_pts*ndof
  int ndof_per_elem_ = 0;                         // target DOFs/element
};

// Builds an OmegaHIntersectionRHSIntegrator for conservative L2 projection.
// Both spaces must use OmegaHLagrangeLayout, scalar, order-0 or order-1 (source
// and target orders may differ), with the same coordinate system.
std::unique_ptr<OmegaHIntersectionRHSIntegrator>
BuildOmegaHConservativeRHSIntegrator(const FunctionSpace& source_space,
                                     const FunctionSpace& target_space);

} // namespace pcms::transfer

#endif // PCMS_TRANSFER_OMEGA_H_INTERSECTION_RHS_INTEGRATOR_HPP
