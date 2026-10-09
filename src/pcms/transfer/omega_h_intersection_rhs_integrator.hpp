#ifndef PCMS_TRANSFER_OMEGA_H_INTERSECTION_RHS_INTEGRATOR_HPP
#define PCMS_TRANSFER_OMEGA_H_INTERSECTION_RHS_INTEGRATOR_HPP

#include "pcms/field/coordinate_system.h"
#include "pcms/field/function_space.h"
#include "pcms/field/layout/omega_h_lagrange.h"
#include "pcms/transfer/linear_form_integrator.hpp"
#include "pcms/transfer/mesh_intersection.hpp"
#include "pcms/transfer/omega_h_intersection_quadrature.hpp"
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
// allocate the owned RHS vector sized to the target DOF count. When built from
// function spaces and the source is an Omega_h Lagrange space, the intersection
// reuses the point search that space already owns; the layout constructor has
// no space to borrow from and builds its own.
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
  /// Reuses a mesh intersection built for the two spaces' discretizations.
  /// Throws pcms_error if it is not between the source and target meshes.
  OmegaHIntersectionRHSIntegrator(
    const FunctionSpace& source_space, const FunctionSpace& target_space,
    std::shared_ptr<const MeshIntersection> intersection);
  OmegaHIntersectionRHSIntegrator(
    std::shared_ptr<const OmegaHLagrangeLayout> source_layout,
    CoordinateSystem source_coordinate_system,
    std::shared_ptr<const OmegaHLagrangeLayout> target_layout,
    CoordinateSystem target_coordinate_system);
  /// Assembles against a quadrature shared with other operators on the same
  /// pair of spaces.
  explicit OmegaHIntersectionRHSIntegrator(
    std::shared_ptr<const OmegaHIntersectionQuadrature> quadrature);
  ~OmegaHIntersectionRHSIntegrator();

  const OmegaHIntersectionQuadrature& GetQuadrature() const noexcept;

  Vec GetVector() const noexcept override;

  CoordinateView<DeviceMemorySpace> GetIntegrationPoints()
    const noexcept override;

  void Assemble(
    Rank2View<const Real, DeviceMemorySpace> sampled_values) override;

private:
  std::shared_ptr<const OmegaHIntersectionQuadrature> quadrature_;
  Vec vec_ = nullptr;
};

// Builds an OmegaHIntersectionRHSIntegrator for conservative L2 projection.
// Both spaces must use OmegaHLagrangeLayout, scalar, Cartesian, order-0 or
// order-1 (source and target orders may differ).
std::unique_ptr<LinearFormIntegrator> BuildOmegaHConservativeRHSIntegrator(
  const FunctionSpace& source_space, const FunctionSpace& target_space);

} // namespace pcms

#endif // PCMS_TRANSFER_OMEGA_H_INTERSECTION_RHS_INTEGRATOR_HPP
