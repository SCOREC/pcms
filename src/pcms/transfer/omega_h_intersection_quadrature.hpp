#ifndef PCMS_TRANSFER_OMEGA_H_INTERSECTION_QUADRATURE_HPP
#define PCMS_TRANSFER_OMEGA_H_INTERSECTION_QUADRATURE_HPP

#include <memory>

#include <Kokkos_Core.hpp>
#include <petscsystypes.h>

#include "pcms/field/coordinate_system.h"
#include "pcms/field/function_space.h"
#include "pcms/field/layout/omega_h_lagrange.h"
#include "pcms/field/point_evaluator.h"
#include "pcms/transfer/mesh_intersection.hpp"
#include "pcms/utility/types.h"

namespace pcms
{

/// Quadrature of the source-target mesh intersection for a pair of Omega_h
/// Lagrange spaces: integration points on the target mesh, the target DOFs and
/// weights each point contributes to, the source element each point lies in,
/// and a source evaluator bound to those points that needs no point search.
/// Built once per space pair and shared by every operator on that pair.
class OmegaHIntersectionQuadrature
{
public:
  /// @param intersection mesh intersection to reuse; computed when null
  OmegaHIntersectionQuadrature(
    const FunctionSpace& source_space, const FunctionSpace& target_space,
    std::shared_ptr<const MeshIntersection> intersection = nullptr);
  /// Layout-only form without a source evaluator (GetSourceEvaluator() is
  /// null); `source_search` and `intersection` may be null.
  OmegaHIntersectionQuadrature(
    std::shared_ptr<const OmegaHLagrangeLayout> source_layout,
    CoordinateSystem source_coordinate_system,
    std::shared_ptr<const OmegaHLagrangeLayout> target_layout,
    CoordinateSystem target_coordinate_system,
    const GridPointSearchVariant* source_search = nullptr,
    const MeshIntersection* intersection = nullptr);

  CoordinateView<DeviceMemorySpace> GetIntegrationPoints() const noexcept;
  LO GetNumPoints() const noexcept;
  Kokkos::View<const LO*, DeviceMemorySpace> GetSourceElementIds()
    const noexcept;
  /// Target DOF (PETSc row) at [point * GetDofsPerElement() + k].
  Kokkos::View<const PetscInt*, DeviceMemorySpace> GetTargetDofs()
    const noexcept;
  /// Weight at [point * GetDofsPerElement() + k], basis times measure.
  Kokkos::View<const Real*, DeviceMemorySpace> GetWeights() const noexcept;
  int GetDofsPerElement() const noexcept;
  PetscInt GetNumTargetDofs() const noexcept;
  const PointEvaluator<Real>* GetSourceEvaluator() const noexcept;
  const std::shared_ptr<const OmegaHLagrangeLayout>& GetSourceLayout()
    const noexcept;
  const std::shared_ptr<const OmegaHLagrangeLayout>& GetTargetLayout()
    const noexcept;

private:
  std::shared_ptr<const OmegaHLagrangeLayout> source_layout_;
  std::shared_ptr<const OmegaHLagrangeLayout> target_layout_;
  Kokkos::View<Real**, DeviceMemorySpace> coords_;
  Kokkos::View<LO*, DeviceMemorySpace> source_elem_ids_;
  Kokkos::View<PetscInt*, DeviceMemorySpace> target_dofs_;
  Kokkos::View<Real*, DeviceMemorySpace> weights_;
  int dofs_per_element_ = 0;
  PetscInt num_target_dofs_ = 0;
  std::unique_ptr<PointEvaluator<Real>> source_evaluator_;
};

} // namespace pcms

#endif // PCMS_TRANSFER_OMEGA_H_INTERSECTION_QUADRATURE_HPP
