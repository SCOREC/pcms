#ifndef PCMS_FIELD_LAYOUT_MFEM_H
#define PCMS_FIELD_LAYOUT_MFEM_H

#include "pcms/field/coordinate_system.h"
#include "pcms/field/field_layout.h"

#include <mfem.hpp>

#include <memory>

namespace pcms
{

// Field layout for an MFEM order-1 H1 (vertex) scalar field.
//
// Each DOF holder is a mesh vertex; DOF-holder ordering follows the MFEM local
// vertex ordering. This is the minimal layout matching the original MFEM
// adapter scope: vertex-only, order-1, single component.
//
// Partitioning relies on DOF-holder coordinates (redev::RCBPtn). The
// classification arrays required by FieldLayout are filled with placeholders
// (entity dim 0, id = local index) because MFEM has no Omega_h-style geometric
// classification; an RCB partition ignores them.
class MFEMLayout : public FieldLayout
{
public:
  MFEMLayout(mfem::ParMesh& pmesh, mfem::ParFiniteElementSpace& pfes,
             CoordinateSystem coordinate_system);

  std::shared_ptr<const Discretization> GetDiscretization()
    const noexcept override;

  int GetNumComponents() const override;
  LO GetNumOwnedDofHolder() const override;
  GO GetNumGlobalDofHolder() const override;

  Rank1View<const bool, HostMemorySpace> GetOwnedHost() const override;
  GlobalIDView<HostMemorySpace> GetGidsHost() const override;
  CoordinateView<DeviceMemorySpace> GetDOFHolderCoordinates() const override;

  [[nodiscard]] bool IsDistributed() const override;
  EntOffsetsArray GetEntOffsets() const override;
  int GetDimension() const override;

  Rank1View<const LO, HostMemorySpace>
  GetDOFHolderClassificationDimensionsHost() const override;

  Rank1View<const LO, HostMemorySpace> GetDOFHolderClassificationIdsHost()
    const override;

  mfem::ParMesh& GetMesh() const noexcept { return pmesh_; }
  mfem::ParFiniteElementSpace& GetFESpace() const noexcept { return pfes_; }

private:
  void AssertVertexScalarSpace() const;

  mfem::ParMesh& pmesh_;
  mfem::ParFiniteElementSpace& pfes_;
  int dim_;
  CoordinateSystem coordinate_system_;

  Kokkos::View<Real**, DeviceMemorySpace> coords_;
  Kokkos::View<bool*, HostMemorySpace> owned_host_;
  Kokkos::View<GO*, HostMemorySpace> gids_host_;
  Kokkos::View<LO*, HostMemorySpace> classification_dims_host_;
  Kokkos::View<LO*, HostMemorySpace> classification_ids_host_;
  std::shared_ptr<const Discretization> discretization_;
};

} // namespace pcms

#endif // PCMS_FIELD_LAYOUT_MFEM_H
