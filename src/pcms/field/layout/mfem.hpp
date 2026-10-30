#ifndef PCMS_FIELD_LAYOUT_MFEM_HPP
#define PCMS_FIELD_LAYOUT_MFEM_HPP

#include "pcms/field/coordinate_system.hpp"
#include "pcms/field/field_layout.h"

#include <mfem.hpp>

#include <memory>

namespace pcms
{

// Field layout for an MFEM order-1 H1 (vertex) scalar field.
class MFEMLayout : public FieldLayout
{
public:
  MFEMLayout(mfem::ParMesh& pmesh, mfem::ParFiniteElementSpace& pfes,
             std::shared_ptr<const CoordinateSystem> coordinate_system);

  std::shared_ptr<const Discretization> GetDiscretization()
    const noexcept override;

  int GetNumComponents() const override;
  LO GetNumLocalDofHolder() const override;
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

  // Build a per-vertex overlap mask from MFEM element attributes
  static Kokkos::View<bool*, HostMemorySpace> OverlapMaskFromAttribute(
    mfem::ParMesh& pmesh, int attribute);

private:
  void AssertVertexScalarSpace() const;

  mfem::ParMesh& pmesh_;
  mfem::ParFiniteElementSpace& pfes_;
  int dim_;

  Kokkos::View<Real**, DeviceMemorySpace> coords_;
  Kokkos::View<bool*, HostMemorySpace> owned_host_;
  Kokkos::View<GO*, HostMemorySpace> gids_host_;
  Kokkos::View<LO*, HostMemorySpace> classification_dims_host_;
  Kokkos::View<LO*, HostMemorySpace> classification_ids_host_;
  std::shared_ptr<const Discretization> discretization_;
};

} // namespace pcms

#endif // PCMS_FIELD_LAYOUT_MFEM_HPP
