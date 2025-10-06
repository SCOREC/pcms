#ifndef PCMS_OMEGA_H_FIELD_LAYOUT_H
#define PCMS_OMEGA_H_FIELD_LAYOUT_H

#include <Omega_h_mesh.hpp>

#include "pcms/arrays.h"
#include "pcms/field_layout.h"
#include "pcms/coordinate_system.h"
#include "pcms/field.h"

#include <array>

namespace pcms
{
class OmegaHFieldLayout : public FieldLayout
{
public:
  OmegaHFieldLayout(Omega_h::Mesh& mesh, std::array<int, 4> nodes_per_dim,
                    int num_components, CoordinateSystem coordinate_system,
                    std::string global_id_name = "global");

  std::unique_ptr<FieldT<Real>> CreateField() const override;

  int GetNumComponents() const override;
  // nodes for standard lagrange FEM
  LO GetNumOwnedDofHolder() const override;
  GO GetNumGlobalDofHolder() const override;

  Rank1View<const bool, HostMemorySpace> GetOwned() const override;
  GlobalIDView<HostMemorySpace> GetGids() const override;
  CoordinateView<HostMemorySpace> GetDOFHolderCoordinates() const override;

  // returns true if the field layout is distributed
  // if the field layout is distributed, the owned and global dofs are the same
  bool IsDistributed() override;

  EntOffsetsArray GetEntOffsets() const override;

  ReversePartitionMap2 GetReversePartitionMap(
    const redev::Partition& partition) const override;

  std::array<int, 4> GetNodesPerDim() const;
  size_t GetNumEnts() const;
  Omega_h::Mesh& GetMesh() const;

private:
  Omega_h::Read<Omega_h::ClassId> GetClassIDs() const;
  Omega_h::Read<Omega_h::I8> GetClassDims() const;

  Omega_h::Mesh& mesh_;
  Omega_h::Write<Omega_h::GO> gids_;
  std::string global_id_name_;
  int num_components_;
  CoordinateSystem coordinate_system_;
  std::array<int, 4> nodes_per_dim_;
  Kokkos::View<Real**> dof_holder_coords_;
  Omega_h::Write<Omega_h::ClassId> class_ids_;
  Omega_h::Write<Omega_h::I8> class_dims_;
  Kokkos::View<bool*> owned_;
};

} // namespace pcms
#endif // PCMS_OMEGA_H_FIELD_LAYOUT_H
