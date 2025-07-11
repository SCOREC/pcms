#ifndef PCMS_OMEGA_H_FIELD_LAYOUT_H
#define PCMS_OMEGA_H_FIELD_LAYOUT_H

#include <Omega_h_mesh.hpp>

#include "pcms/arrays.h"
#include "pcms/field_layout.h"
#include "pcms/coordinate_system.h"
#include "pcms/field.h"

#include <array>

namespace pcms {
class OmegaHFieldLayout : public FieldLayout {
public:
  OmegaHFieldLayout(Omega_h::Mesh& mesh, std::array<int, 4> nodes_per_dim, int num_components, std::string global_id_name = "global");
  int GetNumComponents() const override;
  // nodes for standard lagrange FEM
  LO GetNumOwnedDofHolder() const override;
  GO GetNumGlobalDofHolder() const override;

  GlobalIDView<HostMemorySpace> GetOwnedGids() const override;
  GlobalIDView<HostMemorySpace> GetGids() const override;
  Rank2View<const Real, HostMemorySpace> GetDOFHolderCoordinates() const;

  // returns true if the field layout is distributed
  // if the field layout is distributed, the owned and global dofs are the same
  bool IsDistributed() override;

  size_t GetNumEnts() const;
  std::array<size_t, 4> GetEntOffsets() const override;

  ReversePartitionMap2 GetReversePartitionMap(
    const redev::Partition& partition) const override;

  std::array<int, 4> GetNodesPerDim() const;

private:
  Omega_h::Read<Omega_h::ClassId> GetClassIDs() const;
  Omega_h::Read<Omega_h::I8> GetClassDims() const;

  Omega_h::Mesh& mesh_;
  Omega_h::Write<Omega_h::GO> gids_;
  std::string global_id_name_;
  int num_components_;
  std::array<int, 4> nodes_per_dim_;
  Kokkos::View<Real **> dof_holder_coords_;
  Omega_h::Write<Omega_h::ClassId> class_ids_;
  Omega_h::Write<Omega_h::I8> class_dims_;
};

}
#endif // PCMS_OMEGA_H_FIELD_LAYOUT_H
