#ifndef PCMS_OMEGA_H_FIELD_LAYOUT_H
#define PCMS_OMEGA_H_FIELD_LAYOUT_H

#include <Omega_h_mesh.hpp>

#include "pcms/arrays.h"
#include "pcms/field_layout.h"
#include "pcms/coordinate_system.h"
#include "pcms/field.h"

namespace pcms {
using ReversePartitionMap = std::map<pcms::LO, std::vector<pcms::LO>>;

enum class OmegaHFieldLayoutLocation {
  PieceWise,
  Linear,
};

class OmegaHFieldLayout : public FieldLayout {
public:
  OmegaHFieldLayout(Omega_h::Mesh& mesh, OmegaHFieldLayoutLocation location, int num_components, std::string global_id_name_ = "");
  int GetNumComponents() const override;
  // nodes for standard lagrange FEM
  LO GetNumOwnedDofHolder() const override;
  GO GetNumGlobalDofHolder() const override;

  GlobalIDView<HostMemorySpace> GetOwnedGids() override;
  GlobalIDView<HostMemorySpace> GetGids() const override;
  OmegaHFieldLayoutLocation GetLocation() const;

  // returns true if the field layout is distributed
  // if the field layout is distributed, the owned and global dofs are the same
  bool IsDistributed() override;

  ReversePartitionMap GetReversePartitionMap(
    const redev::Partition& partition) const override;

private:
  Omega_h::Read<Omega_h::ClassId> GetClassIDs() const;
  Omega_h::Read<Omega_h::I8> GetClassDims() const;

  Omega_h::Mesh& mesh_;
  std::string global_id_name_;
  OmegaHFieldLayoutLocation location_;
  int num_components_;
};

}
#endif // PCMS_OMEGA_H_FIELD_LAYOUT_H
