#ifndef PCMS_OMEGA_H_FIELD_LAYOUT_H
#define PCMS_OMEGA_H_FIELD_LAYOUT_H

#include <Omega_h_mesh.hpp>

#include "pcms/arrays.h"
#include "pcms/field_layout.h"
#include "pcms/coordinate_system.h"
#include "pcms/field.h"

namespace pcms {

enum class OmegaHFieldLayoutLocation {
  PieceWise,
  Linear,
};

class OmegaHFieldLayout : public FieldLayout {
public:
  OmegaHFieldLayout(Omega_h::Mesh& mesh, OmegaHFieldLayoutLocation location, int num_components);
  int GetNumComponents() override;
  // nodes for standard lagrange FEM
  LO GetNumOwnedDofHolder() override;
  GO GetNumGlobalDofHolder() override;

  GlobalIDView<HostMemorySpace> GetOwnedGids() override;

  // returns true if the field layout is distributed
  // if the field layout is distributed, the owned and global dofs are the same
  bool IsDistributed() override;

private:
  Omega_h::Mesh& mesh_;
  OmegaHFieldLayoutLocation location_;
  int num_components_;
};

}
#endif // PCMS_OMEGA_H_FIELD_LAYOUT_H
