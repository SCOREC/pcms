#include "omega_h_field.h"
#include "pcms/adapter/omega_h/omega_h_field_layout.h"
#include "omega_h_field_layout.h"
#include "pcms/inclusive_scan.h"
#include "pcms/profile.h"

namespace pcms
{
/*
 * Field Layout
 */

OmegaHFieldLayout::OmegaHFieldLayout(Omega_h::Mesh& mesh,
                                     OmegaHFieldLayoutLocation location,
                                     int num_components,
                                     std::string global_id_name)
  : mesh_(mesh), global_id_name_(global_id_name), location_(location), num_components_(num_components)
{
}

int OmegaHFieldLayout::GetNumComponents() const
{
  return num_components_;
}
// nodes for standard lagrange FEM
LO OmegaHFieldLayout::GetNumOwnedDofHolder() const
{
  switch (location_) {
    case OmegaHFieldLayoutLocation::PieceWise:
      // global ids are faces in 2D and regions in 3D
      return mesh_.nents(mesh_.dim());
      break;
    case OmegaHFieldLayoutLocation::Linear:
      // global ids are vertices
      return mesh_.nents(0);
      break;
  }
  // should never reach here, but return 0 to avoid compiler warning
  return 0;
}

GO OmegaHFieldLayout::GetNumGlobalDofHolder() const
{
  switch (location_) {
    case OmegaHFieldLayoutLocation::PieceWise:
      // global ids are faces in 2D and regions in 3D
      return mesh_.nglobal_ents(mesh_.dim());
      break;
    case OmegaHFieldLayoutLocation::Linear:
      // global ids are vertices
      return mesh_.nglobal_ents(0);
      break;
  }
  // should never reach here, but return 0 to avoid compiler warning
  return 0;
}

GlobalIDView<HostMemorySpace> OmegaHFieldLayout::GetOwnedGids()
{
  PCMS_FUNCTION_TIMER;
  int dim = 0;
  switch (location_) {
    case OmegaHFieldLayoutLocation::PieceWise:
      // global ids are faces in 2D and regions in 3D
      dim = mesh_.dim();
      break;
    case OmegaHFieldLayoutLocation::Linear:
      // ids are just global ids of the vertex
      dim = 0;
      break;
  }
  auto gids = mesh_.globals(dim);
  static_assert(std::is_same_v<HostMemorySpace, DefaultExecutionSpace::memory_space>, "types must match");
  return GlobalIDView<HostMemorySpace>(gids.data(), gids.size());
}

GlobalIDView<HostMemorySpace> OmegaHFieldLayout::GetGids() const
{
  PCMS_FUNCTION_TIMER;
  Omega_h::Read<Omega_h::GO> gid_array;
  if (global_id_name_.empty()) {
    gid_array = mesh_.globals(0);
  } else {
    auto tag = mesh_.get_tagbase(0, global_id_name_);
    if (Omega_h::is<GO>(tag)) {
      gid_array = mesh_.get_array<Omega_h::GO>(0, global_id_name_);
    } else if (Omega_h::is<LO>(tag)) {
      auto array = mesh_.get_array<Omega_h::LO>(0, global_id_name_);
      Omega_h::Write<Omega_h::GO> globals(array.size());
      Omega_h::parallel_for(
        array.size(), OMEGA_H_LAMBDA(int i) { globals[i] = array[i]; });
      gid_array = Omega_h::Read(globals);
    } else {
      std::cerr << "Weird tag type for global arrays.\n";
      std::abort();
    }
  }
  return GlobalIDView<HostMemorySpace>(gid_array.data(), gid_array.size());
}

OmegaHFieldLayoutLocation OmegaHFieldLayout::GetLocation() const
{
  return location_;
}

bool OmegaHFieldLayout::IsDistributed()
{
  return true;
}

Omega_h::Read<Omega_h::ClassId> OmegaHFieldLayout::GetClassIDs() const
{
  PCMS_FUNCTION_TIMER;
  return mesh_.get_array<Omega_h::ClassId>(0, "class_id");
}

Omega_h::Read<Omega_h::I8> OmegaHFieldLayout::GetClassDims() const
{
  PCMS_FUNCTION_TIMER;
  return mesh_.get_array<Omega_h::I8>(0, "class_dim");
}

ReversePartitionMap OmegaHFieldLayout::GetReversePartitionMap(
  const redev::Partition& partition) const
{
  PCMS_FUNCTION_TIMER;
  auto classIds_h = Omega_h::HostRead<Omega_h::ClassId>(GetClassIDs());
  auto classDims_h = Omega_h::HostRead<Omega_h::I8>(GetClassDims());
  const auto coords = Omega_h::HostRead(Omega_h::get_ent_centroids(mesh_, 0));
  auto dim = mesh_.dim();

  // local_index number of vertices going to each destination process by
  // calling getRank - degree array
  std::array<pcms::Real, 3> coord;
  pcms::ReversePartitionMap reverse_partition;
  pcms::LO local_index = 0;
  for (auto i = 0; i < classIds_h.size(); i++) {
    coord[0] = coords[i * dim];
    coord[1] = coords[i * dim + 1];
    coord[2] = (dim == 3) ? coords[i * dim + 2] : 0.0;
    auto dr = std::visit(GetRank{classIds_h[i], classDims_h[i], coord}, partition);
    reverse_partition[dr].emplace_back(local_index++);
  }
  return reverse_partition;
}

} // namespace pcms
