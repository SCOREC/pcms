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

template <typename T>
Omega_h::Write<Omega_h::GO> GetGidsHelper(LO total_ents,
                                          std::array<int, 4> nodes_per_dim,
                                          Omega_h::Mesh& mesh,
                                          const std::string& global_id_name)
{
  PCMS_FUNCTION_TIMER;
  static_assert(
    std::is_same_v<HostMemorySpace, DefaultExecutionSpace::memory_space>,
    "types must match");

  Omega_h::Write<Omega_h::GO> owned_gids(total_ents);
  LO offset = 0;
  for (int i = 0; i <= mesh.dim(); ++i) {
    if (nodes_per_dim[i]) {
      auto dim_gids = mesh.get_array<T>(i, global_id_name);
      Omega_h::parallel_for(
        dim_gids.size(),
        OMEGA_H_LAMBDA(int i) { owned_gids[i + offset] = dim_gids[i]; });
      offset += dim_gids.size();
    }
  }

  return owned_gids;
}

OmegaHFieldLayout::OmegaHFieldLayout(Omega_h::Mesh& mesh,
                                     std::array<int, 4> nodes_per_dim,
                                     int num_components,
                                     std::string global_id_name)
  : mesh_(mesh),
    nodes_per_dim_(nodes_per_dim),
    num_components_(num_components),
    global_id_name_(global_id_name)
{
  LO total_ents = 0;
  for (int i = 0; i <= mesh.dim(); ++i) {
    if (nodes_per_dim[i])
      total_ents += mesh.nglobal_ents(i);
  }

  auto tag = mesh_.get_tagbase(0, global_id_name_);
  if (Omega_h::is<GO>(tag)) {
    gids_ = GetGidsHelper<GO>(total_ents, nodes_per_dim, mesh, global_id_name);
  } else if (Omega_h::is<LO>(tag)) {
    gids_ = GetGidsHelper<LO>(total_ents, nodes_per_dim, mesh, global_id_name);
  } else {
    std::cerr << "Weird tag type for global arrays.\n";
    std::abort();
  }
}

int OmegaHFieldLayout::GetNumComponents() const
{
  return num_components_;
}
// nodes for standard lagrange FEM
LO OmegaHFieldLayout::GetNumOwnedDofHolder() const
{
  LO count = 0;
  for (int i = 0; i <= mesh_.dim(); ++i) {
    count += mesh_.nents(i) * nodes_per_dim_[i];
  }
  return count;
}

GO OmegaHFieldLayout::GetNumGlobalDofHolder() const
{
  LO count = 0;
  for (int i = 0; i <= mesh_.dim(); ++i) {
    count += mesh_.nglobal_ents(i) * nodes_per_dim_[i];
  }
  return count;
}

std::array<int, 4> OmegaHFieldLayout::GetNodesPerDim() const
{
  return nodes_per_dim_;
}

GlobalIDView<HostMemorySpace> OmegaHFieldLayout::GetOwnedGids() const
{
  static_assert(std::is_same_v<HostMemorySpace, DefaultExecutionSpace::memory_space>, "types must match");
  return GlobalIDView<HostMemorySpace>(gids_.data(), gids_.size());
}

GlobalIDView<HostMemorySpace> OmegaHFieldLayout::GetGids() const
{
  return GlobalIDView<HostMemorySpace>(gids_.data(), gids_.size());
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
