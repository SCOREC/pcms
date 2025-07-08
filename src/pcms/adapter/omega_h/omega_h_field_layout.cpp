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
    global_id_name_(global_id_name),
    dof_holder_coords_("", GetNumOwnedDofHolder(), mesh_.dim())
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

  auto coords = mesh_.coords();

  size_t offset = 0;
  for (int i = 0; i <= mesh_.dim(); ++i) {
    if (nodes_per_dim[i] == 1) {
      if (i == 0) {
        Kokkos::parallel_for(
          mesh_.nents(0), KOKKOS_LAMBDA(LO i) {
            dof_holder_coords_(offset + i, 0) = coords[2 * i + 0];
            dof_holder_coords_(offset + i, 1) = coords[2 * i + 1];
          });
      } else if (i == 1) {
        auto edge_verts = mesh_.ask_verts_of(1);
        Kokkos::parallel_for(
          mesh_.nents(1), KOKKOS_LAMBDA(LO i) {
            auto verts = Omega_h::gather_verts<2>(edge_verts, i);
            Real x0 = coords[2 * verts[0] + 0];
            Real y0 = coords[2 * verts[0] + 1];
            Real x1 = coords[2 * verts[1] + 0];
            Real y1 = coords[2 * verts[1] + 1];
            dof_holder_coords_(offset + i, 0) = (x0 + x1) / 2;
            dof_holder_coords_(offset + i, 1) = (y0 + y1) / 2;
          });
      } else {
        std::cerr << "Unsupported" << std::endl;
        std::abort();
      }
    } else if (nodes_per_dim[i] != 0) {
      std::cerr << "Unsupported" << std::endl;
      std::abort();
    }

    offset += mesh.nents(i);
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

Rank2View<const Real, HostMemorySpace>
OmegaHFieldLayout::GetDOFHolderCoordinates() const
{
  return Rank2View<const Real, HostMemorySpace>(
    dof_holder_coords_.data(), dof_holder_coords_.extent(0), 2);
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
