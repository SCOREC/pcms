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

  PCMS_ALWAYS_ASSERT(offset == total_ents);

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
    dof_holder_coords_("", GetNumOwnedDofHolder(), mesh_.dim()),
    class_ids_(GetNumEnts()),
    class_dims_(class_ids_.size()),
    owned_(class_dims_.size())
{
  PCMS_FUNCTION_TIMER;
  LO total_ents = GetNumEnts();

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

  offset = 0;
  for (int i = 0; i <= mesh_.dim(); ++i) {
    if (nodes_per_dim_[i]) {
      auto ids = mesh_.get_array<Omega_h::ClassId>(i, "class_id");
      auto dims = mesh_.get_array<Omega_h::I8>(i, "class_dim");
      auto owned = mesh_.owned(i);
      PCMS_ALWAYS_ASSERT(ids.size() == dims.size() &&
                         dims.size() == mesh_.nents(i));

      Omega_h::parallel_for(
        mesh_.nents(i), OMEGA_H_LAMBDA(LO i) {
          class_ids_[offset + i] = ids[i];
          class_dims_[offset + i] = dims[i];
          owned_[offset + i] = owned[i];
        });
      offset += mesh.nents(i);
    }
  }
}

int OmegaHFieldLayout::GetNumComponents() const
{
  return num_components_;
}

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

Omega_h::Read<Omega_h::I8> OmegaHFieldLayout::GetOwned() const
{
  return owned_;
}

GlobalIDView<HostMemorySpace> OmegaHFieldLayout::GetGids() const
{
  static_assert(std::is_same_v<HostMemorySpace, DefaultExecutionSpace::memory_space>, "types must match");
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
  return Omega_h::Read(class_ids_);
}

Omega_h::Read<Omega_h::I8> OmegaHFieldLayout::GetClassDims() const
{
  PCMS_FUNCTION_TIMER;
  return Omega_h::Read(class_dims_);
}

size_t OmegaHFieldLayout::GetNumEnts() const {
  size_t n = 0;
  for (int i = 0; i <= mesh_.dim(); ++i) {
    if (nodes_per_dim_[i])
      n += mesh_.nents(i);
  }
  return n;
}

std::array<size_t, 5> OmegaHFieldLayout::GetEntOffsets() const
{
  std::array<size_t, 5> offsets{};
  size_t offset = 0;
  for (int i = 0; i < offsets.size(); ++i) {
    offsets[i] = offset;
    if (i <= mesh_.dim() && nodes_per_dim_[i])
      offset += mesh_.nents(i);
  }
  offsets[offsets.size() - 1] = offset;
  return offsets;
}

ReversePartitionMap2 OmegaHFieldLayout::GetReversePartitionMap(
  const redev::Partition& partition) const
{
  PCMS_FUNCTION_TIMER;
  auto classIds_h = Omega_h::HostRead<Omega_h::ClassId>(GetClassIDs());
  auto classDims_h = Omega_h::HostRead<Omega_h::I8>(GetClassDims());
  const auto coords = GetDOFHolderCoordinates();
  auto dim = mesh_.dim();

  PCMS_ALWAYS_ASSERT(classDims_h.size() == classIds_h.size() &&
                     classIds_h.size() == coords.extent(0));

  // local_index number of vertices going to each destination process by
  // calling getRank - degree array
  std::array<pcms::Real, 3> coord;
  ReversePartitionMap2 reverse_partition;
  LO local_index = 0;
  for (int ent_dim = 0; ent_dim <= mesh_.dim(); ++ent_dim) {
    if (nodes_per_dim_[ent_dim] == 0)
      continue;

    for (LO i = 0; i < mesh_.nents(ent_dim); ++i) {
      coord[0] = coords(local_index, 0);
      coord[1] = coords(local_index, 1);
      coord[2] = (dim > 2) ? coords(local_index, 2) : 0.0;

      auto dr = std::visit(
        GetRank{classIds_h[local_index], classDims_h[local_index], coord},
        partition);
      reverse_partition[dr].indices.emplace_back(local_index);
      local_index += 1;

      const auto n = reverse_partition[dr].ent_offsets.size();
      for (int e = ent_dim + 1; e < n; ++e) {
        reverse_partition[dr].ent_offsets[e] += 1;
      }
    }
  }

  return reverse_partition;
}

} // namespace pcms
