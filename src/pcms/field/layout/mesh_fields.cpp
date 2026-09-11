#include "pcms/field/layout/mesh_fields.h"
#include "pcms/utility/assert.h"
#include "pcms/utility/inclusive_scan.h"
#include "pcms/utility/profile.h"
#include <Omega_h_for.hpp>
#include <memory>

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

// this is a workaround to specify the parametric coordinates for MeshFields to
// be replaced when https://github.com/SCOREC/meshFields/issues/70 is resolved
struct ComputeVertexCoordsFunctor
{
  Kokkos::View<Real**> dof_holder_coords_;
  Omega_h::Reals coords_;
  size_t offset_;

  ComputeVertexCoordsFunctor(Kokkos::View<Real**> dof_holder_coords,
                             Omega_h::Reals coords, size_t offset)
    : dof_holder_coords_(dof_holder_coords), coords_(coords), offset_(offset)
  {
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(LO i) const
  {
    dof_holder_coords_(offset_ + i, 0) = coords_[2 * i + 0];
    dof_holder_coords_(offset_ + i, 1) = coords_[2 * i + 1];
  }
};

// this is a workaround to specify the parametric coordinates for MeshFields to
// be replaced when https://github.com/SCOREC/meshFields/issues/70 is resolved
struct ComputeEdgeCoordsFunctor
{
  Kokkos::View<Real**> dof_holder_coords_;
  Omega_h::Reals coords_;
  Omega_h::LOs edge_verts_;
  size_t offset_;

  ComputeEdgeCoordsFunctor(Kokkos::View<Real**> dof_holder_coords,
                           Omega_h::Reals coords, Omega_h::LOs edge_verts,
                           size_t offset)
    : dof_holder_coords_(dof_holder_coords),
      coords_(coords),
      edge_verts_(edge_verts),
      offset_(offset)
  {
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(LO i) const
  {
    auto verts = Omega_h::gather_verts<2>(edge_verts_, i);
    Real x0 = coords_[2 * verts[0] + 0];
    Real y0 = coords_[2 * verts[0] + 1];
    Real x1 = coords_[2 * verts[1] + 0];
    Real y1 = coords_[2 * verts[1] + 1];
    dof_holder_coords_(offset_ + i, 0) = (x0 + x1) / 2;
    dof_holder_coords_(offset_ + i, 1) = (y0 + y1) / 2;
  }
};

struct CopyClassInfoFunctor
{
  Omega_h::Write<Omega_h::ClassId> class_ids_;
  Omega_h::Write<Omega_h::I8> class_dims_;
  Kokkos::View<bool*> owned_;
  Omega_h::Read<Omega_h::ClassId> ids_;
  Omega_h::Read<Omega_h::I8> dims_;
  Omega_h::Read<Omega_h::I8> owned_data_;
  size_t offset_;

  CopyClassInfoFunctor(Omega_h::Write<Omega_h::ClassId> class_ids,
                       Omega_h::Write<Omega_h::I8> class_dims,
                       Kokkos::View<bool*> owned,
                       Omega_h::Read<Omega_h::ClassId> ids,
                       Omega_h::Read<Omega_h::I8> dims,
                       Omega_h::Read<Omega_h::I8> owned_data, size_t offset)
    : class_ids_(class_ids),
      class_dims_(class_dims),
      owned_(owned),
      ids_(ids),
      dims_(dims),
      owned_data_(owned_data),
      offset_(offset)
  {
  }

  OMEGA_H_DEVICE
  void operator()(LO i) const
  {
    class_ids_[offset_ + i] = ids_[i];
    class_dims_[offset_ + i] = dims_[i];
    owned_[offset_ + i] = owned_data_[i];
  }
};

MeshFieldsAdapterLayout::MeshFieldsAdapterLayout(
  Omega_h::Mesh& mesh, std::array<int, 4> nodes_per_dim, int num_components,
  CoordinateSystem coordinate_system, std::string global_id_name)
  : mesh_(mesh),
    global_id_name_(global_id_name),
    num_components_(num_components),
    coordinate_system_(coordinate_system),
    nodes_per_dim_(nodes_per_dim),
    dof_holder_coords_("", GetNumLocalDofHolder(), mesh_.dim()),
    class_ids_(GetNumEnts()),
    class_dims_(class_ids_.size()),
    owned_("", class_dims_.size()),
    owned_host_("", class_dims_.size())
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
        ComputeVertexCoordsFunctor functor(dof_holder_coords_, coords, offset);
        Kokkos::parallel_for(mesh_.nents(0), functor);
      } else if (i == 1) {
        auto edge_verts = mesh_.ask_verts_of(1);
        ComputeEdgeCoordsFunctor functor(dof_holder_coords_, coords, edge_verts,
                                         offset);
        Kokkos::parallel_for(mesh_.nents(1), functor);
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

      CopyClassInfoFunctor functor(class_ids_, class_dims_, owned_, ids, dims,
                                   owned, offset);
      Omega_h::parallel_for(mesh_.nents(i), functor);
      offset += mesh.nents(i);
    }
  }
  gids_host_ = Omega_h::HostWrite<Omega_h::GO>(gids_);
  BuildOwnedViews();

  int n = class_ids_.size();
  classification_dims_host_ =
    Kokkos::View<LO*, HostMemorySpace>("classification_dims", n);
  classification_ids_host_ =
    Kokkos::View<LO*, HostMemorySpace>("classification_ids", n);
  auto class_dims_h = Omega_h::HostRead<Omega_h::I8>(class_dims_);
  auto class_ids_h = Omega_h::HostRead<Omega_h::ClassId>(class_ids_);
  for (int i = 0; i < n; ++i) {
    classification_dims_host_(i) =
      static_cast<LO>(static_cast<unsigned char>(class_dims_h[i]));
    classification_ids_host_(i) = static_cast<LO>(class_ids_h[i]);
  }
  discretization_ = std::make_shared<OmegaHDiscretization>(mesh_);
}

std::shared_ptr<const Discretization>
MeshFieldsAdapterLayout::GetDiscretization() const noexcept
{
  return discretization_;
}

int MeshFieldsAdapterLayout::GetNumComponents() const
{
  return num_components_;
}

LO MeshFieldsAdapterLayout::GetNumLocalDofHolder() const
{
  LO count = 0;
  for (int i = 0; i <= mesh_.dim(); ++i) {
    count += mesh_.nents(i) * nodes_per_dim_[i];
  }
  return count;
}

LO MeshFieldsAdapterLayout::GetNumOwnedDofHolder() const
{
  return num_owned_;
}

GO MeshFieldsAdapterLayout::GetNumGlobalDofHolder() const
{
  LO count = 0;
  for (int i = 0; i <= mesh_.dim(); ++i) {
    count += mesh_.nglobal_ents(i) * nodes_per_dim_[i];
  }
  return count;
}

std::array<int, 4> MeshFieldsAdapterLayout::GetNodesPerDim() const
{
  return nodes_per_dim_;
}

Rank1View<const bool, HostMemorySpace> MeshFieldsAdapterLayout::GetOwnedHost()
  const
{
  Kokkos::deep_copy(owned_host_, owned_);
  return make_const_array_view(owned_host_);
}

GlobalIDView<HostMemorySpace> MeshFieldsAdapterLayout::GetGidsHost() const
{
  return GlobalIDView<HostMemorySpace>(gids_host_.data(), gids_host_.size());
}

GlobalIDView<DeviceMemorySpace> MeshFieldsAdapterLayout::GetGids() const
{
  return GlobalIDView<DeviceMemorySpace>(gids_.data(), gids_.size());
}

CoordinateView<DeviceMemorySpace>
MeshFieldsAdapterLayout::GetDOFHolderCoordinates() const
{
  auto coords_view = MakeConstRank2View(dof_holder_coords_);
  return CoordinateView<DeviceMemorySpace>{coordinate_system_, coords_view};
}

GlobalIDView<HostMemorySpace> MeshFieldsAdapterLayout::GetOwnedGidsHost() const
{
  return GlobalIDView<HostMemorySpace>(owned_gids_host_.data(),
                                       owned_gids_host_.size());
}

GlobalIDView<DeviceMemorySpace> MeshFieldsAdapterLayout::GetOwnedGids() const
{
  return GlobalIDView<DeviceMemorySpace>(owned_gids_.data(),
                                         owned_gids_.size());
}

CoordinateView<DeviceMemorySpace>
MeshFieldsAdapterLayout::GetOwnedDOFHolderCoordinates() const
{
  auto coords_view = MakeConstRank2View(owned_coords_2d_);
  return CoordinateView<DeviceMemorySpace>{coordinate_system_, coords_view};
}

Kokkos::View<const LO*, HostMemorySpace>
MeshFieldsAdapterLayout::GetOwnedToLocalHost() const
{
  return owned_to_local_host_;
}

Kokkos::View<const LO*, DeviceMemorySpace>
MeshFieldsAdapterLayout::GetOwnedToLocal() const
{
  return owned_to_local_;
}

void MeshFieldsAdapterLayout::BuildOwnedViews()
{
  Kokkos::deep_copy(owned_host_, owned_);

  auto owned = BuildOwnedLayoutData(
    owned_host_,
    GlobalIDView<HostMemorySpace>(gids_host_.data(), gids_host_.size()),
    dof_holder_coords_, mesh_.dim());
  num_owned_ = owned.num_owned;
  owned_to_local_host_ = owned.owned_to_local_host;
  owned_gids_host_ = owned.owned_gids_host;
  owned_coords_2d_ = owned.owned_coords_2d;
  owned_to_local_ = owned.owned_to_local;
  owned_gids_ = owned.owned_gids;
}

bool MeshFieldsAdapterLayout::IsDistributed() const
{
  return true;
}

Omega_h::Read<Omega_h::ClassId> MeshFieldsAdapterLayout::GetClassIDs() const
{
  PCMS_FUNCTION_TIMER;
  return Omega_h::Read(class_ids_);
}

Omega_h::Read<Omega_h::I8> MeshFieldsAdapterLayout::GetClassDims() const
{
  PCMS_FUNCTION_TIMER;
  return Omega_h::Read(class_dims_);
}

size_t MeshFieldsAdapterLayout::GetNumEnts() const
{
  size_t n = 0;
  for (int i = 0; i <= mesh_.dim(); ++i) {
    if (nodes_per_dim_[i])
      n += mesh_.nents(i);
  }
  return n;
}

Omega_h::Mesh& MeshFieldsAdapterLayout::GetMesh() const
{
  return mesh_;
}

EntOffsetsArray MeshFieldsAdapterLayout::GetEntOffsets() const
{
  EntOffsetsArray offsets{};
  size_t offset = 0;
  for (size_t i = 0; i < offsets.size(); ++i) {
    offsets[i] = offset;
    if (i <= static_cast<size_t>(mesh_.dim()) && nodes_per_dim_[i])
      offset += mesh_.nents(i);
  }
  offsets[offsets.size() - 1] = offset;
  return offsets;
}

int MeshFieldsAdapterLayout::GetDimension() const
{
  return mesh_.dim();
}

Rank1View<const LO, HostMemorySpace>
MeshFieldsAdapterLayout::GetDOFHolderClassificationDimensionsHost() const
{
  return make_const_array_view(classification_dims_host_);
}

Rank1View<const LO, HostMemorySpace>
MeshFieldsAdapterLayout::GetDOFHolderClassificationIdsHost() const
{
  return make_const_array_view(classification_ids_host_);
}

} // namespace pcms
