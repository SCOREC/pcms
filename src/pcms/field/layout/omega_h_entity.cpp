#include "pcms/field/layout/omega_h_entity.h"

#include "pcms/utility/assert.h"
#include "pcms/utility/mesh_geometry.h"
#include "pcms/utility/omega_h_array_utils.h"

namespace pcms
{

namespace
{

template <typename T>
Omega_h::Write<Omega_h::GO> GetGidsHelper(Omega_h::Mesh& mesh, int entity_dim,
                                          const std::string& global_id_name)
{
  auto dim_gids = mesh.get_array<T>(entity_dim, global_id_name);
  Omega_h::Write<Omega_h::GO> gids(dim_gids.size());
  Omega_h::parallel_for(
    dim_gids.size(), OMEGA_H_LAMBDA(int i) { gids[i] = dim_gids[i]; });
  return gids;
}

Omega_h::Write<Omega_h::GO> BuildGids(Omega_h::Mesh& mesh, int entity_dim,
                                      const std::string& global_id_name)
{
  auto tag = mesh.get_tagbase(entity_dim, global_id_name);
  Omega_h::Write<Omega_h::GO> gids;
  if (Omega_h::is<Omega_h::GO>(tag)) {
    gids = GetGidsHelper<Omega_h::GO>(mesh, entity_dim, global_id_name);
  } else if (Omega_h::is<Omega_h::LO>(tag)) {
    gids = GetGidsHelper<Omega_h::LO>(mesh, entity_dim, global_id_name);
  } else {
    std::cerr << "Weird tag type for global arrays.\n";
    std::abort();
  }
  return gids;
}

Kokkos::View<bool*, DeviceMemorySpace> BuildOwned(Omega_h::Mesh& mesh,
                                                  int entity_dim)
{
  Kokkos::View<bool*, DeviceMemorySpace> owned("owned", mesh.nents(entity_dim));
  auto owned_h = Omega_h::Read<Omega_h::I8>(mesh.owned(entity_dim));
  Kokkos::parallel_for(
    mesh.nents(entity_dim),
    OMEGA_H_LAMBDA(int i) { owned(i) = static_cast<bool>(owned_h[i]); });
  return owned;
}

} // namespace

OmegaHEntityLayout::OmegaHEntityLayout(
  Omega_h::Mesh& mesh, int entity_dim, int num_components,
  std::shared_ptr<const CoordinateSystem> coordinate_system,
  std::string global_id_name)
  // dimension_ is the first member, so this runs before the class_ids_/
  // class_dims_ initializers below read the mesh's classification tags.
  : dimension_((EnsureClassification(mesh), mesh.dim())),
    entity_dim_(entity_dim),
    num_components_(num_components),
    num_global_dof_holder_(mesh.nglobal_ents(entity_dim)),
    gids_(BuildGids(mesh, entity_dim, global_id_name)),
    coords_(get_entity_centroids(mesh, entity_dim)),
    coords_2d_(ConvertCoordsTo2D(coords_, mesh.nents(entity_dim), mesh.dim())),
    class_ids_(mesh.get_array<Omega_h::ClassId>(entity_dim, "class_id")),
    class_dims_(mesh.get_array<Omega_h::I8>(entity_dim, "class_dim")),
    owned_(BuildOwned(mesh, entity_dim)),
    owned_host_("owned_host", owned_.size()),
    discretization_(std::make_shared<OmegaHDiscretization>(mesh))
{
  SetCoordinateSystem(
    ResolveCoordinateSystem(std::move(coordinate_system), dimension_));
  PCMS_ALWAYS_ASSERT(entity_dim_ >= 0 && entity_dim_ <= dimension_);

  Kokkos::deep_copy(owned_host_, owned_);
  BuildOwnedViews();

  class_dims_ = Omega_h::Read<Omega_h::I8>(class_dims_);
  class_ids_ = Omega_h::Read<Omega_h::ClassId>(class_ids_);
  BuildOwnedClassification();
}

std::shared_ptr<const Discretization> OmegaHEntityLayout::GetDiscretization()
  const noexcept
{
  return discretization_;
}

int OmegaHEntityLayout::GetNumComponents() const
{
  return num_components_;
}

LO OmegaHEntityLayout::GetNumDofHolderWithGhosts() const
{
  return static_cast<LO>(coords_.size() / dimension_);
}

LO OmegaHEntityLayout::GetNumDofHolder() const
{
  return num_owned_;
}

bool OmegaHEntityLayout::HasGhosts() const
{
  return num_owned_ != GetNumDofHolderWithGhosts();
}

GO OmegaHEntityLayout::GetNumGlobalDofHolder() const
{
  return num_global_dof_holder_;
}

Rank1View<const bool, HostMemorySpace> OmegaHEntityLayout::GetOwnedHost() const
{
  return make_const_array_view(owned_host_);
}

GlobalIDView<HostMemorySpace> OmegaHEntityLayout::GetGidsHost() const
{
  return GlobalIDView<HostMemorySpace>(owned_gids_host_.data(),
                                       owned_gids_host_.size());
}

CoordinateView<DeviceMemorySpace> OmegaHEntityLayout::GetDOFHolderCoordinates()
  const
{
  using LayoutPolicy =
    detail::default_layout_for_memory_space_t<DeviceMemorySpace>;
  Rank2View<const Real, DeviceMemorySpace, LayoutPolicy> coords_view(
    owned_coords_2d_.data(), num_owned_, dimension_);
  return CoordinateView<DeviceMemorySpace, LayoutPolicy>{GetCoordinateSystem(),
                                                         coords_view};
}

Kokkos::View<const LO*, HostMemorySpace>
OmegaHEntityLayout::GetOwnedIndicesHost() const
{
  return owned_indices_host_;
}

Kokkos::View<const LO*, DeviceMemorySpace> OmegaHEntityLayout::GetOwnedIndices()
  const
{
  return owned_indices_;
}

void OmegaHEntityLayout::BuildOwnedViews()
{
  auto owned = BuildOwnedLayoutData(
    owned_,
    Kokkos::View<const GO*, DeviceMemorySpace>(gids_.data(), gids_.size()),
    coords_2d_, dimension_);
  num_owned_ = owned.num_owned;
  owned_indices_host_ = owned.owned_indices_host;
  owned_gids_host_ = owned.owned_gids_host;
  owned_coords_2d_ = owned.owned_coords_2d;
  owned_indices_ = owned.owned_indices;
}

bool OmegaHEntityLayout::IsDistributed() const
{
  return true;
}

EntOffsetsArray OmegaHEntityLayout::GetEntOffsets() const
{
  return !HasGhosts()
           ? EntOffsetsWithGhosts()
           : OwnedEntOffsets(EntOffsetsWithGhosts(), owned_indices_host_);
}

EntOffsetsArray OmegaHEntityLayout::EntOffsetsWithGhosts() const
{
  EntOffsetsArray offsets{};
  offsets.fill(0);
  const auto n = static_cast<size_t>(GetNumDofHolderWithGhosts());
  for (int i = entity_dim_ + 1; i < ent_offsets_len; ++i)
    offsets[i] = n;
  return offsets;
}

int OmegaHEntityLayout::GetDimension() const
{
  return dimension_;
}

int OmegaHEntityLayout::GetDOFHolderEntityDim() const
{
  return entity_dim_;
}

Rank1View<const LO, HostMemorySpace>
OmegaHEntityLayout::GetDOFHolderClassificationDimensionsHost() const
{
  return make_const_array_view(owned_classification_dims_host_);
}

Rank1View<const LO, HostMemorySpace>
OmegaHEntityLayout::GetDOFHolderClassificationIdsHost() const
{
  return make_const_array_view(owned_classification_ids_host_);
}

void OmegaHEntityLayout::BuildOwnedClassification()
{
  BuildOwnedClassificationHost(class_dims_, class_ids_, owned_indices_,
                               owned_classification_dims_host_,
                               owned_classification_ids_host_);
}

} // namespace pcms
