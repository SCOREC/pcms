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

OmegaHEntityLayout::OmegaHEntityLayout(Omega_h::Mesh& mesh, int entity_dim,
                                       int num_components,
                                       CoordinateSystem coordinate_system,
                                       std::string global_id_name)
  : dimension_(mesh.dim()),
    entity_dim_(entity_dim),
    num_components_(num_components),
    num_global_dof_holder_(mesh.nglobal_ents(entity_dim)),
    coordinate_system_(coordinate_system),
    gids_(BuildGids(mesh, entity_dim, global_id_name)),
    coords_(get_entity_centroids(mesh, entity_dim)),
    coords_2d_(ConvertCoordsTo2D(coords_, mesh.nents(entity_dim), mesh.dim())),
    class_ids_(mesh.get_array<Omega_h::ClassId>(entity_dim, "class_id")),
    class_dims_(mesh.get_array<Omega_h::I8>(entity_dim, "class_dim")),
    owned_(BuildOwned(mesh, entity_dim)),
    owned_host_("owned_host", owned_.size()),
    classification_dims_host_("classification_dims", mesh.nents(entity_dim)),
    classification_ids_host_("classification_ids", mesh.nents(entity_dim)),
    discretization_(std::make_shared<OmegaHDiscretization>(mesh))
{
  PCMS_ALWAYS_ASSERT(entity_dim_ >= 0 && entity_dim_ <= dimension_);

  gids_host_ = Omega_h::HostWrite<Omega_h::GO>(gids_);
  Kokkos::deep_copy(owned_host_, owned_);
  BuildOwnedViews();

  class_dims_ = Omega_h::Read<Omega_h::I8>(class_dims_);
  class_ids_ = Omega_h::Read<Omega_h::ClassId>(class_ids_);
  auto class_dims_host = Omega_h::HostRead<Omega_h::I8>(class_dims_);
  auto class_ids_host = Omega_h::HostRead<Omega_h::ClassId>(class_ids_);
  for (int i = 0; i < mesh.nents(entity_dim_); ++i) {
    classification_dims_host_(i) =
      static_cast<LO>(static_cast<unsigned char>(class_dims_host[i]));
    classification_ids_host_(i) = static_cast<LO>(class_ids_host[i]);
  }
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

LO OmegaHEntityLayout::GetNumLocalDofHolder() const
{
  return static_cast<LO>(coords_.size() / dimension_);
}

LO OmegaHEntityLayout::GetNumOwnedDofHolder() const
{
  return num_owned_;
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
  return GlobalIDView<HostMemorySpace>(gids_host_.data(), gids_host_.size());
}

GlobalIDView<DeviceMemorySpace> OmegaHEntityLayout::GetGids() const
{
  return GlobalIDView<DeviceMemorySpace>(gids_.data(), gids_.size());
}

CoordinateView<DeviceMemorySpace> OmegaHEntityLayout::GetDOFHolderCoordinates()
  const
{
  using LayoutPolicy =
    detail::default_layout_for_memory_space_t<DeviceMemorySpace>;
  Rank2View<const Real, DeviceMemorySpace, LayoutPolicy> coords_view(
    coords_2d_.data(), GetNumLocalDofHolder(), dimension_);
  return CoordinateView<DeviceMemorySpace, LayoutPolicy>{coordinate_system_,
                                                         coords_view};
}

GlobalIDView<HostMemorySpace> OmegaHEntityLayout::GetOwnedGidsHost() const
{
  return GlobalIDView<HostMemorySpace>(owned_gids_host_.data(),
                                       owned_gids_host_.size());
}

GlobalIDView<DeviceMemorySpace> OmegaHEntityLayout::GetOwnedGids() const
{
  return GlobalIDView<DeviceMemorySpace>(owned_gids_.data(),
                                         owned_gids_.size());
}

CoordinateView<DeviceMemorySpace>
OmegaHEntityLayout::GetOwnedDOFHolderCoordinates() const
{
  using LayoutPolicy =
    detail::default_layout_for_memory_space_t<DeviceMemorySpace>;
  Rank2View<const Real, DeviceMemorySpace, LayoutPolicy> coords_view(
    owned_coords_2d_.data(), num_owned_, dimension_);
  return CoordinateView<DeviceMemorySpace, LayoutPolicy>{coordinate_system_,
                                                         coords_view};
}

Kokkos::View<const LO*, HostMemorySpace>
OmegaHEntityLayout::GetOwnedToLocalHost() const
{
  return owned_to_local_host_;
}

Kokkos::View<const LO*, DeviceMemorySpace> OmegaHEntityLayout::GetOwnedToLocal()
  const
{
  return owned_to_local_;
}

void OmegaHEntityLayout::BuildOwnedViews()
{
  auto owned = BuildOwnedLayoutData(
    owned_host_,
    GlobalIDView<HostMemorySpace>(gids_host_.data(), gids_host_.size()),
    coords_2d_, dimension_);
  num_owned_ = owned.num_owned;
  owned_to_local_host_ = owned.owned_to_local_host;
  owned_gids_host_ = owned.owned_gids_host;
  owned_coords_2d_ = owned.owned_coords_2d;
  owned_to_local_ = owned.owned_to_local;
  owned_gids_ = owned.owned_gids;
}

bool OmegaHEntityLayout::IsDistributed() const
{
  return true;
}

EntOffsetsArray OmegaHEntityLayout::GetEntOffsets() const
{
  EntOffsetsArray offsets{};
  offsets.fill(0);
  const auto n = static_cast<size_t>(GetNumLocalDofHolder());
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
  return make_const_array_view(classification_dims_host_);
}

Rank1View<const LO, HostMemorySpace>
OmegaHEntityLayout::GetDOFHolderClassificationIdsHost() const
{
  return make_const_array_view(classification_ids_host_);
}

} // namespace pcms
