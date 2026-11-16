#include "pcms/field/layout/omega_h_lagrange.h"
#include "pcms/utility/assert.h"
#include "pcms/utility/mesh_geometry.h"
#include "pcms/utility/omega_h_array_utils.h"
#include "pcms/utility/profile.h"
#include <stdexcept>

namespace pcms
{

namespace
{
int EntityDimForOrder(int order, int mesh_dim)
{
  switch (order) {
    case 0: return mesh_dim; // one DOF per element
    case 1: return 0;        // one DOF per vertex
    default:
      throw std::invalid_argument(
        "OmegaHLagrangeLayout: only order 0 and 1 are supported");
  }
}

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
  int n = mesh.nents(entity_dim);
  Kokkos::View<bool*, DeviceMemorySpace> owned("owned", n);
  auto src = Omega_h::Read<Omega_h::I8>(mesh.owned(entity_dim));
  Kokkos::parallel_for(
    n, OMEGA_H_LAMBDA(int i) { owned(i) = static_cast<bool>(src[i]); });
  return owned;
}

Kokkos::View<bool*, DeviceMemorySpace> BuildOwned(
  Omega_h::Mesh& mesh, int entity_dim, Omega_h::Read<Omega_h::I8> mask)
{
  auto owned = BuildOwned(mesh, entity_dim);
  auto mask_h = Omega_h::Read<Omega_h::I8>(mask);
  PCMS_ALWAYS_ASSERT(static_cast<int>(mask_h.size()) == mesh.nents(entity_dim));
  Kokkos::parallel_for(
    mesh.nents(entity_dim), OMEGA_H_LAMBDA(int i) {
      owned(i) = owned(i) && static_cast<bool>(mask_h[i]);
    });
  return owned;
}

} // namespace

OmegaHLagrangeLayout::OmegaHLagrangeLayout(
  Omega_h::Mesh& mesh, int order, int num_components,
  std::shared_ptr<const CoordinateSystem> coordinate_system,
  std::string global_id_name)
  : mesh_(mesh),
    order_(order),
    num_components_(num_components),
    global_id_name_(std::move(global_id_name))
{
  PCMS_FUNCTION_TIMER;
  SetCoordinateSystem(
    ResolveCoordinateSystem(std::move(coordinate_system), mesh_.dim()));
  int entity_dim = EntityDimForOrder(order_, mesh_.dim());

  EnsureClassification(mesh_);
  gids_ = BuildGids(mesh_, entity_dim, global_id_name_);
  coords_ = get_entity_centroids(mesh_, entity_dim);
  coords_2d_ = ConvertCoordsTo2D(coords_, mesh_.nents(entity_dim), mesh_.dim());
  owned_ = BuildOwned(mesh_, entity_dim);
  owned_host_ =
    Kokkos::View<bool*, HostMemorySpace>("owned_host", owned_.size());
  Kokkos::deep_copy(owned_host_, owned_);
  BuildOwnedViews();

  class_ids_ = Omega_h::Read<Omega_h::ClassId>(
    mesh_.get_array<Omega_h::ClassId>(entity_dim, "class_id"));
  class_dims_ = Omega_h::Read<Omega_h::I8>(
    mesh_.get_array<Omega_h::I8>(entity_dim, "class_dim"));

  discretization_ = std::make_shared<OmegaHDiscretization>(mesh_);
  BuildOwnedClassification();

  BuildGlobalToLocalPermutation(
    Kokkos::View<const GO*, DeviceMemorySpace>(gids_.data(), gids_.size()),
    EntOffsetsWithGhosts());
}

OmegaHLagrangeLayout::OmegaHLagrangeLayout(
  Omega_h::Mesh& mesh, int order, int num_components,
  std::shared_ptr<const CoordinateSystem> coordinate_system,
  Omega_h::Read<Omega_h::I8> owned_mask, std::string global_id_name)
  : mesh_(mesh),
    order_(order),
    num_components_(num_components),
    global_id_name_(std::move(global_id_name))
{
  PCMS_FUNCTION_TIMER;
  SetCoordinateSystem(
    ResolveCoordinateSystem(std::move(coordinate_system), mesh_.dim()));
  int entity_dim = EntityDimForOrder(order_, mesh_.dim());

  EnsureClassification(mesh_);
  gids_ = BuildGids(mesh_, entity_dim, global_id_name_);
  coords_ = get_entity_centroids(mesh_, entity_dim);
  coords_2d_ = ConvertCoordsTo2D(coords_, mesh_.nents(entity_dim), mesh_.dim());
  owned_ = BuildOwned(mesh_, entity_dim, owned_mask);
  owned_host_ =
    Kokkos::View<bool*, HostMemorySpace>("owned_host", owned_.size());
  Kokkos::deep_copy(owned_host_, owned_);
  BuildOwnedViews();

  class_ids_ = Omega_h::Read<Omega_h::ClassId>(
    mesh_.get_array<Omega_h::ClassId>(entity_dim, "class_id"));
  class_dims_ = Omega_h::Read<Omega_h::I8>(
    mesh_.get_array<Omega_h::I8>(entity_dim, "class_dim"));

  discretization_ = std::make_shared<OmegaHDiscretization>(mesh_);
  BuildOwnedClassification();

  BuildGlobalToLocalPermutation(
    Kokkos::View<const GO*, DeviceMemorySpace>(gids_.data(), gids_.size()),
    EntOffsetsWithGhosts());
}

std::shared_ptr<const Discretization> OmegaHLagrangeLayout::GetDiscretization()
  const noexcept
{
  return discretization_;
}

int OmegaHLagrangeLayout::GetNumComponents() const
{
  return num_components_;
}

LO OmegaHLagrangeLayout::GetNumDofHolderWithGhosts() const
{
  return mesh_.nents(EntityDimForOrder(order_, mesh_.dim()));
}

LO OmegaHLagrangeLayout::GetNumDofHolder() const
{
  return num_owned_;
}

bool OmegaHLagrangeLayout::HasGhosts() const
{
  return num_owned_ != GetNumDofHolderWithGhosts();
}

GO OmegaHLagrangeLayout::GetNumGlobalDofHolder() const
{
  return mesh_.nglobal_ents(EntityDimForOrder(order_, mesh_.dim()));
}

Rank1View<const bool, HostMemorySpace> OmegaHLagrangeLayout::GetOwnedHost()
  const
{
  return make_const_array_view(owned_host_);
}

GlobalIDView<HostMemorySpace> OmegaHLagrangeLayout::GetGidsHost() const
{
  return GlobalIDView<HostMemorySpace>(owned_gids_host_.data(),
                                       owned_gids_host_.size());
}

GlobalIDView<DeviceMemorySpace> OmegaHLagrangeLayout::GetGids() const
{
  return GlobalIDView<DeviceMemorySpace>(gids_.data(), gids_.size());
}

CoordinateView<DeviceMemorySpace>
OmegaHLagrangeLayout::GetDOFHolderCoordinates() const
{
  using LayoutPolicy =
    detail::default_layout_for_memory_space_t<DeviceMemorySpace>;
  Rank2View<const Real, DeviceMemorySpace, LayoutPolicy> coords_view(
    owned_coords_2d_.data(), num_owned_, mesh_.dim());
  return CoordinateView<DeviceMemorySpace, LayoutPolicy>{GetCoordinateSystem(),
                                                         coords_view};
}

Kokkos::View<const LO*, HostMemorySpace>
OmegaHLagrangeLayout::GetOwnedIndicesHost() const
{
  return owned_indices_host_;
}

Kokkos::View<const LO*, DeviceMemorySpace>
OmegaHLagrangeLayout::GetOwnedIndices() const
{
  return owned_indices_;
}

void OmegaHLagrangeLayout::BuildOwnedViews()
{
  auto owned = BuildOwnedLayoutData(
    owned_,
    Kokkos::View<const GO*, DeviceMemorySpace>(gids_.data(), gids_.size()),
    coords_2d_, mesh_.dim());
  num_owned_ = owned.num_owned;
  owned_indices_host_ = owned.owned_indices_host;
  owned_gids_host_ = owned.owned_gids_host;
  owned_coords_2d_ = owned.owned_coords_2d;
  owned_indices_ = owned.owned_indices;
}

bool OmegaHLagrangeLayout::IsDistributed() const
{
  return true;
}

EntOffsetsArray OmegaHLagrangeLayout::GetEntOffsets() const
{
  return !HasGhosts()
           ? EntOffsetsWithGhosts()
           : OwnedEntOffsets(EntOffsetsWithGhosts(), owned_indices_host_);
}

EntOffsetsArray OmegaHLagrangeLayout::EntOffsetsWithGhosts() const
{
  // Slot i holds the starting DOF index for entity dimension i.
  // Slot 4 holds the total DOF count.
  // For order-1 (entity_dim=0): offsets = {0, n, n, n, n}
  // For order-0 (entity_dim=mesh.dim()): offsets = {0,..,0, n, n}
  EntOffsetsArray offsets{};
  offsets.fill(0);
  int entity_dim = EntityDimForOrder(order_, mesh_.dim());
  LO n = mesh_.nents(entity_dim);
  for (int i = entity_dim + 1; i < ent_offsets_len; ++i)
    offsets[i] = n;
  return offsets;
}

int OmegaHLagrangeLayout::GetDimension() const
{
  return mesh_.dim();
}

int OmegaHLagrangeLayout::GetDOFHolderEntityDim() const
{
  return EntityDimForOrder(order_, mesh_.dim());
}

Rank1View<const LO, HostMemorySpace>
OmegaHLagrangeLayout::GetDOFHolderClassificationDimensionsHost() const
{
  return make_const_array_view(owned_classification_dims_host_);
}

Rank1View<const LO, HostMemorySpace>
OmegaHLagrangeLayout::GetDOFHolderClassificationIdsHost() const
{
  return make_const_array_view(owned_classification_ids_host_);
}

void OmegaHLagrangeLayout::BuildOwnedClassification()
{
  BuildOwnedClassificationHost(class_dims_, class_ids_, owned_indices_,
                               owned_classification_dims_host_,
                               owned_classification_ids_host_);
}

int OmegaHLagrangeLayout::GetOrder() const
{
  return order_;
}

Omega_h::Mesh& OmegaHLagrangeLayout::GetMesh() const
{
  return mesh_;
}

} // namespace pcms
