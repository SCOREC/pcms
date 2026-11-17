#include "pcms/field/layout/xgc.h"
#include "pcms/utility/assert.h"
#include "pcms/utility/entity_types.h"
#include "pcms/field/coordinate_systems/cylindrical.hpp"

namespace pcms
{

struct InitilizeXGCMembersFunctor
{
  Kokkos::View<bool*, DeviceMemorySpace> owned_;
  Kokkos::View<GO*, DeviceMemorySpace> gids_;
  Kokkos::View<LO*, DeviceMemorySpace> class_dims_;
  Kokkos::View<LO*, DeviceMemorySpace> class_ids_;
  Rank1View<const ClassificationDimension, DeviceMemorySpace> vert_class_dims_;
  Rank1View<const ClassificationId, DeviceMemorySpace> vert_class_ids_;

  InitilizeXGCMembersFunctor(
    Kokkos::View<bool*, DeviceMemorySpace> owned,
    Kokkos::View<GO*, DeviceMemorySpace> gids,
    Kokkos::View<LO*, DeviceMemorySpace> class_dims,
    Kokkos::View<LO*, DeviceMemorySpace> class_ids,
    Rank1View<const ClassificationDimension, DeviceMemorySpace> vert_class_dims,
    Rank1View<const ClassificationId, DeviceMemorySpace> vert_class_ids)
    : owned_(owned),
      gids_(gids),
      class_dims_(class_dims),
      class_ids_(class_ids),
      vert_class_dims_(vert_class_dims),
      vert_class_ids_(vert_class_ids)
  {
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(LO i) const
  {
    owned_(i) = true;
    gids_(i) = static_cast<GO>(i) + 1;
    // ClassificationDimension is int8_t, a small integer rather than a char
    // NOLINTNEXTLINE(bugprone-signed-char-misuse)
    class_dims_(i) = vert_class_dims_(i);
    class_ids_(i) = vert_class_ids_(i);
  }
};

XGCFieldLayout::XGCFieldLayout(
  const ReverseClassificationVertex& reverse_classification, LO num_plane_nodes)
  : owned_("xgc_owned", num_plane_nodes),
    gids_("xgc_gids", num_plane_nodes),
    class_dims_("xgc_class_dims", num_plane_nodes),
    class_ids_("xgc_class_ids", num_plane_nodes),
    owned_host_("xgc_owned_host", num_plane_nodes),
    gids_host_("xgc_gids_host", num_plane_nodes),
    classification_dims_host_("xgc_classification_dims_host", num_plane_nodes),
    classification_ids_host_("xgc_classification_ids_host", num_plane_nodes),
    coords_("xgc_coords", num_plane_nodes, 2),
    num_plane_nodes_(num_plane_nodes)
{
  // XGC poloidal planes are in the RZ plane of a cylindrical coordinate system
  SetCoordinateSystem(csys::CylindricalRZ::Create());
  auto discretization = std::make_shared<XGCDiscretization>(
    reverse_classification, num_plane_nodes_);
  Kokkos::parallel_for(
    "InitXGCMembers",
    Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0,
                                                            num_plane_nodes_),
    InitilizeXGCMembersFunctor(
      owned_, gids_, class_dims_, class_ids_,
      discretization->GetEntityClassificationDimensions(Vertex),
      discretization->GetEntityClassificationIds(Vertex)));

  Kokkos::deep_copy(owned_host_, owned_);
  Kokkos::deep_copy(gids_host_, gids_);
  Kokkos::deep_copy(classification_dims_host_, class_dims_);
  Kokkos::deep_copy(classification_ids_host_, class_ids_);
  // Copy coordinates to device
  Kokkos::deep_copy(coords_, 0.0);
  discretization_ = std::move(discretization);
}

std::shared_ptr<const Discretization> XGCFieldLayout::GetDiscretization()
  const noexcept
{
  return discretization_;
}

int XGCFieldLayout::GetNumComponents() const
{
  return 1;
}

LO XGCFieldLayout::GetNumOwnedDofHolder() const
{
  return num_plane_nodes_;
}

GO XGCFieldLayout::GetNumGlobalDofHolder() const
{
  return num_plane_nodes_;
}

Rank1View<const bool, HostMemorySpace> XGCFieldLayout::GetOwnedHost() const
{
  return make_const_array_view(owned_host_);
}

GlobalIDView<HostMemorySpace> XGCFieldLayout::GetGidsHost() const
{
  return make_const_array_view(gids_host_);
}

bool XGCFieldLayout::IsDistributed() const
{
  return false;
}

EntOffsetsArray XGCFieldLayout::GetEntOffsets() const
{
  return {0, static_cast<size_t>(num_plane_nodes_),
          static_cast<size_t>(num_plane_nodes_),
          static_cast<size_t>(num_plane_nodes_),
          static_cast<size_t>(num_plane_nodes_)};
}

CoordinateView<DeviceMemorySpace> XGCFieldLayout::GetDOFHolderCoordinates()
  const
{
  return CoordinateView<DeviceMemorySpace>{GetCoordinateSystem(),
                                           MakeConstRank2View(coords_)};
}

int XGCFieldLayout::GetDimension() const
{
  return 2;
}

Rank1View<const LO, HostMemorySpace>
XGCFieldLayout::GetDOFHolderClassificationDimensionsHost() const
{
  return make_const_array_view(classification_dims_host_);
}

Rank1View<const LO, HostMemorySpace>
XGCFieldLayout::GetDOFHolderClassificationIdsHost() const
{
  return make_const_array_view(classification_ids_host_);
}

LO XGCFieldLayout::GetFullDataSize() const noexcept
{
  return num_plane_nodes_;
}

} // namespace pcms
