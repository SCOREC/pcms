#include "pcms/field/layout/xgc.h"
#include "pcms/utility/assert.h"
#include "pcms/field/coordinate_systems/cylindrical.hpp"

namespace pcms
{

struct InitilizeXGCMembersFunctor
{
  Kokkos::View<bool*, DeviceMemorySpace> owned_;
  Kokkos::View<GO*, DeviceMemorySpace> gids_;
  Kokkos::View<LO*, DeviceMemorySpace> class_dims_;
  Kokkos::View<LO*, DeviceMemorySpace> class_ids_;

  InitilizeXGCMembersFunctor(Kokkos::View<bool*, DeviceMemorySpace> owned,
                             Kokkos::View<GO*, DeviceMemorySpace> gids,
                             Kokkos::View<LO*, DeviceMemorySpace> class_dims,
                             Kokkos::View<LO*, DeviceMemorySpace> class_ids)
    : owned_(owned), gids_(gids), class_dims_(class_dims), class_ids_(class_ids)
  {
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(LO i) const
  {
    owned_(i) = false;
    gids_(i) = static_cast<GO>(i) + 1;
    class_dims_(i) = -1;
    class_ids_(i) = -1;
  }
};

struct ClassifyVertsAndOwnedFunctor
{
  pcms::DimID geom;
  Kokkos::View<LO*, DeviceMemorySpace> verts;
  Kokkos::View<bool*, DeviceMemorySpace> owned_;
  Kokkos::View<LO*, DeviceMemorySpace> class_dims_;
  Kokkos::View<LO*, DeviceMemorySpace> class_ids_;

  ClassifyVertsAndOwnedFunctor(pcms::DimID geom,
                               Kokkos::View<LO*, DeviceMemorySpace> verts,
                               Kokkos::View<bool*, DeviceMemorySpace> owned,
                               Kokkos::View<LO*, DeviceMemorySpace> class_dims,
                               Kokkos::View<LO*, DeviceMemorySpace> class_ids)
    : geom(geom),
      verts(verts),
      owned_(owned),
      class_dims_(class_dims),
      class_ids_(class_ids)
  {
  }

  KOKKOS_INLINE_FUNCTION void operator()(const int i) const
  {
    LO vert = verts(i);
    if (vert >= 0 && vert < owned_.extent(0)) {
      owned_(vert) = true;
      class_dims_(vert) = geom.dim;
      class_ids_(vert) = geom.id;
    }
  };
};

XGCFieldLayout::XGCFieldLayout(
  const ReverseClassificationVertex& reverse_classification,
  std::function<int8_t(int, int)> in_overlap, LO num_plane_nodes)
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
  PCMS_ALWAYS_ASSERT(static_cast<bool>(in_overlap));
  Kokkos::parallel_for(
    "InitXGCMembers",
    Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0,
                                                            num_plane_nodes_),
    InitilizeXGCMembersFunctor(owned_, gids_, class_dims_, class_ids_));

  for (const auto& [geom, verts] : reverse_classification) {
    if (!in_overlap(geom.dim, geom.id))
      continue;
    auto verts_host =
      Kokkos::View<LO*, HostMemorySpace>("verts_host", verts.size());
    int idx = 0;
    for (LO vert : verts)
      verts_host(idx++) = vert;
    auto verts_device =
      Kokkos::View<LO*, DeviceMemorySpace>("verts_device", verts.size());
    Kokkos::deep_copy(verts_device, verts_host);
    Kokkos::parallel_for(
      "ClassifyVerts", Kokkos::RangePolicy<>(0, verts.size()),
      ClassifyVertsAndOwnedFunctor(geom, verts_device, owned_, class_dims_,
                                   class_ids_));
    Kokkos::fence(); // Wait for kernel to complete before verts_device is
                     // destroyed, better would be to optimize the
                     // reverse_classification data structure to avoid this copy
                     // and synchronization, but this is simpler for now
  }

  Kokkos::deep_copy(owned_host_, owned_);
  Kokkos::deep_copy(gids_host_, gids_);
  Kokkos::deep_copy(classification_dims_host_, class_dims_);
  Kokkos::deep_copy(classification_ids_host_, class_ids_);
  // Copy coordinates to device
  Kokkos::deep_copy(coords_, 0.0);
  discretization_ = std::make_shared<XGCDiscretization>(reverse_classification,
                                                        num_plane_nodes_);
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
