#include "pcms/backends/xgc/layout.hpp"
#include "pcms/utility/assert.hpp"
#include "pcms/utility/entity_types.hpp"
#include "pcms/utility/coordinate_systems/cylindrical.hpp"

namespace pcms
{

struct InitilizeXGCMembersFunctor
{
  Kokkos::View<bool*, DeviceMemorySpace> owned_;
  Kokkos::View<GO*, DeviceMemorySpace> gids_;

  InitilizeXGCMembersFunctor(Kokkos::View<bool*, DeviceMemorySpace> owned,
                             Kokkos::View<GO*, DeviceMemorySpace> gids)
    : owned_(owned), gids_(gids)
  {
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(LO i) const
  {
    owned_(i) = true;
    gids_(i) = static_cast<GO>(i) + 1;
  }
};

// Sets the classification of the vertices of one geometric entity.
struct ClassifyVertsFunctor
{
  pcms::detail::DimID geom;
  Kokkos::View<LO*, DeviceMemorySpace> verts;
  Kokkos::View<LO*, DeviceMemorySpace> class_dims_;
  Kokkos::View<LO*, DeviceMemorySpace> class_ids_;

  ClassifyVertsFunctor(pcms::detail::DimID geom,
                       Kokkos::View<LO*, DeviceMemorySpace> verts,
                       Kokkos::View<LO*, DeviceMemorySpace> class_dims,
                       Kokkos::View<LO*, DeviceMemorySpace> class_ids)
    : geom(geom), verts(verts), class_dims_(class_dims), class_ids_(class_ids)
  {
  }

  KOKKOS_INLINE_FUNCTION void operator()(const int i) const
  {
    LO vert = verts(i);
    if (vert >= 0 && vert < class_dims_.extent(0)) {
      class_dims_(vert) = geom.dim;
      class_ids_(vert) = geom.id;
    }
  };
};

XGCFieldLayout::XGCFieldLayout(
  const ReverseClassificationVertex& reverse_classification, LO num_plane_nodes,
  std::string name)
  : FieldLayout(std::move(name)),
    owned_("xgc_owned", num_plane_nodes),
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
  Kokkos::parallel_for("InitXGCMembers",
                       Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(
                         0, num_plane_nodes_),
                       InitilizeXGCMembersFunctor(owned_, gids_));
  // vertices that no geometric entity of the reverse classification holds
  // stay unclassified (-1)
  Kokkos::deep_copy(class_dims_, LO{-1});
  Kokkos::deep_copy(class_ids_, LO{-1});
  for (const auto& [geom, verts] : reverse_classification) {
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
      ClassifyVertsFunctor(geom, verts_device, class_dims_, class_ids_));
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
}

int XGCFieldLayout::GetNumComponents() const
{
  return 1;
}

LO XGCFieldLayout::GetNumDofHolder() const
{
  return num_plane_nodes_;
}

GO XGCFieldLayout::GetNumGlobalDofHolder() const
{
  return num_plane_nodes_;
}

Rank1View<const bool, HostMemorySpace> XGCFieldLayout::GetOwnedHost() const
{
  return pcms::detail::make_const_array_view(owned_host_);
}

GlobalIDView<HostMemorySpace> XGCFieldLayout::GetGidsHost() const
{
  return pcms::detail::make_const_array_view(gids_host_);
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
  return pcms::detail::make_const_array_view(classification_dims_host_);
}

Rank1View<const LO, HostMemorySpace>
XGCFieldLayout::GetDOFHolderClassificationIdsHost() const
{
  return pcms::detail::make_const_array_view(classification_ids_host_);
}

LO XGCFieldLayout::GetFullDataSize() const noexcept
{
  return num_plane_nodes_;
}

} // namespace pcms
