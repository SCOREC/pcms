#include "point_cloud_layout.h"

namespace pcms
{
PointCloudLayout::PointCloudLayout(int dim, Kokkos::View<Real**> coords, CoordinateSystem coordinate_system)
  : dim_(dim),
    coordinate_system_(coordinate_system),
    coords_(coords),
    owned_("", coords.extent(0)),
    gids_("", coords.extent(0))
{
  compoents_ = 1;

  Kokkos::parallel_for(
    owned_.size(), KOKKOS_LAMBDA(int i) {
      owned_[i] = 0;
      gids_[i] = 0;
    });
}

int PointCloudLayout::GetNumComponents() const
{
    return compoents_;
}

LO PointCloudLayout::GetNumOwnedDofHolder() const
{
    return coords_.extent(0);
}

GO PointCloudLayout::GetNumGlobalDofHolder() const
{
    return coords_.extent(0);
}

Rank1View<const bool, HostMemorySpace> PointCloudLayout::GetOwned() const
{
    return Rank1View<const bool, HostMemorySpace>{std::data(owned_),
                                                  std::size(owned_)};
}

GlobalIDView<HostMemorySpace> PointCloudLayout::GetGids() const
{
  static_assert(std::is_same_v<HostMemorySpace, DefaultExecutionSpace::memory_space>, "types must match");
  return GlobalIDView<HostMemorySpace>(gids_.data(), gids_.size());
}

CoordinateView<HostMemorySpace>
PointCloudLayout::GetDOFHolderCoordinates() const
{
  Rank2View<const Real, HostMemorySpace> coords_view(coords_.data(),
                                                coords_.extent(0), 2);
  return CoordinateView<HostMemorySpace>{coordinate_system_, coords_view};
}

bool PointCloudLayout::IsDistributed()
{
  return false;
}

size_t PointCloudLayout::GetNumEnts() const
{
  return coords_.extent(0);
}

std::array<size_t, 5> PointCloudLayout::GetEntOffsets() const
{
  std::array<size_t, 5> offsets{};
  for (int i = 0; i < offsets.size(); ++i)
    offsets[i] = coords_.extent(0);
  offsets[0] = 0;
  return offsets;
}

std::array<int, 4> PointCloudLayout::GetNodesPerDim() const
{
  std::array<int, 4> nodes{};
  for (int i = 0; i < nodes.size(); ++i)
    nodes[i] = 0;
  nodes[0] = 1;
  return nodes;
}

ReversePartitionMap2 PointCloudLayout::GetReversePartitionMap(
  const redev::Partition& partition) const
{
  throw std::runtime_error("Unimplemented");
}
} // namespace pcms
