#include "pcms/field/distributed_partition.h"
#include "pcms/field/field_layout.h"
#include "pcms/utility/assert.h"
#include "pcms/utility/memory_spaces.h"
#include <Kokkos_Core.hpp>
#include <algorithm>
#include <limits>

namespace pcms
{

namespace
{
// Copies a (possibly device-resident) coordinate view into a flat, row-major
// host buffer [n][dim], regardless of the view's native Kokkos layout.
std::vector<Real> CopyOwnedCoordsToHost(
  Rank2View<const Real, DeviceMemorySpace> coords, LO n, int dim)
{
  using ViewLayout = std::conditional_t<
    std::is_same_v<typename decltype(coords)::layout_type, Kokkos::layout_left>,
    Kokkos::LayoutLeft, Kokkos::LayoutRight>;
  Kokkos::View<const Real**, ViewLayout, DeviceMemorySpace,
               Kokkos::MemoryUnmanaged>
    device_view(coords.data_handle(), n, dim);
  auto host_view = Kokkos::create_mirror_view_and_copy(HostMemorySpace(),
                                                       device_view);
  std::vector<Real> flat(static_cast<size_t>(n) * static_cast<size_t>(dim));
  for (LO i = 0; i < n; ++i) {
    for (int d = 0; d < dim; ++d) {
      flat[static_cast<size_t>(i) * static_cast<size_t>(dim) +
           static_cast<size_t>(d)] = host_view(i, d);
    }
  }
  return flat;
}
} // namespace

DistributedBoundingBoxPartition::DistributedBoundingBoxPartition(
  int dim, Real padding, std::vector<Real> mins, std::vector<Real> maxs,
  std::vector<char> owns_nothing)
  : num_ranks_(static_cast<int>(owns_nothing.size())),
    dim_(dim),
    padding_(padding),
    mins_(std::move(mins)),
    maxs_(std::move(maxs)),
    owns_nothing_(std::move(owns_nothing))
{
}

DistributedBoundingBoxPartition DistributedBoundingBoxPartition::
  FromGatheredBounds(int dim, Real padding, std::vector<Real> mins,
                     std::vector<Real> maxs, std::vector<char> owns_nothing)
{
  return DistributedBoundingBoxPartition(dim, padding, std::move(mins),
                                         std::move(maxs),
                                         std::move(owns_nothing));
}

void DistributedBoundingBoxPartition::LocalOwnedBounds(
  const FieldLayout& layout, int dim, std::array<Real, 3>& min_out,
  std::array<Real, 3>& max_out, bool& owns_nothing_out)
{
  PCMS_ALWAYS_ASSERT(dim > 0 && dim <= 3);
  min_out.fill(0.0);
  max_out.fill(0.0);

  const LO num_owned = layout.GetNumOwnedDofHolder();
  owns_nothing_out = (num_owned == 0);
  if (owns_nothing_out) {
    return;
  }

  std::vector<Real> local_min(static_cast<size_t>(dim),
                              std::numeric_limits<Real>::max());
  std::vector<Real> local_max(static_cast<size_t>(dim),
                              std::numeric_limits<Real>::lowest());
  auto flat = CopyOwnedCoordsToHost(
    layout.GetOwnedDOFHolderCoordinates().GetValues(), num_owned, dim);
  for (LO i = 0; i < num_owned; ++i) {
    for (int d = 0; d < dim; ++d) {
      const Real v =
        flat[static_cast<size_t>(i) * static_cast<size_t>(dim) +
             static_cast<size_t>(d)];
      local_min[static_cast<size_t>(d)] =
        std::min(local_min[static_cast<size_t>(d)], v);
      local_max[static_cast<size_t>(d)] =
        std::max(local_max[static_cast<size_t>(d)], v);
    }
  }
  for (int d = 0; d < dim; ++d) {
    min_out[static_cast<size_t>(d)] = local_min[static_cast<size_t>(d)];
    max_out[static_cast<size_t>(d)] = local_max[static_cast<size_t>(d)];
  }
}

void DistributedBoundingBoxPartition::GetCandidateRanks(
  const std::array<Real, 3>& point, std::vector<int>& candidate_ranks) const
{
  for (int r = 0; r < num_ranks_; ++r) {
    if (owns_nothing_[static_cast<size_t>(r)]) {
      continue;
    }
    bool inside = true;
    for (int d = 0; d < dim_ && inside; ++d) {
      const size_t idx =
        static_cast<size_t>(r) * static_cast<size_t>(dim_) + static_cast<size_t>(d);
      inside = point[static_cast<size_t>(d)] >= mins_[idx] - padding_ &&
              point[static_cast<size_t>(d)] <= maxs_[idx] + padding_;
    }
    if (inside) {
      candidate_ranks.push_back(r);
    }
  }
}

} // namespace pcms
