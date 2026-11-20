#ifndef PCMS_FIELD_HOLDER_OF_GID_HPP
#define PCMS_FIELD_HOLDER_OF_GID_HPP

// Private: the map from the global id of a DOF holder to the holder.

#include "pcms/field/field_layout.hpp"
#include "pcms/utility/memory_spaces.hpp"
#include "pcms/utility/types.hpp"
#include <Kokkos_Core.hpp>
#include <algorithm>
#include <numeric>
#include <vector>

namespace pcms::detail
{

/// Maps the global id of a DOF holder (of a layout without ghosts) to the
/// holder: a binary search in the sorted gids.
struct HolderOfGid
{
  Kokkos::View<GO*, DeviceMemorySpace> sorted;
  Kokkos::View<LO*, DeviceMemorySpace> holders;

  KOKKOS_INLINE_FUNCTION LO operator()(GO gid) const
  {
    LO lo = 0;
    LO hi = static_cast<LO>(sorted.extent(0));
    while (lo < hi) {
      const LO mid = lo + (hi - lo) / 2;
      if (sorted(mid) < gid) {
        lo = mid + 1;
      } else {
        hi = mid;
      }
    }
    return holders(lo);
  }
};

inline HolderOfGid MakeHolderOfGid(const FieldLayout& layout)
{
  const auto gids = layout.GetGidsHost();
  std::vector<LO> order(gids.extent(0));
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(),
            [&](LO a, LO b) { return gids(a) < gids(b); });
  Kokkos::View<GO*, HostMemorySpace> sorted("sorted_gids", order.size());
  Kokkos::View<LO*, HostMemorySpace> holders("sorted_holders", order.size());
  for (std::size_t i = 0; i < order.size(); ++i) {
    sorted(i) = gids(order[i]);
    holders(i) = order[i];
  }
  return {Kokkos::create_mirror_view_and_copy(DeviceMemorySpace{}, sorted),
          Kokkos::create_mirror_view_and_copy(DeviceMemorySpace{}, holders)};
}

} // namespace pcms::detail

#endif // PCMS_FIELD_HOLDER_OF_GID_HPP
