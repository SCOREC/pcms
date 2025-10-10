#ifndef PCMS_INTERPOLATOR_ALIASES_HPP
#define PCMS_INTERPOLATOR_ALIASES_HPP

#include <Kokkos_Core.hpp>
#include "pcms/arrays.h"

namespace pcms
{

struct Coord
{
  double x, y, z;
};

using range_policy = typename Kokkos::RangePolicy<>;
using team_policy = typename Kokkos::TeamPolicy<>;
using member_type = typename Kokkos::TeamPolicy<>::member_type;

// alias for scratch view

using ScratchSpace =
  typename Kokkos::DefaultExecutionSpace::scratch_memory_space;

using ScratchMatView =
  typename Kokkos::View<double**, Kokkos::LayoutRight, ScratchSpace,
                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
using ScratchVecView =
  typename Kokkos::View<double*, ScratchSpace,
                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

using PointsViewType = Kokkos::View<Coord*>;

struct Points
{
  PointsViewType coordinates;
};

using IntDeviceMatView = Kokkos::View<int**>;
using IntDeviceVecView = Kokkos::View<int*>;
using IntHostMatView = Kokkos::View<int**, Kokkos::HostSpace>;

using RealDefaultRank1View = Rank1View<double, DefaultExecutionSpace>;
using RealConstDefaultRank1View =
  Rank1View<const double, DefaultExecutionSpace>;

} // namespace pcms
#endif
