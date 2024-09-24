#ifndef POINTS_HPP
#define POINTS_HPP

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>

struct Coord {
  double x, y;
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

struct Points {
  PointsViewType coordinates;
};

#endif
