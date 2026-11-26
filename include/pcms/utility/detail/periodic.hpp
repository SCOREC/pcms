#ifndef PCMS_UTILITY_DETAIL_PERIODIC_HPP
#define PCMS_UTILITY_DETAIL_PERIODIC_HPP

#include "pcms/utility/types.hpp"
#include <Kokkos_Core.hpp>

namespace pcms::detail
{

/// @return the image of the difference `d` nearest zero, d - length *
///         round(d / length) (ties away from zero), or `d` if length <= 0
KOKKOS_INLINE_FUNCTION Real NearestImage(Real d, Real length)
{
  return length > 0 ? d - length * Kokkos::round(d / length) : d;
}

/// @return x moved by whole lengths into [lower, lower + length)
/// @param wraps the number of lengths x moved down
KOKKOS_INLINE_FUNCTION Real WrapToPeriod(Real x, Real lower, Real length,
                                         Real& wraps)
{
  Real n = Kokkos::floor((x - lower) / length);
  Real y = x - n * length;
  // rounding can leave y at the upper end or just below the lower end
  if (y >= lower + length) {
    y -= length;
    n += 1;
  }
  if (y < lower) {
    y = lower;
  }
  wraps = n;
  return y;
}

} // namespace pcms::detail

#endif // PCMS_UTILITY_DETAIL_PERIODIC_HPP
