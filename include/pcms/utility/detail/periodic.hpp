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

} // namespace pcms::detail

#endif // PCMS_UTILITY_DETAIL_PERIODIC_HPP
