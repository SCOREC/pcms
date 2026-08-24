#ifndef PCMS_FIELD_POINT_STATUS_HPP
#define PCMS_FIELD_POINT_STATUS_HPP

#include "pcms/utility/memory_spaces.h"
#include <Kokkos_Core.hpp>
#include <cstdint>

namespace pcms
{

// Per-point outcome of a bound coordinate map
// OutsideDomain is handled by OutOfBoundsPolicy and is not necessarily an error
enum class PointStatus : std::int8_t
{
  Valid = 0,
  SingularDifferential,
  OutsideDomain,
  NonConverged,
  Ambiguous
};

using PointStatusView = Kokkos::View<const PointStatus*, DeviceMemorySpace>;

[[nodiscard]] KOKKOS_INLINE_FUNCTION bool PointUsable(PointStatus s) noexcept
{
  return s == PointStatus::Valid || s == PointStatus::SingularDifferential;
}

[[nodiscard]] KOKKOS_INLINE_FUNCTION bool DifferentialUsable(
  PointStatus s) noexcept
{
  return s == PointStatus::Valid;
}

} // namespace pcms

#endif // PCMS_FIELD_POINT_STATUS_HPP
