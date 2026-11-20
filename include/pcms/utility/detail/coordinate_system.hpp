#ifndef PCMS_UTILITY_DETAIL_COORDINATE_SYSTEM_HPP
#define PCMS_UTILITY_DETAIL_COORDINATE_SYSTEM_HPP

#include "pcms/utility/coordinate_system.hpp"
#include <memory>

namespace pcms::detail
{

[[nodiscard]] inline bool SameCoordinateSystem(
  const std::shared_ptr<const CoordinateSystem>& a,
  const std::shared_ptr<const CoordinateSystem>& b) noexcept
{
  if (a == nullptr || b == nullptr) {
    return a == b;
  }
  return *a == *b;
}

/// @return whether the physical components of `system` are its natural ones
///         times per-point scale factors sqrt(g_ii): an orthogonal system
///         whose metric varies
[[nodiscard]] inline bool HasScaleFactors(
  const CoordinateSystem& system) noexcept
{
  return system.IsOrthogonal() && system.HasMetric() && !system.IsNormalized();
}

/// A Cartesian system whose dimension is not yet known ("cartesian
/// (dimension deferred)"): ResolveCoordinateSystem replaces it by the
/// Cartesian system of the data's dimension.
[[nodiscard]] std::shared_ptr<const CoordinateSystem> DeferredCartesian();

/// @return whether `system` is Cartesian, of any dimension (deferred
///         included)
[[nodiscard]] bool IsCartesian(const CoordinateSystem& system) noexcept;

[[nodiscard]] std::shared_ptr<const CoordinateSystem> ResolveCoordinateSystem(
  std::shared_ptr<const CoordinateSystem> system, int dim);

} // namespace pcms::detail

#endif // PCMS_UTILITY_DETAIL_COORDINATE_SYSTEM_HPP
