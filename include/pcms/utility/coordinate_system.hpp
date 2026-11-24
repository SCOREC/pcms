#ifndef PCMS_UTILITY_COORDINATE_SYSTEM_HPP
#define PCMS_UTILITY_COORDINATE_SYSTEM_HPP

#include "pcms/utility/arrays.hpp"
#include "pcms/utility/memory_spaces.hpp"
#include "pcms/utility/types.hpp"
#include <Kokkos_Core.hpp>
#include <memory>
#include <span>
#include <string_view>

namespace pcms
{

/// @brief A coordinate system at a set of points: the metric tensor and the
///        volume factor there.
///
/// Created by CoordinateSystem::Bind. A system whose metric (volume factor)
/// is the same at every point allocates one entry for all of them and
/// returns a view with stride 0 along the points.
class BoundCoordinateSystem
{
public:
  using MetricView =
    Rank3View<const Real, DeviceMemorySpace, Kokkos::layout_stride>;
  using FactorView =
    Rank1View<const Real, DeviceMemorySpace, Kokkos::layout_stride>;

  /// @param num_points number of points
  /// @param metric [points or 1][dimension][dimension]
  /// @param volume_factor [points or 1]
  /// An extent of 1 along the points means the same value at every point.
  BoundCoordinateSystem(
    LO num_points, Kokkos::View<const Real***, DeviceMemorySpace> metric,
    Kokkos::View<const Real*, DeviceMemorySpace> volume_factor)
    : num_points_(num_points),
      metric_(std::move(metric)),
      volume_factor_(std::move(volume_factor))
  {
  }

  /// @return g(p, i, j) = e_i . e_j, the metric tensor in the coordinate basis
  [[nodiscard]] MetricView Metric() const noexcept;

  /// @return the factor of the physical volume element relative to
  ///         dx^1 ... dx^n at each point: sqrt(det g)(p), except for a
  ///         coordinate system reduced by a symmetry, whose factor includes the
  ///         ignored direction (CylindricalRZ: R, the volume per radian of the
  ///         toroidal angle)
  [[nodiscard]] FactorView VolumeFactor() const noexcept;

private:
  LO num_points_;
  Kokkos::View<const Real***, DeviceMemorySpace> metric_;
  Kokkos::View<const Real*, DeviceMemorySpace> volume_factor_;
};

class CoordinateSystem;

/// @brief The transform of a periodic axis that is not a translation
///        (CreateTwistShiftTransform)
class PeriodicTransform
{
public:
  /// @brief Applies the map `wraps(p)` times to point p.
  /// @param points [points][dimension], the identified axis already reduced
  /// @param wraps the number of lengths each point was moved down the
  ///        identified axis
  virtual void Apply(Rank2View<Real, DeviceMemorySpace> points,
                     Rank1View<const Real, DeviceMemorySpace> wraps) const = 0;

  virtual ~PeriodicTransform() = default;
};

/// @brief Points of a coordinate system that are one point: x and x + length
///        along `axis`, after `map` for a twist-shift.
struct PeriodicAxis
{
  int axis = 0;
  /// the fundamental interval [lower, lower + length) of the axis
  Real lower = 0;
  Real length = 0;
  /// what else changes per length crossed; null for a translation
  std::shared_ptr<const PeriodicTransform> transform;

  [[nodiscard]] bool operator==(const PeriodicAxis&) const = default;
};

/// @brief A coordinate system: labels and closed forms, never data.
///
/// A coordinate system says what its coordinates are (kind, dimension,
/// periodic axes) and, when it has one, its closed-form metric (Bind);
/// per-point values come only from bound objects. A parameter system
/// (csys::CreateParameterSystem) has no metric: the metric of a geometry over
/// it belongs to the geometry.
class CoordinateSystem
{
public:
  [[nodiscard]] virtual std::string_view Kind() const noexcept = 0;

  [[nodiscard]] virtual int Dimension() const noexcept = 0;

  /// @return whether the metric is the identity everywhere, so that
  ///         physical and natural components coincide
  [[nodiscard]] virtual bool IsOrthonormal() const noexcept = 0;

  /// @return the periodic axes, applied in order by Geometry::Bind
  [[nodiscard]] virtual std::span<const PeriodicAxis> PeriodicAxes()
    const noexcept
  {
    return {};
  }

  /// @return whether the coordinate system has a closed-form metric, so that
  /// Bind
  ///         gives the metric and volume factor at points
  [[nodiscard]] virtual bool HasMetric() const noexcept { return true; }

  /// @brief The closed-form metric and volume factor at points.
  /// @param points [points][dimension] coordinates in this system
  /// @throws pcms_error if the coordinate system has no metric (HasMetric() is
  /// false)
  [[nodiscard]] virtual std::unique_ptr<BoundCoordinateSystem> Bind(
    Rank2View<const Real, DeviceMemorySpace> points) const = 0;

  [[nodiscard]] virtual bool operator==(
    const CoordinateSystem& other) const noexcept = 0;

  virtual ~CoordinateSystem() = default;
};

} // namespace pcms

#include "pcms/utility/detail/coordinate_system.hpp"

#endif // PCMS_UTILITY_COORDINATE_SYSTEM_HPP
