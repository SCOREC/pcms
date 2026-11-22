#ifndef PCMS_GEOMETRY_PERIODIC_GEOMETRY_HPP
#define PCMS_GEOMETRY_PERIODIC_GEOMETRY_HPP

#include "pcms/geometry/geometry.hpp"
#include "pcms/utility/coordinate_system.hpp"
#include "pcms/utility/types.hpp"
#include <memory>
#include <numbers>
#include <optional>

namespace pcms
{

/// @brief A uniform periodic interval of `num_cells` cells: the cells wrap
///        (the last one is followed by the first) and query points wrap
///        modulo the period at Bind, so every point is found.
///
/// The interval is [origin, origin + period); cell c is
/// [origin + c h, origin + (c + 1) h), h = period / num_cells. The geometric
/// factors of its points are those of its coordinate system.
/// CreateQuadrature is Gauss-Legendre per cell.
/// @param system one-dimensional coordinate system of the interval
/// @param period the period, 2 pi unless given
/// @throws pcms_error unless `system` is one-dimensional, `num_cells` is
///         positive and the period is positive
[[nodiscard]] std::shared_ptr<const Geometry> CreatePeriodicGeometry(
  std::shared_ptr<const CoordinateSystem> system, int num_cells,
  std::optional<Real> period = std::nullopt, Real origin = 0);

/// @brief The equispaced rule of `num_points` points, each weighing
///        period / num_points: exact for the trigonometric polynomials of
///        degree up to num_points - 1.
/// @throws pcms_error if `geometry` is not a periodic geometry or
///         `num_points` is not positive
[[nodiscard]] Quadrature CreateTrapezoidQuadrature(
  const std::shared_ptr<const Geometry>& geometry, int num_points);

} // namespace pcms

#endif // PCMS_GEOMETRY_PERIODIC_GEOMETRY_HPP
