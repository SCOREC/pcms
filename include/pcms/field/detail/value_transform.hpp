// Part of pcms/field/value_transform.hpp, which includes it after the
// declarations it needs.
#include "pcms/field/value_transform.hpp"
#ifndef PCMS_FIELD_DETAIL_VALUE_TRANSFORM_HPP
#define PCMS_FIELD_DETAIL_VALUE_TRANSFORM_HPP

#include "pcms/utility/arrays.hpp"
#include "pcms/utility/coordinate_system.hpp"
#include "pcms/utility/memory_spaces.hpp"
#include "pcms/utility/types.hpp"
#include <Kokkos_Core.hpp>
#include <memory>
#include <optional>

namespace pcms::detail
{

/// A Jacobian at points, J(p, j, i) = d x_j / d xi_i for a map xi -> x; a
/// stride of 0 along the points means the same matrix at every point.
using JacobianView =
  Rank3View<const Real, DeviceMemorySpace, Kokkos::layout_stride>;

/// The metric g(p, a, b) of a coordinate system at points.
using MetricView = BoundCoordinateSystem::MetricView;

/// @brief The value transform of a map whose values flow against it, from
///        `from` (the map's target) to `to` (its source): vectors by J^-1,
///        which needs J square, and covectors by J^T.
/// @param jacobian J(p, j, i) = d from_j / d to_i, [points][dim from][dim to]
/// @throws pcms_error if the extents of `jacobian` do not match the systems
[[nodiscard]] ValueTransform Pullback(
  JacobianView jacobian, std::shared_ptr<const CoordinateSystem> from,
  std::shared_ptr<const CoordinateSystem> to);

/// @brief The value transform along a map, from `from` (its source) to `to`
///        (its target): vectors by J, and covectors by J^-T, which needs J
///        square.
/// @param jacobian J(p, j, i) = d to_j / d from_i, [points][dim to][dim from]
/// @throws pcms_error if the extents of `jacobian` do not match the systems
[[nodiscard]] ValueTransform Pushforward(
  JacobianView jacobian, std::shared_ptr<const CoordinateSystem> from,
  std::shared_ptr<const CoordinateSystem> to);

/// @brief Physical to natural components of `system`: a vector component
///        divided by sqrt(g_ii), a covector component by sqrt(g^ii).
/// @param metric the metric of `system` at the points
[[nodiscard]] ValueTransform ToNatural(
  std::shared_ptr<const CoordinateSystem> system, MetricView metric);

/// @brief Natural to physical components of `system`: a vector component
///        multiplied by sqrt(g_ii), a covector component by sqrt(g^ii).
/// @param metric the metric of `system` at the points
[[nodiscard]] ValueTransform ToPhysical(
  std::shared_ptr<const CoordinateSystem> system, MetricView metric);

/// @brief Vectors to covectors by g_ij and covectors to vectors by g^ij;
///        physical components by g_ij / sqrt(g_ii g_jj) and
///        g^ij / sqrt(g^ii g^jj).
/// @param metric the metric of `system` at the points
/// @param scaling the scaling of the values it takes and writes
[[nodiscard]] ValueTransform LowerIndex(
  std::shared_ptr<const CoordinateSystem> system, MetricView metric,
  ComponentScaling scaling);

/// @brief The metric J^T G J that the map with Jacobian J pulls back from the
///        metric G at the images of the points.
/// @param jacobian J(p, j, i), [points][dim G][dim of the result]
/// @param metric G(p, j, k) at the images
/// @return [points][dim][dim]
[[nodiscard]] Kokkos::View<Real***, DeviceMemorySpace> PullbackMetric(
  JacobianView jacobian, MetricView metric);

/// @brief The transform of values written in `written` to the scaling of
///        `declared` at `points`: ToPhysical or ToNatural from the metric of
///        the declared system there.
/// @return nothing unless the bases differ only in their scaling
/// @throws pcms_error if the system has no metric, or `points` are not in it
[[nodiscard]] std::optional<ValueTransform> ScalingTransform(
  const ValueBasis& written, const ValueBasis& declared,
  const CoordinateView<DeviceMemorySpace>& points);

/// @brief ScalingTransform kept in `cache`, formed again only when the cached
///        one does not take `written` to `declared`.
/// @return the transform, or null if there is none
inline const ValueTransform* CachedScalingTransform(
  std::optional<ValueTransform>& cache, const ValueBasis& written,
  const ValueBasis& declared, const CoordinateView<DeviceMemorySpace>& points)
{
  if (SameValueBasis(written, declared)) {
    return nullptr;
  }
  if (!cache || !SameValueBasis(cache->OutputBasis(written), declared)) {
    cache = ScalingTransform(written, declared, points);
  }
  return cache ? &*cache : nullptr;
}

} // namespace pcms::detail

#endif // PCMS_FIELD_DETAIL_VALUE_TRANSFORM_HPP
