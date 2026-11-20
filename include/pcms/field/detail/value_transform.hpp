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

namespace pcms::detail
{

/// A Jacobian at points, J(p, j, i) = d x_j / d xi_i for a map xi -> x; a
/// stride of 0 along the points means the same matrix at every point.
using JacobianView =
  View<3, const Real, DeviceMemorySpace, Kokkos::layout_stride>;

/// The metric g(p, a, b) of a coordinate system at points.
using MetricView = BoundCoordinateSystem::MetricView;

/// @brief The value transform of a map whose values flow against it: values
///        in the components of `from` (the map's target, where the points
///        were located) are expressed in those of `to` (the map's source,
///        the coordinate system of the query points). Covectors by J^T, with no
///        inverse (out_i = sum_j J(p, j, i) in_j); vectors by J^-1 (the
///        solution of in_j = sum_i J(p, j, i) out_i), which needs J square.
/// @param jacobian J(p, j, i) = d from_j / d to_i, [points][dim from][dim to];
///        it must outlive the transform
/// @param from the map's target coordinate system
/// @param metric the metric of `from` at the mapped points, or empty if it
///        has none; the scale factors of `to` are those of J^T g J, the
///        metric pulled back to the query points
/// @param to the map's source coordinate system
/// @throws pcms_error if the extents of `jacobian` or `metric` do not match
///         the coordinate systems
[[nodiscard]] ValueTransform Pullback(
  JacobianView jacobian, std::shared_ptr<const CoordinateSystem> from,
  MetricView metric, std::shared_ptr<const CoordinateSystem> to);

/// @brief The value transform along a map: components on the coordinates of
///        `from` (the parameter system of a geometry) expressed in those of
///        `to` (its coordinate system). Vectors by J (out_j = sum_i J(p, j, i)
///        in_i, no inverse), covectors by J^-T (the solution of in_i = sum_j
///        J(p, j, i) out_j), which needs J square.
/// @param jacobian J(p, j, i) = d to_j / d from_i, [points][dim to][dim from];
///        it must outlive the transform
/// @param from the parameter system
/// @param to the coordinate system
/// @param metric the metric of `to` at the points, or empty if it has none;
///        the scale factors of `from` are those of J^T g J
/// @throws pcms_error if the extents of `jacobian` or `metric` do not match
///         the coordinate systems
[[nodiscard]] ValueTransform Pushforward(
  JacobianView jacobian, std::shared_ptr<const CoordinateSystem> from,
  std::shared_ptr<const CoordinateSystem> to, MetricView metric);

/// @brief Natural to physical components of an orthogonal coordinate system, or
/// back: a
///        vector component is multiplied by the scale factor sqrt(g_ii) of
///        its axis, a covector component divided by it.
/// @param system the coordinate system of the values
/// @param metric the coordinate system's metric at the points
/// @param to the scaling of the output
[[nodiscard]] ValueTransform Scale(
  std::shared_ptr<const CoordinateSystem> system, MetricView metric,
  ComponentScaling to);

} // namespace pcms::detail

#endif // PCMS_FIELD_DETAIL_VALUE_TRANSFORM_HPP
