#ifndef PCMS_GEOMETRY_PRODUCT_GEOMETRY_HPP
#define PCMS_GEOMETRY_PRODUCT_GEOMETRY_HPP

#include "pcms/geometry/geometry.hpp"
#include "pcms/utility/coordinate_system.hpp"
#include <memory>

namespace pcms
{

/// @brief The product of two geometries: points (a, b) with a in `a` and b in
///        `b`, the coordinates of `a` first (nest products for more factors),
///        in a coordinate system without a metric.
///
/// Bind splits the query coordinates by factor dimension and binds each
/// factor; its BoundPoints keep the factor BoundPoints (ProductBinding in
/// the library's private header) for the spaces of the product, which
/// evaluate factor by factor. A point is found when it is found in both
/// factors. CreateQuadrature is the tensor product of the factors' rules. The
/// product holds no partitioning: it is serial, its factors are not
/// distributed.
/// @param system the product's coordinate system, e.g. a parameter system
///        with a twist-shift identification; null for a parameter system
///        with the factors' identifications
/// @throws pcms_error if a geometry is null, the dimensions add up to more
///         than 6, `system` has a metric or another dimension, or `system`
///         is null and `b` has a twist-shift identification
[[nodiscard]] std::shared_ptr<const Geometry> CreateProductGeometry(
  std::shared_ptr<const Geometry> a, std::shared_ptr<const Geometry> b,
  std::shared_ptr<const CoordinateSystem> system = nullptr);

} // namespace pcms

#endif // PCMS_GEOMETRY_PRODUCT_GEOMETRY_HPP
