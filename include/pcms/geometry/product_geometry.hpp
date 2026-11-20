#ifndef PCMS_GEOMETRY_PRODUCT_GEOMETRY_HPP
#define PCMS_GEOMETRY_PRODUCT_GEOMETRY_HPP

#include "pcms/field/field.hpp"
#include "pcms/geometry/geometry.hpp"
#include "pcms/utility/types.hpp"
#include <memory>
#include <optional>

namespace pcms
{

/// @brief A twist-shift identification of the factors of a product: when the
///        coordinate of `a` along `wrap_axis` wraps around its period, the
///        coordinate of `b` is shifted by a function of the other coordinates
///        of `a`.
///
/// (.., y + period, ..; z) is the point (.., y, ..; z + shift(..)) for the
/// axis y of `a` and the (first) axis z of `b`: a query point whose y is n
/// periods outside [lower, lower + period) is moved to y - n period and
/// z + n shift. `b` wraps its own coordinate (a periodic geometry); `a` need
/// not be periodic, it is only ever asked for points inside the interval.
/// Quadrature points lie in the fundamental domain and are unaffected. A first
/// version for the maintainer to try; the form may change.
struct ProductIdentification
{
  /// axis of the first factor that wraps
  int wrap_axis = 0;
  /// the wrapping interval of that axis: the period is the one the first
  /// factor's coordinate system declares for the axis (CoordinateSystem::
  /// Period), if it does, and may then be left 0; the lower end is a property
  /// of the geometry
  Real lower = 0;
  Real period = 0;
  /// the shift, a function on a geometry whose coordinates are the other axes
  /// of `a` (in order), Cartesian
  std::shared_ptr<const Function<Real>> shift;
};

/// @brief The product of two geometries: points (a, b) with a in `a` and b in
///        `b`, in the Cartesian system of the sum of their dimensions, the
///        coordinates of `a` first (nest products for more factors).
///
/// Bind splits the query coordinates by factor dimension and binds each
/// factor; its BoundPoints keep the factor BoundPoints (ProductBinding in
/// the library's private header) for the spaces of the product, which
/// evaluate factor by factor. A point is found when it is found in both
/// factors. CreateQuadrature is the tensor product of the factors' rules, and
/// the differential is block-diagonal. The product holds no partitioning: it
/// is serial, its factors are not distributed.
/// @param identification the twist-shift of the factors, if any
/// @throws pcms_error if a geometry is null, a factor's system is not
///         Cartesian, the dimensions add up to more than 3, or the
///         identification's axis, period or shift function do not fit `a`
[[nodiscard]] std::shared_ptr<const Geometry> CreateProductGeometry(
  std::shared_ptr<const Geometry> a, std::shared_ptr<const Geometry> b,
  std::optional<ProductIdentification> identification = std::nullopt);

} // namespace pcms

#endif // PCMS_GEOMETRY_PRODUCT_GEOMETRY_HPP
