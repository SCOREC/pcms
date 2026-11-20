#ifndef PCMS_FIELD_EMBEDDING_HPP
#define PCMS_FIELD_EMBEDDING_HPP

#include "pcms/field/coordinate_map.hpp"
#include "pcms/field/field.hpp"
#include "pcms/field/function_space.hpp"
#include <memory>

namespace pcms
{

/// @brief The coordinate map of a position Function: x(xi) = position(xi),
///        from the coordinate system of the position's (parameter) geometry to
///        the coordinate system of its ValueBasis::Coordinates basis.
///
/// Bind evaluates the position and its Gradient at the parameter points (no
/// inverse): the mapped points and the Jacobian J(p, j, i) = d x_j / d xi_i.
/// Values evaluated at the mapped points flow against the map
/// (BoundCoordinateMap::ValueTransform).
/// @param position a Function with a ValueBasis::Coordinates basis, of a space
///        whose evaluator supports Gradient; the parameter geometry and the
///        coordinate system have at most 3 dimensions
/// @throws pcms_error if the position is null, its values are not
///         coordinates, its space has no Gradient or a dimension is above 3
[[nodiscard]] std::shared_ptr<const CoordinateMap> CreateCoordinateMap(
  std::shared_ptr<const Function<Real>> position);

} // namespace pcms

#endif // PCMS_FIELD_EMBEDDING_HPP
