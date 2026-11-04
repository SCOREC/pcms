#ifndef PCMS_SPLINE_FUNCTION_SPACE_H
#define PCMS_SPLINE_FUNCTION_SPACE_H

#include "pcms/field/function_space.hpp"
#include "pcms/geometry/geometry.hpp"

#include <memory>

namespace pcms
{

/// A bicubic spline space on a 2D uniform-grid geometry. Spaces on one
/// geometry accept each other's bound sources.
/// @throws pcms_error if `geometry` is null or not a 2D uniform-grid geometry
[[nodiscard]] std::shared_ptr<const FunctionSpace> CreateSplineSpace(
  std::shared_ptr<const Geometry> geometry);

} // namespace pcms

#endif // PCMS_SPLINE_FUNCTION_SPACE_H
