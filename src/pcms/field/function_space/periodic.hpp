#ifndef PCMS_FIELD_FUNCTION_SPACE_PERIODIC_HPP
#define PCMS_FIELD_FUNCTION_SPACE_PERIODIC_HPP

#include "pcms/field/function_space.hpp"
#include "pcms/geometry/geometry.hpp"
#include <memory>

namespace pcms::detail
{

/// @brief The periodic cubic spline of the node values of a 1D uniform grid
///        identified over its extent.
/// @throws pcms_error for another geometry or fewer than 3 cells
[[nodiscard]] std::shared_ptr<const FunctionSpace> CreatePeriodicSplineSpace(
  std::shared_ptr<const Geometry> geometry, int num_components);

} // namespace pcms::detail

#endif // PCMS_FIELD_FUNCTION_SPACE_PERIODIC_HPP
