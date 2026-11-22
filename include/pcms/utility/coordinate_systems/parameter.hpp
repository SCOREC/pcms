#ifndef PCMS_UTILITY_COORDINATE_SYSTEMS_PARAMETER_HPP
#define PCMS_UTILITY_COORDINATE_SYSTEMS_PARAMETER_HPP

#include "pcms/utility/coordinate_system.hpp"
#include "pcms/utility/types.hpp"
#include <memory>
#include <string>
#include <vector>

namespace pcms::csys
{

/// @brief A parameter system: coordinates that are labels only (the logical
///        indices of a structured grid, flux coordinates), with no metric.
///
/// HasMetric() is false and Bind throws: the metric of a geometry over a
/// parameter system is the geometry's (an embedded geometry's given or
/// induced metric), not the coordinate system's. The basis is not known to be
/// orthogonal, so vector components in it are natural. A parameter system
/// equals only itself.
/// @param name what the coordinates are, e.g. "bout_index"; Kind() returns it
/// @param dimension number of coordinates, 1 to 6
/// @throws pcms_error if the dimension is out of range, or an identification
///         has an axis outside the system or a length that is not positive
[[nodiscard]] std::shared_ptr<const CoordinateSystem> CreateParameterSystem(
  std::string name, int dimension,
  std::vector<Identification> identifications = {});

} // namespace pcms::csys

#endif // PCMS_UTILITY_COORDINATE_SYSTEMS_PARAMETER_HPP
