#ifndef PCMS_UTILITY_COORDINATE_SYSTEMS_CARTESIAN_HPP
#define PCMS_UTILITY_COORDINATE_SYSTEMS_CARTESIAN_HPP

#include "pcms/utility/coordinate_system.hpp"
#include <memory>
#include <vector>

namespace pcms::csys
{

/// Cartesian coordinates of dimension `dim`.
/// @throws pcms_error unless `dim` is 1, 2 or 3, or if an identification has
///         an axis outside the system or a length that is not positive
[[nodiscard]] std::shared_ptr<const CoordinateSystem> CreateCartesian(
  int dim, std::vector<Identification> identifications = {});

} // namespace pcms::csys

#endif // PCMS_UTILITY_COORDINATE_SYSTEMS_CARTESIAN_HPP
