#ifndef PCMS_UTILITY_COORDINATE_SYSTEMS_CARTESIAN_HPP
#define PCMS_UTILITY_COORDINATE_SYSTEMS_CARTESIAN_HPP

#include "pcms/utility/coordinate_system.hpp"
#include <memory>

namespace pcms::csys
{

/// Cartesian coordinates of dimension `dim`.
/// @throws pcms_error unless `dim` is 1, 2 or 3
[[nodiscard]] std::shared_ptr<const CoordinateSystem> CreateCartesian(int dim);

} // namespace pcms::csys

#endif // PCMS_UTILITY_COORDINATE_SYSTEMS_CARTESIAN_HPP
