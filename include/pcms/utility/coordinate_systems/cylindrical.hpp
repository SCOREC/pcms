#ifndef PCMS_UTILITY_COORDINATE_SYSTEMS_CYLINDRICAL_HPP
#define PCMS_UTILITY_COORDINATE_SYSTEMS_CYLINDRICAL_HPP

#include "pcms/utility/coordinate_system.hpp"
#include <memory>

namespace pcms::csys
{

/// Cylindrical (R, theta, Z) coordinates, with theta = atan2(y, x).
[[nodiscard]] std::shared_ptr<const CoordinateSystem>
CreateCylindricalRThetaZ();

/// Cylindrical (R, Z) coordinates of a poloidal plane.
[[nodiscard]] std::shared_ptr<const CoordinateSystem> CreateCylindricalRZ();

} // namespace pcms::csys

#endif // PCMS_UTILITY_COORDINATE_SYSTEMS_CYLINDRICAL_HPP
