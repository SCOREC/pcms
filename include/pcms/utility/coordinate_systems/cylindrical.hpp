#ifndef PCMS_UTILITY_COORDINATE_SYSTEMS_CYLINDRICAL_HPP
#define PCMS_UTILITY_COORDINATE_SYSTEMS_CYLINDRICAL_HPP

#include "pcms/utility/coordinate_system.hpp"
#include "pcms/utility/types.hpp"
#include <memory>
#include <numbers>

namespace pcms::csys
{

/// Cylindrical (R, theta, Z) coordinates, with theta = atan2(y, x), identified
/// along theta over [0, phi_period).
/// @param phi_period 2 pi, or the toroidal wedge of a field period
/// @throws pcms_error if `phi_period` is not positive
[[nodiscard]] std::shared_ptr<const CoordinateSystem> CreateCylindricalRThetaZ(
  Real phi_period = 2 * std::numbers::pi_v<Real>);

/// Cylindrical (R, Z) coordinates of a poloidal plane.
[[nodiscard]] std::shared_ptr<const CoordinateSystem> CreateCylindricalRZ();

} // namespace pcms::csys

#endif // PCMS_UTILITY_COORDINATE_SYSTEMS_CYLINDRICAL_HPP
