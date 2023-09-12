#ifndef PCMS_COUPLING_COORDINATE_TRANSFORM_H
#define PCMS_COUPLING_COORDINATE_TRANSFORM_H
#include "coordinate.h"
#include "coordinate_systems.h"
namespace pcms {
template <typename CoordinateSystemTo, typename CoordinateSystemFrom>
constexpr Coordinate<CoordinateSystemTo> CoordinateTransform(Coordinate<CoordinateSystemFrom> original_coordinate) noexcept{
  static_assert(detail::dependent_always_false<CoordinateSystemTo,CoordinateSystemFrom>::value, "Missing coordinate transform specialization.");
}

template<>
constexpr Coordinate<Cartesian> CoordinateTransform(Coordinate<Cylindrical> coord) noexcept {
  auto r = coord[0];
  auto theta = coord[1];
  auto x = r*cos(theta);
  auto y = r*sin(theta);
  return {x,y,coord[2]};
}
template<>
constexpr Coordinate<Cylindrical> CoordinateTransform(Coordinate<Cartesian> coord) noexcept {
  auto x = coord[0];
  auto y = coord[1];
  auto r = sqrt(x*x+y*y);
  auto theta = atan2(y,x);
  return {r, theta, coord[2]};
}

}

#endif // PCMS_COUPLING_COORDINATE_TRANSFORM_H
