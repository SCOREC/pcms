#ifndef WDM_COUPLING_COORDINATE_TRANSFORM_H
#define WDM_COUPLING_COORDINATE_TRANSFORM_H
#include "coordinate.h"
#include "coordinate_systems.h"
namespace wdmcpl {
template <typename... T> struct dependent_always_false : std::false_type {};

template <typename CoordinateSystemTo, typename CoordinateSystemFrom, typename T, std::enable_if_t<!std::is_same_v<CoordinateSystemTo,CoordinateSystemFrom> && std::is_arithmetic_v<T>,bool> = false >
constexpr Coordinate<CoordinateSystemTo, T> CoordinateTransform(const Coordinate<CoordinateSystemFrom, T> &original_coordinate) noexcept{
  static_assert(dependent_always_false<CoordinateSystemTo,CoordinateSystemFrom>::value, "Missing coordinate transform specialization.");
}
// doing a coordinate transform between the same coordinate systems is a null operation
/*
template<typename T>
constexpr Coordinate<T>& CoordinateTransform(Coordinate<T> & coord) noexcept { return coord; }
 */
// TODO figure out enable_if to disable when T!=Coordinate<>
//template <typename T>
//constexpr T&& CoordinateTransform(T&& coord) noexcept {return std::forward<T>(coord);}
// FIXME Forwarding Null coordinate transform such that no conversions ever necessary...
template <typename CoordinateSystemTo, typename CoordinateSystemFrom, typename T, std::enable_if_t<std::is_same_v<CoordinateSystemTo,CoordinateSystemFrom>,bool> = false >
constexpr const Coordinate<CoordinateSystemTo, T>& CoordinateTransform(const Coordinate<CoordinateSystemFrom, T> &original_coordinate) noexcept{
  return original_coordinate;
}
template <typename CoordinateSystemTo, typename CoordinateSystemFrom, typename T, std::enable_if_t<std::is_same_v<CoordinateSystemTo,CoordinateSystemFrom>,bool> = false >
constexpr Coordinate<CoordinateSystemTo, T>& CoordinateTransform(Coordinate<CoordinateSystemFrom,T> &original_coordinate) noexcept{
  return original_coordinate;
}
template <typename CoordSystem, typename T>
constexpr ScalarData<T>& DataTransform(ScalarData<T> &data) noexcept{ return data; }
template <typename CoordSystem, typename T>
constexpr const ScalarData<T>& DataTransform(const ScalarData<T> &data) noexcept{ return data; }

// if data type not a scalar...
// R CoordinateTransform(T) { }
template<typename CoordinateSystemTo, typename CoordinateSystemFrom, typename T, std::enable_if_t<!std::is_arithmetic_v<T> && !std::is_same_v<CoordinateSystemTo,CoordinateSystemFrom>,bool> = false>
constexpr Coordinate<CoordinateSystemTo, T> CoordinateTransform(const Coordinate<CoordinateSystemFrom,T>& coord) noexcept {
  Coordinate<CoordinateSystemTo, T> new_coordinates;
  // FIXME need loop here

  return new_coordinates;
}

template<typename CoordinateSystemTo, typename T, std::enable_if_t<std::is_arithmetic_v<T> && std::is_same_v<CoordinateSystemTo, Cartesian> ,bool> = false>
constexpr Coordinate<Cartesian, T> CoordinateTransform(const Coordinate<Cylindrical,T>& coord) noexcept {
  auto r = coord[0];
  auto theta = coord[1];
  auto x = r*cos(theta);
  auto y = r*sin(theta);
  return {{x,y,coord[2]}};
}
template<typename CoordinateSystemTo, typename T, std::enable_if_t<std::is_arithmetic_v<T>&& std::is_same_v<CoordinateSystemTo, Cylindrical> ,bool> = false>
constexpr Coordinate<Cylindrical,T> CoordinateTransform(const Coordinate<Cartesian, T>& coord) noexcept {
  auto x = coord[0];
  auto y = coord[1];
  auto r = sqrt(x*x+y*y);
  auto theta = atan2(y,x);
  return {{r, theta, coord[2]}};
}

}

#endif // WDM_COUPLING_COORDINATE_TRANSFORM_H
