#ifndef PCMS_TRANSFER_OMEGA_H_INTERSECTION_QUADRATURE_HPP
#define PCMS_TRANSFER_OMEGA_H_INTERSECTION_QUADRATURE_HPP

#include "pcms/geometry/geometry.hpp"
#include "pcms/geometry/quadrature.hpp"
#include <memory>

namespace pcms::omega_h
{

/// @brief Quadrature of the intersection of two simplex meshes: the reference
///        rule of `order` on the sub-simplices of each target element
///        clipped against the source elements.
/// @param source an Omega_h geometry (omega_h::CreateGeometry)
/// @param target an Omega_h geometry
/// @param order the degree of the rule
/// @return points bound on the target geometry, with `other` naming the
///         source element of each point
/// @throws pcms_error if a geometry is not an Omega_h geometry, the meshes'
///         dimensions differ, there is no rule of the order, or in a build
///         without PCMS_ENABLE_OMEGA_H_MESHFIELDS
[[nodiscard]] Quadrature CreateIntersectionQuadrature(
  std::shared_ptr<const Geometry> source,
  std::shared_ptr<const Geometry> target, int order);

} // namespace pcms::omega_h

#endif // PCMS_TRANSFER_OMEGA_H_INTERSECTION_QUADRATURE_HPP
