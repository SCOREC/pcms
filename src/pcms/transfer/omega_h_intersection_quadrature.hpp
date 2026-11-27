#ifndef PCMS_TRANSFER_OMEGA_H_INTERSECTION_QUADRATURE_HPP
#define PCMS_TRANSFER_OMEGA_H_INTERSECTION_QUADRATURE_HPP

#include "pcms/backends/omega_h/omega_h_geometry.hpp"
#include "pcms/geometry/quadrature.hpp"
#include <memory>

namespace pcms::omega_h
{

/// @brief Quadrature of the intersection of two simplex meshes: the reference
///        rule of `order` on the sub-simplices of each target element
///        clipped against the source elements.
/// @param source the source geometry
/// @param target the target geometry
/// @param order the degree of the rule, at least 1
/// @return points bound on the target geometry, with `other` naming the
///         source element of each point
/// @throws pcms_error if the meshes' dimensions differ or there is no rule
///         of the order
Quadrature CreateIntersectionQuadrature(
  std::shared_ptr<const OmegaHGeometry> source,
  std::shared_ptr<const OmegaHGeometry> target, int order);

} // namespace pcms::omega_h

#endif // PCMS_TRANSFER_OMEGA_H_INTERSECTION_QUADRATURE_HPP
