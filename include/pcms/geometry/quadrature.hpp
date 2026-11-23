#ifndef PCMS_GEOMETRY_QUADRATURE_HPP
#define PCMS_GEOMETRY_QUADRATURE_HPP

#include "pcms/geometry/bound_points.hpp"
#include "pcms/utility/memory_spaces.hpp"
#include "pcms/utility/types.hpp"
#include <Kokkos_Core.hpp>
#include <cstdint>
#include <memory>
#include <vector>

namespace pcms
{

class Geometry;

/// @brief The seed of the Monte Carlo rules when none is given.
inline constexpr std::uint64_t DefaultQuadratureSeed = 8675309;

/// @brief Points with weights such that the integral of f over the cells is
///        sum_p weights(p) f(x_p).
struct Quadrature
{
  /// the points, bound on the geometry that generated them
  std::shared_ptr<const BoundPoints> points;
  /// one weight per point, in the physical measure
  Kokkos::View<const Real*, DeviceMemorySpace> weights;
  /// the same points bound on other geometries
  std::vector<std::shared_ptr<const BoundPoints>> also;
};

/// @brief The Gauss rule of each cell of `geometry` (on simplices, the
///        reference rule), exact for polynomials of degree `order`.
/// @throws pcms_error if the geometry has no such rule of that order
[[nodiscard]] Quadrature CreateGaussQuadrature(
  const std::shared_ptr<const Geometry>& geometry, int order);

/// @brief The equispaced rule of `num_points` points of a 1D uniform grid,
///        each weighing extent / num_points: exact for the trigonometric
///        polynomials of degree up to num_points - 1 when the axis is
///        identified over the extent.
/// @throws pcms_error if `geometry` is not a 1D uniform grid or `num_points`
///         is not positive
[[nodiscard]] Quadrature CreateTrapezoidQuadrature(
  const std::shared_ptr<const Geometry>& geometry, int num_points);

/// @brief Uniform random samples in each cell of `geometry`, each weighing
///        cell measure / samples_per_cell; a cell's samples depend only on
///        `seed` and its global id.
/// @throws pcms_error if the geometry has no Monte Carlo rule or
///         `samples_per_cell` is not positive
[[nodiscard]] Quadrature CreateMonteCarloQuadrature(
  const std::shared_ptr<const Geometry>& geometry, int samples_per_cell,
  std::uint64_t seed = DefaultQuadratureSeed);

/// @brief The points of `quadrature`, a rule of the parameter geometry of
///        `geometry`, bound on `geometry`, with the weights times the volume
///        factor at the points.
/// @throws pcms_error if the quadrature's geometry is not in the geometry's
///         parameter system
[[nodiscard]] Quadrature MapQuadrature(
  const std::shared_ptr<const Geometry>& geometry,
  const Quadrature& quadrature);

} // namespace pcms

#include "pcms/geometry/detail/quadrature.hpp"

#endif // PCMS_GEOMETRY_QUADRATURE_HPP
