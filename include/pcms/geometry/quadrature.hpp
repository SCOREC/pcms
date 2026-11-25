#ifndef PCMS_GEOMETRY_QUADRATURE_HPP
#define PCMS_GEOMETRY_QUADRATURE_HPP

#include "pcms/geometry/bound_points.hpp"
#include "pcms/utility/memory_spaces.hpp"
#include "pcms/utility/types.hpp"
#include <Kokkos_Core.hpp>
#include <cstdint>
#include <memory>

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
  /// cells of a second geometry the points lie in; empty if there is none
  OwnerHints other;
  /// the polynomial degree the rule integrates exactly; -1 if unknown or not
  /// polynomial-exact
  int degree = -1;
};

/// @brief The Gauss rule of each cell of `geometry` (on simplices, the
///        reference rule), exact for polynomials of degree `order`.
/// @throws pcms_error if the geometry has no such rule of that order
[[nodiscard]] Quadrature CreateGaussQuadrature(
  const std::shared_ptr<const Geometry>& geometry, int order);

/// @brief Uniform random samples in each cell of `geometry`, each weighing
///        cell measure / samples_per_cell; a cell's samples depend only on
///        `seed` and its global id.
/// @throws pcms_error if the geometry has no Monte Carlo rule or
///         `samples_per_cell` is not positive
[[nodiscard]] Quadrature CreateMonteCarloQuadrature(
  const std::shared_ptr<const Geometry>& geometry, int samples_per_cell,
  std::uint64_t seed = DefaultQuadratureSeed);

} // namespace pcms

#include "pcms/geometry/detail/quadrature.hpp"

#endif // PCMS_GEOMETRY_QUADRATURE_HPP
