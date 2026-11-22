#ifndef PCMS_GEOMETRY_DETAIL_GEOMETRY_HPP
#define PCMS_GEOMETRY_DETAIL_GEOMETRY_HPP

#include "pcms/geometry/geometry.hpp"
#include <string>
#include <vector>

namespace pcms::detail
{

/// @brief The checks every Neighbourhood makes of its arguments.
/// @param geometry geometry whose Neighbourhood is called
/// @param coords the query points
/// @param options the search options
/// @param caller name of the calling function, for the message
/// @throws pcms_error if `coords` is not in the geometry's coordinate system,
///         that system does not have the identity metric, or
///         `options.max_supports` is neither 0 nor 3 * `options.min_supports`
void CheckNeighbourhoodArguments(
  const Geometry& geometry, const CoordinateView<DeviceMemorySpace>& coords,
  const NeighbourhoodOptions& options, const std::string& caller);

/// @brief Checks that every id of a hint is one of `num_entities` entities.
/// @param ids entity ids, e.g. OwnerHints::entities
/// @param num_entities number of entities the ids index
/// @param caller name of the checking function, for the message
/// @param entity what the ids index, e.g. "element", for the message
/// @throws pcms_error naming the range of the ids if one is outside
///         [0, num_entities)
void CheckEntityIds(Kokkos::View<const LO*, DeviceMemorySpace> ids,
                    LO num_entities, const std::string& caller,
                    const std::string& entity);

/// @brief Multiplies quadrature weights by the volume factor at their bound
///        points, making parameter-coordinate measures physical, or divides
///        them by it.
/// @param weights one weight per point of `bound`; scaled in place
/// @param bound the bound points of the weights
void ScaleByVolumeFactor(Kokkos::View<Real*, DeviceMemorySpace> weights,
                         const BoundPoints& bound, bool divide = false);

/// @brief The quadrature of points whose weights are in the parameter measure.
/// @param weights one weight per point; scaled in place by the volume factor
///        when the geometry's coordinate system has a metric
/// @param owners the owners of the points, as for Geometry::Bind
[[nodiscard]] Quadrature BindQuadrature(
  const Geometry& geometry, const CoordinateView<DeviceMemorySpace>& points,
  Kokkos::View<Real*, DeviceMemorySpace> weights,
  const OwnerHints& owners = {});

/// @brief The Gauss-Legendre rule on [-1, 1] exact for polynomials of degree
///        `order`: order / 2 + 1 nodes.
/// @param nodes the nodes, ascending
/// @param weights the weights of the nodes
void GaussLegendre(int order, std::vector<Real>& nodes,
                   std::vector<Real>& weights);

} // namespace pcms::detail

#endif // PCMS_GEOMETRY_DETAIL_GEOMETRY_HPP
