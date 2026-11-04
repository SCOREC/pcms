#ifndef PCMS_BACKENDS_OMEGA_H_ADJ_SEARCH_HPP
#define PCMS_BACKENDS_OMEGA_H_ADJ_SEARCH_HPP

#include <pcms/geometry/supports.hpp>
#include <pcms/localization/mls_support_helpers.hpp>
#include <pcms/utility/memory_spaces.hpp>
#include <pcms/utility/mesh_geometry.hpp>
#include <pcms/utility/print.hpp>
#include <pcms/utility/types.hpp>

#include <Kokkos_Core.hpp>
#include <Kokkos_MathematicalFunctions.hpp>

namespace pcms::backend::omega_h
{

/// @brief Supports of query points among the sites of `mesh`: a
///        breadth-first search over the site adjacency from the sites of the
///        elements at each point's seed, keeping the sites within the cutoff.
/// @param site_dim 0 for the mesh vertices as sites, mesh.dim() for the
///        element centroids
/// @param target_coords query points, flattened [x0, y0, ..., x1, y1, ...]
/// @param seed_dim mesh.dim() if each point's seed is the element it lies
///        in, 0 if it is a vertex at the point (the search starts from the
///        elements around it)
/// @param seeds one seed per query point; a point with -1 gets no supports
/// @param cutoffDistance initial squared radius of every query point
SupportResults searchNeighbors(
  Omega_h::Mesh& mesh, int site_dim, const Omega_h::Reals& target_coords,
  int seed_dim, Kokkos::View<const LO*, DeviceMemorySpace> seeds,
  Omega_h::Real cutoffDistance, Omega_h::LO min_req_support = 12,
  Omega_h::LO max_allowed_support = 36, bool adapt_radius = true);
} // namespace pcms::backend::omega_h

#endif // PCMS_BACKENDS_OMEGA_H_ADJ_SEARCH_HPP
