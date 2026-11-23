#ifndef PCMS_BACKENDS_OMEGA_H_POINT_SEARCH_HPP
#define PCMS_BACKENDS_OMEGA_H_POINT_SEARCH_HPP

#include <cassert>

#include <utility>
#include <variant>
#include <vector>

#include <Kokkos_Core.hpp>
#include <Omega_h_array.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_bbox.hpp>
#include <Omega_h_shape.hpp>

#include "pcms/localization/point_search.hpp"
#include "pcms/utility/memory_spaces.hpp"
#include "pcms/utility/types.hpp"
#include "pcms/utility/uniform_grid.hpp"
#include "pcms/utility/bounding_box.hpp"
#include "pcms/localization/overlap.hpp"
#include "pcms/localization/shapes.hpp"

namespace pcms::omega_h
{

namespace detail
{
// A simplex from the coordinates of its vertices, one vector per vertex.
template <int Dim>
[[nodiscard]] KOKKOS_INLINE_FUNCTION Simplex<Dim> ToSimplex(
  const Omega_h::Few<Omega_h::Vector<Dim>, Dim + 1>& vertex_coords)
{
  Simplex<Dim> simplex;
  for (int k = 0; k < Dim + 1; ++k) {
    for (int d = 0; d < Dim; ++d) {
      simplex.vertices[k][d] = vertex_coords[k][d];
    }
  }
  return simplex;
}

} // namespace detail

[[nodiscard]] KOKKOS_FUNCTION bool TriangleIntersectsBBox(
  const Omega_h::Matrix<2, 3>& coords, const AABBox<2>& bbox);

template <int Dim>
using PointLocalizationSearch = localization::PointLocalizationSearch<Dim>;
template <int Dim>
using GridPointSearch = localization::GridPointSearch<Dim>;
using localization::GridPointSearch2D;
using localization::GridPointSearch3D;
using localization::PointLocalizationSearch2D;
using localization::PointLocalizationSearch3D;

/// A grid search of either spatial dimension, as owned by components that are
/// dimension-agnostic at compile time (e.g. an evaluator factory over an
/// Omega_h mesh) and handed to consumers that need to locate points on the same
/// mesh without building a second search.
using GridPointSearchVariant =
  std::variant<GridPointSearch2D, GridPointSearch3D>;

/// The grid point search over the elements of `mesh`, of the mesh's dimension,
/// with `divisions` grid cells per axis over the bounding box of the whole mesh
/// (so the call is collective).
GridPointSearchVariant CreateGridPointSearch(Omega_h::Mesh& mesh,
                                             const std::vector<LO>& divisions,
                                             Real tolerance = 1E-12);

} // namespace pcms::omega_h
#endif // PCMS_BACKENDS_OMEGA_H_POINT_SEARCH_HPP
