#ifndef PCMS_BACKENDS_OMEGA_H_POINT_SEARCH_HPP
#define PCMS_BACKENDS_OMEGA_H_POINT_SEARCH_HPP

#include <cassert>

#include <utility>
#include <variant>

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

//
// TODO take a bounding box as we may want a bbox that's bigger than the mesh!
// this function is in the public header for testing, but should not be directly
// used
namespace detail
{
Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>
construct_intersection_map_2d(Omega_h::Mesh& mesh,
                              Kokkos::View<UniformGrid2D[1]> grid,
                              int num_grid_cells, Real inflation = 0);

Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>
construct_intersection_map_3d(Omega_h::Mesh& mesh,
                              Kokkos::View<UniformGrid3D[1]> grid,
                              int num_grid_cells, Real inflation = 0);

/// O(num_cells * nelems) brute-force candidate map. Test oracle for the
/// element-major build above; never use in production.
Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>
construct_intersection_map_reference_2d(Omega_h::Mesh& mesh,
                                        Kokkos::View<UniformGrid2D[1]> grid,
                                        int num_grid_cells);

/// \see construct_intersection_map_reference_2d
Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>
construct_intersection_map_reference_3d(Omega_h::Mesh& mesh,
                                        Kokkos::View<UniformGrid3D[1]> grid,
                                        int num_grid_cells);

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

using localization::detail::GetElementSimplex;
} // namespace detail

[[nodiscard]] KOKKOS_FUNCTION bool triangle_intersects_bbox(
  const Omega_h::Matrix<2, 3>& coords, const AABBox<2>& bbox);

template <int dim>
using PointLocalizationSearch = localization::PointLocalizationSearch<dim>;
using localization::GridPointSearch2D;
using localization::GridPointSearch3D;
using localization::PointLocalizationSearch2D;
using localization::PointLocalizationSearch3D;

/// The grid point search over the elements of `mesh`, with `Nx` x `Ny` grid
/// cells over the bounding box of the whole mesh (so the call is collective).
GridPointSearch2D MakeGridPointSearch2D(Omega_h::Mesh& mesh, LO Nx, LO Ny,
                                        Real tolerance = 1E-12);

/// \see MakeGridPointSearch2D
GridPointSearch3D MakeGridPointSearch3D(Omega_h::Mesh& mesh, LO Nx, LO Ny,
                                        LO Nz, Real tolerance = 1E-12);

/// A grid search of either spatial dimension, as owned by components that are
/// dimension-agnostic at compile time (e.g. an evaluator factory over an
/// Omega_h mesh) and handed to consumers that need to locate points on the same
/// mesh without building a second search.
using GridPointSearchVariant =
  std::variant<GridPointSearch2D, GridPointSearch3D>;

} // namespace pcms::omega_h
#endif // PCMS_BACKENDS_OMEGA_H_POINT_SEARCH_HPP
