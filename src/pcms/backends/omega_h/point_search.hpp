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

#include "pcms/utility/memory_spaces.hpp"
#include "pcms/utility/types.hpp"
#include "pcms/utility/uniform_grid.hpp"
#include "pcms/utility/bounding_box.hpp"
#include "pcms/utility/overlap.hpp"
#include "pcms/utility/shapes.hpp"

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

// The simplex of element `elem` of an Omega_h mesh.
template <int Dim>
[[nodiscard]] KOKKOS_INLINE_FUNCTION Simplex<Dim> GetElementSimplex(
  const Omega_h::Reals& coords, const Omega_h::LOs& elems2verts, LO elem)
{
  const auto verts = Omega_h::gather_verts<Dim + 1>(elems2verts, elem);
  return ToSimplex<Dim>(Omega_h::gather_vectors<Dim + 1, Dim>(coords, verts));
}
} // namespace detail

[[nodiscard]] KOKKOS_FUNCTION bool triangle_intersects_bbox(
  const Omega_h::Matrix<2, 3>& coords, const AABBox<2>& bbox);

/// The search reports the entity the point lies on and the point's parametric
/// coordinates on that entity. It does not choose an element; a consumer that
/// needs one lifts the point into an adjacent element (see OmegaHMeshBinding).
template <int dim>
class PointLocalizationSearch
{
public:
  struct Result
  {
    enum class Dimensionality
    {
      VERTEX = 0,
      EDGE = 1,
      FACE = 2,
      REGION = 3,
      NO_INTERSECT = 4
    };

    Dimensionality dimensionality;
    /// id of the entity of that dimension that contains the point;
    /// OutsideId(nearest entity) for a point outside the mesh; -1 for
    /// NO_INTERSECT
    LO entity;
    /// The point on the entity in reduced coordinates: the entity's
    /// dimension many entries are used (none for a vertex, 1 on an edge, 2 on
    /// a face, 3 in a region), the others are 0. With lambda_0..lambda_d the
    /// barycentric coordinates of the point in the entity, in the order of
    /// the entity's vertices (an edge's from the edge-to-vertex adjacency, a
    /// face's or region's from the element-to-vertex adjacency), the entries
    /// are (lambda_1, ..., lambda_d): the point is
    /// x_0 + sum_j parametric[j] * (x_{j+1} - x_0). For a point outside the
    /// mesh it is approximate: the barycentric coordinates of the nearest
    /// candidate element on the nearest entity, negative ones clamped to 0 and
    /// renormalized, not a projection onto the entity.
    Omega_h::Vector<dim> parametric;

    /// `entity` of a point outside the mesh whose nearest entity is `id`.
    [[nodiscard]] KOKKOS_INLINE_FUNCTION static constexpr LO OutsideId(LO id)
    {
      return -id - 1;
    }

    /// True when the point was not located inside any entity.
    [[nodiscard]] KOKKOS_INLINE_FUNCTION constexpr bool IsOutside() const
    {
      return entity < 0;
    }

    /// Id of the entity of dimension `dimensionality`: the containing entity
    /// of a located point, the nearest entity of an outside point, or -1 for
    /// NO_INTERSECT.
    [[nodiscard]] KOKKOS_INLINE_FUNCTION constexpr LO EntityId() const
    {
      if (dimensionality == Dimensionality::NO_INTERSECT) {
        return -1;
      }
      return entity < 0 ? -entity - 1 : entity;
    }
  };

  static constexpr auto DIM = dim;

  /// \param tolerance Cartesian distance from the facets of a simplex within
  ///        which a point counts as lying on them
  explicit PointLocalizationSearch(Real tolerance) : tolerance_(tolerance) {}

  virtual Kokkos::View<Result*> operator()(
    Kokkos::View<const Real* [dim]> point) const = 0;

  virtual ~PointLocalizationSearch() = default;

protected:
  Real tolerance_;
};

using PointLocalizationSearch2D = PointLocalizationSearch<2>;
using PointLocalizationSearch3D = PointLocalizationSearch<3>;

class GridPointSearch2D : public PointLocalizationSearch2D
{
  using CandidateMapT =
    Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>;

public:
  using Result = PointLocalizationSearch2D::Result;

  GridPointSearch2D(Omega_h::Mesh& mesh, LO Nx, LO Ny);
  GridPointSearch2D(Omega_h::Mesh& mesh, LO Nx, LO Ny, Real tolerance);

  /**
   * Given points in global coordinates, returns for each the entity (vertex,
   * edge, face or region) it lies on and the point's position on that entity,
   * see PointLocalizationSearch::Result. With lambda the barycentric
   * coordinates of the point in a candidate simplex and ||grad lambda_i|| their
   * gradient norms, the point is within Cartesian distance `tolerance` of the
   * plane of facet i iff |lambda_i| <= tolerance ||grad lambda_i||. The facets
   * whose band contains the point decide the entity: the sub-simplex spanned
   * by the vertices of the other facets (a vertex if the point is in the bands
   * of all facets but one, an edge if of all but two, ..., the simplex itself
   * if of none). The entity of lowest dimension over all candidates wins, then
   * the lowest id. Its reduced coordinates are the barycentric coordinates of
   * those vertices renormalized to sum 1, in the entity's vertex order.
   *
   * The tolerance region of an entity is therefore the intersection of facet
   * bands, not a Euclidean ball around it or a segment projection; near a
   * corner of angle alpha it reaches about tolerance / sin(alpha / 2) from the
   * vertex along the bisector. A point beyond the bands of every candidate is
   * outside: the result names, with Result::OutsideId, the entity the nearest
   * candidate (largest signed plane distance) puts it on, treating every facet
   * it is beyond as one it lies on; this is an approximate nearest entity. A
   * grid cell without candidates gives NO_INTERSECT with id -1.
   */
  Kokkos::View<Result*> operator()(
    Kokkos::View<const Real* [DIM]> point) const override;

  /**
   * As operator()(points), into `results` (one per point), starting from a
   * guess per point: the element the point is believed to lie in, or -1 for
   * none. A guess is verified, never trusted. A point whose guessed element,
   * or else an element sharing a facet with it (the mesh's dual graph),
   * contains it within `tolerance` of its facets (the closure band above) is
   * classified in that element alone: its entity and reduced coordinates are
   * the ones the entity-of-simplex rule gives there, without a grid search.
   * Every other point (no guess, a wrong one, a point outside the mesh) is
   * searched for on the grid exactly as operator()(points) does, compacted so
   * that only those points pay for it. `guesses` empty means no guesses.
   *
   * The result equals operator()(points) in the dimensionality, the id and the
   * reduced coordinates (up to rounding) wherever the entity of lowest
   * dimension within the tolerance of the point is a sub-simplex of the
   * guessed element or its neighbour; it can differ only within about
   * tolerance / sin(alpha / 2) of a vertex or edge of sharp angle alpha, where
   * a simplex that does not share a facet with the guess has the point in the
   * band of a lower-dimensional entity.
   */
  void operator()(Kokkos::View<const Real* [DIM]> points,
                  Kokkos::View<const LO*> guesses,
                  Kokkos::View<Result*> results) const;

private:
  /// per entity dimension 1..DIM-1, the element-to-entity and
  /// entity-to-vertex adjacencies
  Omega_h::Few<Omega_h::Adj, DIM> down_;
  Omega_h::Few<Omega_h::Adj, DIM> down_verts_;
  Kokkos::View<UniformGrid2D[1]> grid_{"uniform grid"};
  CandidateMapT candidate_map_;
  Omega_h::LOs elem_verts_;
  Omega_h::Reals coords_;
  /// the elements across each element's facets
  Omega_h::Graph dual_;
};

class GridPointSearch3D : public PointLocalizationSearch3D
{
  using CandidateMapT =
    Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>;

public:
  using Result = PointLocalizationSearch3D::Result;

  GridPointSearch3D(Omega_h::Mesh& mesh, LO Nx, LO Ny, LO Nz);
  GridPointSearch3D(Omega_h::Mesh& mesh, LO Nx, LO Ny, LO Nz, Real tolerance);

  /**
   * Given points in global coordinates, returns for each the entity (vertex,
   * edge, face or region) it lies on and the point's position on that entity,
   * see PointLocalizationSearch::Result. With lambda the barycentric
   * coordinates of the point in a candidate simplex and ||grad lambda_i|| their
   * gradient norms, the point is within Cartesian distance `tolerance` of the
   * plane of facet i iff |lambda_i| <= tolerance ||grad lambda_i||. The facets
   * whose band contains the point decide the entity: the sub-simplex spanned
   * by the vertices of the other facets (a vertex if the point is in the bands
   * of all facets but one, an edge if of all but two, ..., the simplex itself
   * if of none). The entity of lowest dimension over all candidates wins, then
   * the lowest id. Its reduced coordinates are the barycentric coordinates of
   * those vertices renormalized to sum 1, in the entity's vertex order.
   *
   * The tolerance region of an entity is therefore the intersection of facet
   * bands, not a Euclidean ball around it or a segment projection; near a
   * corner of angle alpha it reaches about tolerance / sin(alpha / 2) from the
   * vertex along the bisector. A point beyond the bands of every candidate is
   * outside: the result names, with Result::OutsideId, the entity the nearest
   * candidate (largest signed plane distance) puts it on, treating every facet
   * it is beyond as one it lies on; this is an approximate nearest entity. A
   * grid cell without candidates gives NO_INTERSECT with id -1.
   */
  Kokkos::View<Result*> operator()(
    Kokkos::View<const Real* [DIM]> point) const override;

  /**
   * As operator()(points), into `results` (one per point), starting from a
   * guess per point: the element the point is believed to lie in, or -1 for
   * none. A guess is verified, never trusted. A point whose guessed element,
   * or else an element sharing a facet with it (the mesh's dual graph),
   * contains it within `tolerance` of its facets (the closure band above) is
   * classified in that element alone: its entity and reduced coordinates are
   * the ones the entity-of-simplex rule gives there, without a grid search.
   * Every other point (no guess, a wrong one, a point outside the mesh) is
   * searched for on the grid exactly as operator()(points) does, compacted so
   * that only those points pay for it. `guesses` empty means no guesses.
   *
   * The result equals operator()(points) in the dimensionality, the id and the
   * reduced coordinates (up to rounding) wherever the entity of lowest
   * dimension within the tolerance of the point is a sub-simplex of the
   * guessed element or its neighbour; it can differ only within about
   * tolerance / sin(alpha / 2) of a vertex or edge of sharp angle alpha, where
   * a simplex that does not share a facet with the guess has the point in the
   * band of a lower-dimensional entity.
   */
  void operator()(Kokkos::View<const Real* [DIM]> points,
                  Kokkos::View<const LO*> guesses,
                  Kokkos::View<Result*> results) const;

private:
  /// per entity dimension 1..DIM-1, the element-to-entity and
  /// entity-to-vertex adjacencies
  Omega_h::Few<Omega_h::Adj, DIM> down_;
  Omega_h::Few<Omega_h::Adj, DIM> down_verts_;
  Kokkos::View<UniformGrid<DIM>[1]> grid_{"uniform grid"};
  CandidateMapT candidate_map_;
  Omega_h::LOs elem_verts_;
  Omega_h::Reals coords_;
  /// the elements across each element's facets
  Omega_h::Graph dual_;
};

/// A grid search of either spatial dimension, as owned by components that are
/// dimension-agnostic at compile time (e.g. an evaluator factory over an
/// Omega_h mesh) and handed to consumers that need to locate points on the same
/// mesh without building a second search.
using GridPointSearchVariant =
  std::variant<GridPointSearch2D, GridPointSearch3D>;

} // namespace pcms::omega_h
#endif // PCMS_BACKENDS_OMEGA_H_POINT_SEARCH_HPP
