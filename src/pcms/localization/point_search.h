#ifndef PCMS_COUPLING_POINT_SEARCH_H
#define PCMS_COUPLING_POINT_SEARCH_H
#include <cassert>

#include <variant>

#include <Kokkos_Core.hpp>
#include <Omega_h_array.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_bbox.hpp>
#include <Omega_h_shape.hpp>

#include "pcms/localization/overlap.h"
#include "pcms/utility/types.h"
#include "pcms/utility/uniform_grid.h"
#include "pcms/utility/bounding_box.h"

namespace pcms
{

/**
 * Get the owning element (face for 2D, region for 3D) for a given entity.
 *
 * Given an entity of dimension `entity_dim` with ID `element_id`, returns the
 * smallest ID of the owning element of dimension `mesh_dim` (the mesh's
 * top-level dimension).
 */
LO GetOwningElementId(Omega_h::Mesh& mesh, int mesh_dim, int entity_dim,
                      LO element_id);

template <typename DimEnum>
KOKKOS_INLINE_FUNCTION LO
GetOwningElementIdFromAdj(const Omega_h::Adj& upward_adj, DimEnum entity_dim,
                          DimEnum target_dim, LO element_id)
{
  if (element_id < 0)
    return -1;

  // If entity is already at target dimension, return it directly
  if (static_cast<int>(entity_dim) == static_cast<int>(target_dim))
    return element_id;

  const auto begin = upward_adj.a2ab[element_id];
  const auto end = upward_adj.a2ab[element_id + 1];
  if (begin >= end)
    return -1;

  // Find the smallest owning element ID
  LO owner = upward_adj.ab2b[begin];
  for (auto i = begin + 1; i < end; ++i) {
    const LO candidate = upward_adj.ab2b[i];
    if (candidate < owner)
      owner = candidate;
  }
  return owner;
}

//
// TODO take a bounding box as we may want a bbox that's bigger than the mesh!
// this function is in the public header for testing, but should not be directly
// used
namespace detail
{
Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>
construct_intersection_map_2d(Omega_h::Mesh& mesh,
                              Kokkos::View<Uniform2DGrid[1]> grid,
                              int num_grid_cells);

Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>
construct_intersection_map_3d(Omega_h::Mesh& mesh,
                              Kokkos::View<Uniform3DGrid[1]> grid,
                              int num_grid_cells);

/// O(num_cells * nelems) brute-force candidate map. Test oracle for the
/// element-major build above; never use in production.
Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>
construct_intersection_map_reference_2d(Omega_h::Mesh& mesh,
                                        Kokkos::View<Uniform2DGrid[1]> grid,
                                        int num_grid_cells);

/// \see construct_intersection_map_reference_2d
Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>
construct_intersection_map_reference_3d(Omega_h::Mesh& mesh,
                                        Kokkos::View<Uniform3DGrid[1]> grid,
                                        int num_grid_cells);

/// How the candidate map learns an element's geometry, and a conservative bound
/// on it. The default policy for an Omega_h simplex mesh; defined in the source
/// file alongside the rest of the map machinery. A different element shape means
/// a different policy with the same three members, which is the only thing that
/// has to change to bin that shape against the same grid.
template <int Dim>
struct MeshSimplexElements;

/**
 * Build the grid cell -> element candidate map with `Pred`: the element-major
 * build, and the point search's hot path at construction time.
 */
template <int Dim, typename Pred, typename Elements = MeshSimplexElements<Dim>>
Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>
construct_intersection_map(Omega_h::Mesh& mesh,
                           Kokkos::View<UniformGrid<Dim>[1]> grid,
                           int num_grid_cells);

/// O(num_cells * nelems) brute-force counterpart, for tests. \see
/// construct_intersection_map
template <int Dim, typename Pred, typename Elements = MeshSimplexElements<Dim>>
Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>
construct_intersection_map_reference(Omega_h::Mesh& mesh,
                                     Kokkos::View<UniformGrid<Dim>[1]> grid,
                                     int num_grid_cells);

} // namespace detail

[[nodiscard]] KOKKOS_FUNCTION bool triangle_intersects_bbox(
  const Omega_h::Matrix<2, 3>& coords, const AABBox<2>& bbox);

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
    LO element_id;
    Omega_h::Vector<dim + 1> parametric_coords;

    /// `element_id` of a point outside the mesh whose nearest entity is `id`.
    [[nodiscard]] KOKKOS_INLINE_FUNCTION static constexpr LO OutsideId(LO id)
    {
      return -id - 1;
    }

    /// True when the point was not located inside any entity.
    [[nodiscard]] KOKKOS_INLINE_FUNCTION constexpr bool IsOutside() const
    {
      return element_id < 0;
    }

    /// Id of the entity of dimension `dimensionality`: the containing entity
    /// of a located point, the nearest entity of an outside point, or -1 for
    /// NO_INTERSECT.
    [[nodiscard]] KOKKOS_INLINE_FUNCTION constexpr LO EntityId() const
    {
      if (dimensionality == Dimensionality::NO_INTERSECT) {
        return -1;
      }
      return element_id < 0 ? -element_id - 1 : element_id;
    }
  };

  static constexpr auto DIM = dim;

  using PointSearchTolerances = Kokkos::View<Real[DIM]>;

  explicit PointLocalizationSearch(const PointSearchTolerances& tolerances)
    : tolerances_(tolerances)
  {
    assert(tolerances_.is_allocated());
  }

  virtual Kokkos::View<Result*> operator()(
    Kokkos::View<const Real* [dim]> point) const = 0;

  /**
   * This function provides a temporary solution to retrieve the original
   * behavior of the point search, which previously returned the face ID
   * regardless of which entity the search result belonged to. Many parts of
   * the codebase still rely on this legacy behavior. After updating the point
   * search to return the exact entity, this function can be called to retrieve
   * the old behavior and ensure correctness. Long term, the implementation
   * should be updated to properly handle the new result.
   */
  [[nodiscard]] virtual LO GetOwningElementId(const Result& result) = 0;
  [[nodiscard]] virtual Kokkos::View<LO*> GetOwningElementIds(
    Kokkos::View<const Result*> results) const = 0;
  virtual ~PointLocalizationSearch() = default;

protected:
  PointSearchTolerances tolerances_;
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
  GridPointSearch2D(Omega_h::Mesh& mesh, LO Nx, LO Ny,
                    const PointSearchTolerances& tolerances);

  /**
   * Build the candidate map with a predicate of the caller's choosing:
   *
   *     GridPointSearch2D search{mesh, nx, ny, tolerances, MyPredicate{}};
   *
   * The predicate decides which elements are offered as candidates for a grid
   * cell. The query itself is predicate-independent -- it tests containment of
   * the point in each candidate -- so a narrower predicate costs a pass over
   * more candidates but cannot change which answer is reported.
   */
  template <typename Pred>
  GridPointSearch2D(Omega_h::Mesh& mesh, LO Nx, LO Ny,
                    const PointSearchTolerances& tolerances, Pred);

  /**
   * Given points in global coordinates, returns for each the entity (vertex,
   * edge or face, chosen with the vertex and edge tolerances) containing it
   * and the point's barycentric coordinates in a face adjacent to that entity.
   */
  Kokkos::View<Result*> operator()(
    Kokkos::View<const Real* [DIM]> point) const override;
  [[nodiscard]] LO GetOwningElementId(const Result& result) override;
  [[nodiscard]] Kokkos::View<LO*> GetOwningElementIds(
    Kokkos::View<const Result*> results) const override;

private:
  /// Grid and adjacency: everything the map needs except the predicate.
  void InitializeGeometry(Omega_h::Mesh& mesh, LO Nx, LO Ny);

  /// The candidate map, which is the only predicate-dependent part.
  template <typename Pred>
  void InitializeMap();

  Omega_h::Mesh mesh_;
  Omega_h::Adj tris2edges_adj_;
  Omega_h::Adj tris2verts_adj_;
  Omega_h::Adj edges2verts_adj_;
  Omega_h::Adj edges2faces_up_;
  Omega_h::Adj verts2faces_up_;
  Kokkos::View<Uniform2DGrid[1]> grid_{"uniform grid"};
  LO num_grid_cells_ = 0;
  CandidateMapT candidate_map_;
  Omega_h::LOs tris2verts_;
  Omega_h::Reals coords_;
};

class GridPointSearch3D : public PointLocalizationSearch3D
{
  using CandidateMapT =
    Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>;

public:
  using Result = PointLocalizationSearch3D::Result;

  GridPointSearch3D(Omega_h::Mesh& mesh, LO Nx, LO Ny, LO Nz);
  GridPointSearch3D(Omega_h::Mesh& mesh, LO Nx, LO Ny, LO Nz,
                    const PointSearchTolerances& tolerances);

  /// \see GridPointSearch2D's predicate-taking constructor.
  template <typename Pred>
  GridPointSearch3D(Omega_h::Mesh& mesh, LO Nx, LO Ny, LO Nz,
                    const PointSearchTolerances& tolerances, Pred);

  /**
   * Given points in global coordinates, returns for each the id of a
   * tetrahedron containing it and the point's barycentric coordinates in that
   * tetrahedron. Containment is tested with a Cartesian tolerance band: the
   * point may lie up to `tolerances(2)` outside any face plane (the band is
   * applied per barycentric coordinate, scaled by the gradient norm of that
   * coordinate, so it is the same distance from both sides of a shared face).
   * The result is always labelled REGION. When several tetrahedra contain the
   * point (it lies on a shared face, edge or vertex) the smallest id is
   * returned; this is an interim tie-break that hides a genuine ambiguity,
   * see docs/point_search_shared_face_plan.md.
   *
   * If no candidate contains the point, the result is the candidate whose
   * face is closest, encoded with Result::OutsideId or NO_INTERSECT with
   * id -1 when the grid cell has no candidates at all.
   */
  Kokkos::View<Result*> operator()(
    Kokkos::View<const Real* [DIM]> point) const override;
  [[nodiscard]] LO GetOwningElementId(const Result& result) override;
  [[nodiscard]] Kokkos::View<LO*> GetOwningElementIds(
    Kokkos::View<const Result*> results) const override;

private:
  /// \see GridPointSearch2D::InitializeGeometry
  void InitializeGeometry(Omega_h::Mesh& mesh, LO Nx, LO Ny, LO Nz);

  /// \see GridPointSearch2D::InitializeMap
  template <typename Pred>
  void InitializeMap();

  Omega_h::Mesh mesh_;
  Omega_h::Adj tris2edges_adj_;
  Omega_h::Adj tris2verts_adj_;
  Omega_h::Adj edges2verts_adj_;
  Omega_h::Adj verts2regions_up_;
  Omega_h::Adj edges2regions_up_;
  Omega_h::Adj faces2regions_up_;
  Kokkos::View<UniformGrid<DIM>[1]> grid_{"uniform grid"};
  LO num_grid_cells_ = 0;
  CandidateMapT candidate_map_;
  Omega_h::LOs tris2verts_;
  Omega_h::Reals coords_;
};

template <typename Pred>
GridPointSearch2D::GridPointSearch2D(Omega_h::Mesh& mesh, LO Nx, LO Ny,
                                     const PointSearchTolerances& tolerances,
                                     Pred)
  : PointLocalizationSearch2D(tolerances), mesh_(mesh)
{
  InitializeGeometry(mesh, Nx, Ny);
  InitializeMap<Pred>();
}

template <typename Pred>
void GridPointSearch2D::InitializeMap()
{
  candidate_map_ =
    detail::construct_intersection_map<2, Pred>(mesh_, grid_, num_grid_cells_);
}

template <typename Pred>
GridPointSearch3D::GridPointSearch3D(Omega_h::Mesh& mesh, LO Nx, LO Ny, LO Nz,
                                     const PointSearchTolerances& tolerances,
                                     Pred)
  : PointLocalizationSearch3D(tolerances), mesh_(mesh)
{
  InitializeGeometry(mesh, Nx, Ny, Nz);
  InitializeMap<Pred>();
}

template <typename Pred>
void GridPointSearch3D::InitializeMap()
{
  candidate_map_ =
    detail::construct_intersection_map<3, Pred>(mesh_, grid_, num_grid_cells_);
}

/// A grid search of either spatial dimension, as owned by components that are
/// dimension-agnostic at compile time (e.g. an evaluator factory over an
/// Omega_h mesh) and handed to consumers that need to locate points on the same
/// mesh without building a second search.
using GridPointSearchVariant =
  std::variant<GridPointSearch2D, GridPointSearch3D>;

} // namespace pcms
#endif // PCMS_COUPLING_POINT_SEARCH_H
