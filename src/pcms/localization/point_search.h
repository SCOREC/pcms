#ifndef PCMS_COUPLING_POINT_SEARCH_H
#define PCMS_COUPLING_POINT_SEARCH_H
#include <cassert>

#include <variant>

#include <Kokkos_Core.hpp>
#include <Omega_h_array.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_bbox.hpp>
#include <Omega_h_shape.hpp>

#include <iostream>

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

/**
 * Functor for the element -> grid cell map, the transpose of what the search
 * needs. Row `elem` lists the cells that element overlaps.
 *
 * This is the element-major inverse of the old cell-major functor: rather than
 * testing every element against a cell, each element only visits the cells its
 * own bounding box spans. That makes construction O(nelems * cells_per_elem)
 * instead of O(num_cells * nelems).
 */
template <int Dim, typename Pred>
struct ElemCellOverlapFunctor
{
  ElemCellOverlapFunctor(Omega_h::Mesh& mesh,
                         Kokkos::View<UniformGrid<Dim>[1]> grid)
    : elems2verts_(mesh.ask_elem_verts()), coords_(mesh.coords()), grid_(grid)
  {
    if (mesh.dim() != Dim) {
      std::cerr << "ElemCellOverlapFunctor: mesh dimension does not match the "
                   "functor dimension\n";
      std::terminate();
    }
  }

  /// Two-pass functor. On the first pass we count the cells this element
  /// overlaps; on the second we fill the CSR row with their IDs.
  KOKKOS_INLINE_FUNCTION
  LO operator()(LO elem, LO* fill) const
  {
    const auto elem_verts = Omega_h::gather_verts<Dim + 1>(elems2verts_, elem);
    const auto vertex_coords =
      Omega_h::gather_vectors<Dim + 1, Dim>(coords_, elem_verts);
    const auto& grid = grid_(0);

    const auto range = candidate_cells<Dim, Pred>(vertex_coords, grid);

    Kokkos::Array<LO, Dim> idx = range.lo;
    LO num_overlaps = 0;
    while (true) {
      const LO cell = grid.CellIndexFromAxisIndices(idx);
      if (Overlap<Dim, Pred>::apply(vertex_coords, grid.GetCellBBOX(cell))) {
        if (fill) {
          fill[num_overlaps] = cell;
        }
        ++num_overlaps;
      }
      // odometer over the inclusive range [lo, hi]
      int d = 0;
      for (; d < Dim; ++d) {
        if (++idx[d] <= range.hi[d]) {
          break;
        }
        idx[d] = range.lo[d];
      }
      if (d == Dim) {
        break;
      }
    }
    return num_overlaps;
  }

private:
  Omega_h::LOs elems2verts_;
  Omega_h::Reals coords_;
  Kokkos::View<UniformGrid<Dim>[1]> grid_;
};

/**
 * Build the grid cell -> element candidate map.
 *
 * Constructed as the transpose of the cheap element -> cell map. Kokkos'
 * transpose_crs cannot be used here: it sizes the output row count from
 * in.numRows() (the element count), which both over-allocates the row map and
 * writes out of bounds once num_grid_cells exceeds the element count.
 */
template <int Dim, typename Pred>
Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>
construct_intersection_map(Omega_h::Mesh& mesh,
                           Kokkos::View<UniformGrid<Dim>[1]> grid,
                           int num_grid_cells)
{
  using CrsT = Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>;
  const auto nelems = mesh.nelems();

  CrsT elem2cells{};
  Kokkos::count_and_fill_crs(elem2cells, nelems,
                             ElemCellOverlapFunctor<Dim, Pred>{mesh, grid});

  Kokkos::View<LO*> counts("candidate_counts", num_grid_cells);
  Kokkos::parallel_for(
    "count_cell_candidates", elem2cells.entries.size(), KOKKOS_LAMBDA(int i) {
      Kokkos::atomic_inc(&counts(elem2cells.entries(i)));
    });

  CrsT cell2elems{};
  Kokkos::get_crs_row_map_from_counts(cell2elems.row_map, counts,
                                      "candidate_row_map");
  cell2elems.entries =
    CrsT::entries_type("candidate_entries", elem2cells.entries.size());

  Kokkos::View<LO*> filled("candidate_fill_offsets", num_grid_cells);
  Kokkos::parallel_for(
    "fill_cell_candidates", nelems, KOKKOS_LAMBDA(LO elem) {
      for (auto j = elem2cells.row_map(elem); j < elem2cells.row_map(elem + 1);
           ++j) {
        const LO cell = elem2cells.entries(j);
        const auto slot = Kokkos::atomic_fetch_add(&filled(cell), 1);
        cell2elems.entries(cell2elems.row_map(cell) + slot) = elem;
      }
    });

  // The atomic scatter leaves each row in arbitrary order, where the old
  // cell-major loop emitted ascending element IDs. Restore that ordering: the
  // query takes the first candidate that contains the point, so ordering
  // decides which element wins when a point sits on a shared face, and the
  // out-of-bounds path reports the row's first candidate. Without this the
  // search results would vary run to run. Rows are short, so insertion sort.
  Kokkos::parallel_for(
    "sort_cell_candidates", num_grid_cells, KOKKOS_LAMBDA(LO cell) {
      const auto begin = cell2elems.row_map(cell);
      const auto end = cell2elems.row_map(cell + 1);
      for (auto i = begin + 1; i < end; ++i) {
        const LO value = cell2elems.entries(i);
        auto j = i;
        for (; j > begin && cell2elems.entries(j - 1) > value; --j) {
          cell2elems.entries(j) = cell2elems.entries(j - 1);
        }
        cell2elems.entries(j) = value;
      }
    });
  Kokkos::fence();
  return cell2elems;
}

/**
 * Cell-major brute-force build: for every grid cell, test every element.
 *
 * This is the implementation `construct_intersection_map` replaced. It is
 * O(num_cells * nelems) and must not be used in production, but it is the
 * ground truth for what the candidate map should contain, so it is kept and
 * exposed as an oracle for the equivalence test.
 */
template <int Dim, typename Pred>
struct ReferenceCellMajorFunctor
{
  ReferenceCellMajorFunctor(Omega_h::Mesh& mesh,
                            Kokkos::View<UniformGrid<Dim>[1]> grid)
    : elems2verts_(mesh.ask_elem_verts()),
      coords_(mesh.coords()),
      grid_(grid),
      nelems_(mesh.nelems())
  {
  }

  KOKKOS_INLINE_FUNCTION
  LO operator()(LO row, LO* fill) const
  {
    const auto cell_bbox = grid_(0).GetCellBBOX(row);
    LO num_intersections = 0;
    for (LO elem = 0; elem < nelems_; ++elem) {
      const auto elem_verts =
        Omega_h::gather_verts<Dim + 1>(elems2verts_, elem);
      const auto vertex_coords =
        Omega_h::gather_vectors<Dim + 1, Dim>(coords_, elem_verts);
      if (Overlap<Dim, Pred>::apply(vertex_coords, cell_bbox)) {
        if (fill) {
          fill[num_intersections] = elem;
        }
        ++num_intersections;
      }
    }
    return num_intersections;
  }

private:
  Omega_h::LOs elems2verts_;
  Omega_h::Reals coords_;
  Kokkos::View<UniformGrid<Dim>[1]> grid_;
  LO nelems_;
};

template <int Dim, typename Pred>
Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>
construct_intersection_map_reference(Omega_h::Mesh& mesh,
                                     Kokkos::View<UniformGrid<Dim>[1]> grid,
                                     int num_grid_cells)
{
  Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO> intersection_map{};
  Kokkos::count_and_fill_crs(intersection_map, num_grid_cells,
                             ReferenceCellMajorFunctor<Dim, Pred>{mesh, grid});
  return intersection_map;
}

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
