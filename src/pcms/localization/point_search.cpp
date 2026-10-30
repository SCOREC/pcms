#include "point_search.h"
#include "pcms/localization/overlap.h"
#include <Omega_h_mesh.hpp>
#include <bitset>
#include <cmath>
#include <exception>
#include <iostream>

// From
// https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line#Vector_formulation
template <Omega_h::Int dim>
KOKKOS_INLINE_FUNCTION Omega_h::Real distance_from_line(
  const Omega_h::Vector<dim>& a, const Omega_h::Vector<dim>& b,
  const Omega_h::Vector<dim>& p)
{
  Omega_h::Vector<dim> n = Omega_h::normalize(b - a);
  Omega_h::Vector<dim> ap = a - p;
  return Omega_h::norm(ap - (ap * n) * n);
}

template <Omega_h::Int dim>
KOKKOS_INLINE_FUNCTION bool normal_intersects_segment(
  const Omega_h::Vector<dim> a, const Omega_h::Vector<dim> b,
  const Omega_h::Vector<dim> p)
{
  auto ab = b - a;
  auto ba = a - b;
  auto ap = p - a;
  auto bp = p - b;
  return (ap * ab) * (bp * ba) >= 0;
}

/**
 * Barycentric coordinates of a point in a simplex together with the factors
 * that convert a Cartesian distance from facet i into a band on coordinate i:
 * |lambda[i]| <= tau * grad_norm[i] iff the point is within distance tau of
 * the plane of the facet opposite vertex i.
 */
template <int Dim>
struct BarycentricWithScales
{
  Omega_h::Vector<Dim + 1> lambda;
  Omega_h::Vector<Dim + 1> grad_norm;

  /**
   * Signed Cartesian distance to the most violated facet plane: non-negative
   * inside the simplex, -tau on the boundary of the tolerance band.
   */
  KOKKOS_INLINE_FUNCTION Omega_h::Real SignedPlaneDistance() const
  {
    auto d = lambda[0] / grad_norm[0];
    for (int i = 1; i <= Dim; ++i) {
      d = Kokkos::fmin(d, lambda[i] / grad_norm[i]);
    }
    return d;
  }
};

template <int Dim>
KOKKOS_INLINE_FUNCTION BarycentricWithScales<Dim> barycentric_with_scales(
  const Omega_h::Few<Omega_h::Vector<Dim>, Dim + 1>& vertex_coords,
  const Omega_h::Vector<Dim>& point)
{
  // rows of the inverse basis are the gradients of lambda_1..lambda_Dim;
  // Omega_h matrices are stored by column, so row j is inv[k][j] over k
  const auto inv =
    Omega_h::pseudo_invert(Omega_h::simplex_basis<Dim, Dim>(vertex_coords));
  BarycentricWithScales<Dim> result;
  result.lambda = Omega_h::form_barycentric(inv * (point - vertex_coords[0]));
  auto grad0 = Omega_h::zero_vector<Dim>();
  for (int j = 0; j < Dim; ++j) {
    Omega_h::Vector<Dim> row;
    for (int k = 0; k < Dim; ++k) {
      row[k] = inv[k][j];
    }
    result.grad_norm[j + 1] = Omega_h::norm(row);
    grad0 = grad0 - row;
  }
  result.grad_norm[0] = Omega_h::norm(grad0);
  return result;
}

namespace pcms
{

LO GetOwningElementId(Omega_h::Mesh& mesh, int mesh_dim, int entity_dim,
                      LO element_id)
{
  if (element_id < 0)
    return -1;

  const int target_dim = mesh_dim; // faces for 2D, regions for 3D

  // If entity is already at target dimension, return it directly
  if (entity_dim == target_dim)
    return element_id;

  // Get the upward adjacency from entity_dim to target_dim
  auto upward_adj = mesh.ask_up(entity_dim, target_dim);
  auto a2ab_h = Omega_h::HostRead<LO>(upward_adj.a2ab);
  auto ab2b_h = Omega_h::HostRead<LO>(upward_adj.ab2b);

  const auto begin = a2ab_h[element_id];
  const auto end = a2ab_h[element_id + 1];
  if (begin >= end)
    return -1;

  // Find the smallest owning element ID
  LO owner = ab2b_h[begin];
  for (auto i = begin + 1; i < end; ++i) {
    const LO candidate = ab2b_h[i];
    if (candidate < owner)
      owner = candidate;
  }
  return owner;
}

/**
 * Check if a triangle element represented by 3 coordinates in two dimensions
 * intersects with a bounding box.
 */
[[nodiscard]] KOKKOS_FUNCTION bool triangle_intersects_bbox(
  const Omega_h::Matrix<2, 3>& coords, const AABBox<2>& bbox)
{
  return Overlap<2, TriangleVsBbox>::apply(coords, bbox);
}

namespace detail
{
// The candidate-map machinery. Declared in the header, defined here, and
// explicitly instantiated for the supported predicates at the end of this
// file -- see the note there.

/**
 * The bridge between "an element id" and "the geometry a predicate tests".
 *
 * This is the only shape-specific part of building a candidate map. Everything
 * downstream -- the broad-phase binning and the predicate test -- works on
 * whatever `operands()` returns, so supporting a different element shape is a
 * new policy with these three members and nothing else changes:
 *
 *     struct MyElements {
 *       static MyElements FromMesh(Omega_h::Mesh&);
 *       Operands operands(LO elem) const;              // what Pred::apply
 * takes AABBox<Dim> bound(const Operands&) const;      // a box containing it
 *     };
 *
 * `bound` is allowed to be conservative: the predicate rejects what it
 * over-includes, whereas a bound that is too tight silently drops candidates.
 */
template <int Dim>
struct MeshSimplexElements
{
  Omega_h::LOs elems2verts;
  Omega_h::Reals coords;

  [[nodiscard]] static MeshSimplexElements FromMesh(Omega_h::Mesh& mesh)
  {
    return {mesh.ask_elem_verts(), mesh.coords()};
  }

  [[nodiscard]] KOKKOS_INLINE_FUNCTION Omega_h::Matrix<Dim, Dim + 1> operands(
    LO elem) const
  {
    const auto verts = Omega_h::gather_verts<Dim + 1>(elems2verts, elem);
    return Omega_h::gather_vectors<Dim + 1, Dim>(coords, verts);
  }

  [[nodiscard]] KOKKOS_INLINE_FUNCTION AABBox<Dim> bound(
    const Omega_h::Matrix<Dim, Dim + 1>& ops) const
  {
    return simplex_bbox<Dim>(ops);
  }
};

/**
 * Functor for the element -> grid cell map, the transpose of what the search
 * needs. Row `elem` lists the cells that element overlaps.
 *
 * This is the element-major inverse of the old cell-major functor: rather than
 * testing every element against a cell, each element only visits the cells its
 * own bounding box spans. That makes construction O(nelems * cells_per_elem)
 * instead of O(num_cells * nelems).
 *
 * Operand-agnostic: `Elements` decides what an element *is* and how to bound
 * it, `Pred` decides when it overlaps a cell. Neither is fixed here, so the
 * same functor bins any shape a caller can supply a policy for.
 */
template <int Dim, typename Pred, typename Elements = MeshSimplexElements<Dim>>
struct ElemCellOverlapFunctor
{
  ElemCellOverlapFunctor(Elements elements,
                         Kokkos::View<UniformGrid<Dim>[1]> grid)
    : elements_(elements), grid_(grid)
  {
  }

  /// Two-pass functor. On the first pass we count the cells this element
  /// overlaps; on the second we fill the CSR row with their IDs.
  KOKKOS_INLINE_FUNCTION
  LO operator()(LO elem, LO* fill) const
  {
    const auto& grid = grid_(0);
    const auto ops = elements_.operands(elem);

    const auto range = candidate_cells<Dim, Pred>(elements_.bound(ops), grid);

    Kokkos::Array<LO, Dim> idx = range.lo;
    LO num_overlaps = 0;
    while (true) {
      const LO cell = grid.CellIndexFromAxisIndices(idx);
      if (Overlap<Dim, Pred>::apply(ops, grid.GetCellBBOX(cell))) {
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
  Elements elements_;
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
template <int Dim, typename Pred, typename Elements>
Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>
construct_intersection_map(Omega_h::Mesh& mesh,
                           Kokkos::View<UniformGrid<Dim>[1]> grid,
                           int num_grid_cells)
{
  if (mesh.dim() != Dim) {
    std::cerr
      << "construct_intersection_map: mesh dimension does not match the "
         "map dimension\n";
    std::terminate();
  }

  using CrsT = Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>;
  const auto nelems = mesh.nelems();

  CrsT elem2cells{};
  Kokkos::count_and_fill_crs(elem2cells, nelems,
                             ElemCellOverlapFunctor<Dim, Pred, Elements>{
                               Elements::FromMesh(mesh), grid});

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
template <int Dim, typename Pred, typename Elements = MeshSimplexElements<Dim>>
struct ReferenceCellMajorFunctor
{
  ReferenceCellMajorFunctor(Elements elements,
                            Kokkos::View<UniformGrid<Dim>[1]> grid, LO nelems)
    : elements_(elements), grid_(grid), nelems_(nelems)
  {
  }

  KOKKOS_INLINE_FUNCTION
  LO operator()(LO row, LO* fill) const
  {
    const auto cell_bbox = grid_(0).GetCellBBOX(row);
    LO num_intersections = 0;
    for (LO elem = 0; elem < nelems_; ++elem) {
      const auto ops = elements_.operands(elem);
      if (Overlap<Dim, Pred>::apply(ops, cell_bbox)) {
        if (fill) {
          fill[num_intersections] = elem;
        }
        ++num_intersections;
      }
    }
    return num_intersections;
  }

private:
  Elements elements_;
  Kokkos::View<UniformGrid<Dim>[1]> grid_;
  LO nelems_;
};

template <int Dim, typename Pred, typename Elements>
Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>
construct_intersection_map_reference(Omega_h::Mesh& mesh,
                                     Kokkos::View<UniformGrid<Dim>[1]> grid,
                                     int num_grid_cells)
{
  Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO> intersection_map{};
  Kokkos::count_and_fill_crs(intersection_map, num_grid_cells,
                             ReferenceCellMajorFunctor<Dim, Pred, Elements>{
                               Elements::FromMesh(mesh), grid, mesh.nelems()});
  return intersection_map;
}
Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>
construct_intersection_map_reference_2d(Omega_h::Mesh& mesh,
                                        Kokkos::View<Uniform2DGrid[1]> grid,
                                        int num_grid_cells)
{
  return construct_intersection_map_reference<2, TriangleVsBbox>(
    mesh, grid, num_grid_cells);
}

Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>
construct_intersection_map_reference_3d(Omega_h::Mesh& mesh,
                                        Kokkos::View<Uniform3DGrid[1]> grid,
                                        int num_grid_cells)
{
  return construct_intersection_map_reference<3, SimplexVsBbox>(mesh, grid,
                                                                num_grid_cells);
}

// num_grid_cells should be result of grid.GetNumCells(), take as argument to
// avoid extra copy of grid from gpu to cpu
//
// The 2D and 3D overloads name their predicate here rather than reading it from
// a shared mapping: a "2D refined candidate map" is what this function *is*, so
// the choice belongs to its definition and cannot leak into another call site.
// Callers needing a different predicate use the template directly.
Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>
construct_intersection_map_2d(Omega_h::Mesh& mesh,
                              Kokkos::View<Uniform2DGrid[1]> grid,
                              int num_grid_cells)
{
  return construct_intersection_map<2, TriangleVsBbox>(mesh, grid,
                                                       num_grid_cells);
}

Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>
construct_intersection_map_3d(Omega_h::Mesh& mesh,
                              Kokkos::View<Uniform3DGrid[1]> grid,
                              int num_grid_cells)
{
  return construct_intersection_map<3, SimplexVsBbox>(mesh, grid,
                                                      num_grid_cells);
}
// The supported predicates, instantiated here so that the search classes'
// predicate-taking constructors and the tests can use them from other
// translation units. This is the complete set: naming a predicate that is not
// listed below is a link error, which is the accepted cost of keeping the
// machinery out of the header. Add a line here when adding a predicate that
// callers outside this file need.
template Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>
construct_intersection_map<2, TriangleVsBbox>(Omega_h::Mesh&,
                                              Kokkos::View<UniformGrid<2>[1]>,
                                              int);
template Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>
construct_intersection_map<2, SimplexVsBbox>(Omega_h::Mesh&,
                                             Kokkos::View<UniformGrid<2>[1]>,
                                             int);
template Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>
construct_intersection_map<3, SimplexVsBbox>(Omega_h::Mesh&,
                                             Kokkos::View<UniformGrid<3>[1]>,
                                             int);

template Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>
construct_intersection_map_reference<2, TriangleVsBbox>(
  Omega_h::Mesh&, Kokkos::View<UniformGrid<2>[1]>, int);
template Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>
construct_intersection_map_reference<2, SimplexVsBbox>(
  Omega_h::Mesh&, Kokkos::View<UniformGrid<2>[1]>, int);
template Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>
construct_intersection_map_reference<3, SimplexVsBbox>(
  Omega_h::Mesh&, Kokkos::View<UniformGrid<3>[1]>, int);

} // namespace detail

template <int n, typename Op>
OMEGA_H_INLINE double myreduce(const Omega_h::Vector<n>& x,
                               Op op) OMEGA_H_NOEXCEPT
{
  auto out = x[0];
  for (int i = 1; i < n; ++i)
    out = op(out, x[i]);
  return out;
}

Kokkos::View<GridPointSearch2D::Result*> GridPointSearch2D::operator()(
  Kokkos::View<const Real* [DIM]> points) const
{
  Kokkos::View<GridPointSearch2D::Result*> results("point search result",
                                                   points.extent(0));
  auto num_rows = candidate_map_.numRows();
  // needed so that we don't capture this ptr which will be memory error on cuda
  auto grid = grid_;
  auto candidate_map = candidate_map_;
  auto tris2verts = tris2verts_;
  auto tris2verts_adj = tris2verts_adj_;
  auto tris2edges_adj = tris2edges_adj_;
  auto edges2verts_adj = edges2verts_adj_;
  auto coords = coords_;
  auto tolerances = tolerances_;
  Kokkos::parallel_for(
    points.extent(0), KOKKOS_LAMBDA(int p) {
      Omega_h::Vector<2> point(
        std::initializer_list<double>{points(p, 0), points(p, 1)});
      auto cell_id = grid(0).ClosestCellID(point);
      assert(cell_id < num_rows && cell_id >= 0);
      auto candidates_begin = candidate_map.row_map(cell_id);
      auto candidates_end = candidate_map.row_map(cell_id + 1);

      // Track best entities across all candidates to ensure order invariance
      Omega_h::Real best_vertex_dist = INFINITY;
      LO best_vertex_id = -1;
      int best_vertex_tid = -1; // triangle providing barycentric coords
      Omega_h::Vector<3> best_vertex_bary{0.0, 0.0, 0.0};

      Omega_h::Real best_edge_dist = INFINITY;
      LO best_edge_id = -1;
      int best_edge_tid_min = -1;
      int best_edge_tid_max = -1;
      Omega_h::Vector<3> best_edge_bary_min{0.0, 0.0, 0.0};
      Omega_h::Vector<3> best_edge_bary_max{0.0, 0.0, 0.0};

      bool found_inside = false;
      LO inside_face_id = -1;
      Omega_h::Vector<3> inside_face_bary{0.0, 0.0, 0.0};

      auto begin = candidate_map.row_map(cell_id);
      auto end = candidate_map.row_map(cell_id + 1);
      for (auto ii = begin; ii < end; ++ii) {
        const int triangleID = candidate_map.entries(ii);
        const auto elem_tri2verts =
          Omega_h::gather_verts<3>(tris2verts, triangleID);
        auto vertex_coords =
          Omega_h::gather_vectors<3, 2>(coords, elem_tri2verts);
        auto parametric_coords =
          Omega_h::barycentric_from_global<2, 2>(point, vertex_coords);

        // Check vertices (hierarchy level 1): compute Euclidean distance
        for (int j = 0; j < 3; ++j) {
          const int vertexID = tris2verts_adj.ab2b[triangleID * 3 + j];
          const auto v = Omega_h::get_vector<2>(coords, vertexID);
          const auto dv = Omega_h::norm(point - v);
          if ((dv < best_vertex_dist) ||
              ((dv == best_vertex_dist) && (vertexID < best_vertex_id))) {
            best_vertex_dist = dv;
            best_vertex_id = vertexID;
            best_vertex_tid = triangleID;
            best_vertex_bary = parametric_coords;
          }
        }

        // Check edges (hierarchy level 2): only if projection falls within

        for (int j = 0; j < 3; ++j) {
          const int edgeID = tris2edges_adj.ab2b[triangleID * 3 + j];

          const int va_id = edges2verts_adj.ab2b[edgeID * 2 + 0];
          const int vb_id = edges2verts_adj.ab2b[edgeID * 2 + 1];
          const auto va = Omega_h::get_vector<2>(coords, va_id);
          const auto vb = Omega_h::get_vector<2>(coords, vb_id);

          if (!normal_intersects_segment(va, vb, point))
            continue;

          const auto de = distance_from_line(va, vb, point);
          if ((de < best_edge_dist) ||
              ((de == best_edge_dist) && (edgeID < best_edge_id)) ||
              ((de == best_edge_dist) && (edgeID == best_edge_id) &&
               (triangleID < best_edge_tid_min))) {
            best_edge_dist = de;
            best_edge_id = edgeID;
            best_edge_tid_min = triangleID;
            best_edge_tid_max = triangleID;
            best_edge_bary_min = parametric_coords;
            best_edge_bary_max = parametric_coords;
          } else if ((de == best_edge_dist) && (edgeID == best_edge_id)) {
            // Track the extents of face IDs sharing this nearest edge
            if (triangleID < best_edge_tid_min) {
              best_edge_tid_min = triangleID;
              best_edge_bary_min = parametric_coords;
            }
            if (triangleID > best_edge_tid_max) {
              best_edge_tid_max = triangleID;
              best_edge_bary_max = parametric_coords;
            }
          }
        }

        // Check face interior (hierarchy level 3)
        if (Omega_h::is_barycentric_inside(parametric_coords)) {
          if (!found_inside || (triangleID < inside_face_id)) {
            found_inside = true;
            inside_face_id = triangleID;
            inside_face_bary = parametric_coords;
          }
        }
      }

      const auto vtol = tolerances(0);
      const auto etol = tolerances(1);

      // If we found an inside face, compute its edge distance using
      // barycentric-only helper to check for edge classification while
      // preserving the containing face ID.
      Real inside_edge_dist = INFINITY;
      int inside_edge_argmin = -1;
      Omega_h::Matrix<2, 3> vcoords_in;
      // Track Euclidean-nearest edge of the containing face (if any)
      LO inside_edge_id = -1;
      if (found_inside) {
        const auto elem_tri2verts_in =
          Omega_h::gather_verts<3>(tris2verts, inside_face_id);
        vcoords_in = Omega_h::gather_vectors<3, 2>(coords, elem_tri2verts_in);
        // Euclidean distance to the 3 edges of the containing triangle
        inside_edge_argmin = -1;
        inside_edge_dist = INFINITY;
        inside_edge_id = -1;
        for (int j = 0; j < 3; ++j) {
          const int edgeID = tris2edges_adj.ab2b[inside_face_id * 3 + j];
          const int va_id = edges2verts_adj.ab2b[edgeID * 2 + 0];
          const int vb_id = edges2verts_adj.ab2b[edgeID * 2 + 1];
          const auto va = Omega_h::get_vector<2>(coords, va_id);
          const auto vb = Omega_h::get_vector<2>(coords, vb_id);
          if (!normal_intersects_segment(va, vb, point))
            continue;
          const auto de = distance_from_line(va, vb, point);
          if (de < inside_edge_dist) {
            inside_edge_dist = de;
            inside_edge_argmin = j;
            inside_edge_id = edgeID;
          }
        }
      }

      GridPointSearch2D::Result::Dimensionality dim_out =
        GridPointSearch2D::Result::Dimensionality::NO_INTERSECT;
      LO element_id_out = -1;
      Omega_h::Vector<3> bary_out{0.0, 0.0, 0.0};

      // Apply hierarchy with tolerances
      // Points within tolerance are considered "inside" with positive IDs
      // Only points with no candidates at all get negative IDs
      if (best_vertex_id >= 0 && best_vertex_dist <= vtol) {
        // Point within vertex tolerance - still inside the mesh
        dim_out = GridPointSearch2D::Result::Dimensionality::VERTEX;
        element_id_out = best_vertex_id;
        bary_out = best_vertex_bary;
      } else if ((best_edge_id >= 0 && best_edge_dist <= etol) ||
                 (found_inside && inside_edge_dist <= etol)) {
        // Point within edge tolerance - still inside the mesh
        dim_out = GridPointSearch2D::Result::Dimensionality::EDGE;
        if (found_inside && inside_edge_dist <= etol &&
            inside_edge_argmin >= 0) {
          element_id_out = inside_edge_id;
          bary_out = inside_face_bary;
        } else {
          element_id_out = best_edge_id;
          bary_out = best_edge_bary_min;
        }
      } else if (found_inside) {
        // Point inside face
        dim_out = GridPointSearch2D::Result::Dimensionality::FACE;
        element_id_out = inside_face_id;
        bary_out = inside_face_bary;
      } else {
        // Outside mesh - beyond tolerance of any entity: name the nearest
        // entity, encoded so the id stays negative even for entity 0.
        if (best_vertex_id >= 0 &&
            (best_vertex_dist <= best_edge_dist || best_edge_id < 0)) {
          dim_out = GridPointSearch2D::Result::Dimensionality::VERTEX;
          element_id_out = GridPointSearch2D::Result::OutsideId(best_vertex_id);
          bary_out = best_vertex_bary;
        } else if (best_edge_id >= 0) {
          dim_out = GridPointSearch2D::Result::Dimensionality::EDGE;
          element_id_out = GridPointSearch2D::Result::OutsideId(best_edge_id);
          bary_out = best_edge_bary_min;
        } else {
          // No candidates at all: NO_INTERSECT with id -1
        }
      }

      results(p) = GridPointSearch2D::Result{dim_out, element_id_out, bary_out};
    });

  return results;
}

GridPointSearch2D::GridPointSearch2D(Omega_h::Mesh& mesh, LO Nx, LO Ny)
  : GridPointSearch2D(mesh, Nx, Ny,
                      PointSearchTolerances{"point search 2d tolerances"})
{
  Kokkos::deep_copy(tolerances_, 1E-12);
}

GridPointSearch2D::GridPointSearch2D(Omega_h::Mesh& mesh, LO Nx, LO Ny,
                                     const PointSearchTolerances& tolerances)
  : PointLocalizationSearch(tolerances), mesh_(mesh)
{
  InitializeGeometry(mesh, Nx, Ny);
  // This search's default predicate: the refined triangle-vs-box test.
  InitializeMap<TriangleVsBbox>();
}

void GridPointSearch2D::InitializeGeometry(Omega_h::Mesh& mesh, LO Nx, LO Ny)
{
  auto mesh_bbox = Omega_h::get_bounding_box<2>(&mesh);
  auto grid_h = Kokkos::create_mirror_view(grid_);
  grid_h(0) =
    Uniform2DGrid{.edge_length = {mesh_bbox.max[0] - mesh_bbox.min[0],
                                  mesh_bbox.max[1] - mesh_bbox.min[1]},
                  .bot_left = {mesh_bbox.min[0], mesh_bbox.min[1]},
                  .divisions = {Nx, Ny}};
  Kokkos::deep_copy(grid_, grid_h);
  num_grid_cells_ = grid_h(0).GetNumCells();
  coords_ = mesh.coords();
  tris2verts_ = mesh.ask_elem_verts();
  tris2edges_adj_ = mesh.ask_down(Omega_h::FACE, Omega_h::EDGE);
  tris2verts_adj_ = mesh.ask_down(Omega_h::FACE, Omega_h::VERT);
  edges2verts_adj_ = mesh.ask_down(Omega_h::EDGE, Omega_h::VERT);
  edges2faces_up_ = mesh.ask_up(Omega_h::EDGE, Omega_h::FACE);
  verts2faces_up_ = mesh.ask_up(Omega_h::VERT, Omega_h::FACE);
}

LO GridPointSearch2D::GetOwningElementId(const Result& result)
{
  return pcms::GetOwningElementId(
    mesh_, 2, static_cast<int>(result.dimensionality), result.EntityId());
}

Kokkos::View<LO*> GridPointSearch2D::GetOwningElementIds(
  Kokkos::View<const Result*> results) const
{
  Kokkos::View<LO*> owners("point search owning face ids", results.extent(0));
  auto edges2faces_up = edges2faces_up_;
  auto verts2faces_up = verts2faces_up_;
  constexpr int mesh_dim = 2;
  Kokkos::parallel_for(
    results.extent(0), KOKKOS_LAMBDA(const LO i) {
      const auto result = results(i);
      const LO element_id = result.EntityId();

      LO owner = -1;
      if (element_id >= 0) {
        if (result.dimensionality == Result::Dimensionality::FACE) {
          owner = GetOwningElementIdFromAdj(
            edges2faces_up, Result::Dimensionality::FACE,
            Result::Dimensionality::FACE, element_id);
        } else if (result.dimensionality == Result::Dimensionality::EDGE) {
          owner = GetOwningElementIdFromAdj(
            edges2faces_up, Result::Dimensionality::EDGE,
            Result::Dimensionality::FACE, element_id);
        } else if (result.dimensionality == Result::Dimensionality::VERTEX) {
          owner = GetOwningElementIdFromAdj(
            verts2faces_up, Result::Dimensionality::VERTEX,
            Result::Dimensionality::FACE, element_id);
        }
      }
      owners(i) = owner;
    });
  return owners;
}

Kokkos::View<GridPointSearch3D::Result*> GridPointSearch3D::operator()(
  Kokkos::View<const Real* [DIM]> points) const
{
  Kokkos::View<GridPointSearch3D::Result*> results("point search result",
                                                   points.extent(0));
  auto num_rows = candidate_map_.numRows();
  // needed so that we don't capture this ptr which will be memory error on cuda
  auto grid = grid_;
  auto candidate_map = candidate_map_;
  auto tris2verts = tris2verts_;
  auto coords = coords_;
  auto tolerances = tolerances_;
  Kokkos::parallel_for(
    points.extent(0), KOKKOS_LAMBDA(int p) {
      using Dimensionality = GridPointSearch3D::Result::Dimensionality;
      Omega_h::Vector<DIM> point;
      for (int i = 0; i < DIM; ++i) {
        point[i] = points(p, i);
      }

      auto cell_id = grid(0).ClosestCellID(point);
      assert(cell_id < num_rows && cell_id >= 0);
      const auto candidates_begin = candidate_map.row_map(cell_id);
      const auto candidates_end = candidate_map.row_map(cell_id + 1);
      const auto tau = tolerances(DIM - 1);

      LO inside_elem = -1;
      auto inside_lambda = Omega_h::zero_vector<DIM + 1>();
      LO nearest_elem = -1;
      Real nearest_dist = -INFINITY;
      auto nearest_lambda = Omega_h::zero_vector<DIM + 1>();
      for (auto i = candidates_begin; i < candidates_end; ++i) {
        const LO elem = candidate_map.entries(i);
        const auto elem_verts =
          Omega_h::gather_verts<DIM + 1>(tris2verts, elem);
        const auto vertex_coords =
          Omega_h::gather_vectors<DIM + 1, DIM>(coords, elem_verts);
        const auto bary = barycentric_with_scales<DIM>(vertex_coords, point);
        const auto dist = bary.SignedPlaneDistance();

        if (dist >= -tau) {
          if (inside_elem < 0 || elem < inside_elem) {
            inside_elem = elem;
            inside_lambda = bary.lambda;
          }
        } else if (nearest_elem < 0 || dist > nearest_dist ||
                   (dist == nearest_dist && elem < nearest_elem)) {
          nearest_elem = elem;
          nearest_dist = dist;
          nearest_lambda = bary.lambda;
        }
      }

      if (inside_elem >= 0) {
        results(p) = GridPointSearch3D::Result{Dimensionality::REGION,
                                               inside_elem, inside_lambda};
      } else if (nearest_elem >= 0) {
        results(p) = GridPointSearch3D::Result{
          Dimensionality::REGION,
          GridPointSearch3D::Result::OutsideId(nearest_elem), nearest_lambda};
      } else {
        results(p) = GridPointSearch3D::Result{Dimensionality::NO_INTERSECT, -1,
                                               Omega_h::zero_vector<DIM + 1>()};
      }
    });

  return results;
}

GridPointSearch3D::GridPointSearch3D(Omega_h::Mesh& mesh, LO Nx, LO Ny, LO Nz)
  : GridPointSearch3D(mesh, Nx, Ny, Nz,
                      PointSearchTolerances{"point search 3d tolerances"})
{
  Kokkos::deep_copy(tolerances_, 1E-12);
}

GridPointSearch3D::GridPointSearch3D(Omega_h::Mesh& mesh, LO Nx, LO Ny, LO Nz,
                                     const PointSearchTolerances& tolerances)
  : PointLocalizationSearch(tolerances), mesh_(mesh)
{
  InitializeGeometry(mesh, Nx, Ny, Nz);
  // This search's default predicate: the conservative bounding-box test.
  InitializeMap<SimplexVsBbox>();
}

void GridPointSearch3D::InitializeGeometry(Omega_h::Mesh& mesh, LO Nx, LO Ny,
                                           LO Nz)
{
  auto mesh_bbox = Omega_h::get_bounding_box<3>(&mesh);
  auto grid_h = Kokkos::create_mirror_view(grid_);

  Kokkos::Array<Real, DIM> edge_lengths{};
  Kokkos::Array<Real, DIM> bot_left{};

  for (int i = 0; i < DIM; ++i) {
    edge_lengths[i] = mesh_bbox.max[i] - mesh_bbox.min[i];
    bot_left[i] = mesh_bbox.min[i];
  }

  grid_h(0) = Uniform3DGrid{.edge_length = edge_lengths,
                            .bot_left = bot_left,
                            .divisions = {Nx, Ny, Nz}};

  Kokkos::deep_copy(grid_, grid_h);
  num_grid_cells_ = grid_h(0).GetNumCells();
  coords_ = mesh.coords();
  tris2verts_ = mesh.ask_elem_verts();
  tris2edges_adj_ = mesh.ask_down(Omega_h::FACE, Omega_h::EDGE);
  tris2verts_adj_ = mesh.ask_down(Omega_h::FACE, Omega_h::VERT);
  edges2verts_adj_ = mesh.ask_down(Omega_h::EDGE, Omega_h::VERT);
  verts2regions_up_ = mesh.ask_up(Omega_h::VERT, Omega_h::REGION);
  edges2regions_up_ = mesh.ask_up(Omega_h::EDGE, Omega_h::REGION);
  faces2regions_up_ = mesh.ask_up(Omega_h::FACE, Omega_h::REGION);
}

LO GridPointSearch3D::GetOwningElementId(const Result& result)
{
  return pcms::GetOwningElementId(
    mesh_, 3, static_cast<int>(result.dimensionality), result.EntityId());
}

Kokkos::View<LO*> GridPointSearch3D::GetOwningElementIds(
  Kokkos::View<const Result*> results) const
{
  Kokkos::View<LO*> owners("point search owning region ids", results.extent(0));
  auto verts2regions_up = verts2regions_up_;
  auto edges2regions_up = edges2regions_up_;
  auto faces2regions_up = faces2regions_up_;
  constexpr int mesh_dim = 3;
  Kokkos::parallel_for(
    results.extent(0), KOKKOS_LAMBDA(const LO i) {
      const auto result = results(i);
      const LO element_id = result.EntityId();

      LO owner = -1;
      if (element_id >= 0) {
        if (result.dimensionality == Result::Dimensionality::REGION) {
          owner = GetOwningElementIdFromAdj(
            faces2regions_up, Result::Dimensionality::REGION,
            Result::Dimensionality::REGION, element_id);
        } else if (result.dimensionality == Result::Dimensionality::FACE) {
          owner = GetOwningElementIdFromAdj(
            faces2regions_up, Result::Dimensionality::FACE,
            Result::Dimensionality::REGION, element_id);
        } else if (result.dimensionality == Result::Dimensionality::EDGE) {
          owner = GetOwningElementIdFromAdj(
            edges2regions_up, Result::Dimensionality::EDGE,
            Result::Dimensionality::REGION, element_id);
        } else if (result.dimensionality == Result::Dimensionality::VERTEX) {
          owner = GetOwningElementIdFromAdj(
            verts2regions_up, Result::Dimensionality::VERTEX,
            Result::Dimensionality::REGION, element_id);
        }
      }
      owners(i) = owner;
    });
  return owners;
}
} // namespace pcms
