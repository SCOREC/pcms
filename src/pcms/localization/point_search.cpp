#include "point_search.h"
#include <Omega_h_mesh.hpp>
#include <bitset>
#include <cmath>

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

KOKKOS_INLINE_FUNCTION
AABBox<2> triangle_bbox(const Omega_h::Matrix<2, 3>& coords)
{
  Kokkos::Array<Real, 2> max{coords(0, 0), coords(1, 0)};
  Kokkos::Array<Real, 2> min{coords(0, 0), coords(1, 0)};
  for (int i = 1; i < 3; ++i) {
    max[0] = std::fmax(max[0], coords(0, i));
    max[1] = std::fmax(max[1], coords(1, i));
    min[0] = std::fmin(min[0], coords(0, i));
    min[1] = std::fmin(min[1], coords(1, i));
  }
  return {.center = {(max[0] + min[0]) / 2.0, (max[1] + min[1]) / 2.0},
          .half_width = {(max[0] - min[0]) / 2.0, (max[1] - min[1]) / 2.0}};
}

template <unsigned dim>
KOKKOS_INLINE_FUNCTION AABBox<dim> simplex_bbox(
  const Omega_h::Matrix<dim, dim + 1>& coords)
{
  Kokkos::Array<Real, dim> max;
  Kokkos::Array<Real, dim> min;
  for (int j = 0; j < dim; ++j) {
    max[j] = coords(j, 0);
    min[j] = coords(j, 0);
  }
  for (int i = 1; i < dim + 1; ++i) {
    for (int j = 0; j < dim; ++j) {
      max[j] = std::fmax(max[j], coords(j, i));
      min[j] = std::fmin(min[j], coords(j, i));
    }
  }

  Kokkos::Array<Real, dim> center;
  Kokkos::Array<Real, dim> half_width;

  for (int j = 0; j < dim; ++j) {
    center[j] = (max[j] + min[j]) / 2.0;
    half_width[j] = (max[j] - min[j]) / 2.0;
  }

  return {.center = center, .half_width = half_width};
}

// Liang, You-Dong, and B. A. Barsky. “A New Concept and Method for Line
// Clipping.” ACM Transactions on Graphics 3, no. 1 (January 1984): 1–22.
// https://doi.org/10.1145/357332.357333.
KOKKOS_INLINE_FUNCTION
bool clipt(Real p, Real q, Real& t0, Real& t1)
{
  if (p < 0) {
    auto r = q / p;
    if (r > t1) {
      return false;
    }
    if (r > t0) {
      t0 = r;
    }
  } else if (p > 0) {
    auto r = q / p;
    if (r < t0) {
      return false;
    }
    if (r < t1) {
      t1 = r;
    }
  }
  // p==0
  else {
    if (q < 0) {
      return false;
    }
  }
  return true;
}
// Liang, You-Dong, and B. A. Barsky. “A New Concept and Method for Line
// Clipping.” ACM Transactions on Graphics 3, no. 1 (January 1984): 1–22.
// https://doi.org/10.1145/357332.357333.
KOKKOS_INLINE_FUNCTION
bool line_intersects_bbox(const Omega_h::Vector<2>& p0,
                          const Omega_h::Vector<2>& p1, const AABBox<2>& bbox)
{
  Real t0 = 0;
  Real t1 = 1;
  auto xleft = bbox.center[0] - bbox.half_width[0];
  auto deltax = p1[0] - p0[0];
  if (clipt(-deltax, p0[0] - xleft, t0, t1)) {
    auto xright = bbox.center[0] + bbox.half_width[0];
    if (clipt(deltax, xright - p0[0], t0, t1)) {
      auto deltay = p1[1] - p0[1];
      auto ybottom = bbox.center[1] - bbox.half_width[1];
      if (clipt(-deltay, p0[1] - ybottom, t0, t1)) {
        auto ytop = bbox.center[1] + bbox.half_width[1];
        if (clipt(deltay, ytop - p0[1], t0, t1)) {
          // full liang-barksy algorithm computes new x-y coordinates here, but
          // we just need to check intersections
          return true;
        }
      }
    }
  }
  return false;
}

template <unsigned dim>
[[nodiscard]] KOKKOS_INLINE_FUNCTION bool within_bbox(
  const Omega_h::Vector<dim> coord, const AABBox<dim>& bbox) noexcept
{
  for (int i = 0; i < dim; ++i) {
    if (coord[i] < bbox.center[i] - bbox.half_width[i])
      return false;
    if (coord[i] > bbox.center[i] + bbox.half_width[i])
      return false;
  }
  return true;
}

template <int dim>
[[nodiscard]] KOKKOS_INLINE_FUNCTION bool bbox_verts_within_simplex(
  const AABBox<dim>& bbox, const Omega_h::Matrix<dim, dim + 1>& coords)
{
  // each dimension has a pair of opposing "walls"
  // 2D: { [left, right], [top, bottom] } -> { left, right, top, bottom }
  // 3D: { [left, right], [top, bottom], [front, back] } -> { left, ..., back }
  Kokkos::Array<Real, dim * 2ul> bbox_walls{};
  for (int i = 0; i < dim; i++) {
    bbox_walls[i * 2] = bbox.center[i] - bbox.half_width[i];
    bbox_walls[i * 2 + 1] = bbox.center[i] + bbox.half_width[i];
  }

  // 1 << dim == 2 ** dim == num vertices in ndim bounding box / hypercube
  constexpr unsigned num_verts = 1 << dim;

  // Each vertex is a just a unique combination of walls
  // eg [left, bottom] (2D) or [right, top, front] (3)
  for (unsigned i = 0; i < num_verts; ++i) {
    // conveniently, i acts a bit field representing the current combination
    Omega_h::Vector<dim> vert;
    for (unsigned j = 0; j < dim; ++j) {
      // eg 110 = 6 = [left, top, back]
      // eg 01 = 1 = [right, top]
      vert[j] = (i >> j) & 1 ? bbox_walls[j * 2] : bbox_walls[j * 2 + 1];
    }
    auto xi = Omega_h::barycentric_from_global<dim, dim>(vert, coords);
    if (Omega_h::is_barycentric_inside(xi)) {
      return true;
    }
  }
  return false;
}

/**
 * Check if a triangle element represented by 3 coordinates in two dimensions
 * intersects with a bounding box
 */
[[nodiscard]] KOKKOS_FUNCTION bool triangle_intersects_bbox(
  const Omega_h::Matrix<2, 3>& coords, const AABBox<2>& bbox)
{
  // triangle and grid cell bounding box intersect
  if (intersects(triangle_bbox(coords), bbox)) {
    // if any of the triangle verts inside of bbox
    if (within_bbox<2>(coords[0], bbox) || within_bbox<2>(coords[1], bbox) ||
        within_bbox<2>(coords[2], bbox)) {
      return true;
    }
    // if any of the bbox verts are within the triangle
    if (bbox_verts_within_simplex(bbox, coords)) {
      return true;
    }
    // if any of the triangle's edges intersect with the bounding box
    if (line_intersects_bbox(coords[0], coords[1], bbox) ||
        line_intersects_bbox(coords[1], coords[2], bbox) ||
        line_intersects_bbox(coords[2], coords[0], bbox)) {
      return true;
    }
  }
  return false;
}

template <unsigned dim>
[[nodiscard]] KOKKOS_FUNCTION bool simplex_intersects_bbox(
  const Omega_h::Matrix<dim, dim + 1>& coords, const AABBox<dim>& bbox)
{
  return intersects(simplex_bbox<dim>(coords), bbox);
  // TODO: Add refined cases from triangle_intersects_bbox
}

namespace detail
{
/**
 *  Functor for constructing the mapping from grid cells to intersecting
 * triangles. Each row of the resulting CSR structure represents a grid cell and
 * each row entry corresponds to an ID of an element that intersects that grid
 * cell.
 * \Warning since this works on Omega_h meshes, we currently assume each element
 * is a 2D simplex (triangle)
 * \Warning since this uses Omega_h data which is only available in the
 * "Default" Execution space, the should not be used in an alternative EXE space
 */
/// The bbox-overlap predicate for a simplex against one grid cell. 2D applies
/// the refined edge tests, 3D is a plain bbox overlap; both match what the
/// previous cell-major functors used, so candidate lists are unchanged.
template <int Dim>
[[nodiscard]] KOKKOS_INLINE_FUNCTION bool simplex_overlaps_cell(
  const Omega_h::Matrix<Dim, Dim + 1>& coords, const AABBox<Dim>& bbox)
{
  if constexpr (Dim == 2) {
    return triangle_intersects_bbox(coords, bbox);
  } else {
    return simplex_intersects_bbox<Dim>(coords, bbox);
  }
}

/**
 * Functor for the element -> grid cell map, the transpose of what the search
 * needs. Row `elem` lists the cells that element overlaps.
 *
 * This is the element-major inverse of the old cell-major functor: rather than
 * testing every element against a cell, each element only visits the cells its
 * own bounding box spans. That makes construction O(nelems * cells_per_elem)
 * instead of O(num_cells * nelems).
 */
template <int Dim>
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

    // Per-axis inclusive range of cells the element's bbox spans, widened by
    // one cell on each side. The overlap predicate treats a shared boundary
    // plane as an intersection, so a cell that only touches the bbox counts as
    // a candidate and would be missed by the unwidened range -- which happens
    // constantly on a structured mesh whose vertices land on grid lines. One
    // cell is enough: a cell two or more away cannot touch the bbox at all.
    // The predicate below rejects whatever the widening over-includes.
    Kokkos::Array<LO, Dim> lo, hi, idx;
    for (int d = 0; d < Dim; ++d) {
      Real min = vertex_coords[0][d];
      Real max = vertex_coords[0][d];
      for (int v = 1; v < Dim + 1; ++v) {
        min = Kokkos::min(min, vertex_coords[v][d]);
        max = Kokkos::max(max, vertex_coords[v][d]);
      }
      lo[d] = Kokkos::max(grid.AxisCellIndex(d, min) - 1, 0);
      hi[d] =
        Kokkos::min(grid.AxisCellIndex(d, max) + 1, grid.divisions[d] - 1);
      idx[d] = lo[d];
    }

    LO num_overlaps = 0;
    while (true) {
      const LO cell = grid.CellIndexFromAxisIndices(idx);
      if (simplex_overlaps_cell<Dim>(vertex_coords, grid.GetCellBBOX(cell))) {
        if (fill) {
          fill[num_overlaps] = cell;
        }
        ++num_overlaps;
      }
      // odometer over the inclusive range [lo, hi]
      int d = 0;
      for (; d < Dim; ++d) {
        if (++idx[d] <= hi[d]) {
          break;
        }
        idx[d] = lo[d];
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
template <int Dim>
Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>
construct_intersection_map(Omega_h::Mesh& mesh,
                           Kokkos::View<UniformGrid<Dim>[1]> grid,
                           int num_grid_cells)
{
  using CrsT = Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>;
  const auto nelems = mesh.nelems();

  CrsT elem2cells{};
  Kokkos::count_and_fill_crs(elem2cells, nelems,
                             ElemCellOverlapFunctor<Dim>{mesh, grid});

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
template <int Dim>
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
      if (simplex_overlaps_cell<Dim>(vertex_coords, cell_bbox)) {
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

template <int Dim>
Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>
construct_intersection_map_reference(Omega_h::Mesh& mesh,
                                     Kokkos::View<UniformGrid<Dim>[1]> grid,
                                     int num_grid_cells)
{
  Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO> intersection_map{};
  Kokkos::count_and_fill_crs(intersection_map, num_grid_cells,
                             ReferenceCellMajorFunctor<Dim>{mesh, grid});
  return intersection_map;
}

Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>
construct_intersection_map_reference_2d(Omega_h::Mesh& mesh,
                                        Kokkos::View<Uniform2DGrid[1]> grid,
                                        int num_grid_cells)
{
  return construct_intersection_map_reference<2>(mesh, grid, num_grid_cells);
}

Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>
construct_intersection_map_reference_3d(Omega_h::Mesh& mesh,
                                        Kokkos::View<Uniform3DGrid[1]> grid,
                                        int num_grid_cells)
{
  return construct_intersection_map_reference<3>(mesh, grid, num_grid_cells);
}

// num_grid_cells should be result of grid.GetNumCells(), take as argument to
// avoid extra copy of grid from gpu to cpu
Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>
construct_intersection_map_2d(Omega_h::Mesh& mesh,
                              Kokkos::View<Uniform2DGrid[1]> grid,
                              int num_grid_cells)
{
  return construct_intersection_map<2>(mesh, grid, num_grid_cells);
}

Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>
construct_intersection_map_3d(Omega_h::Mesh& mesh,
                              Kokkos::View<Uniform3DGrid[1]> grid,
                              int num_grid_cells)
{
  return construct_intersection_map<3>(mesh, grid, num_grid_cells);
}
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
        GridPointSearch2D::Result::Dimensionality::REGION;
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
        // Outside mesh - beyond tolerance of any entity
        // Only negate if we truly have no candidates at all.
        if (best_vertex_id >= 0 &&
            (best_vertex_dist <= best_edge_dist || best_edge_id < 0)) {
          dim_out = GridPointSearch2D::Result::Dimensionality::VERTEX;
          element_id_out = best_vertex_id;
          bary_out = best_vertex_bary;
          element_id_out = -element_id_out;
        } else if (best_edge_id >= 0) {
          dim_out = GridPointSearch2D::Result::Dimensionality::EDGE;
          element_id_out = best_edge_id;
          bary_out = best_edge_bary_min;
          element_id_out = -element_id_out;
        } else {
          // No candidates at all: both IDs stay -1
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
  auto mesh_bbox = Omega_h::get_bounding_box<2>(&mesh);
  auto grid_h = Kokkos::create_mirror_view(grid_);
  grid_h(0) =
    Uniform2DGrid{.edge_length = {mesh_bbox.max[0] - mesh_bbox.min[0],
                                  mesh_bbox.max[1] - mesh_bbox.min[1]},
                  .bot_left = {mesh_bbox.min[0], mesh_bbox.min[1]},
                  .divisions = {Nx, Ny}};
  Kokkos::deep_copy(grid_, grid_h);
  // Determine inflation radius from tolerances (max of vertex/edge tol)
  auto tol_h = Kokkos::create_mirror_view(tolerances_);
  Kokkos::deep_copy(tol_h, tolerances_);
  candidate_map_ =
    detail::construct_intersection_map_2d(mesh, grid_, grid_h(0).GetNumCells());
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
  const LO query_id =
    (result.element_id < 0) ? -result.element_id : result.element_id;
  return pcms::GetOwningElementId(
    mesh_, 2, static_cast<int>(result.dimensionality), query_id);
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
      LO element_id = result.element_id;
      if (element_id < 0)
        element_id = -element_id;

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
        results(p) = GridPointSearch3D::Result{Dimensionality::REGION,
                                               -nearest_elem, nearest_lambda};
      } else {
        results(p) = GridPointSearch3D::Result{Dimensionality::REGION, -1,
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
  candidate_map_ =
    detail::construct_intersection_map_3d(mesh, grid_, grid_h(0).GetNumCells());
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
  const LO query_id =
    (result.element_id < 0) ? -result.element_id : result.element_id;
  return pcms::GetOwningElementId(
    mesh_, 3, static_cast<int>(result.dimensionality), query_id);
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
      LO element_id = result.element_id;
      if (element_id < 0)
        element_id = -element_id;

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
