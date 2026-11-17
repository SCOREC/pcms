#include "point_search.h"
#include <Omega_h_mesh.hpp>
#include <bitset>
#include <cmath>
#include <Omega_h_simplex.hpp>

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
                         Kokkos::View<UniformGrid<Dim>[1]> grid, Real inflation)
    : elems2verts_(mesh.ask_elem_verts()),
      coords_(mesh.coords()),
      grid_(grid),
      inflation_(inflation)
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
      lo[d] = Kokkos::max(grid.AxisCellIndex(d, min - inflation_) - 1, 0);
      hi[d] = Kokkos::min(grid.AxisCellIndex(d, max + inflation_) + 1,
                          grid.divisions[d] - 1);
      idx[d] = lo[d];
    }

    LO num_overlaps = 0;
    while (true) {
      const LO cell = grid.CellIndexFromAxisIndices(idx);
      auto cell_bbox = grid.GetCellBBOX(cell);
      for (int d = 0; d < Dim; ++d) {
        cell_bbox.half_width[d] += inflation_;
      }
      if (simplex_overlaps_cell<Dim>(vertex_coords, cell_bbox)) {
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
  Real inflation_;
};

/**
 * Build the grid cell -> element candidate map. Cells are inflated by
 * `inflation` (the search tolerance), so an element within that distance of a
 * point is a candidate of the point's cell.
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
                           int num_grid_cells, Real inflation)
{
  using CrsT = Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>;
  const auto nelems = mesh.nelems();

  CrsT elem2cells{};
  Kokkos::count_and_fill_crs(
    elem2cells, nelems, ElemCellOverlapFunctor<Dim>{mesh, grid, inflation});

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
                              int num_grid_cells, Real inflation)
{
  return construct_intersection_map<2>(mesh, grid, num_grid_cells, inflation);
}

Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>
construct_intersection_map_3d(Omega_h::Mesh& mesh,
                              Kokkos::View<Uniform3DGrid[1]> grid,
                              int num_grid_cells, Real inflation)
{
  return construct_intersection_map<3>(mesh, grid, num_grid_cells, inflation);
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

namespace
{
/**
 * The entity of a simplex a point lies on, by the facet bands of the point's
 * barycentric coordinates: the point is within Cartesian distance tau of the
 * plane of the facet opposite vertex i iff |lambda_i| <= tau ||grad lambda_i||.
 * The facets the point lies on are those whose band contains it; the entity is
 * the sub-simplex spanned by the vertices of the other local indices (a vertex
 * if all facets but one, an edge if all but two, ..., the simplex itself if
 * none). For a point outside, every facet it is beyond counts as one it lies
 * on, which gives the nearest entity approximately.
 *
 * Returns the entity's dimension, sets `id` and the bitmask `mask` of the
 * local vertices spanning it.
 */
template <int Dim>
KOKKOS_INLINE_FUNCTION int entity_of_simplex(
  const BarycentricWithScales<Dim>& b, Omega_h::Real tau, bool outside,
  Omega_h::LO elem, const Omega_h::Few<Omega_h::LO, Dim + 1>& verts,
  const Omega_h::Few<Omega_h::Adj, Dim>& down, unsigned& mask, Omega_h::LO& id)
{
  mask = 0;
  int count = 0;
  int largest = 0;
  for (int i = 0; i <= Dim; ++i) {
    const auto eps = tau * b.grad_norm[i];
    const bool on =
      outside ? b.lambda[i] <= eps : Kokkos::fabs(b.lambda[i]) <= eps;
    if (!on) {
      mask |= 1u << i;
      ++count;
    }
    if (b.lambda[i] > b.lambda[largest]) {
      largest = i;
    }
  }
  if (count == 0) {
    // a simplex smaller than the tolerance: its nearest vertex
    mask = 1u << largest;
    count = 1;
  }
  const int dim = count - 1;
  if (dim == Dim) {
    id = elem;
  } else if (dim == 0) {
    id = verts[largest];
    for (int i = 0; i <= Dim; ++i) {
      if (mask == (1u << i)) {
        id = verts[i];
      }
    }
  } else {
    const int degree = Omega_h::simplex_degree(Dim, dim);
    for (int j = 0; j < degree; ++j) {
      unsigned m = 0;
      for (int v = 0; v <= dim; ++v) {
        m |= 1u << Omega_h::simplex_down_template(Dim, dim, j, v);
      }
      if (m == mask) {
        id = down[dim].ab2b[elem * degree + j];
      }
    }
  }
  return dim;
}

/**
 * The reduced coordinates (lambda_1, ..., lambda_d) of a point on an entity of
 * dimension d >= 1 of a simplex: the barycentric coordinates of the vertices
 * spanning it, renormalized to sum 1 (negative ones clamped to 0), in the
 * entity's own vertex order.
 */
template <int Dim>
KOKKOS_INLINE_FUNCTION Omega_h::Vector<Dim> reduced_coordinates(
  const BarycentricWithScales<Dim>& b, unsigned mask, int dim, Omega_h::LO id,
  const Omega_h::Few<Omega_h::LO, Dim + 1>& verts,
  const Omega_h::Few<Omega_h::Adj, Dim>& down_verts)
{
  auto xi = Omega_h::zero_vector<Dim>();
  if (dim == Dim) {
    for (int j = 0; j < Dim; ++j) {
      xi[j] = b.lambda[j + 1];
    }
  } else if (dim > 0) {
    Omega_h::Real sum = 0;
    for (int i = 0; i <= Dim; ++i) {
      if (mask & (1u << i)) {
        sum += Kokkos::fmax(b.lambda[i], 0.0);
      }
    }
    for (int m = 1; m <= dim; ++m) {
      const auto vertex = down_verts[dim].ab2b[id * (dim + 1) + m];
      for (int i = 0; i <= Dim; ++i) {
        if ((mask & (1u << i)) && verts[i] == vertex) {
          xi[m - 1] = Kokkos::fmax(b.lambda[i], 0.0) / sum;
        }
      }
    }
  }
  return xi;
}

/**
 * The result of a point that lies in the closure band of simplex `elem`
 * (SignedPlaneDistance() >= -tau, the grid search's test for a candidate):
 * the entity of the simplex the point lies on and its reduced coordinates.
 * Returns false, leaving `result` as it was, for a point beyond the band.
 */
template <int Dim, typename Result>
KOKKOS_INLINE_FUNCTION bool classify_in_element(
  const Omega_h::Vector<Dim>& point, Omega_h::LO elem,
  const Omega_h::LOs& elem_verts, const Omega_h::Reals& coords,
  const Omega_h::Few<Omega_h::Adj, Dim>& down,
  const Omega_h::Few<Omega_h::Adj, Dim>& down_verts, Real tau, Result& result)
{
  const auto verts = Omega_h::gather_verts<Dim + 1>(elem_verts, elem);
  const auto b = barycentric_with_scales<Dim>(
    Omega_h::gather_vectors<Dim + 1, Dim>(coords, verts), point);
  if (b.SignedPlaneDistance() < -tau) {
    return false;
  }
  unsigned mask;
  LO id;
  const int dim =
    entity_of_simplex<Dim>(b, tau, false, elem, verts, down, mask, id);
  result = {static_cast<typename Result::Dimensionality>(dim), id,
            reduced_coordinates<Dim>(b, mask, dim, id, verts, down_verts)};
  return true;
}

/**
 * Locates points in a simplex mesh with a uniform grid of candidate simplices.
 * Each candidate of the point's grid cell is classified with
 * entity_of_simplex using the single tolerance `tau`; the entity of lowest
 * dimension wins, then the lowest id (and, for the same entity, the lowest
 * element, whose coordinates give the reduced ones), so the result does not
 * depend on the order of the candidates and neighbouring simplices agree on a
 * shared entity. If no candidate contains the point within `tau` of its facets,
 * the candidate with the largest signed plane distance is classified as
 * outside.
 *
 * Searches the points `misses(0..num_searched)` and writes their results; all
 * points 0..num_searched if `misses` is empty.
 */
template <int Dim>
void grid_search(
  const Kokkos::View<UniformGrid<Dim>[1]>& grid,
  const Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>& candidate_map,
  const Omega_h::LOs& elem_verts, const Omega_h::Reals& coords,
  const Omega_h::Few<Omega_h::Adj, Dim>& down,
  const Omega_h::Few<Omega_h::Adj, Dim>& down_verts, Real tau,
  Kokkos::View<const Real* [Dim]> points, Kokkos::View<const LO*> misses,
  LO num_searched,
  Kokkos::View<typename PointLocalizationSearch<Dim>::Result*> results)
{
  using Search = PointLocalizationSearch<Dim>;
  using Dimensionality = typename Search::Result::Dimensionality;
  const auto num_rows = candidate_map.numRows();
  const bool all = misses.extent(0) == 0;
  Kokkos::parallel_for(
    num_searched, KOKKOS_LAMBDA(int k) {
      const int p = all ? k : misses(k);
      Omega_h::Vector<Dim> point;
      for (int i = 0; i < Dim; ++i) {
        point[i] = points(p, i);
      }
      const auto cell_id = grid(0).ClosestCellID(point);
      assert(cell_id < num_rows && cell_id >= 0);
      (void)num_rows;

      int best_dim = Dim + 1;
      LO best_id = -1;
      LO best_elem = -1;
      auto best_xi = Omega_h::zero_vector<Dim>();
      LO nearest_elem = -1;
      Real nearest_dist = -INFINITY;
      for (auto i = candidate_map.row_map(cell_id);
           i < candidate_map.row_map(cell_id + 1); ++i) {
        const LO elem = candidate_map.entries(i);
        const auto verts = Omega_h::gather_verts<Dim + 1>(elem_verts, elem);
        const auto b = barycentric_with_scales<Dim>(
          Omega_h::gather_vectors<Dim + 1, Dim>(coords, verts), point);
        const auto dist = b.SignedPlaneDistance();
        if (dist < -tau) {
          if (nearest_elem < 0 || dist > nearest_dist ||
              (dist == nearest_dist && elem < nearest_elem)) {
            nearest_elem = elem;
            nearest_dist = dist;
          }
          continue;
        }
        unsigned mask;
        LO id;
        const int dim =
          entity_of_simplex<Dim>(b, tau, false, elem, verts, down, mask, id);
        if (best_id < 0 || dim < best_dim ||
            (dim == best_dim &&
             (id < best_id || (id == best_id && elem < best_elem)))) {
          best_dim = dim;
          best_id = id;
          best_elem = elem;
          best_xi =
            reduced_coordinates<Dim>(b, mask, dim, id, verts, down_verts);
        }
      }

      if (best_id >= 0) {
        results(p) = {static_cast<Dimensionality>(best_dim), best_id, best_xi};
      } else if (nearest_elem >= 0) {
        const auto verts =
          Omega_h::gather_verts<Dim + 1>(elem_verts, nearest_elem);
        const auto b = barycentric_with_scales<Dim>(
          Omega_h::gather_vectors<Dim + 1, Dim>(coords, verts), point);
        unsigned mask;
        LO id;
        const int dim = entity_of_simplex<Dim>(b, tau, true, nearest_elem,
                                               verts, down, mask, id);
        results(p) = {
          static_cast<Dimensionality>(dim), Search::Result::OutsideId(id),
          reduced_coordinates<Dim>(b, mask, dim, id, verts, down_verts)};
      } else {
        results(p) = {Dimensionality::NO_INTERSECT, -1,
                      Omega_h::zero_vector<Dim>()};
      }
    });
}

/**
 * Locates points starting from a guessed element per point. A point with a
 * guess that lies in the closure band of the guessed element, or else of an
 * element across one of its facets (the dual graph), takes the result of that
 * element's classification; the others (no guess, or a wrong one, including
 * points outside the mesh) are collected, counted first and listed only when
 * there are some, and searched for with grid_search. Without guesses every
 * point is searched on the grid.
 */
template <int Dim>
void guessed_search(
  const Kokkos::View<UniformGrid<Dim>[1]>& grid,
  const Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>& candidate_map,
  const Omega_h::LOs& elem_verts, const Omega_h::Reals& coords,
  const Omega_h::Few<Omega_h::Adj, Dim>& down,
  const Omega_h::Few<Omega_h::Adj, Dim>& down_verts, const Omega_h::Graph& dual,
  Real tau, Kokkos::View<const Real* [Dim]> points,
  Kokkos::View<const LO*> guesses,
  Kokkos::View<typename PointLocalizationSearch<Dim>::Result*> results)
{
  using Result = typename PointLocalizationSearch<Dim>::Result;
  constexpr auto Miss = Result::Dimensionality::NO_INTERSECT;
  const auto n = static_cast<LO>(points.extent(0));
  if (guesses.extent(0) == 0) {
    grid_search<Dim>(grid, candidate_map, elem_verts, coords, down, down_verts,
                     tau, points, {}, n, results);
    return;
  }
  const auto offsets = dual.a2ab;
  const auto neighbours = dual.ab2b;
  Kokkos::parallel_for(
    "GridPointSearch::TryGuesses", n, KOKKOS_LAMBDA(const LO i) {
      Omega_h::Vector<Dim> point;
      for (int d = 0; d < Dim; ++d) {
        point[d] = points(i, d);
      }
      Result result{Miss, -1, Omega_h::zero_vector<Dim>()};
      const LO guess = guesses(i);
      if (guess >= 0 &&
          !classify_in_element<Dim>(point, guess, elem_verts, coords, down,
                                    down_verts, tau, result)) {
        for (LO a = offsets[guess]; a < offsets[guess + 1]; ++a) {
          if (classify_in_element<Dim>(point, neighbours[a], elem_verts, coords,
                                       down, down_verts, tau, result)) {
            break;
          }
        }
      }
      results(i) = result;
    });
  LO num_misses = 0;
  Kokkos::parallel_reduce(
    "GridPointSearch::CountMisses", n,
    KOKKOS_LAMBDA(const LO i, LO& count) {
      count += (results(i).dimensionality == Miss) ? 1 : 0;
    },
    num_misses);
  if (num_misses == 0) {
    return;
  }
  Kokkos::View<LO*> misses(
    Kokkos::view_alloc(Kokkos::WithoutInitializing, "point search misses"),
    num_misses);
  Kokkos::parallel_scan(
    "GridPointSearch::CollectMisses", n,
    KOKKOS_LAMBDA(const LO i, LO& before, const bool final) {
      if (results(i).dimensionality == Miss) {
        if (final) {
          misses(before) = i;
        }
        ++before;
      }
    });
  grid_search<Dim>(grid, candidate_map, elem_verts, coords, down, down_verts,
                   tau, points, misses, num_misses, results);
}
} // namespace

Kokkos::View<GridPointSearch2D::Result*> GridPointSearch2D::operator()(
  Kokkos::View<const Real* [DIM]> points) const
{
  Kokkos::View<Result*> results(
    Kokkos::view_alloc(Kokkos::WithoutInitializing, "point search result"),
    points.extent(0));
  (*this)(points, {}, results);
  return results;
}

void GridPointSearch2D::operator()(Kokkos::View<const Real* [DIM]> points,
                                   Kokkos::View<const LO*> guesses,
                                   Kokkos::View<Result*> results) const
{
  guessed_search<DIM>(grid_, candidate_map_, elem_verts_, coords_, down_,
                      down_verts_, dual_, tolerance_, points, guesses, results);
}

GridPointSearch2D::GridPointSearch2D(Omega_h::Mesh& mesh, LO Nx, LO Ny)
  : GridPointSearch2D(mesh, Nx, Ny, 1E-12)
{
}

GridPointSearch2D::GridPointSearch2D(Omega_h::Mesh& mesh, LO Nx, LO Ny,
                                     Real tolerance)
  : PointLocalizationSearch(tolerance), mesh_(mesh)
{
  auto mesh_bbox = Omega_h::get_bounding_box<2>(&mesh);
  auto grid_h = Kokkos::create_mirror_view(grid_);
  grid_h(0) =
    Uniform2DGrid{.edge_length = {mesh_bbox.max[0] - mesh_bbox.min[0],
                                  mesh_bbox.max[1] - mesh_bbox.min[1]},
                  .bot_left = {mesh_bbox.min[0], mesh_bbox.min[1]},
                  .divisions = {Nx, Ny}};
  Kokkos::deep_copy(grid_, grid_h);
  candidate_map_ = detail::construct_intersection_map_2d(
    mesh, grid_, grid_h(0).GetNumCells(), tolerance);
  coords_ = mesh.coords();
  elem_verts_ = mesh.ask_elem_verts();
  dual_ = mesh.ask_dual();
  down_[1] = mesh.ask_down(Omega_h::FACE, Omega_h::EDGE);
  down_verts_[1] = mesh.ask_down(Omega_h::EDGE, Omega_h::VERT);
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
  Kokkos::View<Result*> results(
    Kokkos::view_alloc(Kokkos::WithoutInitializing, "point search result"),
    points.extent(0));
  (*this)(points, {}, results);
  return results;
}

void GridPointSearch3D::operator()(Kokkos::View<const Real* [DIM]> points,
                                   Kokkos::View<const LO*> guesses,
                                   Kokkos::View<Result*> results) const
{
  guessed_search<DIM>(grid_, candidate_map_, elem_verts_, coords_, down_,
                      down_verts_, dual_, tolerance_, points, guesses, results);
}

GridPointSearch3D::GridPointSearch3D(Omega_h::Mesh& mesh, LO Nx, LO Ny, LO Nz)
  : GridPointSearch3D(mesh, Nx, Ny, Nz, 1E-12)
{
}

GridPointSearch3D::GridPointSearch3D(Omega_h::Mesh& mesh, LO Nx, LO Ny, LO Nz,
                                     Real tolerance)
  : PointLocalizationSearch(tolerance), mesh_(mesh)
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
  candidate_map_ = detail::construct_intersection_map_3d(
    mesh, grid_, grid_h(0).GetNumCells(), tolerance);
  coords_ = mesh.coords();
  elem_verts_ = mesh.ask_elem_verts();
  dual_ = mesh.ask_dual();
  down_[1] = mesh.ask_down(Omega_h::REGION, Omega_h::EDGE);
  down_[2] = mesh.ask_down(Omega_h::REGION, Omega_h::FACE);
  down_verts_[1] = mesh.ask_down(Omega_h::EDGE, Omega_h::VERT);
  down_verts_[2] = mesh.ask_down(Omega_h::FACE, Omega_h::VERT);
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
