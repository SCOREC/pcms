#include "pcms/localization/point_search.hpp"
#include "pcms/localization/overlap.hpp"
#include "pcms/utility/memory_spaces.hpp"
#include <bitset>
#include <cmath>
#include <exception>
#include <iostream>

/**
 * Barycentric coordinates of a point in a simplex together with the factors
 * that convert a Cartesian distance from facet i into a band on coordinate i:
 * |lambda[i]| <= tau * grad_norm[i] iff the point is within distance tau of
 * the plane of the facet opposite vertex i.
 */
template <int Dim>
struct BarycentricWithScales
{
  Kokkos::Array<pcms::Real, Dim + 1> lambda;
  Kokkos::Array<pcms::Real, Dim + 1> grad_norm;

  /**
   * Signed Cartesian distance to the most violated facet plane: non-negative
   * inside the simplex, -tau on the boundary of the tolerance band.
   */
  KOKKOS_INLINE_FUNCTION pcms::Real SignedPlaneDistance() const
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
  const pcms::detail::SimplexVertices<Dim>& vertex_coords,
  const Kokkos::Array<pcms::Real, Dim>& point)
{
  // rows of the inverse basis are the gradients of lambda_1..lambda_Dim; the
  // matrices are stored by column, so row j is inv[k][j] over k
  const auto inv =
    pcms::detail::invert(pcms::detail::simplex_basis<Dim>(vertex_coords));
  Kokkos::Array<pcms::Real, Dim> r;
  for (int d = 0; d < Dim; ++d) {
    r[d] = point[d] - vertex_coords[0][d];
  }
  BarycentricWithScales<Dim> result;
  result.lambda =
    pcms::detail::form_barycentric<Dim>(pcms::detail::multiply<Dim>(inv, r));
  Kokkos::Array<pcms::Real, Dim> grad0{};
  for (int j = 0; j < Dim; ++j) {
    pcms::Real norm2 = 0;
    for (int k = 0; k < Dim; ++k) {
      norm2 += inv[k][j] * inv[k][j];
      grad0[k] -= inv[k][j];
    }
    result.grad_norm[j + 1] = Kokkos::sqrt(norm2);
  }
  pcms::Real norm2 = 0;
  for (int k = 0; k < Dim; ++k) {
    norm2 += grad0[k] * grad0[k];
  }
  result.grad_norm[0] = Kokkos::sqrt(norm2);
  return result;
}

namespace pcms::localization
{

namespace detail
{

/**
 * Functor for the shape -> grid cell map, the transpose of what the search
 * needs. Row `i` lists the cells shape `i` overlaps.
 *
 * Any shape works that has a `bounding_box` and an `overlaps` against a box.
 */
template <typename Shape>
struct ShapeCellOverlapFunctor
{
  static constexpr int Dim = dimension_v<Shape>;

  ShapeCellOverlapFunctor(Kokkos::View<const Shape*> shapes,
                          Kokkos::View<UniformGrid<Dim>[1]> grid,
                          Real inflation)
    : shapes_(shapes), grid_(grid), inflation_(inflation)
  {
  }

  KOKKOS_INLINE_FUNCTION
  LO operator()(LO i, LO* fill) const
  {
    const auto& grid = grid_(0);
    const auto& shape = shapes_(i);

    // The shape's bounding box and every cell are widened by the inflation,
    // so a shape within that distance of a cell is a candidate of it.
    auto bound = bounding_box(shape);
    for (int d = 0; d < Dim; ++d) {
      bound.half_width[d] += inflation_;
    }
    const auto range = candidate_cells<Dim>(bound, grid);

    Kokkos::Array<LO, Dim> idx = range.lo;
    LO num_overlaps = 0;
    while (true) {
      const LO cell = grid.CellIndexFromAxisIndices(idx);
      auto cell_bbox = grid.GetCellBBOX(cell);
      for (int d = 0; d < Dim; ++d) {
        cell_bbox.half_width[d] += inflation_;
      }
      if (overlaps(shape, cell_bbox)) {
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
  Kokkos::View<const Shape*> shapes_;
  Kokkos::View<UniformGrid<Dim>[1]> grid_;
  Real inflation_;
};

// helper to get the elements of a mesh as simplices.
template <int Dim>
Kokkos::View<const Simplex<Dim>*> ElementSimplices(
  Kokkos::View<const Real*> coords, Kokkos::View<const LO*> elems2verts)
{
  const auto nelems = static_cast<LO>(elems2verts.extent(0) / (Dim + 1));
  Kokkos::View<Simplex<Dim>*> elements("element_simplices", nelems);
  Kokkos::parallel_for(
    "element_simplices", nelems, KOKKOS_LAMBDA(LO e) {
      elements(e) = detail::GetElementSimplex<Dim>(coords, elems2verts, e);
    });
  return elements;
}

/**
 * Build the grid cell -> shape candidate map. Cells are inflated by
 * `inflation` (the search tolerance), so a shape within that distance of a
 * point is a candidate of the point's cell.
 *
 * Constructed as the transpose of the cheap shape -> cell map. Kokkos'
 * transpose_crs cannot be used here: it sizes the output row count from
 * in.numRows() (the shape count), which both over-allocates the row map and
 * writes out of bounds once num_grid_cells exceeds the shape count.
 */
template <typename Shape>
Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>
construct_intersection_map(
  Kokkos::View<const Shape*> shapes,
  Kokkos::View<UniformGrid<dimension_v<Shape>>[1]> grid, int num_grid_cells,
  Real inflation)
{
  using CrsT = Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>;
  const auto num_shapes = static_cast<LO>(shapes.extent(0));

  CrsT elem2cells{};
  Kokkos::count_and_fill_crs(
    elem2cells, num_shapes,
    ShapeCellOverlapFunctor<Shape>{shapes, grid, inflation});

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
    "fill_cell_candidates", num_shapes, KOKKOS_LAMBDA(LO elem) {
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
 * Cell-major brute-force build: for every grid cell, test every shape.
 *
 * This is the implementation `construct_intersection_map` replaced. It is
 * O(num_cells * nshapes) and must not be used in production, but it is the
 * ground truth for what the candidate map should contain, so it is kept and
 * exposed as an oracle for the equivalence test.
 */
template <typename Shape>
struct ReferenceCellMajorFunctor
{
  static constexpr int Dim = dimension_v<Shape>;

  ReferenceCellMajorFunctor(Kokkos::View<const Shape*> shapes,
                            Kokkos::View<UniformGrid<Dim>[1]> grid)
    : shapes_(shapes), grid_(grid)
  {
  }

  KOKKOS_INLINE_FUNCTION
  LO operator()(LO row, LO* fill) const
  {
    const auto cell_bbox = grid_(0).GetCellBBOX(row);
    LO num_intersections = 0;
    for (LO i = 0; i < static_cast<LO>(shapes_.extent(0)); ++i) {
      if (overlaps(shapes_(i), cell_bbox)) {
        if (fill) {
          fill[num_intersections] = i;
        }
        ++num_intersections;
      }
    }
    return num_intersections;
  }

private:
  Kokkos::View<const Shape*> shapes_;
  Kokkos::View<UniformGrid<Dim>[1]> grid_;
};

template <typename Shape>
Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>
construct_intersection_map_reference(
  Kokkos::View<const Shape*> shapes,
  Kokkos::View<UniformGrid<dimension_v<Shape>>[1]> grid, int num_grid_cells)
{
  Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO> intersection_map{};
  Kokkos::count_and_fill_crs(intersection_map, num_grid_cells,
                             ReferenceCellMajorFunctor<Shape>{shapes, grid});
  return intersection_map;
}

Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>
construct_intersection_map_reference_2d(Kokkos::View<const Real*> coords,
                                        Kokkos::View<const LO*> elem_verts,
                                        Kokkos::View<UniformGrid2D[1]> grid,
                                        int num_grid_cells)
{
  return construct_intersection_map_reference<Simplex<2>>(
    ElementSimplices<2>(coords, elem_verts), grid, num_grid_cells);
}

Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>
construct_intersection_map_reference_3d(Kokkos::View<const Real*> coords,
                                        Kokkos::View<const LO*> elem_verts,
                                        Kokkos::View<UniformGrid3D[1]> grid,
                                        int num_grid_cells)
{
  return construct_intersection_map_reference<Simplex<3>>(
    ElementSimplices<3>(coords, elem_verts), grid, num_grid_cells);
}

// num_grid_cells should be result of grid.GetNumCells(), take as argument to
// avoid extra copy of grid from gpu to cpu
Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>
construct_intersection_map_2d(Kokkos::View<const Real*> coords,
                              Kokkos::View<const LO*> elem_verts,
                              Kokkos::View<UniformGrid2D[1]> grid,
                              int num_grid_cells, Real inflation)
{
  return construct_intersection_map<Simplex<2>>(
    ElementSimplices<2>(coords, elem_verts), grid, num_grid_cells, inflation);
}

Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>
construct_intersection_map_3d(Kokkos::View<const Real*> coords,
                              Kokkos::View<const LO*> elem_verts,
                              Kokkos::View<UniformGrid3D[1]> grid,
                              int num_grid_cells, Real inflation)
{
  return construct_intersection_map<Simplex<3>>(
    ElementSimplices<3>(coords, elem_verts), grid, num_grid_cells, inflation);
}

} // namespace detail

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
  const BarycentricWithScales<Dim>& b, Real tau, bool outside, LO elem,
  const Kokkos::Array<LO, Dim + 1>& verts,
  const Kokkos::Array<Kokkos::View<const LO*>, Dim>& down,
  const Kokkos::Array<Kokkos::Array<unsigned, 6>, Dim>& masks, unsigned& mask,
  LO& id)
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
    // the edges of a triangle, the edges or faces of a tetrahedron
    const int degree = dim == 1 ? (Dim == 2 ? 3 : 6) : 4;
    for (int j = 0; j < degree; ++j) {
      if (masks[dim][j] == mask) {
        id = down[dim](elem * degree + j);
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
KOKKOS_INLINE_FUNCTION Kokkos::Array<Real, Dim> reduced_coordinates(
  const BarycentricWithScales<Dim>& b, unsigned mask, int dim, LO id,
  const Kokkos::Array<LO, Dim + 1>& verts,
  const Kokkos::Array<Kokkos::View<const LO*>, Dim>& down_verts)
{
  Kokkos::Array<Real, Dim> xi{};
  if (dim == Dim) {
    for (int j = 0; j < Dim; ++j) {
      xi[j] = b.lambda[j + 1];
    }
  } else if (dim > 0) {
    Real sum = 0;
    for (int i = 0; i <= Dim; ++i) {
      if (mask & (1u << i)) {
        sum += Kokkos::fmax(b.lambda[i], 0.0);
      }
    }
    for (int m = 1; m <= dim; ++m) {
      const auto vertex = down_verts[dim](id * (dim + 1) + m);
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
  const Kokkos::Array<Real, Dim>& point, LO elem,
  const Kokkos::View<const LO*>& elem_verts,
  const Kokkos::View<const Real*>& coords,
  const Kokkos::Array<Kokkos::View<const LO*>, Dim>& down,
  const Kokkos::Array<Kokkos::Array<unsigned, 6>, Dim>& masks,
  const Kokkos::Array<Kokkos::View<const LO*>, Dim>& down_verts, Real tau,
  Result& result)
{
  Kokkos::Array<LO, Dim + 1> verts;
  for (int k = 0; k < Dim + 1; ++k) {
    verts[k] = elem_verts(elem * (Dim + 1) + k);
  }
  const auto b = barycentric_with_scales<Dim>(
    detail::GetElementSimplex<Dim>(coords, elem_verts, elem).vertices, point);
  if (b.SignedPlaneDistance() < -tau) {
    return false;
  }
  unsigned mask;
  LO id;
  const int dim =
    entity_of_simplex<Dim>(b, tau, false, elem, verts, down, masks, mask, id);
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
  const Kokkos::View<const LO*>& elem_verts,
  const Kokkos::View<const Real*>& coords,
  const Kokkos::Array<Kokkos::View<const LO*>, Dim>& down,
  const Kokkos::Array<Kokkos::Array<unsigned, 6>, Dim>& masks,
  const Kokkos::Array<Kokkos::View<const LO*>, Dim>& down_verts, Real tau,
  Kokkos::View<const Real* [Dim]> points, Kokkos::View<const LO*> misses,
  LO num_searched,
  Kokkos::View<typename PointLocalizationSearch<Dim>::Result*> results)
{
  using Search = PointLocalizationSearch<Dim>;
  using Dimensionality = typename Search::Result::Dimensionality;
  const auto num_rows = candidate_map.numRows();
  const bool all = misses.extent(0) == 0;
  pcms::detail::ForEach(
    "GridPointSearch", num_searched, KOKKOS_LAMBDA(const LO k) {
      const LO p = all ? k : misses(k);
      Kokkos::Array<Real, Dim> point;
      for (int i = 0; i < Dim; ++i) {
        point[i] = points(p, i);
      }
      const auto cell_id = grid(0).ClosestCellID(point);
      assert(cell_id < num_rows && cell_id >= 0);
      (void)num_rows;

      int best_dim = Dim + 1;
      LO best_id = -1;
      LO best_elem = -1;
      Kokkos::Array<Real, Dim> best_xi{};
      LO nearest_elem = -1;
      Real nearest_dist = -INFINITY;
      for (auto i = candidate_map.row_map(cell_id);
           i < candidate_map.row_map(cell_id + 1); ++i) {
        const LO elem = candidate_map.entries(i);
        Kokkos::Array<LO, Dim + 1> verts;
        for (int k = 0; k < Dim + 1; ++k) {
          verts[k] = elem_verts(elem * (Dim + 1) + k);
        }
        const auto b = barycentric_with_scales<Dim>(
          detail::GetElementSimplex<Dim>(coords, elem_verts, elem).vertices,
          point);
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
        const int dim = entity_of_simplex<Dim>(b, tau, false, elem, verts, down,
                                               masks, mask, id);
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
        Kokkos::Array<LO, Dim + 1> verts;
        for (int k = 0; k < Dim + 1; ++k) {
          verts[k] = elem_verts(nearest_elem * (Dim + 1) + k);
        }
        const auto b = barycentric_with_scales<Dim>(
          detail::GetElementSimplex<Dim>(coords, elem_verts, nearest_elem)
            .vertices,
          point);
        unsigned mask;
        LO id;
        const int dim = entity_of_simplex<Dim>(b, tau, true, nearest_elem,
                                               verts, down, masks, mask, id);
        results(p) = {
          static_cast<Dimensionality>(dim), Search::Result::OutsideId(id),
          reduced_coordinates<Dim>(b, mask, dim, id, verts, down_verts)};
      } else {
        results(p) = {Dimensionality::NO_INTERSECT, -1,
                      Kokkos::Array<Real, Dim>{}};
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
  const Kokkos::View<const LO*>& elem_verts,
  const Kokkos::View<const Real*>& coords,
  const Kokkos::Array<Kokkos::View<const LO*>, Dim>& down,
  const Kokkos::Array<Kokkos::Array<unsigned, 6>, Dim>& masks,
  const Kokkos::Array<Kokkos::View<const LO*>, Dim>& down_verts,
  const Kokkos::View<const LO*>& dual_offsets,
  const Kokkos::View<const LO*>& dual_neighbors, Real tau,
  Kokkos::View<const Real* [Dim]> points, Kokkos::View<const LO*> guesses,
  Kokkos::View<typename PointLocalizationSearch<Dim>::Result*> results)
{
  using Result = typename PointLocalizationSearch<Dim>::Result;
  constexpr auto Miss = Result::Dimensionality::NO_INTERSECT;
  const auto n = static_cast<LO>(points.extent(0));
  if (guesses.extent(0) == 0) {
    grid_search<Dim>(grid, candidate_map, elem_verts, coords, down, masks,
                     down_verts, tau, points, {}, n, results);
    return;
  }
  const auto offsets = dual_offsets;
  const auto neighbours = dual_neighbors;
  pcms::detail::ForEach(
    "GridPointSearch::TryGuesses", n, KOKKOS_LAMBDA(const LO i) {
      Kokkos::Array<Real, Dim> point;
      for (int d = 0; d < Dim; ++d) {
        point[d] = points(i, d);
      }
      Result result{Miss, -1, Kokkos::Array<Real, Dim>{}};
      const LO guess = guesses(i);
      if (guess >= 0 &&
          !classify_in_element<Dim>(point, guess, elem_verts, coords, down,
                                    masks, down_verts, tau, result)) {
        for (LO a = offsets(guess); a < offsets(guess + 1); ++a) {
          if (classify_in_element<Dim>(point, neighbours(a), elem_verts, coords,
                                       down, masks, down_verts, tau, result)) {
            break;
          }
        }
      }
      results(i) = result;
    });
  if (pcms::detail::IsSmallHostBatch(static_cast<std::size_t>(n))) {
    LO num_misses = 0;
    for (LO i = 0; i < n; ++i) {
      num_misses += (results(i).dimensionality == Miss) ? 1 : 0;
    }
    if (num_misses == 0) {
      return;
    }
    Kokkos::View<LO*> misses(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "point search misses"),
      num_misses);
    for (LO i = 0, k = 0; i < n; ++i) {
      if (results(i).dimensionality == Miss) {
        misses(k++) = i;
      }
    }
    grid_search<Dim>(grid, candidate_map, elem_verts, coords, down, masks,
                     down_verts, tau, points, misses, num_misses, results);
    return;
  }
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
  grid_search<Dim>(grid, candidate_map, elem_verts, coords, down, masks,
                   down_verts, tau, points, misses, num_misses, results);
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
                      entity_vertex_masks_, down_verts_, dual_offsets_,
                      dual_neighbors_, tolerance_, points, guesses, results);
}

GridPointSearch2D::GridPointSearch2D(
  Kokkos::View<const Real*> coords, Kokkos::View<const LO*> elem_verts,
  Kokkos::Array<Kokkos::View<const LO*>, DIM> down,
  Kokkos::Array<Kokkos::View<const LO*>, DIM> down_verts,
  Kokkos::Array<Kokkos::Array<unsigned, 6>, DIM> entity_vertex_masks,
  Kokkos::View<const LO*> dual_offsets, Kokkos::View<const LO*> dual_neighbors,
  const UniformGrid2D& grid, Real tolerance)
  : PointLocalizationSearch(tolerance),
    down_(down),
    down_verts_(down_verts),
    entity_vertex_masks_(entity_vertex_masks),
    elem_verts_(elem_verts),
    coords_(coords),
    dual_offsets_(dual_offsets),
    dual_neighbors_(dual_neighbors)
{
  auto grid_h = Kokkos::create_mirror_view(grid_);
  grid_h(0) = grid;
  Kokkos::deep_copy(grid_, grid_h);
  candidate_map_ = detail::construct_intersection_map_2d(
    coords, elem_verts, grid_, grid.GetNumCells(), tolerance);
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
                      entity_vertex_masks_, down_verts_, dual_offsets_,
                      dual_neighbors_, tolerance_, points, guesses, results);
}

GridPointSearch3D::GridPointSearch3D(
  Kokkos::View<const Real*> coords, Kokkos::View<const LO*> elem_verts,
  Kokkos::Array<Kokkos::View<const LO*>, DIM> down,
  Kokkos::Array<Kokkos::View<const LO*>, DIM> down_verts,
  Kokkos::Array<Kokkos::Array<unsigned, 6>, DIM> entity_vertex_masks,
  Kokkos::View<const LO*> dual_offsets, Kokkos::View<const LO*> dual_neighbors,
  const UniformGrid3D& grid, Real tolerance)
  : PointLocalizationSearch(tolerance),
    down_(down),
    down_verts_(down_verts),
    entity_vertex_masks_(entity_vertex_masks),
    elem_verts_(elem_verts),
    coords_(coords),
    dual_offsets_(dual_offsets),
    dual_neighbors_(dual_neighbors)
{
  auto grid_h = Kokkos::create_mirror_view(grid_);
  grid_h(0) = grid;
  Kokkos::deep_copy(grid_, grid_h);
  candidate_map_ = detail::construct_intersection_map_3d(
    coords, elem_verts, grid_, grid.GetNumCells(), tolerance);
}
} // namespace pcms::localization
