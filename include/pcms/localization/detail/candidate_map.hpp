#ifndef PCMS_LOCALIZATION_DETAIL_CANDIDATE_MAP_HPP
#define PCMS_LOCALIZATION_DETAIL_CANDIDATE_MAP_HPP

#include <Kokkos_Core.hpp>

#include <type_traits>

#include "pcms/localization/overlap.hpp"
#include "pcms/localization/shapes.hpp"
#include "pcms/utility/bounding_box.hpp"
#include "pcms/utility/types.hpp"
#include "pcms/utility/uniform_grid.hpp"

namespace pcms::detail
{

using CandidateBins = Kokkos::Crs<LO, Kokkos::DefaultExecutionSpace, void, LO>;

template <typename T>
inline constexpr bool IsAABBox = false;
template <int Dim>
inline constexpr bool IsAABBox<AABBox<Dim>> = true;

/**
 * @brief The cells of `grid` overlapped by `box` widened by `pad` cells:
 * clamped to the grid, unclamped along `periodic` axes, over at most a period
 * and a cell
 */
template <int Dim>
KOKKOS_INLINE_FUNCTION CellRange<Dim> PeriodicCells(
  const AABBox<Dim>& box, const UniformGrid<Dim>& grid,
  const Kokkos::Array<bool, Dim>& periodic, LO pad)
{
  auto range = OverlappingCells<Dim>(box, grid, pad);
  for (int d = 0; d < Dim; ++d) {
    if (periodic[d]) {
      const Real n = grid.divisions[d];
      const Real x =
        (box.center[d] - grid.bot_left[d]) * n / grid.edge_length[d];
      const Real h = box.half_width[d] * n / grid.edge_length[d];
      const Real lo = Kokkos::floor(x - h) - pad;
      const Real width = Kokkos::min(Kokkos::floor(x + h) + pad - lo, n);
      // far (or NaN) boxes take the first cell
      range.lo[d] = Kokkos::fabs(lo) < 1e9 ? static_cast<LO>(lo) : 0;
      range.hi[d] = range.lo[d] + (width >= 0 ? static_cast<LO>(width) : 0);
    }
  }
  return range;
}

/** @brief The cell of per-axis indices `idx`, wrapped into the grid */
template <int Dim>
KOKKOS_INLINE_FUNCTION LO WrappedCell(const UniformGrid<Dim>& grid,
                                      Kokkos::Array<LO, Dim> idx)
{
  for (int d = 0; d < Dim; ++d) {
    idx[d] =
      (idx[d] % grid.divisions[d] + grid.divisions[d]) % grid.divisions[d];
  }
  return grid.CellIndexFromAxisIndices(idx);
}

/**
 * Functor for the shape -> grid cell map, the transpose of what the search
 * needs. Row `i` lists the cells shape `i` overlaps.
 *
 * Any shape works that has a `BoundingBox` and an `Overlaps` against a box.
 */
template <typename S>
struct ShapeCellOverlapFunctor
{
  static constexpr int Dim = ShapeDimension<S>;

  ShapeCellOverlapFunctor(Kokkos::View<const S*> shapes,
                          const UniformGrid<Dim>& grid,
                          const Kokkos::Array<bool, Dim>& periodic,
                          Real inflation)
    : shapes_(shapes), grid_(grid), periodic_(periodic), inflation_(inflation)
  {
  }

  KOKKOS_INLINE_FUNCTION
  LO operator()(LO i, LO* fill) const
  {
    const auto& shape = shapes_(i);

    // The shape's bounding box and every cell are widened by the inflation,
    // so a shape within that distance of a cell is a candidate of it.
    auto bound = BoundingBox(shape);
    for (int d = 0; d < Dim; ++d) {
      bound.half_width[d] += inflation_;
    }
    const auto range = PeriodicCells<Dim>(bound, grid_, periodic_, 1);

    Kokkos::Array<LO, Dim> idx = range.lo;
    LO num_overlaps = 0;
    while (true) {
      const LO cell = WrappedCell<Dim>(grid_, idx);
      auto cell_bbox = grid_.GetCellBBOX(cell);
      const auto wrapped = grid_.AxisIndicesFromCellIndex(cell);
      for (int d = 0; d < Dim; ++d) {
        cell_bbox.center[d] +=
          (idx[d] - wrapped[d]) * grid_.edge_length[d] / grid_.divisions[d];
        cell_bbox.half_width[d] += inflation_;
      }
      if (Overlaps(shape, cell_bbox)) {
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
  Kokkos::View<const S*> shapes_;
  UniformGrid<Dim> grid_;
  Kokkos::Array<bool, Dim> periodic_;
  Real inflation_;
};

// Sorts every row of `crs` ascending. Rows are short, so insertion sort.
inline void SortRows(CandidateBins& crs)
{
  const auto row_map = crs.row_map;
  const auto entries = crs.entries;
  Kokkos::parallel_for(
    "sort_candidate_rows", crs.numRows(), KOKKOS_LAMBDA(LO row) {
      const auto begin = row_map(row);
      const auto end = row_map(row + 1);
      for (auto i = begin + 1; i < end; ++i) {
        const LO value = entries(i);
        auto j = i;
        for (; j > begin && entries(j - 1) > value; --j) {
          entries(j) = entries(j - 1);
        }
        entries(j) = value;
      }
    });
}

// `crs` with the repeats in each row removed; rows end up ascending.
inline CandidateBins UniqueRows(CandidateBins crs);

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
template <typename S>
CandidateBins BuildCandidateBins(
  Kokkos::View<const S*> shapes, const UniformGrid<ShapeDimension<S>>& grid,
  const Kokkos::Array<bool, ShapeDimension<S>>& periodic, Real inflation)
{
  const auto num_shapes = static_cast<LO>(shapes.extent(0));
  const auto num_grid_cells = static_cast<LO>(grid.GetNumCells());

  CandidateBins elem2cells{};
  Kokkos::count_and_fill_crs(
    elem2cells, num_shapes,
    ShapeCellOverlapFunctor<S>{shapes, grid, periodic, inflation});

  Kokkos::View<LO*> counts("candidate_counts", num_grid_cells);
  Kokkos::parallel_for(
    "count_cell_candidates", elem2cells.entries.size(), KOKKOS_LAMBDA(int i) {
      Kokkos::atomic_inc(&counts(elem2cells.entries(i)));
    });

  CandidateBins cell2elems{};
  Kokkos::get_crs_row_map_from_counts(cell2elems.row_map, counts,
                                      "candidate_row_map");
  cell2elems.entries =
    CandidateBins::entries_type("candidate_entries", elem2cells.entries.size());

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
  // search results would vary run to run. A shape wider than a period is
  // in a bin once.
  cell2elems = UniqueRows(std::move(cell2elems));
  Kokkos::fence();
  return cell2elems;
}

// Row i: the bin of the cell of point i (the nearest cell for a point outside
// the grid).
template <int Dim>
struct PointCandidatesFunctor
{
  Kokkos::View<const Real* [Dim]> points;
  UniformGrid<Dim> grid;
  Kokkos::Array<bool, Dim> periodic;
  CandidateBins bins;

  KOKKOS_INLINE_FUNCTION
  LO operator()(LO i, LO* fill) const
  {
    AABBox<Dim> point{};
    for (int d = 0; d < Dim; ++d) {
      point.center[d] = points(i, d);
    }
    const LO cell =
      WrappedCell<Dim>(grid, PeriodicCells<Dim>(point, grid, periodic, 0).lo);
    const auto begin = bins.row_map(cell);
    const auto end = bins.row_map(cell + 1);
    if (fill) {
      for (auto k = begin; k < end; ++k) {
        fill[k - begin] = bins.entries(k);
      }
    }
    return static_cast<LO>(end - begin);
  }
};

// Row i: the bins of the cells box i overlaps, with repeats.
template <int Dim>
struct BoxCandidatesFunctor
{
  Kokkos::View<const AABBox<Dim>*> boxes;
  UniformGrid<Dim> grid;
  Kokkos::Array<bool, Dim> periodic;
  CandidateBins bins;

  KOKKOS_INLINE_FUNCTION
  LO operator()(LO i, LO* fill) const
  {
    // a shape touching a cell boundary is in the bins of the cells on both
    // sides, so the box needs no padding
    const auto range = PeriodicCells<Dim>(boxes(i), grid, periodic, 0);
    Kokkos::Array<LO, Dim> idx = range.lo;
    LO count = 0;
    while (true) {
      const LO cell = WrappedCell<Dim>(grid, idx);
      for (auto k = bins.row_map(cell); k < bins.row_map(cell + 1); ++k) {
        if (fill) {
          fill[count] = bins.entries(k);
        }
        ++count;
      }
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
    return count;
  }
};

// `crs` with the repeats in each row removed; rows end up ascending.
inline CandidateBins UniqueRows(CandidateBins crs)
{
  SortRows(crs);
  const auto num_rows = static_cast<LO>(crs.numRows());
  const auto row_map = crs.row_map;
  const auto entries = crs.entries;
  Kokkos::View<LO*> counts("unique_counts", num_rows);
  Kokkos::parallel_for(
    "count_unique_candidates", num_rows, KOKKOS_LAMBDA(LO row) {
      LO count = 0;
      for (auto k = row_map(row); k < row_map(row + 1); ++k) {
        count += (k == row_map(row) || entries(k) != entries(k - 1)) ? 1 : 0;
      }
      counts(row) = count;
    });
  CandidateBins unique{};
  const auto total = Kokkos::get_crs_row_map_from_counts(unique.row_map, counts,
                                                         "unique_row_map");
  unique.entries = CandidateBins::entries_type("unique_entries", total);
  const auto out_row_map = unique.row_map;
  const auto out_entries = unique.entries;
  Kokkos::parallel_for(
    "fill_unique_candidates", num_rows, KOKKOS_LAMBDA(LO row) {
      auto out = out_row_map(row);
      for (auto k = row_map(row); k < row_map(row + 1); ++k) {
        if (k == row_map(row) || entries(k) != entries(k - 1)) {
          out_entries(out++) = entries(k);
        }
      }
    });
  return unique;
}

} // namespace pcms::detail

#endif // PCMS_LOCALIZATION_DETAIL_CANDIDATE_MAP_HPP
