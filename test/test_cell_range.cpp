#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <pcms/localization/overlap.h>
#include <pcms/utility/uniform_grid.h>

// Serial tests of the cell-range layer the routing table is built on: which
// grid cells a bounding box covers (`candidate_cells`, `overlapping_cells`) and
// the invariants of `CellRange`.

using pcms::AABBox;
using pcms::candidate_cells;
using pcms::CellRange;
using pcms::LO;
using pcms::Overlap;
using pcms::overlapping_cells;
using pcms::Real;
using pcms::SimplexVsBbox;
using pcms::TriangleVsBbox;
using pcms::Uniform2DGrid;
using pcms::UniformGrid;

TEST_CASE("bbox vs grid cell range")
{
  Uniform2DGrid grid{
    .edge_length = {10, 10}, .bot_left = {0, 0}, .divisions = {10, 10}};

  SECTION("box inside one cell still covers padded neighbors")
  {
    AABBox<2> box{.center = {0.5, 0.5}, .half_width = {0.1, 0.1}};
    const auto range = overlapping_cells<2>(box, grid);
    // default pad = 1 clamps at the grid border
    REQUIRE(range.lo[0] == 0);
    REQUIRE(range.lo[1] == 0);
    REQUIRE(range.hi[0] == 1);
    REQUIRE(range.hi[1] == 1);
  }

  SECTION("no padding gives the exact span")
  {
    AABBox<2> box{.center = {2.5, 2.5}, .half_width = {0.4, 0.4}};
    const auto range = overlapping_cells<2>(box, grid, 0);
    REQUIRE(range.lo[0] == 2);
    REQUIRE(range.lo[1] == 2);
    REQUIRE(range.hi[0] == 2);
    REQUIRE(range.hi[1] == 2);
  }

  SECTION("box spanning several cells")
  {
    AABBox<2> box{.center = {2.5, 2.5}, .half_width = {1.0, 1.0}};
    const auto range = overlapping_cells<2>(box, grid, 0);
    REQUIRE(range.lo[0] == 1);
    REQUIRE(range.lo[1] == 1);
    REQUIRE(range.hi[0] == 3);
    REQUIRE(range.hi[1] == 3);
  }
}

namespace
{
/// True when grid cell `(cx, cy)` and `box` meet, each treated as closed.
///
/// A shared boundary counts, which is the semantics the overlap predicate has
/// and therefore the set `candidate_cells` has to cover.
[[nodiscard]] bool cell_touches_box(const Uniform2DGrid& grid, LO cx, LO cy,
                                    const AABBox<2>& box)
{
  const Kokkos::Array<LO, 2> cell{cx, cy};
  for (int d = 0; d < 2; ++d) {
    const Real width = grid.edge_length[d] / grid.divisions[d];
    const Real cell_lo = grid.bot_left[d] + cell[d] * width;
    const Real cell_hi = cell_lo + width;
    const Real box_lo = box.center[d] - box.half_width[d];
    const Real box_hi = box.center[d] + box.half_width[d];
    if (cell_hi < box_lo || cell_lo > box_hi) {
      return false;
    }
  }
  return true;
}
} // namespace

/// The property `candidate_cells` exists to guarantee: its range contains every
/// cell a closed-boundary predicate can accept. A range that is too narrow
/// drops candidates silently -- the element is never offered to the query and
/// nothing downstream can tell -- so this is checked as a coverage property
/// rather than against the exact bounds, which are allowed to over-include.
TEST_CASE(
  "candidate_cells covers every cell a closed-boundary predicate accepts")
{
  const Uniform2DGrid grid{
    .edge_length = {10, 10}, .bot_left = {0, 0}, .divisions = {10, 10}};

  // Centres at, near and between grid lines, and extents that land exactly on a
  // grid line: that is when an unwidened range loses the touching cell.
  for (const Real center : {0.5, 2.5, 5.0, 9.5}) {
    for (const Real half : {0.0, 0.25, 0.5, 1.0, 3.0}) {
      const AABBox<2> box{.center = {center, center},
                          .half_width = {half, half}};
      const auto range = candidate_cells<2, TriangleVsBbox>(box, grid);
      for (LO cx = 0; cx < grid.divisions[0]; ++cx) {
        for (LO cy = 0; cy < grid.divisions[1]; ++cy) {
          if (cell_touches_box(grid, cx, cy, box)) {
            REQUIRE(range.Contains(Kokkos::Array<LO, 2>{cx, cy}));
          }
        }
      }
    }
  }
}

TEST_CASE("candidate_cells bins a supplied bound")
{
  const Uniform2DGrid grid{
    .edge_length = {10, 10}, .bot_left = {0, 0}, .divisions = {10, 10}};

  // A degenerate triangle whose bounding box is the single grid line x = y = 3:
  // the cells on either side of that line only touch it, so the padding must
  // come from the predicate rather than from the caller.
  const Omega_h::Matrix<2, 3> degenerate{{3.0, 3.0}, {3.0, 3.0}, {3.0, 3.0}};
  // `candidate_cells` bins a bound it is given; the element side supplies it.
  const auto range = candidate_cells<2, TriangleVsBbox>(
    pcms::detail::simplex_bbox<2>(degenerate), grid);
  REQUIRE(range.Contains(Kokkos::Array<LO, 2>{3, 3}));
  REQUIRE(range.Contains(Kokkos::Array<LO, 2>{2, 2}));
}

TEST_CASE("CellRange states its own invariants")
{
  const auto empty = pcms::CellRange<2>::Empty();
  REQUIRE(empty.IsEmpty());
  REQUIRE(!empty.Contains(Kokkos::Array<LO, 2>{0, 0}));

  const pcms::CellRange<2> unit{.lo = {2, 2}, .hi = {2, 2}};
  REQUIRE(!unit.IsEmpty());
  REQUIRE(unit.Contains(Kokkos::Array<LO, 2>{2, 2}));
  REQUIRE(!unit.Contains(Kokkos::Array<LO, 2>{3, 2}));

  // The empty range must not be what value initialization produces: a
  // zero-initialized CellRange denotes the single cell at axis index 0.
  const pcms::CellRange<2> zeroed{};
  REQUIRE(!zeroed.IsEmpty());
  REQUIRE(zeroed.Contains(Kokkos::Array<LO, 2>{0, 0}));
}

// What dimension does the overlap + grid-cell-range layer actually support?
//
// `Dim` is a template parameter throughout and every loop runs to `Dim`, so
// nothing in these functions is written for 2 or 3 specifically. That is easy
// to lose by accident -- one hardcoded `[2]` or `<3, 3>` and the layer silently
// becomes 2D/3D only -- so it is pinned here at dimensions no mesh can
// represent, which is exactly where a hidden assumption would show up.
//
// Findings, and the boundary they draw:
//
//   generic, and covered here   AABBox, intersects, simplex_bbox, CellRange,
//                               overlapping_cells, candidate_cells,
//                               UniformGrid, Overlap<Dim, SimplexVsBbox>
//   generic over Dim, but the   MeshSimplexElements, ElemCellOverlapFunctor,
//   map builder needs a mesh    construct_intersection_map -- all take `Dim`,
//                               but are fed an Omega_h mesh, which is <= 3D, so
//                               dimensions above 3 cannot be reached through
//                               them
//   NOT generic                 GridPointSearch2D / GridPointSearch3D, which
//   are
//                               concrete classes with hand-written
//                               per-dimension query logic (VERTEX/EDGE/FACE for
//                               2D, REGION for 3D)
//   NOT generic, was silent     Overlap<Dim, SimplexVsSimplex>: nominally
//                               Dim-generic but its face normal is a 2D
//                               perpendicular or a 3D cross product. Now
//                               guarded by a static_assert, so asking for
//                               another dimension is a compile error instead of
//                               reading uninitialized memory.
TEST_CASE("overlap and cell-range layers are dimension-generic")
{
  SECTION("1D")
  {
    const UniformGrid<1> grid{
      .edge_length = {4.0}, .bot_left = {0.0}, .divisions = {4}};
    const AABBox<1> box{.center = {1.5}, .half_width = {0.25}};
    // Cell 1, padded by 1, clamped at the low edge.
    const auto range = candidate_cells<1, SimplexVsBbox>(box, grid);
    REQUIRE(range.lo[0] == 0);
    REQUIRE(range.hi[0] == 2);
    REQUIRE(range.Contains(Kokkos::Array<LO, 1>{1}));
  }

  SECTION("4D")
  {
    const UniformGrid<4> grid{.edge_length = {4.0, 4.0, 4.0, 4.0},
                              .bot_left = {0.0, 0.0, 0.0, 0.0},
                              .divisions = {4, 4, 4, 4}};
    const AABBox<4> box{.center = {1.5, 1.5, 1.5, 1.5},
                        .half_width = {0.25, 0.25, 0.25, 0.25}};
    const auto range = candidate_cells<4, SimplexVsBbox>(box, grid);
    for (int d = 0; d < 4; ++d) {
      REQUIRE(range.lo[d] == 0);
      REQUIRE(range.hi[d] == 2);
    }
    REQUIRE(range.Contains(Kokkos::Array<LO, 4>{1, 1, 1, 1}));
    // Same axis order as 2D: CellIndexFromAxisIndices round-trips.
    const auto indices = Kokkos::Array<LO, 4>{1, 1, 1, 1};
    const auto linear = grid.CellIndexFromAxisIndices(indices);
    const auto back = grid.AxisIndicesFromCellIndex(linear);
    for (int d = 0; d < 4; ++d) {
      REQUIRE(back[d] == indices[d]);
    }
  }

  SECTION("5D range invariants")
  {
    const auto empty = CellRange<5>::Empty();
    REQUIRE(empty.IsEmpty());
    REQUIRE(!empty.Contains(Kokkos::Array<LO, 5>{0, 0, 0, 0, 0}));

    const CellRange<5> unit{.lo = {3, 3, 3, 3, 3}, .hi = {3, 3, 3, 3, 3}};
    REQUIRE(!unit.IsEmpty());
    REQUIRE(unit.Contains(Kokkos::Array<LO, 5>{3, 3, 3, 3, 3}));
    REQUIRE(!unit.Contains(Kokkos::Array<LO, 5>{3, 3, 3, 3, 4}));
  }

  SECTION("SimplexVsBbox predicate in 4D")
  {
    const AABBox<4> box{.center = {1.0, 1.0, 1.0, 1.0},
                        .half_width = {0.5, 0.5, 0.5, 0.5}};

    // A degenerate simplex: all 5 vertices at a point inside the box.
    Omega_h::Matrix<4, 5> inside{};
    for (int v = 0; v < 5; ++v) {
      for (int d = 0; d < 4; ++d) {
        inside(d, v) = 1.0;
      }
    }
    REQUIRE(Overlap<4, SimplexVsBbox>::apply(inside, box));

    // Far outside on every axis.
    Omega_h::Matrix<4, 5> outside{};
    for (int v = 0; v < 5; ++v) {
      for (int d = 0; d < 4; ++d) {
        outside(d, v) = 100.0;
      }
    }
    REQUIRE_FALSE(Overlap<4, SimplexVsBbox>::apply(outside, box));

    // Straddling the box on one axis only: the conservative predicate is true.
    Omega_h::Matrix<4, 5> straddling{};
    for (int v = 0; v < 5; ++v) {
      for (int d = 0; d < 4; ++d) {
        straddling(d, v) = 1.0;
      }
    }
    straddling(3, 0) = 0.5;
    straddling(3, 1) = 1.9;
    REQUIRE(Overlap<4, SimplexVsBbox>::apply(straddling, box));
  }

  SECTION("SimplexVsBbox predicate in 5D")
  {
    const AABBox<5> box{.center = {2.0, 2.0, 2.0, 2.0, 2.0},
                        .half_width = {1.0, 1.0, 1.0, 1.0, 1.0}};
    Omega_h::Matrix<5, 6> inside{};
    for (int v = 0; v < 6; ++v) {
      for (int d = 0; d < 5; ++d) {
        inside(d, v) = 2.0;
      }
    }
    REQUIRE(Overlap<5, SimplexVsBbox>::apply(inside, box));
  }
}
