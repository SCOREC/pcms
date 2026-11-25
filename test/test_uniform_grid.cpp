#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_library.hpp>
#include <pcms/localization/overlap.hpp>
#include <pcms/utility/uniform_grid.hpp>
#include "pcms/backends/omega_h/candidate_rank_table.hpp"
#include "pcms/backends/omega_h/omega_h.hpp"

using pcms::AABBox;
using pcms::CellRange;
using pcms::LO;
using pcms::OverlappingCells;
using pcms::Real;
using pcms::UniformGrid;
using pcms::UniformGrid2D;
using pcms::omega_h::CreateUniformGridFromMesh;

TEST_CASE("uniform grid")
{
  UniformGrid2D uniform_grid{
    .edge_length = {10, 12}, .bot_left = {0, 0}, .divisions = {10, 12}};
  REQUIRE(uniform_grid.GetNumCells() == 120);
  // Omega_h::Vector<2> point{0,0};
  SECTION("Closest Cell ID")
  {
    REQUIRE(0 == uniform_grid.ClosestCellID(Omega_h::Vector<2>{0, 0}));
    REQUIRE(0 == uniform_grid.ClosestCellID(Omega_h::Vector<2>{-5, -5}));
    REQUIRE(1 == uniform_grid.ClosestCellID(Omega_h::Vector<2>{1.5, 0}));
    REQUIRE(119 == uniform_grid.ClosestCellID(Omega_h::Vector<2>{10, 12}));
    REQUIRE(119 == uniform_grid.ClosestCellID(Omega_h::Vector<2>{9.5, 11.5}));
    REQUIRE(119 == uniform_grid.ClosestCellID(Omega_h::Vector<2>{100, 100}));
  }
  SECTION("cell bbox")
  {
    std::array bboxes{uniform_grid.GetCellBBOX(0), uniform_grid.GetCellBBOX(1),
                      uniform_grid.GetCellBBOX(119)};
    REQUIRE(bboxes[0].center[0] == Catch::Approx(0.5));
    REQUIRE(bboxes[0].center[1] == Catch::Approx(0.5));
    REQUIRE(bboxes[1].center[0] == Catch::Approx(1.5));
    REQUIRE(bboxes[1].center[1] == Catch::Approx(0.5));
    REQUIRE(bboxes[2].center[0] == Catch::Approx(9.5));
    REQUIRE(bboxes[2].center[1] == Catch::Approx(11.5));
    for (const auto& bbox : bboxes) {
      for (int i = 0; i < bbox.dim; ++i) {
        REQUIRE(bbox.half_width[i] == Catch::Approx(0.5));
      }
    }
  }
  SECTION("GetTwoDCellIndex")
  {
    auto [i, j] = uniform_grid.GetDimensionedIndex(0);
    REQUIRE(i == 0);
    REQUIRE(j == 0);
    auto [k, l] = uniform_grid.GetDimensionedIndex(119);
    REQUIRE(k == 11);
    REQUIRE(l == 9);
  }
  SECTION("GetCellIndex")
  {
    REQUIRE(0 == uniform_grid.GetCellIndex({0, 0}));
    REQUIRE(119 == uniform_grid.GetCellIndex({11, 9}));
  }
  SECTION("AxisIndicesFromCellIndex")
  {
    // The inverse of CellIndexFromAxisIndices, in coordinate order. Cell 1 is
    // the second cell along x, which the cell-BBOX test below also pins down
    // (its center is (1.5, 0.5)): the grid stores its indices reversed relative
    // to coordinates, so getting this wrong is silent.
    {
      const auto [i, j] = uniform_grid.AxisIndicesFromCellIndex(0);
      REQUIRE(i == 0);
      REQUIRE(j == 0);
    }
    {
      const auto [i, j] = uniform_grid.AxisIndicesFromCellIndex(1);
      REQUIRE(i == 1);
      REQUIRE(j == 0);
    }
    {
      const auto [i, j] = uniform_grid.AxisIndicesFromCellIndex(119);
      REQUIRE(i == 9);
      REQUIRE(j == 11);
    }

    // Round-trip over every cell of the grid, which is the property a caller
    // keyed by linear cell id depends on.
    for (pcms::LO cell = 0; cell < uniform_grid.GetNumCells(); ++cell) {
      const auto axis_indices = uniform_grid.AxisIndicesFromCellIndex(cell);
      REQUIRE(uniform_grid.CellIndexFromAxisIndices(axis_indices) == cell);
      // ...and it must agree with the binning rule used to find a point's cell.
      for (int d = 0; d < 2; ++d) {
        const auto center = uniform_grid.GetCellBBOX(cell).center[d];
        REQUIRE(uniform_grid.AxisCellIndex(d, center) == axis_indices[d]);
      }
    }
  }
}

TEST_CASE("Create uniform grid from omega_h mesh")
{
  // Initialize Omega_h library
  auto lib = Omega_h::Library{};
  auto world = lib.world();

  // Create a simple 2D box mesh: 1.0 x 1.0 domain with 4x4 elements
  auto mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1.0, 1.0, 0.0, 4, 4, 0, false);

  SECTION("Create grid with specified divisions per dimension")
  {
    // Create a 10x12 uniform grid over the mesh
    auto grid = pcms::omega_h::CreateUniformGridFromMesh<2>(mesh, {10, 12});

    // Verify the grid properties
    REQUIRE(grid.GetNumCells() == 120);
    REQUIRE(grid.edge_length[0] == Catch::Approx(1.0));
    REQUIRE(grid.edge_length[1] == Catch::Approx(1.0));
    REQUIRE(grid.bot_left[0] == Catch::Approx(0.0));
    REQUIRE(grid.bot_left[1] == Catch::Approx(0.0));
    REQUIRE(grid.divisions[0] == 10);
    REQUIRE(grid.divisions[1] == 12);

    // Test cell bounding boxes
    auto bbox_first = grid.GetCellBBOX(0);
    REQUIRE(bbox_first.center[0] == Catch::Approx(0.05)); // 1.0/10/2
    REQUIRE(bbox_first.center[1] == Catch::Approx(0.5 / 12.0));

    auto bbox_last = grid.GetCellBBOX(119);
    REQUIRE(bbox_last.center[0] == Catch::Approx(0.95));
    REQUIRE(bbox_last.center[1] == Catch::Approx(11.5 / 12.0));
  }

  SECTION("Create grid with equal divisions")
  {
    // Create a 15x15 uniform grid
    auto grid = pcms::omega_h::CreateUniformGridFromMesh<2>(mesh, 15);

    REQUIRE(grid.GetNumCells() == 225);
    REQUIRE(grid.divisions[0] == 15);
    REQUIRE(grid.divisions[1] == 15);

    // Verify all cells have equal size
    auto cell_width = grid.edge_length[0] / grid.divisions[0];
    auto cell_height = grid.edge_length[1] / grid.divisions[1];
    REQUIRE(cell_width == Catch::Approx(1.0 / 15.0));
    REQUIRE(cell_height == Catch::Approx(1.0 / 15.0));
  }

  SECTION("Test with non-unit domain")
  {
    // Create a larger domain: 5.0 x 3.0
    auto large_mesh = Omega_h::build_box(world, OMEGA_H_SIMPLEX, 5.0, 3.0, 0.0,
                                         10, 6, 0, false);

    auto grid =
      pcms::omega_h::CreateUniformGridFromMesh<2>(large_mesh, {20, 30});

    REQUIRE(grid.GetNumCells() == 600);
    REQUIRE(grid.edge_length[0] == Catch::Approx(5.0));
    REQUIRE(grid.edge_length[1] == Catch::Approx(3.0));
    REQUIRE(grid.bot_left[0] == Catch::Approx(0.0));
    REQUIRE(grid.bot_left[1] == Catch::Approx(0.0));

    // Test cell sizes
    auto bbox = grid.GetCellBBOX(0);
    REQUIRE(bbox.half_width[0] == Catch::Approx(5.0 / 20.0 / 2.0));
    REQUIRE(bbox.half_width[1] == Catch::Approx(3.0 / 30.0 / 2.0));
  }

  SECTION("Test point location in grid cells")
  {
    auto grid = pcms::omega_h::CreateUniformGridFromMesh<2>(mesh, {10, 10});

    // Test corner points
    REQUIRE(grid.ClosestCellID(Omega_h::Vector<2>{0.0, 0.0}) == 0);
    REQUIRE(grid.ClosestCellID(Omega_h::Vector<2>{1.0, 1.0}) == 99);

    // Test middle point
    REQUIRE(grid.ClosestCellID(Omega_h::Vector<2>{0.5, 0.5}) == 55);
  }
}

TEST_CASE("Create uniform grid from 3D mesh")
{
  auto lib = Omega_h::Library{};
  auto world = lib.world();

  // Create a 3D box mesh: 2.0 x 3.0 x 1.5 domain
  auto mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 2.0, 3.0, 1.5, 4, 6, 3, false);

  SECTION("Create 3D grid with specified divisions")
  {
    auto grid = pcms::omega_h::CreateUniformGridFromMesh<3>(mesh, {10, 15, 8});

    REQUIRE(grid.GetNumCells() == 1200); // 10 * 15 * 8
    REQUIRE(grid.edge_length[0] == Catch::Approx(2.0));
    REQUIRE(grid.edge_length[1] == Catch::Approx(3.0));
    REQUIRE(grid.edge_length[2] == Catch::Approx(1.5));
    REQUIRE(grid.divisions[0] == 10);
    REQUIRE(grid.divisions[1] == 15);
    REQUIRE(grid.divisions[2] == 8);
  }

  SECTION("Create 3D grid with equal divisions")
  {
    auto grid = pcms::omega_h::CreateUniformGridFromMesh<3>(mesh, 12);

    REQUIRE(grid.GetNumCells() == 1728); // 12^3
    REQUIRE(grid.divisions[0] == 12);
    REQUIRE(grid.divisions[1] == 12);
    REQUIRE(grid.divisions[2] == 12);
  }
}

TEST_CASE("UniformGrid dimensioned index with unequal 3D divisions")
{
  UniformGrid<3> grid;
  grid.bot_left = {0.0, 0.0, 0.0};
  grid.edge_length = {5.0, 3.0, 4.0};
  grid.divisions = {5, 3, 4};

  for (pcms::LO i = 0; i < grid.GetNumCells(); ++i) {
    REQUIRE(grid.GetCellIndex(grid.GetDimensionedIndex(i)) == i);
  }

  // Entries are in the reverse order of the coordinates: {z, y, x}.
  struct Case
  {
    pcms::LO cell;
    pcms::LO z, y, x;
  };
  const Case cases[] = {{0, 0, 0, 0},  {1, 0, 0, 1},  {5, 0, 1, 0},
                        {15, 1, 0, 0}, {59, 3, 2, 4}, {23, 1, 1, 3}};
  for (const auto& c : cases) {
    INFO("cell " << c.cell);
    const auto index = grid.GetDimensionedIndex(c.cell);
    REQUIRE(index[0] == c.z);
    REQUIRE(index[1] == c.y);
    REQUIRE(index[2] == c.x);
    const auto box = grid.GetCellBBOX(c.cell);
    REQUIRE(box.center[0] == Catch::Approx(c.x + 0.5));
    REQUIRE(box.center[1] == Catch::Approx(c.y + 0.5));
    REQUIRE(box.center[2] == Catch::Approx(c.z + 0.5));
  }
}

// Grid resolution has to be a property of the mesh as a whole, not of one
// partition: search cells are named by a linear id, and an id identifies
// the same region of space on two ranks only if both resolved the grid the
// same way. A resolution derived from `mesh.nelems()` breaks that, because
// `nelems()` counts ghost copies and therefore changes with the partition.
TEST_CASE("global grid resolution")
{
  auto lib = Omega_h::Library{};
  auto world = lib.world();
  auto mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1.0, 1.0, 0.0, 8, 8, 0, false);

  SECTION("the global resolution reproduces the local heuristic serially")
  {
    // Callers used to derive the resolution from mesh.nelems(). On an
    // undecomposed mesh that is already the global count, so switching to the
    // global helper must not change the grid they get.
    REQUIRE(pcms::omega_h::DivisionsPerAxisForGlobalMesh<2>(mesh) ==
            pcms::detail::DivisionsPerAxisForMesh<2>(mesh.nelems()));
    REQUIRE(pcms::omega_h::DivisionsPerAxisForGlobalMesh<3>(mesh) ==
            pcms::detail::DivisionsPerAxisForMesh<3>(mesh.nelems()));
  }

  SECTION("resolution still targets roughly one element per cell")
  {
    const auto global_nelems = mesh.nglobal_ents(mesh.dim());
    const auto n = pcms::omega_h::DivisionsPerAxisForGlobalMesh<2>(mesh);
    REQUIRE(n >= 10); // the min_divisions floor
    REQUIRE(static_cast<Omega_h::GO>(n) * n >= global_nelems);
  }
}

// Linear cell ids and per-axis indices must be inverses on any grid, including
// one whose axes have different division counts. A grid with equal divisions
// on every axis hides an axis mix-up, so this one has three different counts.
TEST_CASE("cell index round trip on a grid with unequal divisions")
{
  const UniformGrid<3> grid{.edge_length = {3.0, 4.0, 5.0},
                            .bot_left = {0.0, 0.0, 0.0},
                            .divisions = {3, 4, 5}};
  for (pcms::LO cell = 0; cell < grid.GetNumCells(); ++cell) {
    CAPTURE(cell);
    const auto axis_indices = grid.AxisIndicesFromCellIndex(cell);
    for (int d = 0; d < 3; ++d) {
      REQUIRE(axis_indices[d] >= 0);
      REQUIRE(axis_indices[d] < grid.divisions[d]);
    }
    REQUIRE(grid.CellIndexFromAxisIndices(axis_indices) == cell);
    // The cell's own box must bin back to the same cell.
    const auto box = grid.GetCellBBOX(cell);
    REQUIRE(grid.ClosestCellID(Omega_h::Vector<3>{box.center[0], box.center[1],
                                                  box.center[2]}) == cell);
  }
}

// Which grid cells a box covers (`OverlappingCells`),
// as used to build the point search's candidate map, and the invariants of
// `CellRange`.

TEST_CASE("bbox vs grid cell range")
{
  UniformGrid2D grid{
    .edge_length = {10, 10}, .bot_left = {0, 0}, .divisions = {10, 10}};

  SECTION("box inside one cell still covers padded neighbors")
  {
    AABBox<2> box{.center = {0.5, 0.5}, .half_width = {0.1, 0.1}};
    const auto range = OverlappingCells<2>(box, grid);
    // default pad = 1 clamps at the grid border
    REQUIRE(range.lo[0] == 0);
    REQUIRE(range.lo[1] == 0);
    REQUIRE(range.hi[0] == 1);
    REQUIRE(range.hi[1] == 1);
  }

  SECTION("no padding gives the exact span")
  {
    AABBox<2> box{.center = {2.5, 2.5}, .half_width = {0.4, 0.4}};
    const auto range = OverlappingCells<2>(box, grid, 0);
    REQUIRE(range.lo[0] == 2);
    REQUIRE(range.lo[1] == 2);
    REQUIRE(range.hi[0] == 2);
    REQUIRE(range.hi[1] == 2);
  }

  SECTION("box spanning several cells")
  {
    AABBox<2> box{.center = {2.5, 2.5}, .half_width = {1.0, 1.0}};
    const auto range = OverlappingCells<2>(box, grid, 0);
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
/// A shared boundary counts, which is the semantics the overlap predicate
/// has and therefore the set `OverlappingCells` has to cover.
[[nodiscard]] bool CellTouchesBox(const UniformGrid2D& grid, LO cx, LO cy,
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

/// The property `OverlappingCells` exists to guarantee: its range contains
/// every cell a closed-boundary predicate can accept. A range that is too
/// narrow drops candidates silently -- the element is never offered to the
/// query and nothing downstream can tell -- so this is checked as a
/// coverage property rather than against the exact bounds, which are
/// allowed to over-include.
TEST_CASE(
  "OverlappingCells covers every cell a closed-boundary predicate accepts")
{
  const UniformGrid2D grid{
    .edge_length = {10, 10}, .bot_left = {0, 0}, .divisions = {10, 10}};

  // Centres at, near and between grid lines, and extents that land exactly
  // on a grid line: that is when an unwidened range loses the touching
  // cell.
  for (const Real center : {0.5, 2.5, 5.0, 9.5}) {
    for (const Real half : {0.0, 0.25, 0.5, 1.0, 3.0}) {
      const AABBox<2> box{.center = {center, center},
                          .half_width = {half, half}};
      const auto range = OverlappingCells<2>(box, grid);
      for (LO cx = 0; cx < grid.divisions[0]; ++cx) {
        for (LO cy = 0; cy < grid.divisions[1]; ++cy) {
          if (CellTouchesBox(grid, cx, cy, box)) {
            REQUIRE(range.Contains(Kokkos::Array<LO, 2>{cx, cy}));
          }
        }
      }
    }
  }
}

TEST_CASE("OverlappingCells bins a supplied bound")
{
  const UniformGrid2D grid{
    .edge_length = {10, 10}, .bot_left = {0, 0}, .divisions = {10, 10}};

  // A degenerate triangle whose bounding box is the single grid line x = y
  // = 3: the cells on either side of that line only touch it, so the
  // padding must come from the predicate rather than from the caller.
  pcms::detail::SimplexVertices<2> degenerate;
  for (int k = 0; k < 3; ++k) {
    degenerate[k] = {3.0, 3.0};
  }
  // `OverlappingCells` bins a bound it is given; the element side supplies
  // it.
  const auto range =
    OverlappingCells<2>(pcms::detail::SimplexBBox<2>(degenerate), grid);
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

// The grid and overlap layer under the point search is written for any
// dimension, though the search classes themselves are 2D and 3D only. These
// checks pin that down above 3D. Every grid has a different division count
// on each axis, and every probe uses different indices per axis, so an axis
// mix-up cannot hide behind a symmetric case.
namespace
{

/// Cells of `grid` whose id, per-axis indices and box disagree, counted on
/// the device: the same functions the search runs in its kernels.
template <unsigned Dim>
pcms::LO CountCellRoundTripFailuresOnDevice(const pcms::UniformGrid<Dim>& grid)
{
  pcms::LO failures = 0;
  Kokkos::parallel_reduce(
    "cell round trip", grid.GetNumCells(),
    KOKKOS_LAMBDA(pcms::LO cell, pcms::LO & n) {
      const auto axis = grid.AxisIndicesFromCellIndex(cell);
      bool ok = grid.CellIndexFromAxisIndices(axis) == cell;
      Omega_h::Vector<Dim> center;
      const auto box = grid.GetCellBBOX(cell);
      for (unsigned d = 0; d < Dim; ++d) {
        ok = ok && axis[d] >= 0 && axis[d] < grid.divisions[d];
        center[d] = box.center[d];
      }
      ok = ok && grid.ClosestCellID(center) == cell;
      n += ok ? 0 : 1;
    },
    failures);
  return failures;
}

} // namespace

// What dimension does the overlap + grid-cell-range layer actually support?
TEST_CASE("uniform grid above 3D", "[nd]")
{
  SECTION("4D grid: cell ids, per-axis indices and cell boxes agree")
  {
    // Cell size 1 on every axis, so a cell's box is easy to state.
    const pcms::UniformGrid<4> grid{.edge_length = {2.0, 3.0, 4.0, 5.0},
                                    .bot_left = {0.0, 0.0, 0.0, 0.0},
                                    .divisions = {2, 3, 4, 5}};
    REQUIRE(grid.GetNumCells() == 120);

    // Ids number cells x-fastest: x + Nx * (y + Ny * (z + Nz * w)).
    const Kokkos::Array<LO, 4> last{1, 2, 3, 4};
    REQUIRE(grid.CellIndexFromAxisIndices(last) == 119);
    const Kokkos::Array<LO, 4> probe{1, 0, 2, 3};
    const LO probe_id = 1 + 2 * (0 + 3 * (2 + 4 * 3));
    REQUIRE(grid.CellIndexFromAxisIndices(probe) == probe_id);
    const auto back = grid.AxisIndicesFromCellIndex(probe_id);
    for (int d = 0; d < 4; ++d) {
      REQUIRE(back[d] == probe[d]);
    }
    const auto box = grid.GetCellBBOX(probe_id);
    for (int d = 0; d < 4; ++d) {
      REQUIRE(box.center[d] == Catch::Approx(probe[d] + 0.5));
      REQUIRE(box.half_width[d] == Catch::Approx(0.5));
    }

    REQUIRE(CountCellRoundTripFailuresOnDevice(grid) == 0);
  }

  SECTION("5D grid: every cell round-trips on the device")
  {
    const pcms::UniformGrid<5> grid{.edge_length = {1.0, 2.0, 3.0, 2.0, 1.0},
                                    .bot_left = {0.0, -1.0, 2.0, 0.5, 0.0},
                                    .divisions = {2, 3, 4, 3, 2}};
    REQUIRE(grid.GetNumCells() == 144);
    REQUIRE(CountCellRoundTripFailuresOnDevice(grid) == 0);
  }

  SECTION("4D candidate cells keep each axis separate")
  {
    const pcms::UniformGrid<4> grid{.edge_length = {2.0, 3.0, 4.0, 5.0},
                                    .bot_left = {0.0, 0.0, 0.0, 0.0},
                                    .divisions = {2, 3, 4, 5}};
    // Inside cell (0, 1, 2, 3); padded by one cell, clamped to the grid.
    const AABBox<4> box{.center = {0.5, 1.5, 2.5, 3.5},
                        .half_width = {0.25, 0.25, 0.25, 0.25}};
    const auto range = pcms::OverlappingCells<4>(box, grid);
    const Kokkos::Array<LO, 4> lo{0, 0, 1, 2};
    const Kokkos::Array<LO, 4> hi{1, 2, 3, 4};
    for (int d = 0; d < 4; ++d) {
      CAPTURE(d);
      REQUIRE(range.lo[d] == lo[d]);
      REQUIRE(range.hi[d] == hi[d]);
    }
    REQUIRE(range.Contains(Kokkos::Array<LO, 4>{0, 1, 2, 3}));
    REQUIRE_FALSE(range.Contains(Kokkos::Array<LO, 4>{0, 1, 2, 1}));
  }

  SECTION("1D")
  {
    const UniformGrid<1> grid{
      .edge_length = {4.0}, .bot_left = {0.0}, .divisions = {4}};
    const AABBox<1> box{.center = {1.5}, .half_width = {0.25}};
    // Cell 1, padded by 1, clamped at the low edge.
    const auto range = OverlappingCells<1>(box, grid);
    REQUIRE(range.lo[0] == 0);
    REQUIRE(range.hi[0] == 2);
    REQUIRE(range.Contains(Kokkos::Array<LO, 1>{1}));
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
}
