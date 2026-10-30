#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_library.hpp>
#include <pcms/utility/uniform_grid.h>

using pcms::CreateUniformGridFromMesh;
using pcms::Uniform2DGrid;
using pcms::UniformGrid;

TEST_CASE("uniform grid")
{
  Uniform2DGrid uniform_grid{
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
    auto grid = CreateUniformGridFromMesh<2>(mesh, {10, 12});

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
    auto grid = CreateUniformGridFromMesh<2>(mesh, 15);

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

    auto grid = CreateUniformGridFromMesh<2>(large_mesh, {20, 30});

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
    auto grid = CreateUniformGridFromMesh<2>(mesh, {10, 10});

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
    auto grid = CreateUniformGridFromMesh<3>(mesh, {10, 15, 8});

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
    auto grid = CreateUniformGridFromMesh<3>(mesh, 12);

    REQUIRE(grid.GetNumCells() == 1728); // 12^3
    REQUIRE(grid.divisions[0] == 12);
    REQUIRE(grid.divisions[1] == 12);
    REQUIRE(grid.divisions[2] == 12);
  }
}

// Grid resolution has to be a property of the mesh as a whole, not of one
// partition: search cells are named by a linear id, and an id identifies the
// same region of space on two ranks only if both resolved the grid the same
// way. A resolution derived from `mesh.nelems()` breaks that, because
// `nelems()` counts ghost copies and therefore changes with the partition.
TEST_CASE("global grid resolution")
{
  auto lib = Omega_h::Library{};
  auto world = lib.world();
  auto mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1.0, 1.0, 0.0, 8, 8, 0, false);

  SECTION("global element count reduces to the local count when nothing is "
          "ghosted")
  {
    REQUIRE(mesh.nghost_layers() == 0);
    REQUIRE(pcms::GlobalElementCount(mesh) == mesh.nelems());
  }

  SECTION("the global resolution reproduces the local heuristic serially")
  {
    // Callers used to derive the resolution from mesh.nelems(). On an
    // undecomposed mesh that is already the global count, so switching to the
    // global helper must not change the grid they get.
    REQUIRE(pcms::DivisionsPerAxisForGlobalMesh<2>(mesh) ==
            pcms::DivisionsPerAxisForMesh<2>(mesh.nelems()));
    REQUIRE(pcms::DivisionsPerAxisForGlobalMesh<3>(mesh) ==
            pcms::DivisionsPerAxisForMesh<3>(mesh.nelems()));
  }

  SECTION("resolution still targets roughly one element per cell")
  {
    const auto global_nelems = pcms::GlobalElementCount(mesh);
    const auto n = pcms::DivisionsPerAxisForGlobalMesh<2>(mesh);
    REQUIRE(n >= 10); // the min_divisions floor
    REQUIRE(n * n >= global_nelems);
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
