#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <pcms/uniform_grid.h>
#include <Omega_h_build.hpp>
#include <Omega_h_library.hpp>

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
