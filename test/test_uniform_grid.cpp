#include <catch2/catch.hpp>
#include <wdmcpl/uniform_grid.h>

using wdmcpl::UniformGrid;

TEST_CASE("uniform grid")
{
  UniformGrid uniform_grid{
    .edge_length = {10, 12}, .bot_left = {0, 0}, .divisions = {10, 12}};
  REQUIRE(uniform_grid.GetNumCells() == 120);
  // Omega_h::Vector<2> point{0,0};
  SECTION("Closest Cell ID")
  {
    REQUIRE(0 == uniform_grid.ClosestCellID({0, 0}));
    REQUIRE(0 == uniform_grid.ClosestCellID({-5, -5}));
    REQUIRE(1 == uniform_grid.ClosestCellID({1.5, 0}));
    REQUIRE(119 == uniform_grid.ClosestCellID({10, 12}));
    REQUIRE(119 == uniform_grid.ClosestCellID({9.5, 11.5}));
    REQUIRE(119 == uniform_grid.ClosestCellID({100, 100}));
  }
  SECTION("cell bbox")
  {
    std::array bboxes{uniform_grid.GetCellBBOX(0), uniform_grid.GetCellBBOX(1),
                      uniform_grid.GetCellBBOX(119)};
    REQUIRE(bboxes[0].center[0] == Approx(0.5));
    REQUIRE(bboxes[0].center[1] == Approx(0.5));
    REQUIRE(bboxes[1].center[0] == Approx(1.5));
    REQUIRE(bboxes[1].center[1] == Approx(0.5));
    REQUIRE(bboxes[2].center[0] == Approx(9.5));
    REQUIRE(bboxes[2].center[1] == Approx(11.5));
    for (const auto& bbox : bboxes) {
      for (int i = 0; i < bbox.dim; ++i) {
        REQUIRE(bbox.half_width[i] == Approx(0.5));
      }
    }
  }
  SECTION("GetTwoDCellIndex") {
    auto [i,j] = uniform_grid.GetTwoDCellIndex(0);
    REQUIRE(i == 0);
    REQUIRE(j == 0);
    auto [k,l] = uniform_grid.GetTwoDCellIndex(119);
    REQUIRE(k == 11);
    REQUIRE(l == 9);
  }
  SECTION("GetCellIndex") {
    REQUIRE(0 == uniform_grid.GetCellIndex(0,0));
    REQUIRE(119 == uniform_grid.GetCellIndex(11,9));
  }
}