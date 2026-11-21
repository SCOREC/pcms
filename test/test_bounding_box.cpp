#include <catch2/catch_test_macros.hpp>
#include <pcms/localization/overlap.hpp>

using pcms::AABBox;
using pcms::overlaps;
#include <iostream>

TEST_CASE("test intersection")
{
  const AABBox<2> unit_square{.center = {0, 0}, .half_width = {0.5, 0.5}};
  REQUIRE(overlaps(unit_square, unit_square));
  AABBox<2> shifted_unit_square{.center = {0, 0.4}, .half_width = {0.5, 0.5}};
  REQUIRE(overlaps(unit_square, shifted_unit_square));
  shifted_unit_square.center[1] *= 2;
  REQUIRE(overlaps(unit_square, shifted_unit_square));
  shifted_unit_square.center[1] *= 2;
  REQUIRE(!overlaps(unit_square, shifted_unit_square));
  // test on negative side
  shifted_unit_square.center[1] *= -1;
  REQUIRE(!overlaps(unit_square, shifted_unit_square));
  shifted_unit_square.center[1] /= 2;
  REQUIRE(overlaps(unit_square, shifted_unit_square));
}
