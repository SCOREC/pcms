#include <catch2/catch_test_macros.hpp>
#include <wdmcpl/bounding_box.h>

using wdmcpl::AABBox;
using wdmcpl::intersects;
#include <iostream>

TEST_CASE("test intersection") {
  const AABBox<2> unit_square{.center={0,0}, .half_width={0.5,0.5}};
  REQUIRE(intersects(unit_square,unit_square));
  AABBox<2> shifted_unit_square{.center={0,0.4}, .half_width={0.5,0.5}};
  REQUIRE(intersects(unit_square,shifted_unit_square));
  shifted_unit_square.center[1] *= 2;
  REQUIRE(intersects(unit_square,shifted_unit_square));
  shifted_unit_square.center[1] *= 2;
  REQUIRE(!intersects(unit_square,shifted_unit_square));
  // test on negative side
  shifted_unit_square.center[1] *= -1;
  REQUIRE(!intersects(unit_square,shifted_unit_square));
  shifted_unit_square.center[1] /= 2;
  REQUIRE(intersects(unit_square,shifted_unit_square));

}