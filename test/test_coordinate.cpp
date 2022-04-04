#include <catch2/catch.hpp>
#include <wdmcpl/coordinate.h>
#include <wdmcpl/coordinate_systems.h>

using wdmcpl::Coordinate;
using wdmcpl::Cartesian;
using wdmcpl::Cylindrical;

// coordinates from different coordinate systems should not be equality comparable
static_assert(!std::equality_comparable_with<Coordinate<Cartesian>,Coordinate<Cylindrical>>);

TEST_CASE( "Coordinate Strong type works", "[coordinate]" ) {
  Coordinate<Cartesian> cart{1.0,2.0,3.0};
  Coordinate<Cylindrical> cyl{4.0,5.0,6.0};
  SECTION("Values()")
  {
    auto values = cart.Values();
    REQUIRE(values[0] == Approx(1.0));
    REQUIRE(values[1] == Approx(2.0));
    REQUIRE(values[2] == Approx(3.0));
    values = cyl.Values();
    REQUIRE(values[0] == Approx(4.0));
    REQUIRE(values[1] == Approx(5.0));
    REQUIRE(values[2] == Approx(6.0));
  }
  SECTION("operator[]"){
    REQUIRE(cart[0] == Approx(1.0));
    REQUIRE(cart[1] == Approx(2.0));
    REQUIRE(cart[2] == Approx(3.0));
    REQUIRE(cyl[0] == Approx(4.0));
    REQUIRE(cyl[1] == Approx(5.0));
    REQUIRE(cyl[2] == Approx(6.0));
  }
}