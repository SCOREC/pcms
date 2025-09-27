#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <pcms/coordinate.h>
#include <pcms/coordinate_systems.h>
#include "mdspan/mdspan.hpp"
#include <pcms/arrays.h>

using pcms::Cartesian;
using pcms::Coordinate;
using pcms::CoordinateElement;
using pcms::Cylindrical;
using pcms::Real;

#if __cplusplus >= 202002L
// coordinates from different coordinate systems should not be equality
// comparable
static_assert(!std::equality_comparable_with<Coordinate<Cartesian>,
                                             Coordinate<Cylindrical>>);
#endif

TEST_CASE("Coordinate Strong type works", "[coordinate]")
{
  Coordinate<Cartesian> cart{1.0, 2.0, 3.0};
  Coordinate<Cylindrical> cyl{4.0, 5.0, 6.0};
  SECTION("Values()")
  {
    auto values = cart.Values();
    REQUIRE(values[0] == Catch::Approx(1.0));
    REQUIRE(values[1] == Catch::Approx(2.0));
    REQUIRE(values[2] == Catch::Approx(3.0));
    values = cyl.Values();
    REQUIRE(values[0] == Catch::Approx(4.0));
    REQUIRE(values[1] == Catch::Approx(5.0));
    REQUIRE(values[2] == Catch::Approx(6.0));
  }
  SECTION("operator[]")
  {
    REQUIRE(cart[0] == Catch::Approx(1.0));
    REQUIRE(cart[1] == Catch::Approx(2.0));
    REQUIRE(cart[2] == Catch::Approx(3.0));
    REQUIRE(cyl[0] == Catch::Approx(4.0));
    REQUIRE(cyl[1] == Catch::Approx(5.0));
    REQUIRE(cyl[2] == Catch::Approx(6.0));
  }
}
/*
TEST_CASE( "Coordinate value type works", "[coordinate]" ) {
  using CartVal = CoordinateElement<Cartesian> ;
  namespace stdex = std::experimental;
  pcms::CoordinateMDArray<double,Cartesian> a(2);

  Kokkos::mdarray<
    double,
    Kokkos::extents<int, Kokkos::dynamic_extent, 3>,
    Kokkos::layout_right, Omega_h::Write<double> > bla(2);

  double count = 1.0;
  for(int i=0; i<a.extent(0); ++i) {
    for(int j=0; j<a.extent(1); ++j) {
      bla(i,j) = count;
      a(i,j) = CartVal{count++};
    }
  }
  REQUIRE(a(0,0).underlying() == Catch::Approx(1.0));
  REQUIRE(a(0,1).underlying() == Catch::Approx(2.0));
  REQUIRE(a(0,2).underlying() == Catch::Approx(3.0));
  REQUIRE(a(1,0).underlying() == Catch::Approx(4.0));
  REQUIRE(a(1,1).underlying() == Catch::Approx(5.0));
  REQUIRE(a(1,2).underlying() == Catch::Approx(6.0));

  REQUIRE(bla(0,0) == Catch::Approx(1.0));
  REQUIRE(bla(0,1) == Catch::Approx(2.0));
  REQUIRE(bla(0,2) == Catch::Approx(3.0));
  REQUIRE(bla(1,0) == Catch::Approx(4.0));
  REQUIRE(bla(1,1) == Catch::Approx(5.0));
  REQUIRE(bla(1,2) == Catch::Approx(6.0));
}
 */
