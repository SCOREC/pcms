#include <catch2/catch.hpp>
#include <wdmcpl/coordinate.h>
#include <wdmcpl/coordinate_systems.h>
#include <wdmcpl/external/mdspan.hpp>
#include <wdmcpl/arrays.h>
#include <Omega_h_array.hpp>

using wdmcpl::Cartesian;
using wdmcpl::Coordinate;
using wdmcpl::Cartesian;
using wdmcpl::Cylindrical;
using wdmcpl::CoordinateElement;
using wdmcpl::Real;

// Version required for equality_comparable
#if __cplusplus >= 202002L
// coordinates from different coordinate systems should not be equality comparable
static_assert(!std::equality_comparable_with<Coordinate<Cartesian>,Coordinate<Cylindrical>>);
static_assert(!std::equality_comparable_with<CoordinateElement<Cartesian>,
#endif

// CoordinateValueType should be memcpyable this allows ("bitcast") between
// underlying type when needed
static_assert(std::is_trivially_copyable_v<CoordinateElement<Cartesian,Real>>);
static_assert(std::is_trivially_copyable_v<CoordinateElement<Cylindrical,Real>>);
static_assert(sizeof(Real) == sizeof(CoordinateElement<Cartesian,Real>));
static_assert(sizeof(Real) == sizeof(CoordinateElement<Cylindrical,Real>));

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
/*
TEST_CASE( "Coordinate value type works", "[coordinate]" ) {
  using CartVal = CoordinateElement<Cartesian> ;
  namespace stdex = std::experimental;
  wdmcpl::CoordinateMDArray<double,Cartesian> a(2);

  std::experimental::mdarray<
    double,
    std::experimental::extents<int, std::experimental::dynamic_extent, 3>,
    std::experimental::layout_right, Omega_h::Write<double> > bla(2);

  double count = 1.0;
  for(int i=0; i<a.extent(0); ++i) {
    for(int j=0; j<a.extent(1); ++j) {
      bla(i,j) = count;
      a(i,j) = CartVal{count++};
    }
  }
  REQUIRE(a(0,0).underlying() == Approx(1.0));
  REQUIRE(a(0,1).underlying() == Approx(2.0));
  REQUIRE(a(0,2).underlying() == Approx(3.0));
  REQUIRE(a(1,0).underlying() == Approx(4.0));
  REQUIRE(a(1,1).underlying() == Approx(5.0));
  REQUIRE(a(1,2).underlying() == Approx(6.0));

  REQUIRE(bla(0,0) == Approx(1.0));
  REQUIRE(bla(0,1) == Approx(2.0));
  REQUIRE(bla(0,2) == Approx(3.0));
  REQUIRE(bla(1,0) == Approx(4.0));
  REQUIRE(bla(1,1) == Approx(5.0));
  REQUIRE(bla(1,2) == Approx(6.0));
}
 */
