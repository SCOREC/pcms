#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <pcms/localization/shapes.hpp>
#include <pcms/localization/overlap.hpp>

#include <type_traits>

// Tests of the shapes and of the overlap operations on them: every operation is
// checked on its own, and each case of the exact triangle-box test is reached
// through a configuration only that case can detect.

using pcms::AABBox;
using pcms::Point;
using pcms::Real;
using pcms::Simplex;

namespace
{
Simplex<2> Triangle(Real x0, Real y0, Real x1, Real y1, Real x2, Real y2)
{
  Simplex<2> s;
  s.vertices[0] = {x0, y0};
  s.vertices[1] = {x1, y1};
  s.vertices[2] = {x2, y2};
  return s;
}

AABBox<2> Box(Real lo_x, Real hi_x, Real lo_y, Real hi_y)
{
  return {.center = {(lo_x + hi_x) / 2, (lo_y + hi_y) / 2},
          .half_width = {(hi_x - lo_x) / 2, (hi_y - lo_y) / 2}};
}
} // namespace

TEST_CASE("shape traits give the tag and the dimension")
{
  static_assert(std::is_same_v<pcms::tag_t<Point<2>>, pcms::PointTag>);
  static_assert(std::is_same_v<pcms::tag_t<AABBox<3>>, pcms::BoxTag>);
  static_assert(std::is_same_v<pcms::tag_t<Simplex<2>>, pcms::SimplexTag>);
  static_assert(
    std::is_same_v<pcms::tag_t<pcms::Segment<2>>, pcms::SegmentTag>);
  static_assert(pcms::dimension_v<Point<4>> == 4);
  static_assert(pcms::dimension_v<AABBox<3>> == 3);
  static_assert(pcms::dimension_v<Simplex<2>> == 2);
  SUCCEED();
}

TEST_CASE("point against box and simplex")
{
  const auto box = Box(0, 1, 0, 1);
  const auto tri = Triangle(0, 0, 1, 0, 0, 1);

  REQUIRE(pcms::overlaps(Point<2>{{0.5, 0.5}}, box));
  REQUIRE(pcms::overlaps(Point<2>{{1.0, 0.5}}, box)); // closed: on the edge
  REQUIRE_FALSE(pcms::overlaps(Point<2>{{1.5, 0.5}}, box));

  REQUIRE(pcms::overlaps(Point<2>{{0.2, 0.2}}, tri));
  REQUIRE_FALSE(pcms::overlaps(Point<2>{{0.8, 0.8}}, tri));

  // Either operand order.
  REQUIRE(pcms::overlaps(box, Point<2>{{0.5, 0.5}}));
  REQUIRE(pcms::overlaps(tri, Point<2>{{0.2, 0.2}}));
  REQUIRE_FALSE(pcms::overlaps(tri, Point<2>{{0.8, 0.8}}));
}

TEST_CASE("segment against box")
{
  using pcms::Segment;
  const auto box = Box(0, 1, 0, 1);
  // Crossing the box without either end inside it.
  REQUIRE(pcms::overlaps(Segment<2>{{-1.0, 0.5}, {2.0, 0.5}}, box));
  // Entirely inside.
  REQUIRE(pcms::overlaps(Segment<2>{{0.2, 0.2}, {0.8, 0.8}}, box));
  // Passing beside the box.
  REQUIRE_FALSE(pcms::overlaps(Segment<2>{{-1.0, 1.5}, {2.0, 1.5}}, box));
  // Closed: touching an edge counts.
  REQUIRE(pcms::overlaps(Segment<2>{{1.0, -1.0}, {1.0, 2.0}}, box));
  // Either operand order.
  REQUIRE(pcms::overlaps(box, Segment<2>{{-1.0, 0.5}, {2.0, 0.5}}));
  REQUIRE_FALSE(pcms::overlaps(box, Segment<2>{{-1.0, 1.5}, {2.0, 1.5}}));
}

TEST_CASE("exact triangle-box test, one case at a time")
{
  const auto tri = Triangle(0, 0, 4, 0, 0, 4);

  SECTION("a triangle vertex inside the box")
  {
    REQUIRE(pcms::overlaps(tri, Box(3.5, 5, -1, 0.5)));
  }
  SECTION("a box corner inside the triangle")
  {
    REQUIRE(pcms::overlaps(tri, Box(1, 2, 1, 2)));
  }
  SECTION("only an edge crosses the box")
  {
    // A vertical strip across the base: no triangle vertex is in it, and none
    // of its corners is in the triangle.
    const auto flat = Triangle(0, 0, 4, 0, 0, 1);
    REQUIRE(pcms::overlaps(flat, Box(1.9, 2.1, -1, 3)));
  }
  SECTION("disjoint although the bounding boxes overlap")
  {
    const auto corner = Box(3, 4, 3, 4);
    REQUIRE_FALSE(pcms::overlaps(tri, corner));
    // Their bounding boxes do overlap: in 2D the triangle itself is tested.
    REQUIRE(
      pcms::overlaps(pcms::detail::simplex_bbox<2>(tri.vertices), corner));
  }
  SECTION("either operand order")
  {
    REQUIRE(pcms::overlaps(Box(1, 2, 1, 2), tri));
    REQUIRE_FALSE(pcms::overlaps(Box(3, 4, 3, 4), tri));
  }
}

TEST_CASE("simplex-box overlap outside 2D tests the bounding box")
{
  // The tetrahedron (0,0,0), (1,0,0), (0,1,0), (0,0,1) and a box in its
  // bounding box's far corner, which the tetrahedron itself does not reach.
  Simplex<3> tet{};
  for (int k = 1; k < 4; ++k) {
    tet.vertices[k][k - 1] = 1.0;
  }
  const AABBox<3> far_corner{.center = {0.9, 0.9, 0.9},
                             .half_width = {0.05, 0.05, 0.05}};
  REQUIRE(pcms::overlaps(tet, far_corner));
  const AABBox<3> outside{.center = {2.0, 2.0, 2.0},
                          .half_width = {0.5, 0.5, 0.5}};
  REQUIRE_FALSE(pcms::overlaps(tet, outside));
}

TEST_CASE("simplex-box overlap above 3D", "[nd]")
{
  SECTION("4D simplex-vs-box overlap needs every axis to overlap")
  {
    const AABBox<4> box{.center = {1.0, 1.0, 1.0, 1.0},
                        .half_width = {0.5, 0.5, 0.5, 0.5}};
    // A simplex spanning [0.8, 1.2] on axes 0-2. On axis 3 it either reaches
    // into the box or stops short of it.
    const auto simplex = [](pcms::Real axis3_lo, pcms::Real axis3_hi) {
      Simplex<4> s{};
      for (int v = 0; v < 5; ++v) {
        for (int d = 0; d < 3; ++d) {
          s.vertices[v][d] = (v % 2 == 0) ? 0.8 : 1.2;
        }
        s.vertices[v][3] = (v == 0) ? axis3_lo : axis3_hi;
      }
      return s;
    };
    REQUIRE(pcms::overlaps(simplex(1.4, 2.0), box));
    // Overlapping on axes 0-2 but entirely past the box on axis 3.
    REQUIRE_FALSE(pcms::overlaps(simplex(1.6, 2.0), box));
  }
  SECTION("simplex-box overlap in 4D")
  {
    const AABBox<4> box{.center = {1.0, 1.0, 1.0, 1.0},
                        .half_width = {0.5, 0.5, 0.5, 0.5}};

    // A degenerate simplex: all 5 vertices at a point inside the box.
    Simplex<4> inside{};
    for (int v = 0; v < 5; ++v) {
      for (int d = 0; d < 4; ++d) {
        inside.vertices[v][d] = 1.0;
      }
    }
    REQUIRE(overlaps(inside, box));

    // Far outside on every axis.
    Simplex<4> outside{};
    for (int v = 0; v < 5; ++v) {
      for (int d = 0; d < 4; ++d) {
        outside.vertices[v][d] = 100.0;
      }
    }
    REQUIRE_FALSE(overlaps(outside, box));

    // Straddling the box on one axis only: the bounding-box test is true.
    Simplex<4> straddling{};
    for (int v = 0; v < 5; ++v) {
      for (int d = 0; d < 4; ++d) {
        straddling.vertices[v][d] = 1.0;
      }
    }
    straddling.vertices[0][3] = 0.5;
    straddling.vertices[1][3] = 1.9;
    REQUIRE(overlaps(straddling, box));
  }

  SECTION("simplex-box overlap in 5D")
  {
    const AABBox<5> box{.center = {2.0, 2.0, 2.0, 2.0, 2.0},
                        .half_width = {1.0, 1.0, 1.0, 1.0, 1.0}};
    Simplex<5> inside{};
    for (int v = 0; v < 6; ++v) {
      for (int d = 0; d < 5; ++d) {
        inside.vertices[v][d] = 2.0;
      }
    }
    REQUIRE(overlaps(inside, box));
  }
}

namespace
{
// m * inverse(m) is the identity: column j of it is m applied to column j of
// the inverse.
template <int Dim>
void CheckInverse(const pcms::detail::SquareMatrix<Dim>& m)
{
  const auto inverse = pcms::detail::invert(m);
  for (int j = 0; j < Dim; ++j) {
    const auto column = pcms::detail::multiply<Dim>(m, inverse[j]);
    for (int i = 0; i < Dim; ++i) {
      REQUIRE(column[i] == Catch::Approx(i == j ? 1.0 : 0.0).margin(1e-14));
    }
  }
}
} // namespace

TEST_CASE("invert gives the inverse of a 1x1 to 3x3 matrix")
{
  using pcms::detail::SquareMatrix;
  CheckInverse<1>(SquareMatrix<1>{{{4.0}}});
  CheckInverse<2>(SquareMatrix<2>{{{2.0, 1.0}, {-1.0, 3.0}}});
  CheckInverse<3>(
    SquareMatrix<3>{{{2.0, 0.5, 0.0}, {1.0, 3.0, -1.0}, {0.0, 1.0, 4.0}}});
}

TEST_CASE("barycentric_from_global")
{
  using pcms::detail::barycentric_from_global;

  SECTION("1D")
  {
    const pcms::detail::SimplexVertices<1> seg{{{{1.0}}, {{3.0}}}};
    const auto xi = barycentric_from_global<1>({2.5}, seg);
    REQUIRE(xi[0] == Catch::Approx(0.25));
    REQUIRE(xi[1] == Catch::Approx(0.75));
  }
  SECTION("2D")
  {
    const auto tri = Triangle(0, 0, 1, 0, 0.5, 1);
    // each vertex has weight 1 on itself and 0 on the others
    for (int k = 0; k < 3; ++k) {
      const auto xi = barycentric_from_global<2>(tri.vertices[k], tri.vertices);
      for (int j = 0; j < 3; ++j) {
        REQUIRE(xi[j] == Catch::Approx(j == k ? 1.0 : 0.0).margin(1e-14));
      }
    }
    const auto xi = barycentric_from_global<2>({0.5, 0.5}, tri.vertices);
    REQUIRE(xi[0] == Catch::Approx(0.25));
    REQUIRE(xi[1] == Catch::Approx(0.25));
    REQUIRE(xi[2] == Catch::Approx(0.5));
  }
  SECTION("3D")
  {
    Simplex<3> tet{};
    for (int k = 1; k < 4; ++k) {
      tet.vertices[k][k - 1] = 1.0;
    }
    // the centroid has equal weights
    const auto xi =
      barycentric_from_global<3>({0.25, 0.25, 0.25}, tet.vertices);
    for (int k = 0; k < 4; ++k) {
      REQUIRE(xi[k] == Catch::Approx(0.25));
    }
    // on the unit tetrahedron the weights of vertices 1-3 are the coordinates
    const auto p = barycentric_from_global<3>({0.1, 0.2, 0.3}, tet.vertices);
    REQUIRE(p[0] == Catch::Approx(0.4));
    REQUIRE(p[1] == Catch::Approx(0.1));
    REQUIRE(p[2] == Catch::Approx(0.2));
    REQUIRE(p[3] == Catch::Approx(0.3));
  }
}

TEST_CASE("is_barycentric_inside")
{
  using pcms::detail::is_barycentric_inside;
  using Xi = Kokkos::Array<Real, 3>;

  REQUIRE(is_barycentric_inside(Xi{0.2, 0.3, 0.5}));
  // closed: a point on an edge or at a vertex is inside
  REQUIRE(is_barycentric_inside(Xi{0.0, 0.5, 0.5}));
  REQUIRE(is_barycentric_inside(Xi{1.0, 0.0, 0.0}));
  REQUIRE_FALSE(is_barycentric_inside(Xi{-0.1, 0.6, 0.5}));
  REQUIRE_FALSE(is_barycentric_inside(Xi{1.1, -0.05, -0.05}));
  // the fuzz widens the band on both sides
  REQUIRE(is_barycentric_inside(Xi{-0.1, 0.6, 0.5}, 0.2));
  REQUIRE_FALSE(is_barycentric_inside(Xi{-0.3, 0.8, 0.5}, 0.2));
}
