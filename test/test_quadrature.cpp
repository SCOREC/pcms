#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include <pcms/field/function_spaces.hpp>
#include <pcms/geometry/quadrature.hpp>
#include <pcms/backends/omega_h/omega_h.hpp>
#include <pcms/backends/omega_h/omega_h_geometry.hpp>
#include <pcms/transfer/omega_h_intersection_quadrature.hpp>
#include <pcms/utility/coordinate_systems/cylindrical.hpp>
#include "field_test_utils.hpp"
#include <cmath>
#include <functional>

namespace
{

// sum_p w_p f(x_p) over a quadrature, on the host
double Integrate(const pcms::Quadrature& quadrature,
                 const std::function<double(const double*)>& f)
{
  const auto points = pcms::test::CopyCoordinatesToHost(
    quadrature.points->LocatedPoints().GetValues());
  const auto weights = Kokkos::create_mirror_view_and_copy(
    pcms::HostMemorySpace{}, quadrature.weights);
  REQUIRE(points.extent(0) == weights.extent(0));
  double sum = 0;
  for (std::size_t p = 0; p < points.extent(0); ++p) {
    double x[3] = {points(p, 0), points(p, 1),
                   points.extent(1) > 2 ? points(p, 2) : 0.0};
    sum += weights(p) * f(x);
  }
  return sum;
}

double WeightSum(const pcms::Quadrature& quadrature)
{
  return Integrate(quadrature, [](const double*) { return 1.0; });
}

// The integral of x^a y^b z^c over the unit square or cube
double ExactMonomial(int a, int b, int c, int dim)
{
  const double exact = 1.0 / ((a + 1) * (b + 1));
  return dim == 3 ? exact / (c + 1) : exact;
}

// Quadratures over an Omega_h geometry
template <typename Make>
void CheckPolynomialExactness(Omega_h::Mesh& mesh, Make make)
{
  const int dim = mesh.dim();
  auto geometry =
    pcms::omega_h::CreateGeometry(mesh, pcms::detail::DeferredCartesian());
  const int max_order = GENERATE(1, 2);
  CAPTURE(dim, max_order);
  const auto quadrature = make(geometry, max_order);
  CHECK(quadrature.points->GetGeometry() == geometry);
  CHECK(quadrature.points->MissingPoints().extent(0) == 0);
  CHECK(WeightSum(quadrature) == Catch::Approx(1.0).margin(1e-12));
  for (int a = 0; a <= max_order; ++a) {
    for (int b = 0; a + b <= max_order; ++b) {
      for (int c = 0; a + b + c <= max_order && (dim == 3 || c == 0); ++c) {
        CAPTURE(a, b, c);
        const double integral = Integrate(quadrature, [&](const double* x) {
          return std::pow(x[0], a) * std::pow(x[1], b) * std::pow(x[2], c);
        });
        CHECK(integral ==
              Catch::Approx(ExactMonomial(a, b, c, dim)).margin(1e-12));
      }
    }
  }
}

} // namespace

TEST_CASE("CreateQuadrature: the reference rule integrates polynomials up to "
          "its order exactly (2D)",
          "[quadrature]")
{
  Omega_h::Library lib;
  auto mesh = pcms::test::BuildSquareMesh(lib, 3);
  CheckPolynomialExactness(mesh, [](const auto& geometry, int order) {
    return pcms::CreateGaussQuadrature(geometry, order);
  });
}

TEST_CASE("CreateQuadrature: the reference rule integrates polynomials up to "
          "its order exactly (3D)",
          "[quadrature]")
{
  Omega_h::Library lib;
  auto mesh = pcms::test::BuildUnitCube(lib, 2);
  CheckPolynomialExactness(mesh, [](const auto& geometry, int order) {
    return pcms::CreateGaussQuadrature(geometry, order);
  });
}

TEST_CASE("CreateQuadrature: Monte Carlo weights sum to the area and the "
          "estimate converges; the seed fixes the points",
          "[quadrature]")
{
  Omega_h::Library lib;
  auto mesh = pcms::test::BuildSquareMesh(lib, 3);
  auto geometry =
    pcms::omega_h::CreateGeometry(mesh, pcms::detail::DeferredCartesian());
  const auto square = [](const double* x) { return x[0] * x[0]; };
  double previous_error = 1.0;
  for (const int samples : {4, 64, 1024}) {
    const auto quadrature = pcms::CreateMonteCarloQuadrature(geometry, samples);
    CAPTURE(samples);
    REQUIRE(quadrature.points->NumPoints() == mesh.nelems() * samples);
    CHECK(quadrature.points->MissingPoints().extent(0) == 0);
    CHECK(WeightSum(quadrature) == Catch::Approx(1.0).margin(1e-12));
    const double error = std::abs(Integrate(quadrature, square) - 1.0 / 3.0);
    CAPTURE(error);
    previous_error = error;
  }
  CHECK(previous_error < 0.01);

  const auto first = pcms::CreateMonteCarloQuadrature(geometry, 8, 7);
  const auto second = pcms::CreateMonteCarloQuadrature(geometry, 8, 7);
  CHECK(Integrate(first, square) == Integrate(second, square));
  const auto third = pcms::CreateMonteCarloQuadrature(geometry, 8, 8);
  CHECK(Integrate(first, square) != Integrate(third, square));
  // the rule of a seed does not change between builds
  CHECK(Integrate(first, square) == 0.32930288482711945);
}

TEST_CASE("CreateQuadrature: unsupported rules throw", "[quadrature]")
{
  Omega_h::Library lib;
  auto mesh = pcms::test::BuildSquareMesh(lib, 1);
  auto geometry =
    pcms::omega_h::CreateGeometry(mesh, pcms::detail::DeferredCartesian());
  CHECK_THROWS_AS(pcms::CreateGaussQuadrature(geometry, 0), pcms::pcms_error);
  // MeshFields' simplex rules go up to order 2
  CHECK_THROWS_AS(pcms::CreateGaussQuadrature(geometry, 3), pcms::pcms_error);
}

TEST_CASE("CreateIntersectionQuadrature: weights sum to the domain, the rule "
          "is exact on the pieces, and both meshes are bound without a miss",
          "[quadrature][intersection]")
{
  Omega_h::Library lib;
  const int dim = GENERATE(2, 3);
  auto source_mesh = dim == 2 ? pcms::test::BuildSquareMesh(lib, 3)
                              : pcms::test::BuildUnitCube(lib, 2);
  auto target_mesh = dim == 2 ? pcms::test::BuildUnitSquare(lib, 0)
                              : pcms::test::BuildUnitCube(lib, 1);
  const auto source =
    std::dynamic_pointer_cast<const pcms::omega_h::OmegaHGeometry>(
      pcms::omega_h::CreateGeometry(source_mesh,
                                    pcms::detail::DeferredCartesian()));
  const auto target =
    std::dynamic_pointer_cast<const pcms::omega_h::OmegaHGeometry>(
      pcms::omega_h::CreateGeometry(target_mesh,
                                    pcms::detail::DeferredCartesian()));
  const int order = 2;
  const auto quadrature =
    pcms::omega_h::CreateIntersectionQuadrature(source, target, order);
  CHECK(quadrature.points->GetGeometry() == target);
  CHECK(quadrature.points->MissingPoints().extent(0) == 0);
  REQUIRE(quadrature.also.size() == 1);
  CHECK(quadrature.also[0]->GetGeometry() == source);
  CHECK(quadrature.also[0]->MissingPoints().extent(0) == 0);
  CHECK(WeightSum(quadrature) == Catch::Approx(1.0).margin(1e-12));
  for (int a = 0; a <= order; ++a) {
    for (int b = 0; a + b <= order; ++b) {
      const int c = dim == 3 ? std::max(0, order - a - b) : 0;
      const double integral = Integrate(quadrature, [&](const double* x) {
        return std::pow(x[0], a) * std::pow(x[1], b) * std::pow(x[2], c);
      });
      CAPTURE(a, b, c);
      CHECK(integral ==
            Catch::Approx(ExactMonomial(a, b, c, dim)).margin(1e-12));
    }
  }
}

namespace
{

std::shared_ptr<const pcms::omega_h::OmegaHGeometry> CartesianGeometry(
  Omega_h::Mesh& mesh)
{
  return std::dynamic_pointer_cast<const pcms::omega_h::OmegaHGeometry>(
    pcms::omega_h::CreateGeometry(mesh, pcms::detail::DeferredCartesian()));
}

} // namespace

TEST_CASE("CreateIntersectionQuadrature: a tet inside a tet",
          "[quadrature][3d]")
{
  // Source is the reference tet (volume 1/6). Target is the similar tet at the
  // (1,0,0) corner, cut by the plane x=0.5, which sits entirely inside the
  // source (volume (1/6)*(1/2)^3 = 1/48).
  Omega_h::Library lib;
  auto source_mesh = pcms::test::BuildReferenceTet(lib);
  auto target_mesh =
    pcms::test::BuildTet(lib, Omega_h::Reals({1.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.5,
                                              0.5, 0.0, 0.5, 0.0, 0.5}));
  const auto quadrature = pcms::omega_h::CreateIntersectionQuadrature(
    CartesianGeometry(source_mesh), CartesianGeometry(target_mesh), 2);
  CHECK(WeightSum(quadrature) == Catch::Approx(1.0 / 48.0).margin(1e-14));
  const auto points = pcms::test::CopyCoordinatesToHost(
    quadrature.points->LocatedPoints().GetValues());
  REQUIRE(points.extent(0) > 0);
  for (std::size_t i = 0; i < points.extent(0); ++i) {
    CHECK(points(i, 0) >= 0.5 - 1e-12);
    CHECK(points(i, 1) >= -1e-12);
    CHECK(points(i, 2) >= -1e-12);
    CHECK(points(i, 0) + points(i, 1) + points(i, 2) <= 1.0 + 1e-12);
  }
}

TEST_CASE("CreateIntersectionQuadrature: dual tets overlap in an octahedron",
          "[quadrature][3d]")
{
  // Dual tets in the cube [-1,1]^3. Intersection = {|x|+|y|+|z| <= 1},
  // volume 4/3.
  Omega_h::Library lib;
  auto source_mesh = pcms::test::BuildTet(
    lib, Omega_h::Reals(
           {1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 1.0, -1.0, 1.0, -1.0}));
  auto target_mesh =
    pcms::test::BuildTet(lib, Omega_h::Reals({-1.0, -1.0, -1.0, -1.0, 1.0, 1.0,
                                              1.0, -1.0, 1.0, 1.0, 1.0, -1.0}));
  const auto quadrature = pcms::omega_h::CreateIntersectionQuadrature(
    CartesianGeometry(source_mesh), CartesianGeometry(target_mesh), 2);
  CHECK(WeightSum(quadrature) == Catch::Approx(4.0 / 3.0).margin(1e-10));
  CHECK(Integrate(quadrature, [](const double* x) {
          return x[0] + x[1] + x[2];
        }) == Catch::Approx(0.0).margin(1e-10));
  const auto points = pcms::test::CopyCoordinatesToHost(
    quadrature.points->LocatedPoints().GetValues());
  for (std::size_t i = 0; i < points.extent(0); ++i) {
    CHECK(std::abs(points(i, 0)) + std::abs(points(i, 1)) +
            std::abs(points(i, 2)) <=
          1.0 + 1e-10);
  }
}

TEST_CASE("quadratures integrate in the physical measure of a cylindrical "
          "system",
          "[quadrature]")
{
  Omega_h::Library lib;
  auto mesh = pcms::test::BuildUnitSquare(lib, 1);
  // the (R, Z) plane over [0, 1]^2 shifted to R in [1, 2]: the volume
  // element per radian is R dR dZ
  auto coords = mesh.coords();
  Omega_h::Write<pcms::Real> shifted(coords.size());
  Omega_h::parallel_for(
    mesh.nverts(), OMEGA_H_LAMBDA(Omega_h::LO v) {
      shifted[2 * v] = coords[2 * v] + 1.0;
      shifted[2 * v + 1] = coords[2 * v + 1];
    });
  mesh.set_coords(Omega_h::Reals(shifted));
  auto cylindrical =
    pcms::omega_h::CreateGeometry(mesh, pcms::csys::CreateCylindricalRZ());
  const auto quadrature = pcms::CreateGaussQuadrature(cylindrical, 2);
  const auto weights = Kokkos::create_mirror_view_and_copy(
    pcms::HostMemorySpace(), quadrature.weights);
  const auto points = pcms::test::CopyCoordinatesToHost(
    quadrature.points->LocatedPoints().GetValues());
  pcms::Real volume = 0;
  pcms::Real moment = 0;
  for (std::size_t i = 0; i < weights.size(); ++i) {
    volume += weights[i];
    moment += weights[i] * points(i, 0);
  }
  // int R dR dZ = 3/2, int R^2 dR dZ = 7/3
  CHECK_THAT(volume, Catch::Matchers::WithinAbs(1.5, pcms::test::ExactTol));
  CHECK_THAT(moment,
             Catch::Matchers::WithinAbs(7.0 / 3.0, pcms::test::ExactTol));
}
