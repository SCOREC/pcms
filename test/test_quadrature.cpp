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
