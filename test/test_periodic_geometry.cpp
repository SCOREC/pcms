#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <Kokkos_Core.hpp>
#include <pcms/utility/arrays.hpp>
#include <pcms/utility/detail/arrays.hpp>
#include <pcms/utility/detail/periodic.hpp>
#include <pcms/utility/coordinate_systems/cartesian.hpp>
#include "field_test_utils.hpp"
#include "periodic_test_utils.hpp"
#include "pcms/backends/uniform_grid/uniform_grid_geometry.hpp"
#include <cmath>
#include <numbers>

namespace
{

using pcms::Real;
constexpr Real TwoPi = 2 * std::numbers::pi_v<Real>;

// sum_p w_p f(x_p)
template <typename F>
Real Integrate(const pcms::Quadrature& quadrature, F f)
{
  const auto points = pcms::test::CopyCoordinatesToHost(
    quadrature.points->LocatedPoints().GetValues());
  const auto weights = Kokkos::create_mirror_view_and_copy(
    pcms::HostMemorySpace{}, quadrature.weights);
  Real sum = 0;
  for (std::size_t p = 0; p < points.extent(0); ++p) {
    sum += weights(p) * f(points(p, 0));
  }
  return sum;
}

} // namespace

TEST_CASE("Periodic geometry: points wrap modulo the period at Bind",
          "[periodic]")
{
  const auto geometry = pcms::test::PeriodicGrid(4, 2.0, 1.0);
  Kokkos::View<Real**, pcms::HostMemorySpace> host("p", 5, 1);
  const Real x[5] = {1.0, 1.6, 0.9, 7.25, -0.4};
  for (int i = 0; i < 5; ++i) {
    host(i, 0) = x[i];
  }
  const auto points = pcms::detail::ToDeviceView(host);
  const auto bound =
    geometry->Bind(pcms::CoordinateView<pcms::DeviceMemorySpace>(
      geometry->GetCoordinateSystem(), pcms::MakeConstRank2View(points)));
  CHECK(bound->MissingPoints().extent(0) == 0);
  const auto& binding =
    dynamic_cast<const pcms::detail::UniformGridBinding<1>&>(*bound);
  const auto located =
    pcms::test::CopyCoordinatesToHost(bound->LocatedPoints().GetValues());
  const auto cell = Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace{},
                                                        binding.CellIndices());
  const auto xi = Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace{},
                                                      binding.Parametric());
  // [1, 3), h = 0.5
  const int expected_cell[5] = {0, 1, 3, 0, 1};
  const Real expected_xi[5] = {0.0, 0.2, 0.8, 0.5, 0.2};
  const Real expected_phase[5] = {0.0, 0.3, 0.95, 0.125, 0.3};
  for (int i = 0; i < 5; ++i) {
    CHECK(cell(i) == expected_cell[i]);
    CHECK(xi(i, 0) == Catch::Approx(expected_xi[i]).margin(1e-13));
    CHECK(located(i, 0) ==
          Catch::Approx(1.0 + 2.0 * expected_phase[i]).margin(1e-13));
  }
}

TEST_CASE("Periodic geometry: points at a multiple of the period reduce to the "
          "lower end",
          "[periodic]")
{
  for (const auto& [period, x] :
       {std::pair{TwoPi, 33 * TwoPi},
        std::pair{TwoPi, std::nextafter(33 * TwoPi, 0.0)},
        std::pair{0.1, 3.9}}) {
    CAPTURE(period, x);
    const auto geometry = pcms::test::PeriodicGrid(4, period);
    const auto points = pcms::test::DevicePoints({{x}});
    const auto bound =
      geometry->Bind(pcms::CoordinateView<pcms::DeviceMemorySpace>(
        geometry->GetCoordinateSystem(), pcms::MakeConstRank2View(points)));
    const auto located =
      pcms::test::CopyCoordinatesToHost(bound->LocatedPoints().GetValues());
    CHECK(located(0, 0) >= 0);
    CHECK(located(0, 0) < period);
  }
}

TEST_CASE("Periodic helpers: nearest images and reduction to the period",
          "[periodic]")
{
  using pcms::detail::NearestImage;
  using pcms::detail::WrapToPeriod;
  const Real length = 2.5;
  CHECK(NearestImage(length / 2, length) == -length / 2);
  CHECK(NearestImage(-length / 2, length) == length / 2);
  CHECK(NearestImage(0.3 + 4 * length, length) ==
        Catch::Approx(0.3).margin(1e-14));
  CHECK(NearestImage(7.0, 0.0) == 7.0);
  Real wraps = 0;
  CHECK(WrapToPeriod(1.0 + length, 1.0, length, wraps) == 1.0);
  CHECK(wraps == 1.0);
  CHECK(WrapToPeriod(0.5, 1.0, length, wraps) ==
        Catch::Approx(0.5 + length).margin(1e-14));
  CHECK(wraps == -1.0);
  for (const Real x : {33 * TwoPi, std::nextafter(33 * TwoPi, 0.0)}) {
    CAPTURE(x);
    const Real y = WrapToPeriod(x, 0.0, TwoPi, wraps);
    CHECK(y >= 0.0);
    CHECK(y < TwoPi);
  }
}

TEST_CASE("Periodic geometry: a grid outside the fundamental interval of its "
          "periodic axis is refused",
          "[periodic]")
{
  pcms::UniformGrid<1> grid;
  grid.bot_left = {-TwoPi / 2};
  grid.edge_length = {TwoPi};
  grid.divisions = {4};
  CHECK_THROWS_AS(pcms::uniform_grid::CreateGeometry(
                    grid, pcms::csys::CreateCartesian(1, {{0, 0.0, TwoPi}})),
                  pcms::pcms_error);
}

TEST_CASE("Periodic geometry: the geometric factors are its system's",
          "[periodic]")
{
  const auto geometry = pcms::test::PeriodicGrid(5, 2.0);
  const auto quadrature = pcms::CreateGaussQuadrature(geometry, 2);
  const auto jacobian = quadrature.points->Jacobian();
  const auto volume = quadrature.points->VolumeFactor();
  REQUIRE(jacobian.extent(0) == quadrature.points->NumPoints());
  Kokkos::View<Real*, pcms::DeviceMemorySpace> out("out", 2);
  Kokkos::parallel_for(
    "copy", Kokkos::RangePolicy<pcms::DeviceMemorySpace::execution_space>(0, 1),
    KOKKOS_LAMBDA(const int) {
      out(0) = jacobian(3, 0, 0);
      out(1) = volume(3);
    });
  const auto host =
    Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace{}, out);
  CHECK(host(0) == 1.0);
  CHECK(host(1) == 1.0);
}

TEST_CASE("Periodic geometry: Gauss per cell integrates polynomials of its "
          "order exactly",
          "[periodic][quadrature]")
{
  const auto geometry = pcms::test::PeriodicGrid(3, 3.0, -1.0);
  for (int order = 0; order <= 9; ++order) {
    CAPTURE(order);
    const auto quadrature = pcms::CreateGaussQuadrature(geometry, order);
    // integral of x^order over [-1, 2)
    const Real exact =
      (std::pow(2.0, order + 1) - std::pow(-1.0, order + 1)) / (order + 1);
    CHECK(Integrate(quadrature, [&](Real x) { return std::pow(x, order); }) ==
          Catch::Approx(exact).margin(1e-12));
  }
}

TEST_CASE("Periodic geometry: the trapezoid rule is exact for trigonometric "
          "polynomials",
          "[periodic][quadrature]")
{
  const auto geometry = pcms::test::PeriodicGrid(4);
  const int n = 9; // exact up to degree 8
  const auto quadrature = pcms::CreateTrapezoidQuadrature(geometry, n);
  for (int j = 0; j <= 4; ++j) {
    for (int k = 0; k <= 4; ++k) {
      CAPTURE(j, k);
      const Real cc = Integrate(
        quadrature, [&](Real x) { return std::cos(j * x) * std::cos(k * x); });
      const Real expected = j != k ? 0.0 : (j == 0 ? TwoPi : TwoPi / 2);
      CHECK(cc == Catch::Approx(expected).margin(1e-12));
    }
  }
  // the rule of degree n - 1 is not exact beyond it
  CHECK(Integrate(quadrature, [](Real x) { return std::cos(9 * x); }) ==
        Catch::Approx(TwoPi).margin(1e-12));
}
