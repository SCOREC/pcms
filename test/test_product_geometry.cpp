#include "periodic_test_utils.hpp"
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "embedding_test_utils.hpp"
#include "field_test_utils.hpp"
#include <Omega_h_library.hpp>
#include <pcms/geometry/product_geometry.hpp>
#include <pcms/utility/coordinate_systems/cartesian.hpp>
#include "pcms/field/function_space/product_binding.hpp"
#include "pcms/backends/uniform_grid/uniform_grid_geometry.hpp"
#include <numbers>

namespace
{

using pcms::Real;
constexpr Real TwoPi = 2 * std::numbers::pi_v<Real>;

std::shared_ptr<const pcms::Geometry> Periodic(int cells, Real period,
                                               Real origin = 0)
{
  return pcms::test::PeriodicGrid(cells, period, origin);
}

template <typename F>
Real Integrate(const pcms::Quadrature& quadrature, F f)
{
  const auto points = pcms::test::CopyCoordinatesToHost(
    quadrature.points->LocatedPoints().GetValues());
  const auto weights = Kokkos::create_mirror_view_and_copy(
    pcms::HostMemorySpace{}, quadrature.weights);
  Real sum = 0;
  for (std::size_t p = 0; p < points.extent(0); ++p) {
    std::vector<Real> x;
    for (std::size_t d = 0; d < points.extent(1); ++d) {
      x.push_back(points(p, d));
    }
    sum += weights(p) * f(x);
  }
  return sum;
}

// The geometric factors at the points of `bound`, on the host: the Jacobian
// and the metric flattened [point][i * columns + j] (rows: the system,
// columns: the parameters), and the volume factor.
struct HostFactors
{
  int rows = 0;
  int columns = 0;
  Kokkos::View<Real**, pcms::DeviceMemorySpace>::host_mirror_type jacobian;
  Kokkos::View<Real**, pcms::DeviceMemorySpace>::host_mirror_type metric;
  Kokkos::View<Real*, pcms::DeviceMemorySpace>::host_mirror_type volume;
  [[nodiscard]] Real J(int p, int i, int j) const
  {
    return jacobian(p, i * columns + j);
  }
};

HostFactors FactorsOnHost(const pcms::BoundPoints& bound)
{
  const auto j = bound.Jacobian();
  const auto g = bound.Metric();
  const auto v = bound.VolumeFactor();
  const int rows = static_cast<int>(j.extent(1));
  const int columns = static_cast<int>(j.extent(2));
  const auto n = bound.NumPoints();
  Kokkos::View<Real**, pcms::DeviceMemorySpace> jacobian(
    "j", n, static_cast<std::size_t>(rows) * columns);
  Kokkos::View<Real**, pcms::DeviceMemorySpace> metric(
    "g", n, static_cast<std::size_t>(columns) * columns);
  Kokkos::View<Real*, pcms::DeviceMemorySpace> volume("v", n);
  Kokkos::parallel_for(
    "copy", Kokkos::RangePolicy<pcms::DeviceMemorySpace::execution_space>(0, n),
    KOKKOS_LAMBDA(const pcms::LO p) {
      for (int a = 0; a < rows; ++a) {
        for (int b = 0; b < columns; ++b) {
          jacobian(p, a * columns + b) = j(p, a, b);
        }
      }
      for (int a = 0; a < columns; ++a) {
        for (int b = 0; b < columns; ++b) {
          metric(p, a * columns + b) = g(p, a, b);
        }
      }
      volume(p) = v(p);
    });
  return {
    rows, columns,
    Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace{}, jacobian),
    Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace{}, metric),
    Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace{}, volume)};
}

} // namespace

TEST_CASE("Product geometry: Bind splits the coordinates by factor and keeps "
          "the factors' bound points",
          "[product][geometry]")
{
  const auto a = Periodic(4, 2.0);
  const auto b = Periodic(3, 3.0, 1.0);
  const auto product = pcms::CreateProductGeometry(a, b);
  CHECK(product->GetCoordinateSystem()->Dimension() == 2);
  const auto points = pcms::test::DevicePoints({{0.6, 1.5}, {-0.1, 4.2}});
  const auto bound =
    product->Bind(pcms::CoordinateView<pcms::DeviceMemorySpace>(
      product->GetCoordinateSystem(), pcms::MakeConstRank2View(points)));
  CHECK(bound->MissingPoints().extent(0) == 0);
  const auto& binding = dynamic_cast<const pcms::ProductBinding&>(*bound);
  CHECK(binding.FactorA()->GetGeometry() == a);
  CHECK(binding.FactorB()->GetGeometry() == b);
  const auto cell_a = Kokkos::create_mirror_view_and_copy(
    pcms::HostMemorySpace{},
    dynamic_cast<const pcms::detail::UniformGridBinding<1>&>(*binding.FactorA())
      .CellIndices());
  const auto cell_b = Kokkos::create_mirror_view_and_copy(
    pcms::HostMemorySpace{},
    dynamic_cast<const pcms::detail::UniformGridBinding<1>&>(*binding.FactorB())
      .CellIndices());
  // a: h = 0.5, 0.6 -> 1, -0.1 -> 1.9 -> 3; b: [1, 4), h = 1, 1.5 -> 0,
  // 4.2 -> 1.2 -> 0
  CHECK(cell_a(0) == 1);
  CHECK(cell_a(1) == 3);
  CHECK(cell_b(0) == 0);
  CHECK(cell_b(1) == 0);
  CHECK_THROWS(pcms::CreateProductGeometry(a, nullptr));
}

TEST_CASE("Product geometry: tensor quadrature and the factors' geometric "
          "factors",
          "[product][geometry]")
{
  const auto product =
    pcms::CreateProductGeometry(Periodic(4, 2.0), Periodic(3, 3.0));
  const auto quadrature = pcms::CreateGaussQuadrature(product, 3);
  CHECK(Integrate(quadrature, [](const auto&) { return 1.0; }) ==
        Catch::Approx(6.0).margin(1e-12));
  // x^3 y^2 over [0, 2) x [0, 3) = 4 * 9
  CHECK(Integrate(quadrature, [](const auto& x) {
          return x[0] * x[0] * x[0] * x[1] * x[1];
        }) == Catch::Approx(36.0).margin(1e-11));
  // the factors' coordinate systems are Cartesian: identities
  const auto factors = FactorsOnHost(*quadrature.points);
  CHECK(factors.J(5, 0, 0) == 1.0);
  CHECK(factors.J(5, 1, 1) == 1.0);
  CHECK(factors.J(5, 0, 1) == 0.0);
  CHECK(factors.volume(5) == 1.0);
}

TEST_CASE("Product geometry: the geometric factors are the block composition "
          "of the factors'",
          "[product][geometry]")
{
  // an embedded plane x = A xi + b (J(i, j) = A[i][j], volume |det A| = 2.85)
  // times a periodic interval
  const Real a[2][2] = {{2.0, 0.5}, {0.3, 1.5}};
  const auto grid = pcms::test::MakeUnitSquareGrid(4);
  const auto plane = pcms::CreateEmbeddedGeometry(
    grid, pcms::test::MakeSplinePosition(
            grid,
            [&](Real xi1, Real xi2) {
              return std::array<Real, 2>{a[0][0] * xi1 + a[0][1] * xi2,
                                         a[1][0] * xi1 + a[1][1] * xi2};
            },
            pcms::csys::CreateCartesian(2)));
  const auto product = pcms::CreateProductGeometry(plane, Periodic(3, 3.0));
  // the images of xi = (0.3, 0.6) and (0.8, 0.2), at two angles
  const auto points =
    pcms::test::DevicePoints({{0.9, 0.99, 0.5}, {1.7, 0.54, 2.5}});
  const auto bound =
    product->Bind(pcms::CoordinateView<pcms::DeviceMemorySpace>(
      product->GetCoordinateSystem(), pcms::MakeConstRank2View(points)));
  const auto factors = FactorsOnHost(*bound);
  REQUIRE(factors.rows == 3);
  REQUIRE(factors.columns == 3);
  for (const int p : {0, 1}) {
    CAPTURE(p);
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 2; ++j) {
        CHECK_THAT(factors.J(p, i, j),
                   Catch::Matchers::WithinAbs(a[i][j], 1e-10));
      }
      CHECK(factors.J(p, i, 2) == 0.0);
      CHECK(factors.J(p, 2, i) == 0.0);
    }
    CHECK(factors.J(p, 2, 2) == 1.0);
    CHECK(factors.metric(p, 2 * 3 + 2) == 1.0);
    CHECK_THAT(factors.volume(p), Catch::Matchers::WithinAbs(2.85, 1e-10));
  }
}

TEST_CASE("Product geometry: Omega_h mesh times a periodic interval",
          "[product][geometry]")
{
  Omega_h::Library lib;
  auto mesh = pcms::test::BuildSquareMesh(lib, 3);
  const auto plane =
    pcms::omega_h::CreateGeometry(mesh, pcms::detail::DeferredCartesian());
  const auto product = pcms::CreateProductGeometry(plane, Periodic(5, TwoPi));
  CHECK(product->GetCoordinateSystem()->Dimension() == 3);
#ifdef PCMS_ENABLE_OMEGA_H_MESHFIELDS
  const auto quadrature = pcms::CreateGaussQuadrature(product, 2);
  CHECK(Integrate(quadrature, [](const auto&) { return 1.0; }) ==
        Catch::Approx(TwoPi).margin(1e-12));
  CHECK(Integrate(quadrature, [](const auto& x) { return x[0] * x[1]; }) ==
        Catch::Approx(0.25 * TwoPi).margin(1e-11));
#endif
  // a point outside the plane's mesh is missing; ERROR refuses it
  const auto points =
    pcms::test::DevicePoints({{0.3, 0.4, 9.0}, {1.5, 0.4, 1.0}});
  const pcms::CoordinateView<pcms::DeviceMemorySpace> view(
    product->GetCoordinateSystem(), pcms::MakeConstRank2View(points));
  CHECK_THROWS(product->Bind(view));
  const auto bound = product->Bind(
    view, pcms::OutOfBoundsPolicy{pcms::OutOfBoundsMode::FILL, 0.0});
  const auto missing = Kokkos::create_mirror_view_and_copy(
    pcms::HostMemorySpace{}, bound->MissingPoints());
  REQUIRE(missing.extent(0) == 1);
  CHECK(missing(0) == 1);
  // the factors of the Omega_h factor are its coordinate system's
  CHECK(FactorsOnHost(*quadrature.points).volume(0) == 1.0);
}

TEST_CASE("Product geometry: Rebind re-locates in place of the previous bind",
          "[product][geometry]")
{
  const auto product =
    pcms::CreateProductGeometry(Periodic(4, 2.0), Periodic(3, 3.0));
  auto points = pcms::test::DevicePoints({{0.1, 0.2}, {0.3, 0.4}});
  auto bound = product->Bind(pcms::CoordinateView<pcms::DeviceMemorySpace>(
    product->GetCoordinateSystem(), pcms::MakeConstRank2View(points)));
  points = pcms::test::DevicePoints({{1.9, 2.9}, {0.3, 0.4}});
  bound = product->Rebind(
    std::move(bound),
    pcms::CoordinateView<pcms::DeviceMemorySpace>(
      product->GetCoordinateSystem(), pcms::MakeConstRank2View(points)));
  const auto& binding = dynamic_cast<const pcms::ProductBinding&>(*bound);
  const auto cell = Kokkos::create_mirror_view_and_copy(
    pcms::HostMemorySpace{},
    dynamic_cast<const pcms::detail::UniformGridBinding<1>&>(*binding.FactorA())
      .CellIndices());
  CHECK(cell(0) == 3);
}
