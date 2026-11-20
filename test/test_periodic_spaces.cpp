#include "periodic_test_utils.hpp"
#include <pcms/geometry/periodic_geometry.hpp>
#include <pcms/transfer/forms.hpp>
#include <pcms/transfer/projection.hpp>
#include <pcms/transfer/linear_solver.hpp>
#include <pcms/utility/coordinate_systems/cartesian.hpp>
#include <chrono>
#include <numbers>

namespace
{

using pcms::Real;
using pcms::test::EvaluateAt;
using pcms::test::SetDofs;
constexpr Real TwoPi = 2 * std::numbers::pi_v<Real>;

std::shared_ptr<const pcms::Geometry> Periodic(int cells, Real period = TwoPi,
                                               Real origin = 0)
{
  return pcms::CreatePeriodicGeometry(pcms::csys::CreateCartesian(1), cells,
                                      period, origin);
}

std::vector<std::vector<Real>> Points(int n, Real lo, Real hi)
{
  std::vector<std::vector<Real>> points;
  for (int i = 0; i < n; ++i) {
    points.push_back({lo + (hi - lo) * (i + 0.37) / n});
  }
  return points;
}

} // namespace

TEST_CASE("Fourier space: evaluation and derivative of a mode list with "
          "stride, including wrapped points",
          "[periodic][fourier]")
{
  const Real period = 4.0;
  const Real origin = -1.0;
  const auto geometry = Periodic(8, period, origin);
  const auto space = pcms::CreateFourierSpace(geometry, {0, 1, 3}, 2);
  // DOFs: 1; cos(2 t), sin(2 t); cos(6 t), sin(6 t)
  REQUIRE(space->GetLayout()->GetNumDofHolder() == 5);
  const std::vector<Real> c = {0.5, 2.0, -1.0, 0.25, 3.0};
  auto u = space->CreateFunction<Real>();
  SetDofs(u, c);
  const Real w = TwoPi / period;
  const auto f = [&](Real x) {
    const Real t = w * (x - origin);
    return c[0] + c[1] * std::cos(2 * t) + c[2] * std::sin(2 * t) +
           c[3] * std::cos(6 * t) + c[4] * std::sin(6 * t);
  };
  const auto df = [&](Real x) {
    const Real t = w * (x - origin);
    return w * (-2 * c[1] * std::sin(2 * t) + 2 * c[2] * std::cos(2 * t) -
                6 * c[3] * std::sin(6 * t) + 6 * c[4] * std::cos(6 * t));
  };
  // inside, below and above the interval
  const auto points = Points(24, -9.0, 11.0);
  const auto value = EvaluateAt(*space, u, points);
  const auto slope =
    EvaluateAt(*space, u, points, pcms::DifferentialOperator::Gradient);
  for (std::size_t p = 0; p < points.size(); ++p) {
    CHECK(value[p][0] == Catch::Approx(f(points[p][0])).margin(1e-11));
    CHECK(slope[p][0] == Catch::Approx(df(points[p][0])).margin(1e-10));
  }
  // the wrap: a point one period away has the same value
  const auto a = EvaluateAt(*space, u, {{0.3}});
  const auto b = EvaluateAt(*space, u, {{0.3 + period}});
  CHECK(a[0][0] == Catch::Approx(b[0][0]).margin(1e-11));
  CHECK_THROWS(pcms::CreateFourierSpace(geometry, {-1}));
  CHECK_THROWS(pcms::CreateFourierSpace(geometry, {1}, 0));
}

TEST_CASE("Periodic Lagrange order 1: nodes, midpoints and slopes",
          "[periodic][lagrange]")
{
  const int n = 6;
  const auto geometry = Periodic(n, 3.0, 1.0);
  const Real h = 0.5;
  const auto space = pcms::CreatePeriodicLagrangeSpace(geometry);
  auto u = space->CreateFunction<Real>();
  const std::vector<Real> g = {1, 4, -2, 0.5, 3, 7};
  SetDofs(u, g);
  for (int i = 0; i < n; ++i) {
    const Real x = 1.0 + h * i;
    const Real next = g[(i + 1) % n];
    const auto v = EvaluateAt(*space, u, {{x}, {x + 0.5 * h}});
    CHECK(v[0][0] == Catch::Approx(g[i]).margin(1e-13));
    CHECK(v[1][0] == Catch::Approx(0.5 * (g[i] + next)).margin(1e-13));
    const auto d = EvaluateAt(*space, u, {{x + 0.25 * h}},
                              pcms::DifferentialOperator::Gradient);
    CHECK(d[0][0] == Catch::Approx((next - g[i]) / h).margin(1e-12));
  }
  // the last cell wraps to the first node
  const auto wrapped = EvaluateAt(*space, u, {{1.0 + 3.0 - 0.25 * h}});
  CHECK(wrapped[0][0] ==
        Catch::Approx(0.25 * g[5] + 0.75 * g[0]).margin(1e-12));
}

TEST_CASE("Periodic Hermite: cubics exactly away from the wrap, sine to "
          "fourth order",
          "[periodic][hermite]")
{
  const int n = 8;
  const Real period = 2.0;
  const Real h = period / n;
  const auto geometry = Periodic(n, period);
  const auto space = pcms::CreatePeriodicHermiteSpace(geometry);
  REQUIRE(space->GetLayout()->GetNumDofHolder() == 2 * n);
  // a cubic's values and derivatives at the nodes: the cells whose two nodes
  // are not across the wrap reproduce it
  const auto p = [](Real x) { return 1 + x - 2 * x * x + 0.5 * x * x * x; };
  const auto dp = [](Real x) { return 1 - 4 * x + 1.5 * x * x; };
  std::vector<Real> dofs;
  for (int i = 0; i < n; ++i) {
    dofs.push_back(p(h * i));
    dofs.push_back(dp(h * i));
  }
  auto u = space->CreateFunction<Real>();
  SetDofs(u, dofs);
  const auto points = Points(14, 0.0, period - h);
  const auto value = EvaluateAt(*space, u, points);
  const auto slope =
    EvaluateAt(*space, u, points, pcms::DifferentialOperator::Gradient);
  for (std::size_t q = 0; q < points.size(); ++q) {
    CHECK(value[q][0] == Catch::Approx(p(points[q][0])).margin(1e-12));
    CHECK(slope[q][0] == Catch::Approx(dp(points[q][0])).margin(1e-11));
  }
  // sin on the whole circle
  const auto error = [&](int cells) {
    const auto fine = pcms::CreatePeriodicHermiteSpace(Periodic(cells, TwoPi));
    std::vector<Real> d;
    for (int i = 0; i < cells; ++i) {
      d.push_back(std::sin(TwoPi * i / cells));
      d.push_back(std::cos(TwoPi * i / cells));
    }
    auto s = fine->CreateFunction<Real>();
    SetDofs(s, d);
    Real worst = 0;
    const auto pts = Points(50, -1.0, 8.0);
    const auto v = EvaluateAt(*fine, s, pts);
    for (std::size_t q = 0; q < pts.size(); ++q) {
      worst = std::max(worst, std::abs(v[q][0] - std::sin(pts[q][0])));
    }
    return worst;
  };
  const Real coarse = error(8);
  const Real fine = error(16);
  CHECK(coarse < 2e-3);
  CHECK(coarse / fine > 12.0); // fourth order: 16
}

TEST_CASE("Periodic cubic spline: interpolates the nodes, reproduces "
          "constants, approximates a sine",
          "[periodic][spline]")
{
  const auto error = [](int cells, bool slope) {
    const auto space = pcms::CreatePeriodicSplineSpace(Periodic(cells, TwoPi));
    std::vector<Real> d;
    for (int i = 0; i < cells; ++i) {
      d.push_back(std::sin(TwoPi * i / cells));
    }
    auto s = space->CreateFunction<Real>();
    SetDofs(s, d);
    // interpolation of the node values
    std::vector<std::vector<Real>> nodes;
    for (int i = 0; i < cells; ++i) {
      nodes.push_back({TwoPi * i / cells});
    }
    const auto at_nodes = EvaluateAt(*space, s, nodes);
    for (int i = 0; i < cells; ++i) {
      CHECK(at_nodes[i][0] == Catch::Approx(d[i]).margin(1e-12));
    }
    Real worst = 0;
    const auto pts = Points(60, -1.0, 8.0);
    const auto v = EvaluateAt(*space, s, pts,
                              slope ? pcms::DifferentialOperator::Gradient
                                    : pcms::DifferentialOperator::Value);
    for (std::size_t q = 0; q < pts.size(); ++q) {
      const Real exact = slope ? std::cos(pts[q][0]) : std::sin(pts[q][0]);
      worst = std::max(worst, std::abs(v[q][0] - exact));
    }
    return worst;
  };
  CHECK(error(16, false) < 5e-4);
  CHECK(error(16, false) / error(32, false) > 12.0);
  CHECK(error(32, true) < 2e-3);
  const auto space = pcms::CreatePeriodicSplineSpace(Periodic(7, 2.0, 1.0));
  auto one = space->CreateFunction<Real>();
  SetDofs(one, std::vector<Real>(7, 3.5));
  const auto v = EvaluateAt(*space, one, Points(9, 0.0, 5.0));
  const auto d = EvaluateAt(*space, one, Points(9, 0.0, 5.0),
                            pcms::DifferentialOperator::Gradient);
  for (std::size_t q = 0; q < v.size(); ++q) {
    CHECK(v[q][0] == Catch::Approx(3.5).margin(1e-12));
    CHECK(d[q][0] == Catch::Approx(0.0).margin(1e-11));
  }
  CHECK_THROWS(pcms::CreatePeriodicSplineSpace(Periodic(2)));
}

TEST_CASE("Periodic cubic spline: the circulant moments match the dense "
          "solve, and a fine grid sets up fast",
          "[periodic][spline]")
{
  for (const int n : {8, 64}) {
    CAPTURE(n);
    const Real h = TwoPi / n;
    const auto f = pcms::test::Pseudo(static_cast<std::size_t>(n), 3);
    // M = S f with the dense S
    std::vector<Real> m(n, 0.0);
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        Real sum = 0;
        for (int k = 0; k < n; ++k) {
          const Real phi = TwoPi * k / n;
          sum += (2 * std::cos(phi) - 2) / h /
                 (h / 6 * (4 + 2 * std::cos(phi))) * std::cos(phi * (i - j));
        }
        m[i] += sum / n * f[j];
      }
    }
    const auto space = pcms::CreatePeriodicSplineSpace(Periodic(n));
    auto s = space->CreateFunction<Real>();
    SetDofs(s, f);
    const auto pts = Points(3 * n, 0.0, TwoPi);
    const auto v = EvaluateAt(*space, s, pts);
    for (std::size_t q = 0; q < pts.size(); ++q) {
      const int c = std::min(static_cast<int>(pts[q][0] / h), n - 1);
      const int c1 = (c + 1) % n;
      const Real t = pts[q][0] / h - c;
      const Real u = 1 - t;
      const Real expected =
        u * f[c] + t * f[c1] +
        h * h / 6 * ((u * u * u - u) * m[c] + (t * t * t - t) * m[c1]);
      CHECK(v[q][0] == Catch::Approx(expected).margin(1e-12));
    }
  }
  const auto start = std::chrono::steady_clock::now();
  const auto fine = pcms::CreatePeriodicSplineSpace(Periodic(1024));
  CHECK(std::chrono::steady_clock::now() - start < std::chrono::seconds(2));
  CHECK(fine->GetLayout()->GetNumDofHolder() == 1024);
}

TEST_CASE("Periodic spaces: <E x, y> = <x, E^T y> and the matrix is E",
          "[periodic][adjoint]")
{
  const auto geometry = Periodic(9, 3.0, -1.0);
  const auto points = Points(31, -5.0, 6.0);
  pcms::test::CheckAdjointAndMatrix(
    pcms::CreateFourierSpace(geometry, {0, 1, 2, 5}, 2), points);
  pcms::test::CheckAdjointAndMatrix(pcms::CreatePeriodicLagrangeSpace(geometry),
                                    points);
  pcms::test::CheckAdjointAndMatrix(pcms::CreatePeriodicHermiteSpace(geometry),
                                    points);
  pcms::test::CheckAdjointAndMatrix(pcms::CreatePeriodicSplineSpace(geometry),
                                    points);
}

TEST_CASE("Periodic spaces: natural quadrature integrates products of two "
          "functions exactly",
          "[periodic][quadrature]")
{
  const auto geometry = Periodic(6, 2.0, 0.5);
  for (const auto& space : {pcms::CreateFourierSpace(geometry, {0, 1, 2, 3}),
                            pcms::CreatePeriodicLagrangeSpace(geometry)}) {
    const auto quadrature = space->CreateQuadrature();
    const auto weights = Kokkos::create_mirror_view_and_copy(
      pcms::HostMemorySpace{}, quadrature.weights);
    Real total = 0;
    for (std::size_t p = 0; p < weights.extent(0); ++p) {
      total += weights(p);
    }
    CHECK(total == Catch::Approx(2.0).margin(1e-13));
  }
  // the mass matrix of Fourier modes is diagonal: 2 pi/period... and CG takes
  // two iterations (two distinct eigenvalues)
  const auto fourier = pcms::CreateFourierSpace(Periodic(4), {0, 1, 2, 3, 4});
  const pcms::transfer::detail::BilinearForm mass(fourier,
                                                  fourier->CreateQuadrature());
  pcms::transfer::detail::LinearSolver solver(mass, pcms::LinearSolveOptions{});
  auto b = fourier->CreateFunction<Real>();
  SetDofs(b, {1, 2, -3, 4, 5, -6, 7, 8, -9});
  auto x = fourier->CreateFunction<Real>();
  CHECK(solver.Solve(b, x) <= 2);
  const auto sol = pcms::test::GetDofs(x);
  CHECK(sol[0] == Catch::Approx(1 / TwoPi).epsilon(1e-12));
  CHECK(sol[1] == Catch::Approx(2 / (TwoPi / 2)).epsilon(1e-12));
}

TEST_CASE("Periodic spaces: Fourier projection reproduces the analytic "
          "coefficients; spline and Hermite projections are accurate",
          "[periodic][projection]")
{
  const auto geometry = Periodic(6);
  const auto source_space = pcms::CreateFourierSpace(geometry, {0, 1, 2});
  auto source = source_space->CreateFunction<Real>();
  const std::vector<Real> c = {0.5, 2.0, -1.0, 0.25, 3.0};
  SetDofs(source, c);
  for (const int harmonics : {1, 2, 4}) {
    CAPTURE(harmonics);
    std::vector<int> modes;
    for (int m = 0; m <= harmonics; ++m) {
      modes.push_back(m);
    }
    const auto target_space = pcms::CreateFourierSpace(geometry, modes);
    pcms::ProjectionOptions options;
    options.quadrature =
      pcms::CreateFourierSpace(geometry, {0, 1, 2, 3, 4, 5, 6})
        ->CreateQuadrature();
    auto target = target_space->CreateFunction<Real>();
    pcms::CreateProjection(source_space, target_space, options)
      ->Apply(source, target);
    const auto values = pcms::test::GetDofs(target);
    for (std::size_t m = 0; m < values.size(); ++m) {
      CHECK(values[m] ==
            Catch::Approx(m < c.size() ? c[m] : 0.0).margin(1e-12));
    }
  }
  // sin(x) + 0.3 cos(2 x) onto the spline and Hermite spaces of 16 cells
  const auto f = [](Real x) { return std::sin(x) + 0.3 * std::cos(2 * x); };
  const auto fine_space = pcms::CreateFourierSpace(geometry, {1, 2});
  auto smooth = fine_space->CreateFunction<Real>();
  SetDofs(smooth, {0.0, 1.0, 0.3, 0.0});
  for (const auto& space : {pcms::CreatePeriodicSplineSpace(Periodic(16)),
                            pcms::CreatePeriodicHermiteSpace(Periodic(16))}) {
    pcms::ProjectionOptions options;
    options.quadrature = pcms::CreateGaussQuadrature(space->GetGeometry(), 10);
    auto target = space->CreateFunction<Real>();
    pcms::CreateProjection(fine_space, space, options)->Apply(smooth, target);
    const auto pts = Points(40, 0.0, TwoPi);
    const auto v = EvaluateAt(*space, target, pts);
    for (std::size_t q = 0; q < pts.size(); ++q) {
      CHECK(v[q][0] == Catch::Approx(f(pts[q][0])).margin(2e-3));
    }
  }
}
