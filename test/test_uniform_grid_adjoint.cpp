#include "evaluator_test_utils.hpp"
#include "field_test_utils.hpp"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <pcms/backends/uniform_grid/geometry.hpp>
#include <pcms/transfer/l2_inner_product.hpp>
#include <pcms/transfer/projection.hpp>
#include <pcms/utility/coordinate_systems/cartesian.hpp>
#include <pcms/utility/coordinate_systems/cylindrical.hpp>
#include <pcms/utility/coordinate_systems/parameter.hpp>
#include <cmath>
#include <numbers>
#include "pcms/field/point_rows.hpp"
#include <pcms/transfer/detail/projection.hpp>

namespace
{

using pcms::DeviceMemorySpace;
using pcms::HostMemorySpace;
using pcms::Real;
using pcms::test::GetDofs;
using pcms::test::Pseudo;
using pcms::test::SetDofs;
constexpr Real TwoPi = 2 * std::numbers::pi_v<Real>;

std::shared_ptr<const pcms::Geometry> Grid2(
  std::shared_ptr<const pcms::CoordinateSystem> system =
    pcms::csys::CreateCartesian(2),
  int nx = 5, int ny = 4)
{
  pcms::UniformGrid<2> grid;
  grid.bot_left = {-1.0, 0.5};
  grid.edge_length = {3.0, 2.0};
  grid.divisions = {nx, ny};
  return pcms::uniform_grid::CreateGeometry(grid, std::move(system));
}

std::shared_ptr<const pcms::Geometry> Grid3(
  std::shared_ptr<const pcms::CoordinateSystem> system =
    pcms::csys::CreateCartesian(3))
{
  pcms::UniformGrid<3> grid;
  grid.bot_left = {-1.0, 0.5, 0.0};
  grid.edge_length = {3.0, 2.0, 1.5};
  grid.divisions = {4, 3, 5};
  return pcms::uniform_grid::CreateGeometry(grid, std::move(system));
}

/// @brief n points spread over the box [lo, lo + length).
std::vector<std::vector<Real>> Points(int n, const std::vector<Real>& lo,
                                      const std::vector<Real>& length)
{
  std::vector<std::vector<Real>> points;
  for (int i = 0; i < n; ++i) {
    std::vector<Real> x;
    for (std::size_t d = 0; d < lo.size(); ++d) {
      x.push_back(
        lo[d] +
        length[d] *
          std::fmod(0.1234 + 0.6180339 * (i + 1) * (d + 1) + 0.31 * d, 1.0));
    }
    points.push_back(x);
  }
  return points;
}

std::shared_ptr<const pcms::BoundPoints> BindPoints(
  const pcms::Geometry& geometry, const std::vector<std::vector<Real>>& points,
  pcms::OutOfBoundsPolicy policy = {})
{
  const auto device = pcms::test::DevicePoints(points);
  return geometry.Bind(
    pcms::CoordinateView<DeviceMemorySpace>(geometry.GetCoordinateSystem(),
                                            pcms::MakeConstRank2View(device)),
    policy);
}

Real Dot(const std::vector<Real>& a, const std::vector<Real>& b)
{
  REQUIRE(a.size() == b.size());
  Real sum = 0;
  for (std::size_t i = 0; i < a.size(); ++i) {
    sum += a[i] * b[i];
  }
  return sum;
}

/// @return E^T y for the evaluator of `space` at `bound`
std::vector<Real> Transpose(const pcms::FunctionSpace& space,
                            std::shared_ptr<const pcms::BoundPoints> bound,
                            const std::vector<Real>& y)
{
  const auto evaluator = space.CreatePointEvaluator<Real>(std::move(bound));
  auto functional = space.CreateFunction<Real>("functional");
  const auto device_y =
    pcms::test::CreateDeviceRank2View(y, space.GetLayout()->GetNumComponents());
  evaluator->EvaluateTranspose(
    pcms::ValueView<const Real, DeviceMemorySpace>(
      functional.GetData().GetValueBasis(), pcms::MakeConstRank2View(device_y)),
    functional);
  return GetDofs(functional);
}

} // namespace

TEST_CASE("Uniform-grid Lagrange spaces: <E x, y> = <x, E^T y> and the "
          "rows are E",
          "[uniform_grid][adjoint]")
{
  for (const int order : {0, 1}) {
    CAPTURE(order);
    pcms::test::CheckAdjointAndRows(
      pcms::CreateLagrangeSpace(Grid2(), order, 1),
      Points(23, {-1.0, 0.5}, {3.0, 2.0}));
    pcms::test::CheckAdjointAndRows(
      pcms::CreateLagrangeSpace(Grid3(), order, 1),
      Points(41, {-1.0, 0.5, 0.0}, {3.0, 2.0, 1.5}));
  }
}

TEST_CASE("Uniform-grid Lagrange spaces: the transpose of several "
          "components, and of points outside the grid",
          "[uniform_grid][adjoint]")
{
  // the last two points are outside the grid
  auto points = Points(17, {-1.0, 0.5}, {3.0, 2.0});
  points.push_back({-1.5, 0.7});
  points.push_back({1.0, 3.1});
  for (const int order : {0, 1}) {
    const auto space = pcms::CreateLagrangeSpace(Grid2(), order, 2);
    const auto num_values =
      static_cast<std::size_t>(space->GetLayout()->GetNumDofHolder()) * 2;
    const auto x = Pseudo(num_values, 3);
    const auto y = Pseudo(points.size() * 2, 4);
    for (const auto mode : {pcms::OutOfBoundsMode::FILL,
                            pcms::OutOfBoundsMode::NEAREST_BOUNDARY}) {
      CAPTURE(order, static_cast<int>(mode));
      const auto bound =
        BindPoints(*space->GetGeometry(), points, {.mode = mode});
      auto u = space->CreateFunction<Real>("u");
      SetDofs(u, x);
      std::vector<Real> ex;
      for (const auto& row : pcms::test::EvaluateAt(
             *space, u, points, pcms::DifferentialOperator::Value,
             {.mode = mode})) {
        ex.insert(ex.end(), row.begin(), row.end());
      }
      if (mode == pcms::OutOfBoundsMode::FILL) {
        for (std::size_t k = ex.size() - 4; k < ex.size(); ++k) {
          ex[k] = 0;
        }
      }
      const auto ety = Transpose(*space, bound, y);
      CHECK(Dot(ex, y) ==
            Catch::Approx(Dot(x, ety)).epsilon(1e-12).margin(1e-12));
    }
  }
}

TEST_CASE("Uniform-grid Lagrange spaces: the rows of the gradient are those "
          "of the partials",
          "[uniform_grid][adjoint]")
{
  const auto space = pcms::CreateLagrangeSpace(Grid3(), 1, 1);
  const auto points = Points(13, {-1.0, 0.5, 0.0}, {3.0, 2.0, 1.5});
  const auto bound = BindPoints(*space->GetGeometry(), points);
  const auto evaluator = space->CreatePointEvaluator<Real>(
    bound, pcms::DifferentialOperator::Gradient);
  const auto* source = pcms::detail::AsRowsSource(*evaluator);
  REQUIRE(source != nullptr);
  const auto x =
    Pseudo(static_cast<std::size_t>(space->GetLayout()->GetNumDofHolder()), 5);
  auto u = space->CreateFunction<Real>("u");
  SetDofs(u, x);
  const auto gradient = pcms::test::EvaluateAt(
    *space, u, points, pcms::DifferentialOperator::Gradient);
  const auto& rows = source->Rows();
  const auto offsets =
    Kokkos::create_mirror_view_and_copy(HostMemorySpace{}, rows.offsets);
  const auto columns =
    Kokkos::create_mirror_view_and_copy(HostMemorySpace{}, rows.holders);
  const auto values =
    Kokkos::create_mirror_view_and_copy(HostMemorySpace{}, rows.values);
  REQUIRE(offsets.extent(0) == points.size() * 3 + 1);
  for (std::size_t p = 0; p < points.size(); ++p) {
    for (int j = 0; j < 3; ++j) {
      Real sum = 0;
      for (pcms::LO k = offsets(p * 3 + j); k < offsets(p * 3 + j + 1); ++k) {
        sum += values(k) * x[columns(k)];
      }
      CHECK(sum == Catch::Approx(gradient[p][j]).epsilon(1e-12).margin(1e-12));
    }
  }
}

namespace
{

/// @return the integral of x^a over [lo, hi]
Real MonomialIntegral(int a, Real lo, Real hi)
{
  return (std::pow(hi, a + 1) - std::pow(lo, a + 1)) / (a + 1);
}

/// @return the sum of the weights times f at the points of `quadrature`
template <typename F>
Real Integrate(const pcms::Quadrature& quadrature, F f)
{
  const auto x = pcms::test::CopyCoordinatesToHost(
    quadrature.points->LocatedPoints().GetValues());
  const auto w =
    Kokkos::create_mirror_view_and_copy(HostMemorySpace{}, quadrature.weights);
  Real sum = 0;
  for (std::size_t p = 0; p < w.extent(0); ++p) {
    sum += w(p) * f(&x(p, 0));
  }
  return sum;
}

} // namespace

TEST_CASE("Uniform-grid quadrature is exact for polynomials of its degree in "
          "each coordinate",
          "[uniform_grid][quadrature]")
{
  const std::vector<Real> lo = {-1.0, 0.5, 0.0};
  const std::vector<Real> hi = {2.0, 2.5, 1.5};
  for (const auto& system : {pcms::csys::CreateCartesian(3),
                             pcms::csys::CreateParameterSystem("xi", 3)}) {
    const auto grid = Grid3(system);
    for (int order = 1; order <= 5; ++order) {
      CAPTURE(order, system->Kind());
      const auto quadrature = pcms::CreateGaussQuadrature(grid, order);
      const int per_axis = order / 2 + 1;
      CHECK(quadrature.points->NumPoints() ==
            4 * 3 * 5 * per_axis * per_axis * per_axis);
      CHECK(quadrature.points->GetGeometry() == grid);
      for (int a = 0; a <= order; ++a) {
        for (int c = 0; c <= order; c += 2) {
          const int b = order - a;
          const Real expected = MonomialIntegral(a, lo[0], hi[0]) *
                                MonomialIntegral(b, lo[1], hi[1]) *
                                MonomialIntegral(c, lo[2], hi[2]);
          const Real sum = Integrate(quadrature, [&](const Real* x) {
            return std::pow(x[0], a) * std::pow(x[1], b) * std::pow(x[2], c);
          });
          CHECK(sum == Catch::Approx(expected).epsilon(1e-12).margin(1e-12));
        }
      }
    }
  }
  const auto rz = Grid2(pcms::csys::CreateCylindricalRZ());
  const auto quadrature = pcms::CreateGaussQuadrature(rz, 3);
  const Real sum = Integrate(
    quadrature, [](const Real* x) { return x[0] * x[0] * x[1] * x[1] * x[1]; });
  CHECK(sum == Catch::Approx(MonomialIntegral(3, -1.0, 2.0) *
                             MonomialIntegral(3, 0.5, 2.5))
                 .epsilon(1e-12));
  CHECK_THROWS(pcms::CreateMonteCarloQuadrature(rz, 2));
  CHECK_THROWS(pcms::CreateGaussQuadrature(rz, -1));
  CHECK(pcms::CreateLagrangeSpace(Grid2(), 0, 1)
          ->GetQuadrature()
          .points->NumPoints() == 20);
  CHECK(pcms::CreateLagrangeSpace(Grid2(), 1, 1)
          ->GetQuadrature()
          .points->NumPoints() == 80);
}

TEST_CASE("Galerkin projection into a uniform grid keeps the integral",
          "[uniform_grid][projection]")
{
  // nested cells: the target rules integrate source x target basis exactly
  const auto system = pcms::csys::CreateCylindricalRZ();
  pcms::UniformGrid<2> grid;
  grid.bot_left = {1.0, -0.5};
  grid.edge_length = {2.0, 1.0};
  grid.divisions = {6, 4};
  const auto source_space = pcms::CreateLagrangeSpace(
    pcms::uniform_grid::CreateGeometry(grid, system), 1, 1);
  grid.divisions = {12, 8};
  const auto target_geometry = pcms::uniform_grid::CreateGeometry(grid, system);
  auto source = source_space->CreateFunction<Real>("n");
  SetDofs(source, Pseudo(static_cast<std::size_t>(
                           source_space->GetLayout()->GetNumDofHolder()),
                         7));
  const Real expected = pcms::L2Integral(source_space).Apply(source);
  for (const int order : {0, 1}) {
    const auto target_space =
      pcms::CreateLagrangeSpace(target_geometry, order, 1);
    const pcms::L2Integral integral(target_space);
    for (const auto mass :
         {pcms::MassMatrixType::Consistent, pcms::MassMatrixType::Lumped}) {
      CAPTURE(order, static_cast<int>(mass));
      pcms::detail::ProjectionOptions options;
      options.mass_matrix_type = mass;
      options.solve.relative_tolerance = 1e-14;
      auto target = target_space->CreateFunction<Real>("n");
      pcms::detail::CreateProjection(source_space, target_space, options)
        ->Apply(source, target);
      CHECK(integral.Apply(target) == Catch::Approx(expected).epsilon(1e-12));
    }
  }
}
