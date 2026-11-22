#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <Kokkos_Core.hpp>
#include <pcms/backends/uniform_grid/geometry.hpp>
#include <pcms/field/function_spaces.hpp>
#include <pcms/field/gradient_source.hpp>
#include <pcms/field/query_points.hpp>
#include <pcms/transfer/copy.hpp>
#include <pcms/transfer/interpolator.hpp>
#include <pcms/utility/coordinate_systems/cartesian.hpp>
#include "field_test_utils.hpp"
#include <vector>

using Catch::Matchers::WithinAbs;
using pcms::DifferentialOperator;
using pcms::OutOfBoundsMode;
using pcms::OutOfBoundsPolicy;
using pcms::Real;

namespace
{

constexpr Real Tol = pcms::test::EvalTol;

KOKKOS_INLINE_FUNCTION Real Linear(Real x)
{
  return 3.0 * x - 1.0;
}

KOKKOS_INLINE_FUNCTION Real Cubic(Real x)
{
  return x * x * x - 2.0 * x * x + x + 3.0;
}

KOKKOS_INLINE_FUNCTION Real CubicDerivative(Real x)
{
  return 3.0 * x * x - 4.0 * x + 1.0;
}

std::shared_ptr<const pcms::Geometry> MakeGeometry(Real bot_left, Real length,
                                                   pcms::LO divisions)
{
  pcms::UniformGrid1D grid;
  grid.bot_left = {bot_left};
  grid.edge_length = {length};
  grid.divisions = {divisions};
  return pcms::uniform_grid::CreateGeometry(grid,
                                            pcms::csys::CreateCartesian(1));
}

/// `op` of `f` at `pts`, one component, on a host view.
Kokkos::View<Real**, pcms::HostMemorySpace> EvaluateAt(
  const pcms::FunctionSpace& space, const pcms::Function<Real>& f,
  const std::vector<Real>& pts,
  DifferentialOperator op = DifferentialOperator::Value,
  OutOfBoundsPolicy policy = {})
{
  const auto bound = pcms::test::Bind(*space.GetGeometry(), pts, policy);
  const auto evaluator = space.CreatePointEvaluator<Real>(bound, op);
  return pcms::test::EvaluateToHost(
    *evaluator, f, evaluator->OutputBasis(f.GetData().GetValueBasis()));
}

} // namespace

TEST_CASE("1D uniform-grid layouts hold the vertices or the cell centers")
{
  const auto geometry = MakeGeometry(-1.0, 4.0, 8);
  const auto p1 = pcms::CreateLagrangeSpace(geometry, 1, 1);
  const auto p0 = pcms::CreateLagrangeSpace(geometry, 0, 1);
  REQUIRE(p1->GetLayout()->GetNumDofHolder() == 9);
  REQUIRE(p0->GetLayout()->GetNumDofHolder() == 8);
  REQUIRE(p1->GetLayout()->GetDimension() == 1);
  const auto coords = pcms::test::CopyCoordinatesToHost(
    p1->GetLayout()->GetDOFHolderCoordinates().GetValues());
  const auto centers = pcms::test::CopyCoordinatesToHost(
    p0->GetLayout()->GetDOFHolderCoordinates().GetValues());
  REQUIRE(coords.extent(1) == 1);
  REQUIRE_THAT(coords(0, 0), WithinAbs(-1.0, Tol));
  REQUIRE_THAT(coords(8, 0), WithinAbs(3.0, Tol));
  REQUIRE_THAT(centers(0, 0), WithinAbs(-0.75, Tol));
  REQUIRE_THAT(centers(7, 0), WithinAbs(2.75, Tol));
}

TEST_CASE("1D linear Lagrange reproduces linear functions")
{
  const auto space =
    pcms::CreateLagrangeSpace(MakeGeometry(-1.0, 4.0, 5), 1, 1);
  auto f = space->CreateFunction<Real>("f");
  pcms::test::SetField(f, KOKKOS_LAMBDA(Real x) { return Linear(x); });
  // interior, cell boundaries and both ends of the grid
  const std::vector<Real> pts{-1.0, -0.9, -0.2, 0.0, 1.3, 2.9, 3.0};
  const auto values = EvaluateAt(*space, f, pts);
  for (size_t i = 0; i < pts.size(); ++i) {
    CAPTURE(pts[i]);
    REQUIRE_THAT(values(i, 0), WithinAbs(Linear(pts[i]), Tol));
  }
}

TEST_CASE("1D piecewise-constant Lagrange evaluates the cell value")
{
  const auto space = pcms::CreateLagrangeSpace(MakeGeometry(0.0, 4.0, 4), 0, 1);
  auto f = space->CreateFunction<Real>("f");
  pcms::test::SetField(f, KOKKOS_LAMBDA(Real x) { return Linear(x); });
  const auto values = EvaluateAt(*space, f, {0.2, 1.7, 3.9});
  REQUIRE_THAT(values(0, 0), WithinAbs(Linear(0.5), Tol));
  REQUIRE_THAT(values(1, 0), WithinAbs(Linear(1.5), Tol));
  REQUIRE_THAT(values(2, 0), WithinAbs(Linear(3.5), Tol));
  REQUIRE_THROWS(space->CreatePointEvaluator<Real>(
    pcms::test::Bind(*space->GetGeometry(), {0.5}),
    DifferentialOperator::Gradient));
}

TEST_CASE("1D multi-component Lagrange")
{
  const auto space = pcms::CreateLagrangeSpace(MakeGeometry(0.0, 1.0, 4), 1, 2);
  auto f = space->CreateFunction<Real>("f");
  pcms::test::SetFieldComponents(f, [](Real x, Real* out) {
    out[0] = x;
    out[1] = 2.0 - x;
  });
  const auto bound = pcms::test::Bind(*space->GetGeometry(), {0.3, 0.8});
  const auto evaluator = space->CreatePointEvaluator<Real>(bound);
  const auto values =
    pcms::test::EvaluateToHost(*evaluator, f, f.GetData().GetValueBasis(), 2);
  REQUIRE_THAT(values(0, 0), WithinAbs(0.3, Tol));
  REQUIRE_THAT(values(0, 1), WithinAbs(1.7, Tol));
  REQUIRE_THAT(values(1, 0), WithinAbs(0.8, Tol));
}

TEST_CASE("1D cubic spline reproduces cubics, ends included")
{
  const auto space = pcms::CreateSplineSpace(MakeGeometry(-1.0, 4.0, 8));
  auto f = space->CreateFunction<Real>("f");
  pcms::test::SetField(f, KOKKOS_LAMBDA(Real x) { return Cubic(x); });
  // not-a-knot boundary conditions: exact on the whole grid, ends included
  const std::vector<Real> pts{-1.0, -0.93, -0.3, 0.0, 0.77, 1.9, 2.99, 3.0};
  const auto values = EvaluateAt(*space, f, pts);
  const auto gradient =
    EvaluateAt(*space, f, pts, DifferentialOperator::Gradient);
  for (size_t i = 0; i < pts.size(); ++i) {
    CAPTURE(pts[i]);
    REQUIRE_THAT(values(i, 0), WithinAbs(Cubic(pts[i]), Tol));
    REQUIRE_THAT(gradient(i, 0), WithinAbs(CubicDerivative(pts[i]), Tol));
  }
  REQUIRE_THROWS(space->CreatePointEvaluator<Real>(
    pcms::test::Bind(*space->GetGeometry(), {0.5}),
    DifferentialOperator::Hessian));
}

TEST_CASE("1D spline of a smooth function converges")
{
  const auto space = pcms::CreateSplineSpace(MakeGeometry(0.0, 3.0, 32));
  auto f = space->CreateFunction<Real>("f");
  pcms::test::SetField(f, KOKKOS_LAMBDA(Real x) { return Kokkos::sin(x); });
  const std::vector<Real> pts{0.05, 0.9, 1.57, 2.4, 2.95};
  const auto values = EvaluateAt(*space, f, pts);
  const auto gradient =
    EvaluateAt(*space, f, pts, DifferentialOperator::Gradient);
  for (size_t i = 0; i < pts.size(); ++i) {
    CAPTURE(pts[i]);
    REQUIRE_THAT(values(i, 0), WithinAbs(std::sin(pts[i]), 1e-6));
    REQUIRE_THAT(gradient(i, 0), WithinAbs(std::cos(pts[i]), 1e-4));
  }
}

TEST_CASE("1D evaluation fills points outside the grid")
{
  const OutOfBoundsPolicy fill{OutOfBoundsMode::FILL, -7.0};
  for (const bool spline : {false, true}) {
    CAPTURE(spline);
    const auto geometry = MakeGeometry(0.0, 2.0, 8);
    const auto space = spline ? pcms::CreateSplineSpace(geometry)
                              : pcms::CreateLagrangeSpace(geometry, 1, 1);
    auto f = space->CreateFunction<Real>("f");
    pcms::test::SetField(f, KOKKOS_LAMBDA(Real x) { return Linear(x); });
    const auto values = EvaluateAt(*space, f, {-0.5, 1.0, 2.5},
                                   DifferentialOperator::Value, fill);
    REQUIRE(values(0, 0) == -7.0);
    REQUIRE_THAT(values(1, 0), WithinAbs(Linear(1.0), Tol));
    REQUIRE(values(2, 0) == -7.0);
    const auto bound =
      pcms::test::Bind(*space->GetGeometry(), {-0.5, 1.0, 2.5}, fill);
    REQUIRE(bound->MissingPoints().extent(0) == 2);
  }
  // NEAREST_BOUNDARY extrapolates linearly from the closest cell
  const auto space = pcms::CreateLagrangeSpace(MakeGeometry(0.0, 2.0, 8), 1, 1);
  auto f = space->CreateFunction<Real>("f");
  pcms::test::SetField(f, KOKKOS_LAMBDA(Real x) { return Linear(x); });
  const auto values =
    EvaluateAt(*space, f, {-0.5, 2.5}, DifferentialOperator::Value,
               OutOfBoundsPolicy{OutOfBoundsMode::NEAREST_BOUNDARY, 0.0});
  REQUIRE_THAT(values(0, 0), WithinAbs(Linear(-0.5), Tol));
  REQUIRE_THAT(values(1, 0), WithinAbs(Linear(2.5), Tol));
}

TEST_CASE("1D spaces interpolate between grids")
{
  const auto source = pcms::CreateSplineSpace(MakeGeometry(-1.0, 4.0, 8));
  auto f = source->CreateFunction<Real>("f");
  pcms::test::SetField(f, KOKKOS_LAMBDA(Real x) { return Cubic(x); });

  // the target grid is a different, finer one inside the source grid
  const auto target_geometry = MakeGeometry(-0.5, 3.0, 20);
  const auto value_space = pcms::CreateLagrangeSpace(target_geometry, 1, 1);
  const auto gradient_space = pcms::CreateLagrangeSpace(target_geometry, 1, 1);
  auto value = value_space->CreateFunction<Real>("f");
  auto gradient =
    gradient_space->CreateFunction<Real>("df", pcms::values::Covector);
  pcms::CreateInterpolator<Real>(source, value_space)->Apply(f, value);
  pcms::CreateInterpolator<Real>(pcms::Gradient(source), gradient_space)
    ->Apply(f, gradient);

  const auto coords = pcms::test::CopyCoordinatesToHost(
    value_space->GetLayout()->GetDOFHolderCoordinates().GetValues());
  const auto v = value.GetDOFHolderDataHost().GetValues();
  const auto g = gradient.GetDOFHolderDataHost().GetValues();
  REQUIRE(coords.extent(0) == 21);
  for (size_t i = 0; i < coords.extent(0); ++i) {
    CAPTURE(i, coords(i, 0));
    REQUIRE_THAT(v(i, 0), WithinAbs(Cubic(coords(i, 0)), Tol));
    REQUIRE_THAT(g(i, 0), WithinAbs(CubicDerivative(coords(i, 0)), Tol));
  }

  // linear Lagrange to linear Lagrange is exact for a linear function
  const auto lagrange =
    pcms::CreateLagrangeSpace(MakeGeometry(0.0, 1.0, 4), 1, 1);
  const auto fine = pcms::CreateLagrangeSpace(MakeGeometry(0.1, 0.8, 13), 1, 1);
  auto lf = lagrange->CreateFunction<Real>("f");
  auto ff = fine->CreateFunction<Real>("f");
  pcms::test::SetField(lf, KOKKOS_LAMBDA(Real x) { return Linear(x); });
  pcms::CreateInterpolator<Real>(lagrange, fine)->Apply(lf, ff);
  const auto fine_coords = pcms::test::CopyCoordinatesToHost(
    fine->GetLayout()->GetDOFHolderCoordinates().GetValues());
  const auto fv = ff.GetDOFHolderDataHost().GetValues();
  for (size_t i = 0; i < fine_coords.extent(0); ++i) {
    REQUIRE_THAT(fv(i, 0), WithinAbs(Linear(fine_coords(i, 0)), Tol));
  }

  // copy between functions of one space
  auto copy = lagrange->CreateFunction<Real>("copy");
  pcms::CreateCopy<Real>(lagrange, lagrange)->Apply(lf, copy);
  REQUIRE_THAT(copy.GetDOFHolderDataHost().GetValues()(3, 0),
               WithinAbs(Linear(0.75), Tol));
}

TEST_CASE("1D QueryPoints locate twice")
{
  const auto space = pcms::CreateSplineSpace(MakeGeometry(-1.0, 4.0, 8));
  auto f = space->CreateFunction<Real>("f");
  pcms::test::SetField(f, KOKKOS_LAMBDA(Real x) { return Cubic(x); });
  const auto basis = f.GetData().GetValueBasis();

  pcms::QueryPoints query(space->GetGeometry(),
                          OutOfBoundsPolicy{OutOfBoundsMode::FILL, -1.0});
  for (const auto& pts :
       {std::vector<Real>{0.1, 5.0, 2.2}, std::vector<Real>{-0.4, 0.9}}) {
    const auto coords = pcms::test::CreateDeviceCoordinateView(
      pts, space->GetCoordinateSystem(), 1);
    query.Locate(coords.coordinate_view);
    Kokkos::View<Real**, pcms::DeviceMemorySpace> out("out", pts.size(), 1);
    query.Evaluate(f, pcms::ValueView<Real, pcms::DeviceMemorySpace>(
                        basis, pcms::MakeRank2View(out)));
    const auto host =
      Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace(), out);
    for (size_t i = 0; i < pts.size(); ++i) {
      CAPTURE(pts[i]);
      if (pts[i] > 3.0) {
        REQUIRE(host(i, 0) == -1.0);
      } else {
        REQUIRE_THAT(host(i, 0), WithinAbs(Cubic(pts[i]), Tol));
      }
    }
  }
}
