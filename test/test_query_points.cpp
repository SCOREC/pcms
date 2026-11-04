#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <Kokkos_Core.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_library.hpp>
#include <pcms/utility/coordinate_systems/cartesian.hpp>
#include <pcms/backends/uniform_grid/geometry.hpp>
#include <pcms/field/function_space/lagrange.hpp>
#include <pcms/field/function_space/spline.hpp>
#include <pcms/field/query_points.hpp>
#include "field_test_utils.hpp"
#include <vector>

using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;
using pcms::DeviceMemorySpace;
using pcms::DifferentialOperator;
using pcms::HostMemorySpace;
using pcms::OutOfBoundsMode;
using pcms::OutOfBoundsPolicy;
using pcms::Real;
using pcms::ValueBasis;

namespace
{

KOKKOS_INLINE_FUNCTION Real Cubic(Real x, Real y)
{
  return x * x + 3.0 * x * y - y * y * y;
}

/// Evaluates `op` of `f` at the located points of `query` into a new
/// [NumPoints()][width] view tagged `basis` and returns a host copy (in the
/// layout of the device view, which a GPU build does not share with the
/// default host view).
auto EvaluateToHost(pcms::QueryPoints& query, const pcms::Function<Real>& f,
                    const ValueBasis& basis, int width = 1,
                    DifferentialOperator op = DifferentialOperator::Value)
{
  Kokkos::View<Real**, DeviceMemorySpace> out("out", query.Bound()->NumPoints(),
                                              width);
  query.Evaluate(
    f,
    pcms::ValueView<Real, DeviceMemorySpace>(basis, pcms::MakeRank2View(out)),
    op);
  return Kokkos::create_mirror_view_and_copy(HostMemorySpace(), out);
}

} // namespace

TEST_CASE("QueryPoints evaluates at the points of the last Locate")
{
  auto lib = Omega_h::Library{};
  auto mesh = pcms::test::BuildSquareMesh(lib, 8);
  auto space = pcms::test::MakeP1Space(mesh);
  auto f = space->CreateFunction<Real>("f");
  pcms::test::SetField(
    f, OMEGA_H_LAMBDA(Real x, Real y) { return pcms::test::linear_f(x, y); });
  const auto basis = f.GetData().GetValueBasis();
  // the second point of `first` is outside the unit square
  const std::vector<Real> first{0.25, 0.5, 2.0, 2.0, 0.75, 0.1};
  const std::vector<Real> second{0.6, 0.6, 0.3, 0.2};
  auto first_coords =
    pcms::test::CreateDeviceCoordinateView(first, space->GetCoordinateSystem());
  auto second_coords = pcms::test::CreateDeviceCoordinateView(
    second, space->GetCoordinateSystem());

  pcms::QueryPoints query(space->GetGeometry(),
                          OutOfBoundsPolicy{OutOfBoundsMode::FILL, -1.0});
  REQUIRE_THROWS_WITH(
    query.Evaluate(f,
                   pcms::test::TagLike(
                     f, Kokkos::View<Real**, DeviceMemorySpace>("out", 3, 1))),
    ContainsSubstring("call Locate first"));
  query.Locate(first_coords.coordinate_view);
  REQUIRE(query.Bound()->MissingPoints().extent(0) == 1);
  auto values = EvaluateToHost(query, f, basis);
  REQUIRE_THAT(values(0, 0),
               WithinAbs(pcms::test::linear_f(0.25, 0.5), pcms::test::EvalTol));
  REQUIRE(values(1, 0) == -1.0);
  REQUIRE_THAT(values(2, 0),
               WithinAbs(pcms::test::linear_f(0.75, 0.1), pcms::test::EvalTol));

  // Locate drops the evaluators, so the bound source is reused, and the
  // evaluator created after it sees the new points: the explicit path's values
  const auto* address = query.Bound().get();
  query.Locate(second_coords.coordinate_view);
  REQUIRE(query.Bound().get() == address);
  values = EvaluateToHost(query, f, basis);
  const auto explicit_values = pcms::test::EvaluateToHost(
    *space->CreatePointEvaluator<Real>(
      space->GetGeometry()->Bind(second_coords.coordinate_view)),
    f, basis);
  for (int i = 0; i < 2; ++i) {
    REQUIRE_THAT(values(i, 0),
                 WithinAbs(explicit_values(i, 0), pcms::test::ExactTol));
  }

  // a function on another geometry, an operator the space does not support
  auto other = pcms::test::MakeP1Space(mesh);
  auto g = other->CreateFunction<Real>("g");
  REQUIRE_THROWS_WITH(EvaluateToHost(query, g, basis),
                      ContainsSubstring("another geometry"));
  REQUIRE_THROWS_WITH(
    EvaluateToHost(query, f, basis, 2, DifferentialOperator::Gradient),
    ContainsSubstring("Gradient is not supported"));
}

TEST_CASE("QueryPoints and the one-shot evaluation apply the operator")
{
  pcms::UniformGrid<2> grid;
  grid.bot_left = {0.0, -1.0};
  grid.edge_length = {2.0, 2.0};
  grid.divisions = {16, 16};
  auto space = pcms::CreateSplineSpace(
    pcms::CreateUniformGridGeometry(grid, pcms::csys::Cartesian::Create(2)));
  auto f = space->CreateFunction<Real>("f");
  pcms::test::SetField(
    f, OMEGA_H_LAMBDA(Real x, Real y) { return Cubic(x, y); });
  const std::vector<Real> pts{0.3, -0.7, 1.1, 0.2, 0.55, 0.45};
  auto coords =
    pcms::test::CreateDeviceCoordinateView(pts, space->GetCoordinateSystem());
  const auto bound = space->GetGeometry()->Bind(coords.coordinate_view);
  const auto gradient =
    space->CreatePointEvaluator<Real>(bound, DifferentialOperator::Gradient);
  const ValueBasis basis = gradient->OutputBasis(f.GetData().GetValueBasis());
  const auto expected = pcms::test::EvaluateToHost(*gradient, f, basis, 2);

  pcms::QueryPoints query(space->GetGeometry());
  query.Locate(coords.coordinate_view);
  const auto through_query =
    EvaluateToHost(query, f, basis, 2, DifferentialOperator::Gradient);
  Kokkos::View<Real**, DeviceMemorySpace> out("out", pts.size() / 2, 2);
  pcms::Evaluate(
    f, coords.coordinate_view,
    pcms::ValueView<Real, DeviceMemorySpace>(basis, pcms::MakeRank2View(out)),
    DifferentialOperator::Gradient);
  const auto one_shot =
    Kokkos::create_mirror_view_and_copy(HostMemorySpace(), out);
  for (size_t i = 0; i < pts.size() / 2; ++i) {
    for (int d = 0; d < 2; ++d) {
      CAPTURE(i, d);
      REQUIRE_THAT(through_query(i, d),
                   WithinAbs(expected(i, d), pcms::test::ExactTol));
      REQUIRE_THAT(one_shot(i, d),
                   WithinAbs(expected(i, d), pcms::test::ExactTol));
    }
  }
}
