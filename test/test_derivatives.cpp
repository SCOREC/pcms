#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <Kokkos_Core.hpp>
#include <pcms/utility/coordinate_systems/cartesian.hpp>
#include <pcms/utility/coordinate_systems/cylindrical.hpp>
#include <pcms/field/function_spaces.hpp>
#include <pcms/backends/uniform_grid/geometry.hpp>
#include <pcms/field/gradient_source.hpp>
#include <pcms/transfer/interpolator.hpp>
#include "field_test_utils.hpp"
#include <vector>

using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;
using pcms::DeviceMemorySpace;
using pcms::Real;
using pcms::ValueBasis;

namespace
{

// A bicubic not-a-knot spline reproduces this cubic exactly, and so its
// gradient.
KOKKOS_INLINE_FUNCTION Real Cubic(Real x, Real y)
{
  return x * x + 3.0 * x * y - y * y * y;
}

std::shared_ptr<const pcms::FunctionSpace> MakeSplineSpace()
{
  pcms::UniformGrid<2> grid;
  grid.bot_left = {0.0, -1.0};
  grid.edge_length = {2.0, 2.0};
  grid.divisions = {16, 16};
  return pcms::CreateSplineSpace(
    pcms::CreateUniformGridGeometry(grid, pcms::csys::Cartesian::Create(2)));
}

} // namespace

TEST_CASE("Spline gradients")
{
  auto space = MakeSplineSpace();
  auto f = space->CreateFunction<Real>("f");
  pcms::test::SetField(
    f, OMEGA_H_LAMBDA(Real x, Real y) { return Cubic(x, y); });

  const std::vector<Real> pts = {0.3, -0.7, 1.1, 0.2, 1.7, 0.9, 0.55, 0.45};
  const size_t n = pts.size() / 2;
  const auto bound = pcms::test::Bind(*space->GetGeometry(), pts);
  auto values = space->CreatePointEvaluator<Real>(bound);
  auto gradient = space->CreatePointEvaluator<Real>(
    bound, pcms::DifferentialOperator::Gradient);

  // The gradient is a physical covector on the space's coordinate system.
  const ValueBasis basis = gradient->OutputBasis(f.GetData().GetValueBasis());
  REQUIRE(basis.variance == pcms::values::Covector);
  REQUIRE(basis.scaling == pcms::ComponentScaling::Physical);
  REQUIRE(
    pcms::SameCoordinateSystem(basis.system, space->GetCoordinateSystem()));

  const auto grad = pcms::test::EvaluateToHost(*gradient, f, basis, 2);
  const auto value =
    pcms::test::EvaluateToHost(*values, f, f.GetData().GetValueBasis());
  for (size_t i = 0; i < n; ++i) {
    const Real x = pts[2 * i], y = pts[2 * i + 1];
    CAPTURE(i, x, y);
    REQUIRE_THAT(value(i, 0), WithinAbs(Cubic(x, y), pcms::test::EvalTol));
    REQUIRE_THAT(grad(i, 0), WithinAbs(2.0 * x + 3.0 * y, pcms::test::EvalTol));
    REQUIRE_THAT(grad(i, 1),
                 WithinAbs(3.0 * x - 3.0 * y * y, pcms::test::EvalTol));
  }

  // A gradient written into a view of its shape tagged as values is rejected.
  Kokkos::View<Real**, DeviceMemorySpace> wrong("wrong", n, 2);
  REQUIRE_THROWS_WITH(gradient->Evaluate(f, pcms::test::TagLike(f, wrong)),
                      ContainsSubstring("does not match the basis"));
  REQUIRE_THROWS_WITH(space->CreatePointEvaluator<Real>(
                        bound, pcms::DifferentialOperator::Hessian),
                      ContainsSubstring("Hessian is not supported"));
}

TEST_CASE("Evaluators reject fields of another function space")
{
  // Two spaces over the same grid still have distinct layouts.
  auto space = MakeSplineSpace();
  auto other = MakeSplineSpace();
  auto f = space->CreateFunction<Real>("f");
  auto g = other->CreateFunction<Real>("g");
  pcms::test::SetField(
    f, OMEGA_H_LAMBDA(Real x, Real y) { return Cubic(x, y); });
  pcms::test::SetField(
    g, OMEGA_H_LAMBDA(Real x, Real y) { return Cubic(x, y); });

  const std::vector<Real> pts = {0.3, -0.7, 1.1, 0.2};
  const auto bound = pcms::test::Bind(*space->GetGeometry(), pts);
  auto values = space->CreatePointEvaluator<Real>(bound);
  auto gradient = pcms::Gradient(space)->CreatePointEvaluator<Real>(bound);

  const auto value_basis = f.GetData().GetValueBasis();
  const auto gradient_basis = gradient->OutputBasis(value_basis);
  REQUIRE_NOTHROW(pcms::test::EvaluateToHost(*values, f, value_basis));
  REQUIRE_NOTHROW(pcms::test::EvaluateToHost(*gradient, f, gradient_basis, 2));
  REQUIRE_THROWS_WITH(pcms::test::EvaluateToHost(*values, g, value_basis),
                      ContainsSubstring("not a field of the function space"));
  REQUIRE_THROWS_WITH(
    pcms::test::EvaluateToHost(*gradient, g, gradient_basis, 2),
    ContainsSubstring("not a field of the function space"));
}

TEST_CASE("Gradients fill the out-of-bounds rows")
{
  auto space = MakeSplineSpace();
  auto f = space->CreateFunction<Real>("f");
  pcms::test::SetField(
    f, OMEGA_H_LAMBDA(Real x, Real y) { return Cubic(x, y); });
  const std::vector<Real> pts = {0.5, 0.5, 3.0, 0.0};
  const auto bound = pcms::test::Bind(
    *space->GetGeometry(), pts,
    pcms::OutOfBoundsPolicy{pcms::OutOfBoundsMode::FILL, -7.0});
  auto gradient = space->CreatePointEvaluator<Real>(
    bound, pcms::DifferentialOperator::Gradient);
  const auto grad = pcms::test::EvaluateToHost(
    *gradient, f, gradient->OutputBasis(f.GetData().GetValueBasis()), 2);
  REQUIRE_THAT(grad(0, 0),
               WithinAbs(2.0 * 0.5 + 3.0 * 0.5, pcms::test::EvalTol));
  REQUIRE(grad(1, 0) == -7.0);
  REQUIRE(grad(1, 1) == -7.0);
}

TEST_CASE("Interpolating a field and its gradient over one localization")
{
  std::shared_ptr<const pcms::FunctionSpace> source = MakeSplineSpace();
  auto f = source->CreateFunction<Real>("f");
  pcms::test::SetField(
    f, OMEGA_H_LAMBDA(Real x, Real y) { return Cubic(x, y); });

  // Targets on a coarser grid strictly inside the source grid.
  pcms::UniformGrid<2> grid;
  grid.bot_left = {0.2, -0.8};
  grid.edge_length = {1.5, 1.4};
  grid.divisions = {5, 7};
  auto value_space = pcms::CreateLagrangeSpace(
    pcms::CreateUniformGridGeometry(grid, source->GetCoordinateSystem()), 1, 1);
  auto gradient_space = pcms::CreateLagrangeSpace(
    pcms::CreateUniformGridGeometry(grid, source->GetCoordinateSystem()), 1, 2);
  auto value = value_space->CreateFunction<Real>("f");
  auto gradient =
    gradient_space->CreateFunction<Real>("grad f", pcms::values::Covector);

  const auto value_layout = value_space->GetLayout();
  const auto bound =
    source->GetGeometry()->Bind(value_layout->GetDOFHolderCoordinates());
  const auto grad_source = pcms::Gradient(source);
  pcms::Interpolator<Real> values(*source, *value_space, bound);
  pcms::Interpolator<Real> gradients(*grad_source, *gradient_space, bound);
  values.Apply(f, value);
  gradients.Apply(f, gradient);

  const auto coords = pcms::test::CopyCoordinatesToHost(
    value_space->GetLayout()->GetDOFHolderCoordinates().GetValues());
  const auto v = value.GetDOFHolderDataHost().GetValues();
  const auto g = gradient.GetDOFHolderDataHost().GetValues();
  for (size_t i = 0; i < coords.extent(0); ++i) {
    const Real x = coords(i, 0), y = coords(i, 1);
    CAPTURE(i, x, y);
    REQUIRE_THAT(v(i, 0), WithinAbs(Cubic(x, y), pcms::test::EvalTol));
    REQUIRE_THAT(g(i, 0), WithinAbs(2.0 * x + 3.0 * y, pcms::test::EvalTol));
    REQUIRE_THAT(g(i, 1),
                 WithinAbs(3.0 * x - 3.0 * y * y, pcms::test::EvalTol));
  }

  // A gradient source binds through the source, so it can also bind itself.
  pcms::Interpolator<Real> self_bound(*grad_source, *gradient_space);
  self_bound.Apply(f, gradient);
  REQUIRE_THAT(gradient.GetDOFHolderDataHost().GetValues()(0, 0),
               WithinAbs(g(0, 0), pcms::test::ExactTol));

  // Values are not written into a covector target.
  REQUIRE_THROWS_WITH(
    pcms::Interpolator<Real>(*source, *gradient_space, bound)
      .Apply(f, gradient),
    ContainsSubstring("differs from the target's declared basis"));
  // The bound points must be the target's DOF holders.
  pcms::UniformGrid<2> other_grid = grid;
  other_grid.divisions = {4, 7};
  auto other_space = pcms::CreateLagrangeSpace(
    pcms::CreateUniformGridGeometry(other_grid, source->GetCoordinateSystem()),
    1, 1);
  REQUIRE_THROWS_WITH(pcms::Interpolator<Real>(*source, *other_space, bound),
                      ContainsSubstring("target's DOF-holder sites"));
  // Second derivatives are not available from the spline yet.
  REQUIRE_THROWS_WITH(
    pcms::Gradient(grad_source),
    ContainsSubstring(
      "Apply: Gradient is not supported by the wrapped source"));
}

TEST_CASE("Unsupported operators are refused by name")
{
  auto space = MakeSplineSpace();
  const auto bound = pcms::test::Bind(*space->GetGeometry(), {0.3, -0.7});
  REQUIRE_FALSE(space->Supports(pcms::DifferentialOperator::Curl));
  REQUIRE_THROWS_WITH(
    space->CreatePointEvaluator<Real>(bound, pcms::DifferentialOperator::Curl),
    ContainsSubstring("Curl is not supported"));
  // an operator source only supports Value
  auto grad = pcms::Gradient(space);
  REQUIRE(grad->Supports(pcms::DifferentialOperator::Value));
  REQUIRE_FALSE(grad->Supports(pcms::DifferentialOperator::Gradient));
  // an operator the wrapped source lacks is refused by name when applied
  REQUIRE_THROWS_WITH(
    pcms::Curl(space),
    ContainsSubstring("Apply: Curl is not supported by the wrapped source"));
}

TEST_CASE("Gradients of vectors need Cartesian coordinates")
{
  // GradientBasis directly: rank-1 stored basis
  auto cyl = pcms::csys::CylindricalRThetaZ::Create();
  pcms::ValueBasis v{cyl, pcms::ComponentScaling::Physical,
                     pcms::values::Vector};
  REQUIRE_THROWS_WITH(pcms::GradientBasis(v, cyl),
                      ContainsSubstring("connection"));
  auto cart = pcms::csys::Cartesian::Create(3);
  pcms::ValueBasis w{cart, pcms::ComponentScaling::Physical,
                     pcms::values::Vector};
  REQUIRE(pcms::GradientBasis(w, cart).Rank() == 2);
  // the parity of the stored basis is kept
  w.parity = pcms::Parity::Odd;
  REQUIRE(pcms::GradientBasis(w, cart).parity == pcms::Parity::Odd);
  pcms::ValueBasis s{cart, pcms::ComponentScaling::Physical,
                     pcms::values::Scalar, pcms::Parity::Odd};
  REQUIRE(pcms::GradientBasis(s, cart).parity == pcms::Parity::Odd);
}
