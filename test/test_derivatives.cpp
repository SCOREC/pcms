#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <Kokkos_Core.hpp>
#include <pcms/utility/coordinate_systems/cartesian.hpp>
#include <pcms/utility/coordinate_systems/cylindrical.hpp>
#include <pcms/field/function_spaces.hpp>
#include <pcms/backends/uniform_grid/geometry.hpp>
#include <pcms/field/gradient_source.hpp>
#include <pcms/transfer/copy.hpp>
#include <pcms/transfer/interpolator.hpp>
#include <pcms/utility/coordinate_systems/parameter.hpp>
#include "field_test_utils.hpp"
#include <random>
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
    pcms::uniform_grid::CreateGeometry(grid, pcms::csys::CreateCartesian(2)));
}

} // namespace

TEST_CASE("Spline gradients")
{
  auto space = MakeSplineSpace();
  auto f = space->CreateFunction<Real>("f");
  pcms::test::SetField(
    f, KOKKOS_LAMBDA(Real x, Real y) { return Cubic(x, y); });

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
  REQUIRE(pcms::detail::SameCoordinateSystem(basis.system,
                                             space->GetCoordinateSystem()));

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
    f, KOKKOS_LAMBDA(Real x, Real y) { return Cubic(x, y); });
  pcms::test::SetField(
    g, KOKKOS_LAMBDA(Real x, Real y) { return Cubic(x, y); });

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
    f, KOKKOS_LAMBDA(Real x, Real y) { return Cubic(x, y); });
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
    f, KOKKOS_LAMBDA(Real x, Real y) { return Cubic(x, y); });

  // Targets on a coarser grid strictly inside the source grid.
  pcms::UniformGrid<2> grid;
  grid.bot_left = {0.2, -0.8};
  grid.edge_length = {1.5, 1.4};
  grid.divisions = {5, 7};
  auto value_space = pcms::CreateLagrangeSpace(
    pcms::uniform_grid::CreateGeometry(grid, source->GetCoordinateSystem()), 1,
    1);
  auto gradient_space = pcms::CreateLagrangeSpace(
    pcms::uniform_grid::CreateGeometry(grid, source->GetCoordinateSystem()), 1,
    2);
  auto value = value_space->CreateFunction<Real>("f");
  auto gradient = gradient_space->CreateFunction<Real>(
    "grad f", {.variance = pcms::values::Covector});

  const auto value_layout = value_space->GetLayout();
  const auto bound =
    source->GetGeometry()->Bind(value_layout->GetDOFHolderCoordinates());
  const auto grad_source = pcms::Gradient(source);
  const auto values =
    pcms::CreateInterpolator<Real>(source, value_space, bound);
  const auto gradients =
    pcms::CreateInterpolator<Real>(grad_source, gradient_space, bound);
  values->Apply(f, value);
  gradients->Apply(f, gradient);

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
  const auto self_bound =
    pcms::CreateInterpolator<Real>(grad_source, gradient_space);
  self_bound->Apply(f, gradient);
  REQUIRE_THAT(gradient.GetDOFHolderDataHost().GetValues()(0, 0),
               WithinAbs(g(0, 0), pcms::test::ExactTol));

  // Values are not written into a covector target.
  REQUIRE_THROWS_WITH(
    pcms::CreateInterpolator<Real>(source, gradient_space, bound)
      ->Apply(f, gradient),
    ContainsSubstring("differs from the target's declared basis"));
  // The bound points must be the target's DOF holders.
  pcms::UniformGrid<2> other_grid = grid;
  other_grid.divisions = {4, 7};
  auto other_space =
    pcms::CreateLagrangeSpace(pcms::uniform_grid::CreateGeometry(
                                other_grid, source->GetCoordinateSystem()),
                              1, 1);
  REQUIRE_THROWS_WITH(
    pcms::CreateInterpolator<Real>(source, other_space, bound),
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
  auto cyl = pcms::csys::CreateCylindricalRThetaZ();
  pcms::ValueBasis v{cyl, pcms::ComponentScaling::Physical,
                     pcms::values::Vector};
  REQUIRE_THROWS_WITH(pcms::detail::GradientBasis(v, cyl),
                      ContainsSubstring("connection"));
  auto cart = pcms::csys::CreateCartesian(3);
  pcms::ValueBasis w{cart, pcms::ComponentScaling::Physical,
                     pcms::values::Vector};
  REQUIRE(pcms::detail::GradientBasis(w, cart).Rank() == 2);
  // the parity of the stored basis is kept
  w.parity = pcms::Parity::Odd;
  REQUIRE(pcms::detail::GradientBasis(w, cart).parity == pcms::Parity::Odd);
  pcms::ValueBasis s{cart, pcms::ComponentScaling::Physical,
                     pcms::values::Scalar, pcms::Parity::Odd};
  REQUIRE(pcms::detail::GradientBasis(s, cart).parity == pcms::Parity::Odd);
}

namespace
{

// P1 Lagrange spaces in (R, theta, Z): a source grid, and a target grid
// strictly inside it.
std::shared_ptr<const pcms::FunctionSpace> CylindricalSpace(bool target,
                                                            int components)
{
  pcms::UniformGrid<3> grid;
  grid.bot_left = {1.0, 0.0, -1.0};
  grid.edge_length = {1.0, 1.0, 2.0};
  grid.divisions = {4, 4, 4};
  if (target) {
    grid.bot_left = {1.2, 0.1, -0.8};
    grid.edge_length = {0.6, 0.8, 1.6};
    grid.divisions = {3, 3, 3};
  }
  return pcms::CreateLagrangeSpace(
    pcms::uniform_grid::CreateGeometry(grid,
                                       pcms::csys::CreateCylindricalRThetaZ()),
    1, components);
}

void SetRandom(pcms::Function<Real>& f, std::mt19937& rng)
{
  std::uniform_real_distribution<Real> uniform(-1.0, 1.0);
  const auto host = f.GetDOFHolderDataHost();
  Kokkos::View<Real**, pcms::HostMemorySpace> values("random", host.extent(0),
                                                     host.extent(1));
  for (size_t i = 0; i < values.extent(0); ++i) {
    for (size_t c = 0; c < values.extent(1); ++c) {
      values(i, c) = uniform(rng);
    }
  }
  f.SetDOFHolderDataHost(pcms::ValueView<const Real, pcms::HostMemorySpace>(
    host.GetBasis(), pcms::MakeConstRank2View(values)));
}

Real Dot(const pcms::Function<Real>& a, const pcms::Function<Real>& b)
{
  const auto x = a.GetDOFHolderDataHost().GetValues();
  const auto y = b.GetDOFHolderDataHost().GetValues();
  Real sum = 0;
  for (size_t i = 0; i < x.extent(0); ++i) {
    for (size_t c = 0; c < x.extent(1); ++c) {
      sum += x(i, c) * y(i, c);
    }
  }
  return sum;
}

} // namespace

TEST_CASE("A cylindrical gradient interpolated into a physical covector",
          "[gradient][scaling]")
{
  const auto source = CylindricalSpace(false, 1);
  const auto target = CylindricalSpace(true, 3);
  auto f = source->CreateFunction<Real>("f");
  pcms::test::SetFieldComponents(f, [](Real r, Real theta, Real z, Real* out) {
    out[0] = 2.0 * r + 3.0 * theta - z;
  });
  auto gradient = target->CreateFunction<Real>(
    "grad f", {.scaling = pcms::ComponentScaling::Physical,
               .variance = pcms::values::Covector});
  pcms::CreateInterpolator<Real>(pcms::Gradient(source), target)
    ->Apply(f, gradient);
  const auto coords = pcms::test::CopyCoordinatesToHost(
    target->GetLayout()->GetDOFHolderCoordinates().GetValues());
  const auto g = gradient.GetDOFHolderDataHost().GetValues();
  for (size_t i = 0; i < coords.extent(0); ++i) {
    CAPTURE(i);
    REQUIRE_THAT(g(i, 0), WithinAbs(2.0, pcms::test::EvalTol));
    REQUIRE_THAT(g(i, 1), WithinAbs(3.0 / coords(i, 0), pcms::test::EvalTol));
    REQUIRE_THAT(g(i, 2), WithinAbs(-1.0, pcms::test::EvalTol));
  }
}

TEST_CASE("An interpolator that converts the scaling keeps the fill value",
          "[gradient][scaling]")
{
  const auto source = CylindricalSpace(false, 1);
  pcms::UniformGrid<3> grid;
  grid.bot_left = {1.2, 0.1, -0.8};
  grid.edge_length = {1.8, 0.8, 1.6};
  grid.divisions = {3, 3, 3};
  const auto target =
    pcms::CreateLagrangeSpace(pcms::uniform_grid::CreateGeometry(
                                grid, pcms::csys::CreateCylindricalRThetaZ()),
                              1, 3);
  auto f = source->CreateFunction<Real>("f");
  pcms::test::SetFieldComponents(f, [](Real r, Real theta, Real z, Real* out) {
    out[0] = 2.0 * r + 3.0 * theta - z;
  });
  auto gradient = target->CreateFunction<Real>(
    "grad f", {.scaling = pcms::ComponentScaling::Physical,
               .variance = pcms::values::Covector});
  pcms::CreateInterpolator<Real>(
    pcms::Gradient(source), target,
    pcms::OutOfBoundsPolicy{pcms::OutOfBoundsMode::FILL, -7.0})
    ->Apply(f, gradient);
  const auto coords = pcms::test::CopyCoordinatesToHost(
    target->GetLayout()->GetDOFHolderCoordinates().GetValues());
  const auto g = gradient.GetDOFHolderDataHost().GetValues();
  int filled = 0;
  for (size_t i = 0; i < coords.extent(0); ++i) {
    if (coords(i, 0) > 2.0) {
      ++filled;
      for (int c = 0; c < 3; ++c) {
        CHECK(g(i, c) == -7.0);
      }
    }
  }
  CHECK(filled > 0);
}

TEST_CASE("An interpolator that converts the scaling has the transpose of its "
          "evaluation",
          "[scaling][adjoint]")
{
  const auto source = CylindricalSpace(false, 3);
  const auto target = CylindricalSpace(true, 3);
  const pcms::ValueBasis natural{.scaling = pcms::ComponentScaling::Natural,
                                 .variance = pcms::values::Vector};
  const pcms::ValueBasis physical{.scaling = pcms::ComponentScaling::Physical,
                                  .variance = pcms::values::Vector};
  auto x = source->CreateFunction<Real>("x", natural);
  auto functional = source->CreateFunction<Real>("functional", natural);
  auto ax = target->CreateFunction<Real>("Ax", physical);
  auto y = target->CreateFunction<Real>("y", physical);
  std::mt19937 rng(11);
  SetRandom(x, rng);
  SetRandom(y, rng);
  const auto op = pcms::CreateInterpolator<Real>(source, target);
  REQUIRE(op->SupportsTranspose());
  op->Apply(x, ax);
  op->ApplyTranspose(y, functional);
  const Real lhs = Dot(ax, y);
  const Real rhs = Dot(x, functional);
  CHECK_THAT(lhs, WithinAbs(rhs, 1e-12 * (1 + std::abs(rhs))));
}

TEST_CASE("A copy converts natural components to physical ones and back",
          "[scaling]")
{
  const auto space = CylindricalSpace(true, 3);
  auto natural = space->CreateFunction<Real>(
    "n", {.scaling = pcms::ComponentScaling::Natural,
          .variance = pcms::values::Vector});
  auto physical = space->CreateFunction<Real>(
    "p", {.scaling = pcms::ComponentScaling::Physical,
          .variance = pcms::values::Vector});
  auto back = space->CreateFunction<Real>(
    "b", {.scaling = pcms::ComponentScaling::Natural,
          .variance = pcms::values::Vector});
  pcms::test::SetFieldComponents(natural, [](Real, Real, Real, Real* out) {
    out[0] = 1.0;
    out[1] = 1.0;
    out[2] = 1.0;
  });
  const auto copy = pcms::CreateCopy<Real>(space, space);
  copy->Apply(natural, physical);
  copy->Apply(physical, back);
  const auto coords = pcms::test::CopyCoordinatesToHost(
    space->GetLayout()->GetDOFHolderCoordinates().GetValues());
  const auto p = physical.GetDOFHolderDataHost().GetValues();
  const auto b = back.GetDOFHolderDataHost().GetValues();
  for (size_t i = 0; i < coords.extent(0); ++i) {
    CAPTURE(i);
    CHECK_THAT(p(i, 0), WithinAbs(1.0, pcms::test::ExactTol));
    CHECK_THAT(p(i, 1), WithinAbs(coords(i, 0), pcms::test::ExactTol));
    CHECK_THAT(p(i, 2), WithinAbs(1.0, pcms::test::ExactTol));
    for (int c = 0; c < 3; ++c) {
      CHECK_THAT(b(i, c), WithinAbs(1.0, pcms::test::ExactTol));
    }
  }
}

TEST_CASE("Converting the scaling on a system without a metric throws",
          "[scaling]")
{
  pcms::UniformGrid<2> grid;
  grid.bot_left = {0.0, 0.0};
  grid.edge_length = {1.0, 1.0};
  grid.divisions = {2, 2};
  const auto space = pcms::CreateLagrangeSpace(
    pcms::uniform_grid::CreateGeometry(
      grid, pcms::csys::CreateParameterSystem("xi", 2)),
    1, 2);
  const auto natural = space->CreateFunction<Real>(
    "n", {.scaling = pcms::ComponentScaling::Natural,
          .variance = pcms::values::Vector});
  auto physical = space->CreateFunction<Real>(
    "p", {.scaling = pcms::ComponentScaling::Physical,
          .variance = pcms::values::Vector});
  CHECK_THROWS_WITH(
    pcms::CreateCopy<Real>(space, space)->Apply(natural, physical),
    ContainsSubstring("which it does not have"));
}
