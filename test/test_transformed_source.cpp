#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <Kokkos_Core.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_library.hpp>
#include <pcms/field/coordinate_map.hpp>
#include <pcms/field/coordinate_systems/cartesian.hpp>
#include <pcms/field/coordinate_systems/cylindrical.hpp>
#include <pcms/field/function_space/lagrange.h>
#include <pcms/field/transformed_source.hpp>
#include <pcms/transfer/interpolator.h>
#include "field_test_utils.h"
#include <cmath>
#include <vector>

using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;
using pcms::ComponentScaling;
using pcms::DeviceMemorySpace;
using pcms::EvaluationRequest;
using pcms::GO;
using pcms::HostMemorySpace;
using pcms::LagrangeFunctionSpace;
using pcms::OutOfBoundsMode;
using pcms::OutOfBoundsPolicy;
using pcms::Real;
using pcms::SameValueBasis;
using pcms::TransformedSource;
using pcms::ValueBasis;
using pcms::ValueView;
namespace values = pcms::values;

namespace
{

constexpr double kTol = 1e-10;

std::shared_ptr<LagrangeFunctionSpace> BuildCylindricalSource(
  Omega_h::Library& lib, Omega_h::Mesh& mesh, int num_components)
{
  mesh = Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 2.0, 1.6, 1.1, 8, 8,
                            6, false);
  return LagrangeFunctionSpace::FromMesh(
    mesh, 1, num_components, pcms::csys::CylindricalRThetaZ::Create(), "global",
    LagrangeFunctionSpace::Backend::OmegaH);
}

std::shared_ptr<LagrangeFunctionSpace> BuildCartesianTarget(
  Omega_h::Library& lib, Omega_h::Mesh& mesh, int num_components,
  double x_extent = 1.0)
{
  mesh = Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, x_extent, 1.0, 1.0, 4,
                            4, 4, false);
  return LagrangeFunctionSpace::FromMesh(
    mesh, 1, num_components, pcms::csys::Cartesian::Create(3), "global",
    LagrangeFunctionSpace::Backend::OmegaH);
}

ValueBasis CylindricalVector()
{
  return ValueBasis{pcms::csys::CylindricalRThetaZ::Create(),
                    ComponentScaling::Physical, values::Vector};
}

ValueBasis CartesianVector()
{
  return ValueBasis{pcms::csys::Cartesian::Create(3),
                    ComponentScaling::Physical, values::Vector};
}

Kokkos::View<Real**, HostMemorySpace> EvaluateToHost(
  const pcms::PointEvaluator<Real>& evaluator, const pcms::Field<Real>& field,
  const ValueBasis& tag, int n, int width)
{
  Kokkos::View<Real**, DeviceMemorySpace> out("out", n, width);
  evaluator.Evaluate(
    field, ValueView<Real, DeviceMemorySpace>(tag, pcms::MakeRank2View(out)));
  Kokkos::View<Real**, HostMemorySpace> host("host", n, width);
  pcms::DeepCopyMismatchLayouts(host, out);
  return host;
}

} // namespace

TEST_CASE("TransformedSource: scalar evaluation maps the query points")
{
  auto lib = Omega_h::Library{};
  Omega_h::Mesh src_mesh(&lib), tgt_mesh(&lib);
  auto src = BuildCylindricalSource(lib, src_mesh, 1);
  auto tgt = BuildCartesianTarget(lib, tgt_mesh, 1);
  auto f = src->CreateFunction<Real>("f");
  pcms::test::SetField(
    f, OMEGA_H_LAMBDA(Real r, Real theta, Real z) {
      return 2.0 * r + 3.0 * theta - z;
    });

  TransformedSource src_as_cartesian(
    src, std::make_shared<pcms::CartesianToCylindrical>());
  REQUIRE(pcms::SameCoordinateSystem(src_as_cartesian.GetCoordinateSystem(),
                                     tgt->GetCoordinateSystem()));

  const auto pts = tgt->GetLayout()->GetDOFHolderCoordinates();
  const int n = static_cast<int>(pts.GetValues().extent(0));
  auto evaluator = src_as_cartesian.CreatePointEvaluator<Real>(
    EvaluationRequest::FromCoordinates(pts));
  REQUIRE(SameValueBasis(evaluator->OutputBasis(ValueBasis{}), ValueBasis{}));

  const auto coords = pcms::test::CopyCoordinatesToHost(pts.GetValues());
  const auto result = EvaluateToHost(*evaluator, f, ValueBasis{}, n, 1);
  for (int i = 0; i < n; ++i) {
    const Real x = coords(i, 0), y = coords(i, 1), z = coords(i, 2);
    CAPTURE(i);
    REQUIRE_THAT(result(i, 0), WithinAbs(2.0 * std::sqrt(x * x + y * y) +
                                           3.0 * std::atan2(y, x) - z,
                                         kTol));
  }
}

TEST_CASE("TransformedSource: vector evaluation rotates into the target basis")
{
  auto lib = Omega_h::Library{};
  Omega_h::Mesh src_mesh(&lib), tgt_mesh(&lib);
  auto src = BuildCylindricalSource(lib, src_mesh, 3);
  auto tgt = BuildCartesianTarget(lib, tgt_mesh, 3);
  auto b =
    src->CreateFunction<Real>("b", values::Vector, ComponentScaling::Physical);
  pcms::test::SetFieldComponents(b, [](Real r, Real, Real z, Real* out) {
    out[0] = r;
    out[1] = 0.0;
    out[2] = z;
  });

  TransformedSource src_as_cartesian(
    src, std::make_shared<pcms::CartesianToCylindrical>());
  const auto pts = tgt->GetLayout()->GetDOFHolderCoordinates();
  const int n = static_cast<int>(pts.GetValues().extent(0));
  auto evaluator = src_as_cartesian.CreatePointEvaluator<Real>(
    EvaluationRequest::FromCoordinates(pts));
  REQUIRE(SameValueBasis(evaluator->OutputBasis(CylindricalVector()),
                         CartesianVector()));
  REQUIRE_THROWS_WITH(evaluator->OutputBasis(CartesianVector()),
                      ContainsSubstring("cannot re-express"));

  const auto coords = pcms::test::CopyCoordinatesToHost(pts.GetValues());
  const auto result = EvaluateToHost(*evaluator, b, CartesianVector(), n, 3);
  for (int i = 0; i < n; ++i) {
    CAPTURE(i);
    REQUIRE_THAT(result(i, 0), WithinAbs(coords(i, 0), kTol));
    REQUIRE_THAT(result(i, 1), WithinAbs(coords(i, 1), kTol));
    REQUIRE_THAT(result(i, 2), WithinAbs(coords(i, 2), kTol));
  }

  SECTION("evaluate gates the field's stored basis and the output tag")
  {
    auto borrowed = src->CreateFunction<Real>("c", values::Vector,
                                              pcms::csys::Cartesian::Create(3));
    Kokkos::View<Real**, DeviceMemorySpace> out("out", n, 3);
    REQUIRE_THROWS_WITH(
      evaluator->Evaluate(
        borrowed, ValueView<Real, DeviceMemorySpace>(CartesianVector(),
                                                     pcms::MakeRank2View(out))),
      ContainsSubstring("cannot re-express"));
    REQUIRE_THROWS_WITH(
      evaluator->Evaluate(b, ValueView<Real, DeviceMemorySpace>(
                               CylindricalVector(), pcms::MakeRank2View(out))),
      ContainsSubstring("does not match"));
  }
}

TEST_CASE("TransformedSource: FILL rows survive the rotation untouched")
{
  auto lib = Omega_h::Library{};
  Omega_h::Mesh src_mesh(&lib), tgt_mesh(&lib);
  auto src = BuildCylindricalSource(lib, src_mesh, 3);
  auto tgt = BuildCartesianTarget(lib, tgt_mesh, 3, 3.0);
  auto b =
    src->CreateFunction<Real>("b", values::Vector, ComponentScaling::Physical);
  pcms::test::SetFieldComponents(b, [](Real r, Real, Real z, Real* out) {
    out[0] = r;
    out[1] = 0.0;
    out[2] = z;
  });

  const Real fill = -7.0;
  TransformedSource src_as_cartesian(
    src, std::make_shared<pcms::CartesianToCylindrical>());
  const auto pts = tgt->GetLayout()->GetDOFHolderCoordinates();
  const int n = static_cast<int>(pts.GetValues().extent(0));
  auto evaluator = src_as_cartesian.CreatePointEvaluator<Real>(
    EvaluationRequest::FromCoordinates(
      pts, OutOfBoundsPolicy{OutOfBoundsMode::FILL, fill}));

  const auto coords = pcms::test::CopyCoordinatesToHost(pts.GetValues());
  const auto result = EvaluateToHost(*evaluator, b, CartesianVector(), n, 3);
  int num_filled = 0;
  for (int i = 0; i < n; ++i) {
    const Real x = coords(i, 0), y = coords(i, 1);
    CAPTURE(i);
    if (std::sqrt(x * x + y * y) > 2.0) {
      ++num_filled;
      REQUIRE(result(i, 0) == fill);
      REQUIRE(result(i, 1) == fill);
      REQUIRE(result(i, 2) == fill);
    } else {
      REQUIRE_THAT(result(i, 0), WithinAbs(x, kTol));
      REQUIRE_THAT(result(i, 1), WithinAbs(y, kTol));
      REQUIRE_THAT(result(i, 2), WithinAbs(coords(i, 2), kTol));
    }
  }
  REQUIRE(num_filled > 0);
}

TEST_CASE("TransformedSource: construction errors")
{
  auto lib = Omega_h::Library{};
  Omega_h::Mesh src_mesh(&lib);
  auto src = BuildCylindricalSource(lib, src_mesh, 3);

  REQUIRE_THROWS_WITH(TransformedSource(src, nullptr),
                      ContainsSubstring("must not be null"));
  REQUIRE_THROWS_WITH(
    TransformedSource(src, std::make_shared<pcms::CylindricalToCartesian>()),
    ContainsSubstring("source's coordinate system"));
}

TEST_CASE("TransformedSource: evaluator creation and evaluation errors")
{
  auto lib = Omega_h::Library{};
  Omega_h::Mesh src_mesh(&lib);
  auto src = BuildCylindricalSource(lib, src_mesh, 3);
  auto cart_to_cyl = std::make_shared<pcms::CartesianToCylindrical>();
  const std::vector<Real> cart_pts = {0.5, 0.5, 0.5, 0.2, 0.9, 0.1};
  auto query = pcms::test::CreateDeviceCoordinateView(
    cart_pts, pcms::csys::Cartesian::Create(3), 3);
  TransformedSource src_as_cartesian(src, cart_to_cyl);

  SECTION("query points must be in the transformed system")
  {
    auto cyl_query = pcms::test::CreateDeviceCoordinateView(
      cart_pts, pcms::csys::CylindricalRThetaZ::Create(), 3);
    REQUIRE_THROWS_WITH(
      src_as_cartesian.CreatePointEvaluator<Real>(
        EvaluationRequest::FromCoordinates(cyl_query.coordinate_view)),
      ContainsSubstring("not the factory's coordinate system"));
  }

  SECTION("non-Real evaluation surfaces the source backend's own limit")
  {
    // No in-tree backend evaluates integer fields at points, so a
    // TransformedSource simply forwards the inner factory's refusal; the
    // T != Real guard in TransformedPointEvaluator::Evaluate is defensive.
    REQUIRE_THROWS_WITH(
      src_as_cartesian.CreatePointEvaluator<GO>(
        EvaluationRequest::FromCoordinates(query.coordinate_view)),
      ContainsSubstring("only supports double"));
  }
}

TEST_CASE("TransformedSource: one source view serves many fields and targets")
{
  auto lib = Omega_h::Library{};
  Omega_h::Mesh src_mesh(&lib), tgt_a_mesh(&lib), tgt_b_mesh(&lib);
  auto src = BuildCylindricalSource(lib, src_mesh, 3);
  auto tgt_a = BuildCartesianTarget(lib, tgt_a_mesh, 3);
  auto tgt_b = BuildCartesianTarget(lib, tgt_b_mesh, 3, 0.8);
  auto b =
    src->CreateFunction<Real>("b", values::Vector, ComponentScaling::Physical);
  auto e =
    src->CreateFunction<Real>("e", values::Vector, ComponentScaling::Physical);
  pcms::test::SetFieldComponents(b, [](Real r, Real, Real z, Real* out) {
    out[0] = r;
    out[1] = 0.0;
    out[2] = z;
  });
  pcms::test::SetFieldComponents(e, [](Real r, Real, Real z, Real* out) {
    out[0] = 2.0 * r;
    out[1] = 0.0;
    out[2] = -z;
  });

  // Built from the source and the map alone, so it is reused for
  // every field on that space and every target space.
  TransformedSource src_as_cartesian(
    src, std::make_shared<pcms::CartesianToCylindrical>());
  pcms::Interpolator<Real> to_a(src_as_cartesian, *tgt_a);
  pcms::Interpolator<Real> to_b(src_as_cartesian, *tgt_b);

  auto b_a = tgt_a->CreateFunction<Real>("b", values::Vector);
  auto e_a = tgt_a->CreateFunction<Real>("e", values::Vector);
  auto b_b = tgt_b->CreateFunction<Real>("b", values::Vector);
  to_a.Apply(b, b_a);
  to_a.Apply(e, e_a);
  to_b.Apply(b, b_b);

  const auto coords = pcms::test::CopyCoordinatesToHost(
    tgt_a->GetLayout()->GetDOFHolderCoordinates().GetValues());
  const auto rb = b_a.GetDOFHolderDataHost();
  const auto re = e_a.GetDOFHolderDataHost();
  for (int i = 0; i < static_cast<int>(coords.extent(0)); ++i) {
    CAPTURE(i);
    REQUIRE_THAT(rb(i, 0), WithinAbs(coords(i, 0), kTol));
    REQUIRE_THAT(re(i, 0), WithinAbs(2.0 * coords(i, 0), kTol));
    REQUIRE_THAT(re(i, 2), WithinAbs(-coords(i, 2), kTol));
  }
}

TEST_CASE("TransformedSource: nesting composes maps and rotations")
{
  auto lib = Omega_h::Library{};
  Omega_h::Mesh src_mesh(&lib);
  auto src = BuildCylindricalSource(lib, src_mesh, 3);
  auto b =
    src->CreateFunction<Real>("b", values::Vector, ComponentScaling::Physical);
  pcms::test::SetFieldComponents(b, [](Real r, Real, Real z, Real* out) {
    out[0] = r;
    out[1] = 0.0;
    out[2] = z;
  });

  auto src_as_cartesian = std::make_shared<TransformedSource>(
    src, std::make_shared<pcms::CartesianToCylindrical>());
  TransformedSource src_as_cylindrical(
    src_as_cartesian, std::make_shared<pcms::CylindricalToCartesian>());

  const std::vector<Real> cyl_pts = {0.7, 0.3, 0.4, 1.0, 0.9,
                                     0.9, 1.5, 1.2, 0.2};
  auto query = pcms::test::CreateDeviceCoordinateView(
    cyl_pts, pcms::csys::CylindricalRThetaZ::Create(), 3);
  const int n = 3;
  auto evaluator = src_as_cylindrical.CreatePointEvaluator<Real>(
    EvaluationRequest::FromCoordinates(query.coordinate_view));
  const auto result = EvaluateToHost(*evaluator, b, CylindricalVector(), n, 3);
  for (int i = 0; i < n; ++i) {
    CAPTURE(i);
    REQUIRE_THAT(result(i, 0), WithinAbs(cyl_pts[3 * i], kTol));
    REQUIRE_THAT(result(i, 1), WithinAbs(0.0, kTol));
    REQUIRE_THAT(result(i, 2), WithinAbs(cyl_pts[3 * i + 2], kTol));
  }
}

TEST_CASE("Interpolator over a TransformedSource transfers across systems")
{
  auto lib = Omega_h::Library{};
  Omega_h::Mesh src_mesh(&lib), tgt_mesh(&lib);
  auto src = BuildCylindricalSource(lib, src_mesh, 3);
  auto tgt = BuildCartesianTarget(lib, tgt_mesh, 3);
  auto b_src =
    src->CreateFunction<Real>("b", values::Vector, ComponentScaling::Physical);
  auto b_tgt = tgt->CreateFunction<Real>("b", values::Vector);
  pcms::test::SetFieldComponents(b_src, [](Real r, Real, Real z, Real* out) {
    out[0] = r;
    out[1] = 0.0;
    out[2] = z;
  });

  TransformedSource src_as_cartesian(
    src, std::make_shared<pcms::CartesianToCylindrical>());
  pcms::Interpolator<Real> op(src_as_cartesian, *tgt);
  op.Apply(b_src, b_tgt);

  const auto coords = pcms::test::CopyCoordinatesToHost(
    tgt->GetLayout()->GetDOFHolderCoordinates().GetValues());
  const auto result = b_tgt.GetDOFHolderDataHost();
  const int n = static_cast<int>(coords.extent(0));
  for (int i = 0; i < n; ++i) {
    CAPTURE(i);
    REQUIRE_THAT(result(i, 0), WithinAbs(coords(i, 0), kTol));
    REQUIRE_THAT(result(i, 1), WithinAbs(coords(i, 1), kTol));
    REQUIRE_THAT(result(i, 2), WithinAbs(coords(i, 2), kTol));
  }

  SECTION("the target must declare the basis the evaluator writes")
  {
    auto wrong = tgt->CreateFunction<Real>(
      "w", values::Vector, pcms::csys::CylindricalRThetaZ::Create(),
      ComponentScaling::Physical);
    REQUIRE_THROWS_WITH(op.Apply(b_src, wrong),
                        ContainsSubstring("target's declared basis"));
  }

  SECTION("a plain Interpolator refuses to bridge coordinate systems")
  {
    REQUIRE_THROWS_WITH(pcms::Interpolator<Real>(*src, *tgt),
                        ContainsSubstring("coordinate system"));
  }
}
