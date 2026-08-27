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

  TransformedSource view(*src, std::make_shared<pcms::CartesianToCylindrical>(),
                         ValueBasis{}, ValueBasis{});
  REQUIRE(pcms::SameCoordinateSystem(view.GetCoordinateSystem(),
                                     tgt->GetCoordinateSystem()));

  const auto pts = tgt->GetLayout()->GetDOFHolderCoordinates();
  const int n = static_cast<int>(pts.GetValues().extent(0));
  auto evaluator =
    view.CreatePointEvaluator<Real>(EvaluationRequest::FromCoordinates(pts));
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

  TransformedSource view(*src, std::make_shared<pcms::CartesianToCylindrical>(),
                         CylindricalVector(), CartesianVector());
  const auto pts = tgt->GetLayout()->GetDOFHolderCoordinates();
  const int n = static_cast<int>(pts.GetValues().extent(0));
  auto evaluator =
    view.CreatePointEvaluator<Real>(EvaluationRequest::FromCoordinates(pts));
  REQUIRE(SameValueBasis(evaluator->OutputBasis(CylindricalVector()),
                         CartesianVector()));
  REQUIRE_THROWS_WITH(evaluator->OutputBasis(CartesianVector()),
                      ContainsSubstring("not the source basis"));

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
      ContainsSubstring("not the source basis"));
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
  TransformedSource view(*src, std::make_shared<pcms::CartesianToCylindrical>(),
                         CylindricalVector(), CartesianVector());
  const auto pts = tgt->GetLayout()->GetDOFHolderCoordinates();
  const int n = static_cast<int>(pts.GetValues().extent(0));
  auto evaluator =
    view.CreatePointEvaluator<Real>(EvaluationRequest::FromCoordinates(
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
  auto cart_to_cyl = std::make_shared<pcms::CartesianToCylindrical>();

  REQUIRE_THROWS_WITH(
    TransformedSource(*src, nullptr, ValueBasis{}, ValueBasis{}),
    ContainsSubstring("must not be null"));
  REQUIRE_THROWS_WITH(
    TransformedSource(*src, std::make_shared<pcms::CylindricalToCartesian>(),
                      ValueBasis{}, ValueBasis{}),
    ContainsSubstring("source's coordinate system"));
  REQUIRE_THROWS_WITH(
    TransformedSource(*src, cart_to_cyl, CylindricalVector(), ValueBasis{}),
    ContainsSubstring("different ranks"));
  REQUIRE_THROWS_WITH(
    TransformedSource(*src, cart_to_cyl, CylindricalVector(),
                      ValueBasis{pcms::csys::Cartesian::Create(3),
                                 ComponentScaling::Physical, values::Covector}),
    ContainsSubstring("different variances"));
  REQUIRE_THROWS_WITH(
    TransformedSource(*src, cart_to_cyl,
                      ValueBasis{pcms::csys::CylindricalRThetaZ::Create(),
                                 ComponentScaling::Physical, values::Tensor},
                      ValueBasis{pcms::csys::Cartesian::Create(3),
                                 ComponentScaling::Physical, values::Tensor}),
    ContainsSubstring("not implemented"));
  REQUIRE_THROWS_WITH(
    TransformedSource(*src, cart_to_cyl, CartesianVector(), CartesianVector()),
    ContainsSubstring("source_basis must be"));
  REQUIRE_THROWS_WITH(TransformedSource(*src, cart_to_cyl, CylindricalVector(),
                                        CylindricalVector()),
                      ContainsSubstring("target_basis must be"));
}

TEST_CASE("TransformedSource: evaluator creation errors")
{
  auto lib = Omega_h::Library{};
  Omega_h::Mesh src_mesh(&lib);
  auto src = BuildCylindricalSource(lib, src_mesh, 3);
  auto cart_to_cyl = std::make_shared<pcms::CartesianToCylindrical>();
  const std::vector<Real> cart_pts = {0.5, 0.5, 0.5, 0.2, 0.9, 0.1};
  auto query = pcms::test::CreateDeviceCoordinateView(
    cart_pts, pcms::csys::Cartesian::Create(3), 3);

  SECTION("query points must be in the transformed system")
  {
    TransformedSource view(*src, cart_to_cyl, ValueBasis{}, ValueBasis{});
    auto cyl_query = pcms::test::CreateDeviceCoordinateView(
      cart_pts, pcms::csys::CylindricalRThetaZ::Create(), 3);
    REQUIRE_THROWS_WITH(
      view.CreatePointEvaluator<Real>(
        EvaluationRequest::FromCoordinates(cyl_query.coordinate_view)),
      ContainsSubstring("not the factory's coordinate system"));
  }

  SECTION("basis transformations require Real")
  {
    TransformedSource view(*src, cart_to_cyl, CylindricalVector(),
                           CartesianVector());
    REQUIRE_THROWS_WITH(
      view.CreatePointEvaluator<GO>(
        EvaluationRequest::FromCoordinates(query.coordinate_view)),
      ContainsSubstring("require T == Real"));
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

  TransformedSource as_cartesian(
    *src, std::make_shared<pcms::CartesianToCylindrical>(), CylindricalVector(),
    CartesianVector());
  TransformedSource back_to_cylindrical(
    as_cartesian, std::make_shared<pcms::CylindricalToCartesian>(),
    CartesianVector(), CylindricalVector());

  const std::vector<Real> cyl_pts = {0.7, 0.3, 0.4, 1.0, 0.9,
                                     0.9, 1.5, 1.2, 0.2};
  auto query = pcms::test::CreateDeviceCoordinateView(
    cyl_pts, pcms::csys::CylindricalRThetaZ::Create(), 3);
  const int n = 3;
  auto evaluator = back_to_cylindrical.CreatePointEvaluator<Real>(
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

  TransformedSource view(*src, std::make_shared<pcms::CartesianToCylindrical>(),
                         b_src.GetData().GetValueBasis(),
                         b_tgt.GetData().GetValueBasis());
  pcms::Interpolator<Real> op(view, *tgt);
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
