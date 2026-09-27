#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <Kokkos_Core.hpp>
#include <pcms/field/coordinate_system.hpp>
#include <pcms/field/coordinate_map.hpp>
#include <pcms/field/coordinate_view.hpp>
#include <pcms/field/value_view.hpp>
#include <cmath>
#include <vector>
#include "pcms/field/coordinate_systems/cartesian.hpp"
#include "pcms/field/coordinate_systems/cylindrical.hpp"
#include "field_test_utils.h"

using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;
using pcms::BoundCoordinateMap;
using pcms::ComponentScaling;
using pcms::CoordinateView;
using pcms::DeviceMemorySpace;
using pcms::FieldValueType;
using pcms::HostMemorySpace;
using pcms::Real;
using pcms::SameCoordinateSystem;
using pcms::ValueBasis;
using pcms::ValueView;
using pcms::Variance;

namespace
{

ValueBasis CylindricalVector()
{
  return ValueBasis{pcms::csys::CylindricalRThetaZ::Create(),
                    ComponentScaling::Physical, pcms::values::Vector};
}

ValueBasis CartesianVector()
{
  return ValueBasis{pcms::csys::Cartesian::Create(3),
                    ComponentScaling::Physical, pcms::values::Vector};
}

ValueBasis NaturalCylindricalVector()
{
  return ValueBasis{pcms::csys::CylindricalRThetaZ::Create(),
                    ComponentScaling::Natural, pcms::values::Vector};
}

} // namespace

TEST_CASE("coordinate systems: structural equality")
{
  // Independently constructed systems of the same family and dimension are
  // the same system: equality is structural, never pointer identity.
  REQUIRE(SameCoordinateSystem(pcms::csys::Cartesian::Create(3),
                               pcms::csys::Cartesian::Create(3)));
  REQUIRE(pcms::csys::Cartesian::Create(3) != pcms::csys::Cartesian::Create(3));
  REQUIRE(SameCoordinateSystem(pcms::csys::CylindricalRThetaZ::Create(),
                               pcms::csys::CylindricalRThetaZ::Create()));
  // Distinct dimensions and families are distinct systems.
  REQUIRE_FALSE(SameCoordinateSystem(pcms::csys::Cartesian::Create(2),
                                     pcms::csys::Cartesian::Create(3)));
  REQUIRE_FALSE(SameCoordinateSystem(pcms::csys::CylindricalRZ::Create(),
                                     pcms::csys::CylindricalRThetaZ::Create()));
  REQUIRE(pcms::csys::Cartesian::Create(2)->Dimension() == 2);
  REQUIRE(pcms::csys::CylindricalRZ::Create()->Dimension() == 2);
  REQUIRE(pcms::csys::CylindricalRThetaZ::Create()->Dimension() == 3);

  // The analytic orthogonal coordinate systems are orthogonal
  REQUIRE(pcms::csys::Cartesian::Create(3)->HasOrthogonalBasis());
  REQUIRE(pcms::csys::CylindricalRThetaZ::Create()->HasOrthogonalBasis());
}

TEST_CASE("coordinate systems: ResolveCoordinateSystem resolves the deferred "
          "Cartesian placeholder from the data dimension")
{
  const auto placeholder = pcms::csys::Cartesian::Deferred();
  REQUIRE(SameCoordinateSystem(pcms::ResolveCoordinateSystem(placeholder, 2),
                               pcms::csys::Cartesian::Create(2)));
  REQUIRE(SameCoordinateSystem(pcms::ResolveCoordinateSystem(placeholder, 3),
                               pcms::csys::Cartesian::Create(3)));
  // A concrete coordinate system validates against the data dimension instead.
  const auto concrete = pcms::csys::CylindricalRThetaZ::Create();
  REQUIRE(SameCoordinateSystem(pcms::ResolveCoordinateSystem(concrete, 3),
                               pcms::csys::CylindricalRThetaZ::Create()));
  REQUIRE_THROWS_WITH(pcms::ResolveCoordinateSystem(concrete, 2),
                      ContainsSubstring("coordinate columns"));
  REQUIRE_THROWS_WITH(pcms::ResolveCoordinateSystem(nullptr, 2),
                      ContainsSubstring("null"));
  // A Cartesian that already states its dimension is validated like any
  // other system, not silently rewritten to the data's dimension.
  REQUIRE_THROWS_WITH(
    pcms::ResolveCoordinateSystem(pcms::csys::Cartesian::Create(3), 2),
    ContainsSubstring("coordinate columns"));
  // The unresolved placeholder's Kind self-explains in error messages.
  REQUIRE(placeholder->Kind() == "cartesian (dimension deferred)");
  // Deferred() is the only route to a dimension-0 Cartesian.
  REQUIRE_THROWS_WITH(pcms::csys::Cartesian::Create(0),
                      ContainsSubstring("1..3"));
  REQUIRE_THROWS_WITH(pcms::csys::Cartesian::Create(4),
                      ContainsSubstring("1..3"));
}

TEST_CASE("CoordinateView validates the dimension; family form derives it")
{
  auto pts2 = pcms::test::CreateDeviceCoordinateView(
    {0.0, 1.0, 2.0, 3.0}, pcms::csys::Cartesian::Deferred(), 2);
  // Family form: dimension from extent(1).
  CoordinateView<DeviceMemorySpace> v(pcms::csys::Cartesian::Deferred(),
                                      pcms::MakeConstRank2View(pts2.view));
  REQUIRE(SameCoordinateSystem(v.GetCoordinateSystem(),
                               pcms::csys::Cartesian::Create(2)));
  // Mis-declared buffer fails at its entry point.
  REQUIRE_THROWS_WITH(
    CoordinateView<DeviceMemorySpace>(pcms::csys::CylindricalRThetaZ::Create(),
                                      pcms::MakeConstRank2View(pts2.view)),
    ContainsSubstring("coordinate columns"));
}

TEST_CASE("CartesianToCylindrical: Bind produces known values")
{
  const std::vector<Real> xyz = {1.0,  0.0, 0.5,  //
                                 0.0,  2.0, -1.0, //
                                 -3.0, 0.0, 2.0,  //
                                 1.0,  1.0, 0.0};
  auto data = pcms::test::CreateDeviceRank2View(xyz, 3);
  const pcms::CartesianToCylindrical map;
  REQUIRE(SameCoordinateSystem(map.GetSourceCoordinateSystem(),
                               pcms::csys::Cartesian::Create(3)));
  REQUIRE(SameCoordinateSystem(map.GetTargetCoordinateSystem(),
                               pcms::csys::CylindricalRThetaZ::Create()));

  const auto bound =
    map.Bind(pcms::test::MakeCoords(data, pcms::csys::Cartesian::Create(3)));
  REQUIRE(SameCoordinateSystem(bound->MappedPoints().GetCoordinateSystem(),
                               pcms::csys::CylindricalRThetaZ::Create()));
  REQUIRE(bound->NumPoints() == 4);

  auto out =
    pcms::test::CopyCoordinatesToHost(bound->MappedPoints().GetValues());
  for (int i = 0; i < 4; ++i) {
    const Real x = xyz[3 * i], y = xyz[3 * i + 1], z = xyz[3 * i + 2];
    CAPTURE(i);
    REQUIRE_THAT(out(i, 0),
                 WithinAbs(std::sqrt(x * x + y * y), pcms::test::ExactTol));
    REQUIRE_THAT(out(i, 1), WithinAbs(std::atan2(y, x), pcms::test::ExactTol));
    REQUIRE_THAT(out(i, 2), WithinAbs(z, pcms::test::ExactTol));
  }

  // Values run opposite to points: a cylindrical field becomes Cartesian.
  REQUIRE(pcms::SameValueBasis(bound->OutputBasis(CylindricalVector()),
                               CartesianVector()));
}

TEST_CASE("round trip via sequential Bind calls")
{
  const std::vector<Real> xyz = {0.3, -0.4, 1.0, 2.0, 1.0, -0.5, 0.0, 0.0, 3.0};
  auto data = pcms::test::CreateDeviceRank2View(xyz, 3);
  const auto stage1 = pcms::CartesianToCylindrical{}.Bind(
    pcms::test::MakeCoords(data, pcms::csys::Cartesian::Create(3)));
  const auto stage2 =
    pcms::CylindricalToCartesian{}.Bind(stage1->MappedPoints());
  auto out =
    pcms::test::CopyCoordinatesToHost(stage2->MappedPoints().GetValues());
  for (size_t i = 0; i < xyz.size() / 3; ++i) {
    for (int d = 0; d < 3; ++d) {
      REQUIRE_THAT(out(i, d), WithinAbs(xyz[3 * i + d], pcms::test::ExactTol));
    }
  }
}

TEST_CASE("Bind returns owned output and completes before returning")
{
  const std::vector<Real> xyz = {1.0, 0.0, 3.0};
  auto data = pcms::test::CreateDeviceRank2View(xyz, 3);
  const auto bound = pcms::CartesianToCylindrical{}.Bind(
    pcms::test::MakeCoords(data, pcms::csys::Cartesian::Create(3)));
  // Clobber the caller's buffer after Bind; the owned output stands.
  Kokkos::deep_copy(data, -99.0);
  auto out =
    pcms::test::CopyCoordinatesToHost(bound->MappedPoints().GetValues());
  REQUIRE_THAT(out(0, 0), WithinAbs(1.0, pcms::test::ExactTol));
  REQUIRE_THAT(out(0, 1), WithinAbs(0.0, pcms::test::ExactTol));
  REQUIRE_THAT(out(0, 2), WithinAbs(3.0, pcms::test::ExactTol));
}

TEST_CASE("one map instance serves multiple point sets")
{
  const pcms::CartesianToCylindrical map;
  auto a = pcms::test::CreateDeviceRank2View({1.0, 0.0, 0.0}, 3);
  auto b =
    pcms::test::CreateDeviceRank2View({0.0, 2.0, 1.0, -3.0, 0.0, 2.0}, 3);
  const auto bound_a =
    map.Bind(pcms::test::MakeCoords(a, pcms::csys::Cartesian::Create(3)));
  const auto bound_b =
    map.Bind(pcms::test::MakeCoords(b, pcms::csys::Cartesian::Create(3)));
  REQUIRE(bound_a->NumPoints() == 1);
  REQUIRE(bound_b->NumPoints() == 2);
  auto out_a =
    pcms::test::CopyCoordinatesToHost(bound_a->MappedPoints().GetValues());
  auto out_b =
    pcms::test::CopyCoordinatesToHost(bound_b->MappedPoints().GetValues());
  REQUIRE_THAT(out_a(0, 1), WithinAbs(0.0, pcms::test::ExactTol));
  REQUIRE_THAT(out_b(0, 1), WithinAbs(M_PI / 2.0, pcms::test::ExactTol));
  REQUIRE_THAT(out_b(1, 1), WithinAbs(M_PI, pcms::test::ExactTol));
}

TEST_CASE("Bind rejects a source system that is not the map's own")
{
  const std::vector<Real> pts = {1.0, 0.0, 0.0};
  auto data = pcms::test::CreateDeviceRank2View(pts, 3);
  REQUIRE_THROWS_WITH(
    pcms::CartesianToCylindrical{}.Bind(
      pcms::test::MakeCoords(data, pcms::csys::CylindricalRThetaZ::Create())),
    ContainsSubstring("source coordinate system"));
}

TEST_CASE("bound map: manufactured pushforward rotates vector components")
{
  // Radial unit vectors in cylindrical components must become
  // (cos theta, sin theta, 0) in Cartesian components.
  const std::vector<Real> rtz = {1.0, 0.0,        0.0, //
                                 2.0, M_PI / 2.0, 1.0, //
                                 0.5, M_PI / 4.0, -1.0};
  std::vector<Real> xyz(rtz.size());
  for (int i = 0; i < 3; ++i) {
    xyz[3 * i] = rtz[3 * i] * std::cos(rtz[3 * i + 1]);
    xyz[3 * i + 1] = rtz[3 * i] * std::sin(rtz[3 * i + 1]);
    xyz[3 * i + 2] = rtz[3 * i + 2];
  }
  auto data = pcms::test::CreateDeviceRank2View(xyz, 3);
  const auto bound = pcms::CartesianToCylindrical{}.Bind(
    pcms::test::MakeCoords(data, pcms::csys::Cartesian::Create(3)));
  REQUIRE(bound->NumPoints() == 3);

  auto values = pcms::test::CreateDeviceRank2View(
    {1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0}, 3);
  const auto result = bound->TransformValues(ValueView<Real, DeviceMemorySpace>(
    CylindricalVector(), pcms::MakeRank2View(values)));
  REQUIRE(pcms::SameValueBasis(result.GetBasis(), CartesianVector()));
  REQUIRE(result.GetValues().data_handle() == values.data());
  auto out =
    pcms::test::CopyCoordinatesToHost(pcms::MakeConstRank2View(values));
  for (int i = 0; i < 3; ++i) {
    const Real theta = rtz[3 * i + 1];
    CAPTURE(i);
    REQUIRE_THAT(out(i, 0), WithinAbs(std::cos(theta), pcms::test::ExactTol));
    REQUIRE_THAT(out(i, 1), WithinAbs(std::sin(theta), pcms::test::ExactTol));
    REQUIRE_THAT(out(i, 2), WithinAbs(0.0, pcms::test::ExactTol));
  }
}

TEST_CASE("bound map: the two cylindrical maps invert each other")
{
  const Real theta = 0.7;
  const std::vector<Real> cart = {2.0 * std::cos(theta), 2.0 * std::sin(theta),
                                  0.3};
  const std::vector<Real> cyl = {2.0, theta, 0.3};
  auto cart_data = pcms::test::CreateDeviceRank2View(cart, 3);
  auto cyl_data = pcms::test::CreateDeviceRank2View(cyl, 3);
  const auto to_cart = pcms::CartesianToCylindrical{}.Bind(
    pcms::test::MakeCoords(cart_data, pcms::csys::Cartesian::Create(3)));
  const auto to_cyl = pcms::CylindricalToCartesian{}.Bind(
    pcms::test::MakeCoords(cyl_data, pcms::csys::CylindricalRThetaZ::Create()));

  const std::vector<Real> v = {0.0, 1.0, 0.0}; // theta-direction unit vector
  auto values = pcms::test::CreateDeviceRank2View(v, 3);
  const auto as_cart =
    to_cart->TransformValues(ValueView<Real, DeviceMemorySpace>(
      CylindricalVector(), pcms::MakeRank2View(values)));
  auto a = pcms::test::CopyCoordinatesToHost(pcms::MakeConstRank2View(values));
  REQUIRE_THAT(a(0, 0), WithinAbs(-std::sin(theta), pcms::test::ExactTol));
  REQUIRE_THAT(a(0, 1), WithinAbs(std::cos(theta), pcms::test::ExactTol));

  const auto round_trip = to_cyl->TransformValues(as_cart);
  REQUIRE(pcms::SameValueBasis(round_trip.GetBasis(), CylindricalVector()));
  auto b = pcms::test::CopyCoordinatesToHost(pcms::MakeConstRank2View(values));
  for (int d = 0; d < 3; ++d) {
    REQUIRE_THAT(b(0, d), WithinAbs(v[d], pcms::test::ExactTol));
  }
}

TEST_CASE("bound map: TransformValues validates the values it is given")
{
  const std::vector<Real> xyz = {1.0, 0.5, 0.0};
  auto data = pcms::test::CreateDeviceRank2View(xyz, 3);
  const auto bound = pcms::CartesianToCylindrical{}.Bind(
    pcms::test::MakeCoords(data, pcms::csys::Cartesian::Create(3)));

  SECTION("a basis the map cannot re-express")
  {
    auto values = pcms::test::CreateDeviceRank2View({1.0, 0.0, 0.0}, 3);
    REQUIRE_THROWS_WITH(
      bound->TransformValues(ValueView<Real, DeviceMemorySpace>(
        NaturalCylindricalVector(), pcms::MakeRank2View(values))),
      ContainsSubstring("cannot re-express"));
    REQUIRE_THROWS_WITH(bound->OutputBasis(NaturalCylindricalVector()),
                        ContainsSubstring("cannot re-express"));
  }
  SECTION("values already in the output basis pass through unchanged")
  {
    auto values = pcms::test::CreateDeviceRank2View({0.3, -0.4, 0.5}, 3);
    REQUIRE(pcms::SameValueBasis(bound->OutputBasis(CartesianVector()),
                                 CartesianVector()));
    const auto result =
      bound->TransformValues(ValueView<Real, DeviceMemorySpace>(
        CartesianVector(), pcms::MakeRank2View(values)));
    REQUIRE(pcms::SameValueBasis(result.GetBasis(), CartesianVector()));
    auto out =
      pcms::test::CopyCoordinatesToHost(pcms::MakeConstRank2View(values));
    REQUIRE(out(0, 0) == 0.3);
    REQUIRE(out(0, 1) == -0.4);
    REQUIRE(out(0, 2) == 0.5);
  }
  SECTION("a value count that is not the bound point count")
  {
    auto values =
      pcms::test::CreateDeviceRank2View({1.0, 0.0, 0.0, 0.0, 1.0, 0.0}, 3);
    REQUIRE_THROWS_WITH(
      bound->TransformValues(ValueView<Real, DeviceMemorySpace>(
        CylindricalVector(), pcms::MakeRank2View(values))),
      ContainsSubstring("value count does not match"));
  }
}

TEST_CASE("a points-only bound map cannot re-express component values")
{
  auto data = pcms::test::CreateDeviceRank2View({1.0, 2.0, 3.0}, 3);
  Kokkos::View<Real**, DeviceMemorySpace> mapped("mapped", 1, 3);
  const BoundCoordinateMap bound(pcms::csys::Cartesian::Create(3), mapped);
  REQUIRE(bound.NumPoints() == 1);
  REQUIRE_THROWS_WITH(bound.OutputBasis(CartesianVector()),
                      ContainsSubstring("maps points only"));
}

TEST_CASE("ValueView validates component counts and basis")
{
  auto data = pcms::test::CreateDeviceRank2View({1.0, 2.0}, 2);
  // A rank-1 value on a 3-dimensional system needs 3 components.
  REQUIRE_THROWS_WITH(
    (ValueView<const Real, DeviceMemorySpace>(
      ValueBasis{pcms::csys::Cartesian::Create(3), ComponentScaling::Physical,
                 pcms::values::Vector},
      pcms::MakeConstRank2View(data))),
    ContainsSubstring("require 3 components"));
  // Component count comes from the system's dimension: a 2-dimensional system
  // takes 2.
  REQUIRE_NOTHROW(ValueView<const Real, DeviceMemorySpace>(
    ValueBasis{pcms::csys::CylindricalRZ::Create(), ComponentScaling::Physical,
               pcms::values::Vector},
    pcms::MakeConstRank2View(data)));
  // Non-scalar values require a basis coordinate system.
  auto data3 = pcms::test::CreateDeviceRank2View({1.0, 2.0, 3.0}, 3);
  REQUIRE_THROWS_WITH(
    (ValueView<const Real, DeviceMemorySpace>(
      ValueBasis{nullptr, ComponentScaling::Physical, pcms::values::Vector},
      pcms::MakeConstRank2View(data3))),
    ContainsSubstring("basis coordinate system"));
  // Physical components are undefined on a non-orthogonal basis, and rank 0
  // (an empty signature) ignores the basis at any component count.
  REQUIRE_NOTHROW(ValueView<const Real, DeviceMemorySpace>(
    ValueBasis{}, pcms::MakeConstRank2View(data)));
}

TEST_CASE("value declarations reject rank > 2")
{
  std::vector<Real> flat(27, 0.0);
  auto data = pcms::test::CreateDeviceRank2View(flat, 27);
  REQUIRE_THROWS_WITH(
    (ValueView<const Real, DeviceMemorySpace>(
      ValueBasis{
        pcms::csys::Cartesian::Create(3), ComponentScaling::Physical,
        pcms::values::Of({Variance::Contravariant, Variance::Contravariant,
                          Variance::Contravariant})},
      pcms::MakeConstRank2View(data))),
    ContainsSubstring("rank-3"));
}
