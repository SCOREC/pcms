#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <Kokkos_Core.hpp>
#include <pcms/utility/coordinate_system.hpp>
#include <pcms/field/coordinate_map.hpp>
#include <pcms/utility/coordinate_view.hpp>
#include <pcms/field/value_transform.hpp>
#include <pcms/field/value_view.hpp>
#include <cmath>
#include <random>
#include <vector>
#include "pcms/utility/coordinate_systems/cartesian.hpp"
#include "pcms/utility/coordinate_systems/cylindrical.hpp"
#include "field_test_utils.hpp"

using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;
using pcms::BoundCoordinateMap;
using pcms::ComponentScaling;
using pcms::CoordinateView;
using pcms::DeviceMemorySpace;
using pcms::FieldValueType;
using pcms::HostMemorySpace;
using pcms::Real;
using pcms::ValueBasis;
using pcms::ValueView;
using pcms::Variance;
using pcms::detail::SameCoordinateSystem;

namespace
{

ValueBasis CylindricalVector()
{
  return ValueBasis{pcms::csys::CreateCylindricalRThetaZ(),
                    ComponentScaling::Physical, pcms::values::Vector};
}

ValueBasis CartesianVector()
{
  return ValueBasis{pcms::csys::CreateCartesian(3), ComponentScaling::Physical,
                    pcms::values::Vector};
}

ValueBasis NaturalCylindricalVector()
{
  return ValueBasis{pcms::csys::CreateCylindricalRThetaZ(),
                    ComponentScaling::Natural, pcms::values::Vector};
}

// Re-expresses `values` in place through the bound map's value transform and
// returns the basis they are in afterwards.
ValueBasis TransformInPlace(const BoundCoordinateMap& bound,
                            const ValueBasis& basis,
                            Kokkos::View<Real**, DeviceMemorySpace> values)
{
  const auto& transform = bound.ValueTransform();
  const ValueBasis out = transform.OutputBasis(basis);
  transform.Apply(
    ValueView<Real, DeviceMemorySpace>(basis, pcms::MakeRank2View(values)),
    ValueView<Real, DeviceMemorySpace>(out, pcms::MakeRank2View(values)));
  return out;
}

} // namespace

TEST_CASE("coordinate systems: structural equality")
{
  // Independently constructed systems of the same family and dimension are
  // the same system: equality is structural, never pointer identity.
  REQUIRE(pcms::detail::SameCoordinateSystem(pcms::csys::CreateCartesian(3),
                                             pcms::csys::CreateCartesian(3)));
  REQUIRE(pcms::csys::CreateCartesian(3) != pcms::csys::CreateCartesian(3));
  REQUIRE(
    pcms::detail::SameCoordinateSystem(pcms::csys::CreateCylindricalRThetaZ(),
                                       pcms::csys::CreateCylindricalRThetaZ()));
  // Distinct dimensions and families are distinct systems.
  REQUIRE_FALSE(pcms::detail::SameCoordinateSystem(
    pcms::csys::CreateCartesian(2), pcms::csys::CreateCartesian(3)));
  REQUIRE_FALSE(pcms::detail::SameCoordinateSystem(
    pcms::csys::CreateCylindricalRZ(), pcms::csys::CreateCylindricalRThetaZ()));
  REQUIRE(pcms::csys::CreateCartesian(2)->Dimension() == 2);
  REQUIRE(pcms::csys::CreateCylindricalRZ()->Dimension() == 2);
  REQUIRE(pcms::csys::CreateCylindricalRThetaZ()->Dimension() == 3);

  // The analytic orthogonal coordinate systems are orthogonal
  REQUIRE(pcms::csys::CreateCartesian(3)->IsOrthogonal());
  REQUIRE(pcms::csys::CreateCylindricalRThetaZ()->IsOrthogonal());
}

TEST_CASE("coordinate systems: ResolveCoordinateSystem resolves the deferred "
          "Cartesian placeholder from the data dimension")
{
  const auto placeholder = pcms::detail::DeferredCartesian();
  REQUIRE(pcms::detail::SameCoordinateSystem(
    pcms::detail::ResolveCoordinateSystem(placeholder, 2),
    pcms::csys::CreateCartesian(2)));
  REQUIRE(pcms::detail::SameCoordinateSystem(
    pcms::detail::ResolveCoordinateSystem(placeholder, 3),
    pcms::csys::CreateCartesian(3)));
  // A concrete coordinate system validates against the data dimension instead.
  const auto concrete = pcms::csys::CreateCylindricalRThetaZ();
  REQUIRE(pcms::detail::SameCoordinateSystem(
    pcms::detail::ResolveCoordinateSystem(concrete, 3),
    pcms::csys::CreateCylindricalRThetaZ()));
  REQUIRE_THROWS_WITH(pcms::detail::ResolveCoordinateSystem(concrete, 2),
                      ContainsSubstring("coordinate columns"));
  REQUIRE_THROWS_WITH(pcms::detail::ResolveCoordinateSystem(nullptr, 2),
                      ContainsSubstring("null"));
  // A Cartesian that already states its dimension is validated like any
  // other system, not silently rewritten to the data's dimension.
  REQUIRE_THROWS_WITH(
    pcms::detail::ResolveCoordinateSystem(pcms::csys::CreateCartesian(3), 2),
    ContainsSubstring("coordinate columns"));
  // The unresolved placeholder's Kind self-explains in error messages.
  REQUIRE(placeholder->Kind() == "cartesian (dimension deferred)");
  // Deferred() is the only route to a dimension-0 Cartesian.
  REQUIRE_THROWS_WITH(pcms::csys::CreateCartesian(0),
                      ContainsSubstring("1..3"));
  REQUIRE_THROWS_WITH(pcms::csys::CreateCartesian(4),
                      ContainsSubstring("1..3"));
}

TEST_CASE("CoordinateView validates the dimension; family form derives it")
{
  auto pts2 = pcms::test::CreateDeviceCoordinateView(
    {0.0, 1.0, 2.0, 3.0}, pcms::detail::DeferredCartesian(), 2);
  // Family form: dimension from extent(1).
  CoordinateView<DeviceMemorySpace> v(pcms::detail::DeferredCartesian(),
                                      pcms::MakeConstRank2View(pts2.view));
  REQUIRE(pcms::detail::SameCoordinateSystem(v.GetCoordinateSystem(),
                                             pcms::csys::CreateCartesian(2)));
  // Mis-declared buffer fails at its entry point.
  REQUIRE_THROWS_WITH(
    CoordinateView<DeviceMemorySpace>(pcms::csys::CreateCylindricalRThetaZ(),
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
  const auto map = pcms::CreateCartesianToCylindrical();
  REQUIRE(pcms::detail::SameCoordinateSystem(map->GetSourceCoordinateSystem(),
                                             pcms::csys::CreateCartesian(3)));
  REQUIRE(pcms::detail::SameCoordinateSystem(
    map->GetTargetCoordinateSystem(), pcms::csys::CreateCylindricalRThetaZ()));

  const auto bound =
    map->Bind(pcms::test::MakeCoords(data, pcms::csys::CreateCartesian(3)));
  REQUIRE(pcms::detail::SameCoordinateSystem(
    bound->MappedPoints().GetCoordinateSystem(),
    pcms::csys::CreateCylindricalRThetaZ()));
  REQUIRE(bound->NumPoints() == 4);

  auto out =
    pcms::test::CopyCoordinatesToHost(bound->MappedPoints().GetValues());
  for (size_t i = 0; i < 4; ++i) {
    const Real x = xyz[3 * i], y = xyz[3 * i + 1], z = xyz[3 * i + 2];
    CAPTURE(i);
    REQUIRE_THAT(out(i, 0),
                 WithinAbs(std::sqrt(x * x + y * y), pcms::test::ExactTol));
    REQUIRE_THAT(out(i, 1), WithinAbs(std::atan2(y, x), pcms::test::ExactTol));
    REQUIRE_THAT(out(i, 2), WithinAbs(z, pcms::test::ExactTol));
  }

  // Values run opposite to points: a cylindrical field becomes Cartesian.
  REQUIRE(pcms::detail::SameValueBasis(
    bound->ValueTransform().OutputBasis(CylindricalVector()),
    CartesianVector()));
}

TEST_CASE("round trip via sequential Bind calls")
{
  const std::vector<Real> xyz = {0.3, -0.4, 1.0, 2.0, 1.0, -0.5, 0.0, 0.0, 3.0};
  auto data = pcms::test::CreateDeviceRank2View(xyz, 3);
  const auto stage1 = pcms::CreateCartesianToCylindrical()->Bind(
    pcms::test::MakeCoords(data, pcms::csys::CreateCartesian(3)));
  const auto stage2 =
    pcms::CreateCylindricalToCartesian()->Bind(stage1->MappedPoints());
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
  const auto bound = pcms::CreateCartesianToCylindrical()->Bind(
    pcms::test::MakeCoords(data, pcms::csys::CreateCartesian(3)));
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
  const auto map = pcms::CreateCartesianToCylindrical();
  auto a = pcms::test::CreateDeviceRank2View({1.0, 0.0, 0.0}, 3);
  auto b =
    pcms::test::CreateDeviceRank2View({0.0, 2.0, 1.0, -3.0, 0.0, 2.0}, 3);
  const auto bound_a =
    map->Bind(pcms::test::MakeCoords(a, pcms::csys::CreateCartesian(3)));
  const auto bound_b =
    map->Bind(pcms::test::MakeCoords(b, pcms::csys::CreateCartesian(3)));
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
    pcms::CreateCartesianToCylindrical()->Bind(
      pcms::test::MakeCoords(data, pcms::csys::CreateCylindricalRThetaZ())),
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
  for (size_t i = 0; i < 3; ++i) {
    xyz[3 * i] = rtz[3 * i] * std::cos(rtz[3 * i + 1]);
    xyz[3 * i + 1] = rtz[3 * i] * std::sin(rtz[3 * i + 1]);
    xyz[3 * i + 2] = rtz[3 * i + 2];
  }
  auto data = pcms::test::CreateDeviceRank2View(xyz, 3);
  const auto bound = pcms::CreateCartesianToCylindrical()->Bind(
    pcms::test::MakeCoords(data, pcms::csys::CreateCartesian(3)));
  REQUIRE(bound->NumPoints() == 3);

  auto values = pcms::test::CreateDeviceRank2View(
    {1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0}, 3);
  REQUIRE(pcms::detail::SameValueBasis(
    TransformInPlace(*bound, CylindricalVector(), values), CartesianVector()));
  auto out =
    pcms::test::CopyCoordinatesToHost(pcms::MakeConstRank2View(values));
  for (size_t i = 0; i < 3; ++i) {
    const Real theta = rtz[3 * i + 1];
    CAPTURE(i);
    REQUIRE_THAT(out(i, 0), WithinAbs(std::cos(theta), pcms::test::ExactTol));
    REQUIRE_THAT(out(i, 1), WithinAbs(std::sin(theta), pcms::test::ExactTol));
    REQUIRE_THAT(out(i, 2), WithinAbs(0.0, pcms::test::ExactTol));
  }
}

TEST_CASE("bound map: points on the axis rotate by the angle they map to")
{
  auto data = pcms::test::CreateDeviceRank2View({0.0, 0.0, 2.0}, 3);
  const auto bound = pcms::CreateCartesianToCylindrical()->Bind(
    pcms::test::MakeCoords(data, pcms::csys::CreateCartesian(3)));
  auto mapped =
    pcms::test::CopyCoordinatesToHost(bound->MappedPoints().GetValues());
  REQUIRE_THAT(mapped(0, 0), WithinAbs(0.0, pcms::test::ExactTol));
  const Real theta = mapped(0, 1);

  const std::vector<Real> c = {0.3, -1.2, 0.7};
  auto values = pcms::test::CreateDeviceRank2View(
    {c[0] * std::cos(theta) + c[1] * std::sin(theta),
     -c[0] * std::sin(theta) + c[1] * std::cos(theta), c[2]},
    3);
  REQUIRE(pcms::detail::SameValueBasis(
    TransformInPlace(*bound, CylindricalVector(), values), CartesianVector()));
  auto out =
    pcms::test::CopyCoordinatesToHost(pcms::MakeConstRank2View(values));
  for (size_t d = 0; d < 3; ++d) {
    REQUIRE_THAT(out(0, d), WithinAbs(c[d], pcms::test::ExactTol));
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
  const auto to_cart = pcms::CreateCartesianToCylindrical()->Bind(
    pcms::test::MakeCoords(cart_data, pcms::csys::CreateCartesian(3)));
  const auto to_cyl = pcms::CreateCylindricalToCartesian()->Bind(
    pcms::test::MakeCoords(cyl_data, pcms::csys::CreateCylindricalRThetaZ()));

  const std::vector<Real> v = {0.0, 1.0, 0.0}; // theta-direction unit vector
  auto values = pcms::test::CreateDeviceRank2View(v, 3);
  const auto as_cart = TransformInPlace(*to_cart, CylindricalVector(), values);
  auto a = pcms::test::CopyCoordinatesToHost(pcms::MakeConstRank2View(values));
  REQUIRE_THAT(a(0, 0), WithinAbs(-std::sin(theta), pcms::test::ExactTol));
  REQUIRE_THAT(a(0, 1), WithinAbs(std::cos(theta), pcms::test::ExactTol));

  REQUIRE(pcms::detail::SameValueBasis(
    TransformInPlace(*to_cyl, as_cart, values), CylindricalVector()));
  auto b = pcms::test::CopyCoordinatesToHost(pcms::MakeConstRank2View(values));
  for (int d = 0; d < 3; ++d) {
    REQUIRE_THAT(b(0, d), WithinAbs(v[d], pcms::test::ExactTol));
  }
}

TEST_CASE("bound map: the value transform validates the values it is given")
{
  const std::vector<Real> xyz = {1.0, 0.5, 0.0};
  auto data = pcms::test::CreateDeviceRank2View(xyz, 3);
  const auto bound = pcms::CreateCartesianToCylindrical()->Bind(
    pcms::test::MakeCoords(data, pcms::csys::CreateCartesian(3)));

  SECTION("natural components are re-expressed too")
  {
    // d/dtheta = R (-sin, cos, 0)
    auto values = pcms::test::CreateDeviceRank2View({0.0, 1.0, 0.0}, 3);
    REQUIRE(pcms::detail::SameValueBasis(
      TransformInPlace(*bound, NaturalCylindricalVector(), values),
      CartesianVector()));
    auto out =
      pcms::test::CopyCoordinatesToHost(pcms::MakeConstRank2View(values));
    REQUIRE_THAT(out(0, 0), WithinAbs(-0.5, pcms::test::ExactTol));
    REQUIRE_THAT(out(0, 1), WithinAbs(1.0, pcms::test::ExactTol));
    REQUIRE_THAT(out(0, 2), WithinAbs(0.0, pcms::test::ExactTol));
  }
  SECTION("a basis of another system")
  {
    REQUIRE_THROWS_WITH(bound->ValueTransform().OutputBasis(ValueBasis{
                          pcms::csys::CreateCylindricalRZ(),
                          ComponentScaling::Physical, pcms::values::Vector}),
                        ContainsSubstring("takes components in"));
  }
  SECTION("values already in the output basis pass through unchanged")
  {
    auto values = pcms::test::CreateDeviceRank2View({0.3, -0.4, 0.5}, 3);
    REQUIRE(pcms::detail::SameValueBasis(
      TransformInPlace(*bound, CartesianVector(), values), CartesianVector()));
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
    REQUIRE_THROWS_WITH(TransformInPlace(*bound, CylindricalVector(), values),
                        ContainsSubstring("one row per point"));
  }
}

TEST_CASE("bound map: physical covectors rotate like physical vectors")
{
  const Real theta = 0.7;
  auto data = pcms::test::CreateDeviceRank2View(
    {2.0 * std::cos(theta), 2.0 * std::sin(theta), 0.3}, 3);
  const auto bound = pcms::CreateCartesianToCylindrical()->Bind(
    pcms::test::MakeCoords(data, pcms::csys::CreateCartesian(3)));
  auto values = pcms::test::CreateDeviceRank2View({0.0, 1.0, 0.0}, 3);
  const ValueBasis covector{pcms::csys::CreateCylindricalRThetaZ(),
                            ComponentScaling::Physical, pcms::values::Covector};
  const ValueBasis out = TransformInPlace(*bound, covector, values);
  REQUIRE(out.scaling == ComponentScaling::Physical);
  REQUIRE(out.variance == pcms::values::Covector);
  auto host =
    pcms::test::CopyCoordinatesToHost(pcms::MakeConstRank2View(values));
  REQUIRE_THAT(host(0, 0), WithinAbs(-std::sin(theta), pcms::test::ExactTol));
  REQUIRE_THAT(host(0, 1), WithinAbs(std::cos(theta), pcms::test::ExactTol));
  REQUIRE_THAT(host(0, 2), WithinAbs(0.0, pcms::test::ExactTol));
}

TEST_CASE("the plane insertion passes covectors and refuses vectors")
{
  auto data = pcms::test::CreateDeviceRank2View({2.0, 0.5, 3.0, -1.0}, 2);
  const auto bound =
    pcms::CreateCoordinateEmbedding(pcms::csys::CreateCylindricalRZ(),
                                    pcms::csys::CreateCylindricalRThetaZ(), 1,
                                    0.3)
      ->Bind(pcms::test::MakeCoords(data, pcms::csys::CreateCylindricalRZ()));
  const auto& transform = bound->ValueTransform();
  const ValueBasis covector{pcms::csys::CreateCylindricalRThetaZ(),
                            ComponentScaling::Natural, pcms::values::Covector};
  const ValueBasis out = transform.OutputBasis(covector);
  REQUIRE(SameCoordinateSystem(out.system, pcms::csys::CreateCylindricalRZ()));
  // (dR, dtheta, dZ) components of a covector, two points
  auto in =
    pcms::test::CreateDeviceRank2View({1.0, 5.0, 2.0, 3.0, 7.0, 4.0}, 3);
  Kokkos::View<Real**, DeviceMemorySpace> result("result", 2, 2);
  transform.Apply(
    ValueView<Real, DeviceMemorySpace>(covector, pcms::MakeRank2View(in)),
    ValueView<Real, DeviceMemorySpace>(out, pcms::MakeRank2View(result)));
  auto host =
    pcms::test::CopyCoordinatesToHost(pcms::MakeConstRank2View(result));
  REQUIRE(host(0, 0) == 1.0);
  REQUIRE(host(0, 1) == 2.0);
  REQUIRE(host(1, 0) == 3.0);
  REQUIRE(host(1, 1) == 4.0);
  REQUIRE_THROWS_WITH(transform.OutputBasis(NaturalCylindricalVector()),
                      ContainsSubstring("not square"));
}

TEST_CASE("ValueView validates component counts and basis")
{
  auto data = pcms::test::CreateDeviceRank2View({1.0, 2.0}, 2);
  // A rank-1 value on a 3-dimensional system needs 3 components (per block:
  // values at points may hold several vectors per point).
  REQUIRE_THROWS_WITH(
    (ValueView<const Real, DeviceMemorySpace>(
      ValueBasis{pcms::csys::CreateCartesian(3), ComponentScaling::Physical,
                 pcms::values::Vector},
      pcms::MakeConstRank2View(data))),
    ContainsSubstring("a multiple of 3 components"));
  // Component count comes from the system's dimension: a 2-dimensional system
  // takes 2.
  REQUIRE_NOTHROW(ValueView<const Real, DeviceMemorySpace>(
    ValueBasis{pcms::csys::CreateCylindricalRZ(), ComponentScaling::Physical,
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
      ValueBasis{pcms::csys::CreateCartesian(3), ComponentScaling::Physical,
                 pcms::VarianceSignature{Variance::Contravariant,
                                         Variance::Contravariant,
                                         Variance::Contravariant}},
      pcms::MakeConstRank2View(data))),
    ContainsSubstring("rank-3"));
}

namespace
{

using DeviceValues = Kokkos::View<Real**, DeviceMemorySpace>;

pcms::Rank2View<Real, DeviceMemorySpace> MakeView(DeviceValues values)
{
  return pcms::MakeRank2View(values);
}

pcms::Rank2View<const Real, DeviceMemorySpace> MakeConstView(
  const DeviceValues& values)
{
  return pcms::MakeConstRank2View(values);
}

DeviceValues RandomValues(int rows, int columns, std::mt19937& engine)
{
  std::uniform_real_distribution<Real> uniform(-1.0, 1.0);
  DeviceValues values("values", rows, columns);
  const auto host = Kokkos::create_mirror_view(values);
  for (int p = 0; p < rows; ++p) {
    for (int c = 0; c < columns; ++c) {
      host(p, c) = uniform(engine);
    }
  }
  Kokkos::deep_copy(values, host);
  return values;
}

// A Jacobian [points][rows][columns]: the identity (where square) plus a
// perturbation, so that a square one is invertible.
Kokkos::View<Real***, DeviceMemorySpace> RandomJacobian(int n, int rows,
                                                        int columns,
                                                        std::mt19937& engine)
{
  std::uniform_real_distribution<Real> uniform(-0.4, 0.4);
  Kokkos::View<Real***, DeviceMemorySpace> jacobian("jacobian", n, rows,
                                                    columns);
  const auto host = Kokkos::create_mirror_view(jacobian);
  for (int p = 0; p < n; ++p) {
    for (int i = 0; i < rows; ++i) {
      for (int j = 0; j < columns; ++j) {
        host(p, i, j) = (i == j ? 1.0 : 0.0) + uniform(engine);
      }
    }
  }
  Kokkos::deep_copy(jacobian, host);
  return jacobian;
}

Real Pairing(const DeviceValues& a, const DeviceValues& b)
{
  const auto ha = Kokkos::create_mirror_view_and_copy(HostMemorySpace(), a);
  const auto hb = Kokkos::create_mirror_view_and_copy(HostMemorySpace(), b);
  Real sum = 0;
  for (size_t p = 0; p < ha.extent(0); ++p) {
    for (size_t c = 0; c < ha.extent(1); ++c) {
      sum += ha(p, c) * hb(p, c);
    }
  }
  return sum;
}

// <T v, w> = <v, T^T w> for random v and w; `blocks` values per point.
void CheckTranspose(const pcms::ValueTransform& transform,
                    const ValueBasis& in_basis, int blocks,
                    std::mt19937& engine)
{
  const ValueBasis out_basis = transform.OutputBasis(in_basis);
  const int n = transform.NumPoints();
  const int din = blocks * in_basis.system->Dimension();
  const int dout = blocks * out_basis.system->Dimension();
  const auto v = RandomValues(n, din, engine);
  const auto w = RandomValues(n, dout, engine);
  DeviceValues tv("tv", n, dout);
  DeviceValues tw("tw", n, din);
  transform.Apply(ValueView<Real, DeviceMemorySpace>(in_basis, MakeView(v)),
                  ValueView<Real, DeviceMemorySpace>(out_basis, MakeView(tv)));
  transform.ApplyTranspose(
    ValueView<Real, DeviceMemorySpace>(out_basis, MakeView(w)),
    ValueView<Real, DeviceMemorySpace>(in_basis, MakeView(tw)));
  const Real lhs = Pairing(tv, w);
  const Real rhs = Pairing(v, tw);
  CHECK_THAT(lhs, WithinAbs(rhs, 1e-12 * (1 + std::abs(rhs))));
}

ValueBasis Basis(std::shared_ptr<const pcms::CoordinateSystem> system,
                 ComponentScaling scaling, pcms::VarianceSignature variance)
{
  return ValueBasis{std::move(system), scaling, variance};
}

} // namespace

TEST_CASE("ValueTransform: ApplyTranspose is the transpose of Apply",
          "[value_transform]")
{
  std::mt19937 engine(7);
  const int n = 5;
  const auto cart2 = pcms::csys::CreateCartesian(2);
  const auto cart3 = pcms::csys::CreateCartesian(3);
  const auto cyl = pcms::csys::CreateCylindricalRThetaZ();
  // points of the coordinate systems: (R, theta, Z) with R > 0
  const auto points3 = RandomValues(n, 3, engine);
  const auto cyl_points =
    pcms::test::CreateDeviceRank2View({1.0, 0.1, 0.0, 2.0, 0.2, 1.0, 0.5, 3.0,
                                       -1.0, 3.0, -2.0, 0.5, 1.5, 1.0, 0.0},
                                      3);
  // the metrics of the coordinate systems at the points
  const auto g_cart3 = cart3->Bind(MakeConstView(points3));
  const auto g_cyl = cyl->Bind(MakeConstView(cyl_points));
  const auto square = RandomJacobian(n, 3, 3, engine);
  const auto wide = RandomJacobian(n, 3, 2, engine);
  const auto j_square =
    pcms::View<3, const Real, DeviceMemorySpace>(square.data(), n, 3, 3);
  const auto j_wide =
    pcms::View<3, const Real, DeviceMemorySpace>(wide.data(), n, 3, 2);
  const auto natural = ComponentScaling::Natural;
  const auto physical = ComponentScaling::Physical;
  using pcms::values::Covector;
  using pcms::values::Vector;

  SECTION("Pushforward, square")
  {
    const auto t =
      pcms::detail::Pushforward(j_square, cart3, cart3, g_cart3->Metric());
    CheckTranspose(t, Basis(cart3, physical, Vector), 1, engine);
    CheckTranspose(t, Basis(cart3, physical, Covector), 1, engine);
    CheckTranspose(t, Basis(cart3, physical, Covector), 2, engine);
  }
  SECTION("Pushforward, rectangular: vectors only")
  {
    const auto t =
      pcms::detail::Pushforward(j_wide, cart2, cart3, g_cart3->Metric());
    CheckTranspose(t, Basis(cart2, physical, Vector), 1, engine);
    REQUIRE_THROWS_WITH(t.OutputBasis(Basis(cart2, physical, Covector)),
                        ContainsSubstring("not square"));
  }
  SECTION("Pullback, square, into a coordinate system with scale factors")
  {
    const auto t =
      pcms::detail::Pullback(j_square, cart3, g_cart3->Metric(), cyl);
    CheckTranspose(t, Basis(cart3, physical, Vector), 1, engine);
    CheckTranspose(t, Basis(cart3, physical, Covector), 1, engine);
    REQUIRE(t.OutputBasis(Basis(cart3, physical, Vector)).scaling == physical);
  }
  SECTION("Pullback, from a coordinate system with scale factors")
  {
    const auto t =
      pcms::detail::Pullback(j_square, cyl, g_cyl->Metric(), cart3);
    CheckTranspose(t, Basis(cyl, physical, Vector), 1, engine);
    CheckTranspose(t, Basis(cyl, natural, Vector), 1, engine);
    CheckTranspose(t, Basis(cyl, physical, Covector), 1, engine);
    CheckTranspose(t, Basis(cyl, natural, Covector), 1, engine);
  }
  SECTION("Pullback, rectangular: covectors only")
  {
    const auto t =
      pcms::detail::Pullback(j_wide, cart3, g_cart3->Metric(), cart2);
    CheckTranspose(t, Basis(cart3, physical, Covector), 1, engine);
    REQUIRE_THROWS_WITH(t.OutputBasis(Basis(cart3, physical, Vector)),
                        ContainsSubstring("not square"));
  }
  SECTION("Scale")
  {
    const auto to_physical =
      pcms::detail::Scale(cyl, g_cyl->Metric(), physical);
    const auto to_natural = pcms::detail::Scale(cyl, g_cyl->Metric(), natural);
    CheckTranspose(to_physical, Basis(cyl, natural, Vector), 1, engine);
    CheckTranspose(to_physical, Basis(cyl, natural, Covector), 1, engine);
    CheckTranspose(to_natural, Basis(cyl, physical, Vector), 1, engine);
    CheckTranspose(to_natural, Basis(cyl, physical, Covector), 3, engine);
    // a natural covector at R = 2: the theta component is divided by R
    auto values =
      pcms::test::CreateDeviceRank2View({1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0,
                                         3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0},
                                        3);
    const auto in = Basis(cyl, natural, Covector);
    to_physical.Apply(ValueView<Real, DeviceMemorySpace>(in, MakeView(values)),
                      ValueView<Real, DeviceMemorySpace>(
                        to_physical.OutputBasis(in), MakeView(values)));
    auto host =
      pcms::test::CopyCoordinatesToHost(pcms::MakeConstRank2View(values));
    REQUIRE_THAT(host(1, 0), WithinAbs(1.0, pcms::test::ExactTol));
    REQUIRE_THAT(host(1, 1), WithinAbs(1.0, pcms::test::ExactTol));
    REQUIRE_THAT(host(1, 2), WithinAbs(3.0, pcms::test::ExactTol));
  }
}
