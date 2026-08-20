#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <Kokkos_Core.hpp>
#include <pcms/field/coordinate_system.hpp>
#include <pcms/field/coordinate_view.hpp>
#include <vector>
#include "pcms/field/coordinate_systems/cartesian.hpp"
#include "pcms/field/coordinate_systems/cylindrical.hpp"
#include "field_test_utils.h"

using Catch::Matchers::ContainsSubstring;
using pcms::CoordinateView;
using pcms::DeviceMemorySpace;
using pcms::HostMemorySpace;
using pcms::Real;
using pcms::SameCoordinateSystem;

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
  REQUIRE(
    pcms::csys::CylindricalRThetaZ::Create()->HasOrthogonalBasis());
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
