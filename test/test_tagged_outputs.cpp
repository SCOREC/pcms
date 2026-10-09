#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <Kokkos_Core.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_library.hpp>
#include <pcms/field/function_space/lagrange.h>
#include <pcms/transfer/interpolator.h>
#include "field_test_utils.h"
#include <memory>
#include <vector>
#include "pcms/field/coordinate_systems/cartesian.hpp"
#include "pcms/field/coordinate_systems/cylindrical.hpp"

using Catch::Matchers::ContainsSubstring;
using pcms::DeviceMemorySpace;
using pcms::Field;
using pcms::LagrangeFunctionSpace;
using pcms::Real;
using pcms::ValueView;
namespace values = pcms::values;
using Catch::Matchers::WithinAbs;
using pcms::ComponentScaling;
using pcms::HostMemorySpace;
using pcms::ValueBasis;
using pcms::Variance;
using pcms::VarianceSignature;

namespace
{

// A 3D simplex mesh whose coordinates are interpreted as (r, theta, z).
Omega_h::Mesh BuildCylindricalMesh(Omega_h::Library& lib, int divisions)
{
  return Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 2.0, 1.6, 1.1,
                            divisions, divisions, divisions, false);
}

std::shared_ptr<LagrangeFunctionSpace> BuildCylindricalSpace(
  Omega_h::Mesh& mesh, int num_components)
{
  return LagrangeFunctionSpace::FromMesh(
    mesh, 1, num_components, pcms::csys::CylindricalRThetaZ::Create(), "global",
    LagrangeFunctionSpace::Backend::OmegaH);
}

int NumDOFHolders(const LagrangeFunctionSpace& space)
{
  return static_cast<int>(
    space.GetLayout()->GetDOFHolderCoordinates().GetValues().extent(0));
}

} // namespace

TEST_CASE("PointEvaluator::Evaluate gates on the source field's stored basis")
{
  auto lib = Omega_h::Library{};
  auto mesh = BuildCylindricalMesh(lib, 6);
  auto space = BuildCylindricalSpace(mesh, 3);
  // Borrowed-basis storage: Cartesian components on a cylindrical-system space.
  auto cartesian = space->CreateFunction<Real>(
    "b", values::Vector, pcms::csys::Cartesian::Create(3));
  auto cylindrical = space->CreateFunction<Real>("b_native", values::Vector,
                                                 ComponentScaling::Physical);
  pcms::test::SetFieldComponents(cartesian, [](Real, Real, Real, Real* out) {
    out[0] = 1.0;
    out[1] = 0.0;
    out[2] = 0.0;
  });

  const std::vector<Real> pts = {0.7, 0.3, 0.4, 1.0, 0.9, 0.9};
  auto query = pcms::test::CreateDeviceCoordinateView(
    pts, pcms::csys::CylindricalRThetaZ::Create(), 3);
  auto evaluator = space->CreatePointEvaluator<Real>(
    pcms::EvaluationRequest::FromCoordinates(query.coordinate_view));
  const int n = static_cast<int>(pts.size()) / 3;
  Kokkos::View<Real**, DeviceMemorySpace> out("out", n, 3);

  SECTION("a matching tag is accepted")
  {
    REQUIRE_NOTHROW(evaluator->Evaluate(
      cartesian,
      ValueView<Real, DeviceMemorySpace>(cartesian.GetData().GetValueBasis(),
                                         pcms::MakeRank2View(out))));
  }

  SECTION("a mismatched basis claim is rejected")
  {
    REQUIRE_THROWS_WITH(
      evaluator->Evaluate(cartesian, ValueView<Real, DeviceMemorySpace>(
                                       cylindrical.GetData().GetValueBasis(),
                                       pcms::MakeRank2View(out))),
      ContainsSubstring("does not match the basis this call writes"));
  }
}

TEST_CASE("tagged writes gate on the field's declaration")
{
  auto lib = Omega_h::Library{};
  auto mesh = BuildCylindricalMesh(lib, 6);
  auto space = BuildCylindricalSpace(mesh, 3);
  auto b = space->CreateFunction<Real>("b", values::Vector,
                                       ComponentScaling::Physical);
  const int n = static_cast<int>(
    space->GetLayout()->GetDOFHolderCoordinates().GetValues().extent(0));
  std::vector<Real> data(static_cast<size_t>(n) * 3, 1.0);
  const pcms::Rank2View<const Real, HostMemorySpace> raw(data.data(), n, 3);

  SECTION("a matching tag is accepted")
  {
    b.SetDOFHolderDataHost(ValueView<const Real, HostMemorySpace>(
      ValueBasis{pcms::csys::CylindricalRThetaZ::Create(),
                 ComponentScaling::Physical,
                 VarianceSignature{Variance::Contravariant}},
      raw));
    REQUIRE_THAT(b.GetDOFHolderDataHost()(0, 0),
                 WithinAbs(1.0, 1e-14)); // tagged getter forwards indexing
  }

  SECTION("a mismatched basis claim is rejected")
  {
    REQUIRE_THROWS_WITH(
      b.SetDOFHolderDataHost(ValueView<const Real, HostMemorySpace>(
        ValueBasis{pcms::csys::Cartesian::Create(3), ComponentScaling::Physical,
                   VarianceSignature{Variance::Contravariant}},
        raw)),
      ContainsSubstring("does not match this field's declaration"));
  }

  SECTION("a mismatched value type is rejected")
  {
    REQUIRE_THROWS_WITH(
      b.SetDOFHolderDataHost(
        ValueView<const Real, HostMemorySpace>(ValueBasis{}, raw)),
      ContainsSubstring("does not match this field's declaration"));
  }

  SECTION("the unchecked path inherits the declaration unverified")
  {
    REQUIRE_NOTHROW(b.SetDOFHolderDataUncheckedHost(raw));
  }
}
