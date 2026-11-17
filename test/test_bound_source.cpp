#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <Kokkos_Core.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_library.hpp>
#include <pcms/geometry/bound_source.hpp>
#include <pcms/field/coordinate_systems/cylindrical.hpp>
#include <pcms/field/function_space/lagrange.h>
#include "field_test_utils.h"
#include <vector>

using Catch::Matchers::ContainsSubstring;
using pcms::HostMemorySpace;
using pcms::OutOfBoundsMode;
using pcms::OutOfBoundsPolicy;
using pcms::Real;

TEST_CASE("Evaluators created from one bound source share its points")
{
  auto lib = Omega_h::Library{};
  auto mesh = pcms::test::BuildSquareMesh(lib, 8);
  auto space = pcms::test::MakeP1Space(mesh);
  auto field = space->CreateFunction<Real>("f");
  const auto linear = OMEGA_H_LAMBDA(Real x, Real y)
  {
    return pcms::test::linear_f(x, y);
  };
  pcms::test::SetField(field, linear);

  const auto pts = pcms::test::StandardEvalCoords2D();
  const size_t n = pts.size() / 2;
  auto bound = [&] {
    // The bound source keeps what it needs, so the caller's points may go.
    return pcms::test::Bind(*space->GetGeometry(), pts);
  }();
  REQUIRE(bound->NumPoints() == static_cast<pcms::LO>(n));
  REQUIRE(bound->Policy().mode == OutOfBoundsMode::ERROR);

  auto first = space->CreatePointEvaluator<Real>(bound);
  auto second = space->CreatePointEvaluator<Real>(bound);
  pcms::test::CheckEvaluation(*first, field, pts, linear);
  pcms::test::CheckEvaluation(*second, field, pts, linear);
}

TEST_CASE("A bound source is accepted only by sources on its geometry")
{
  auto lib = Omega_h::Library{};
  auto mesh = pcms::test::BuildSquareMesh(lib, 8);
  auto space = pcms::test::MakeP1Space(mesh);
  auto other = pcms::test::MakeP1Space(mesh);
  auto bound =
    pcms::test::Bind(*space->GetGeometry(), pcms::test::StandardEvalCoords2D());

  REQUIRE_THROWS_WITH(other->CreatePointEvaluator<Real>(bound),
                      ContainsSubstring("another geometry"));
  REQUIRE_THROWS_WITH(space->CreatePointEvaluator<Real>(nullptr),
                      ContainsSubstring("must not be null"));
  REQUIRE_THROWS_WITH(space->CreatePointEvaluator<Real>(
                        bound, pcms::DifferentialOperator::Divergence),
                      ContainsSubstring("Divergence is not supported"));
}

TEST_CASE("Binding checks the query coordinate system")
{
  auto lib = Omega_h::Library{};
  auto mesh = pcms::test::BuildSquareMesh(lib, 8);
  auto space = pcms::test::MakeP1Space(mesh);
  auto coords = pcms::test::CreateDeviceCoordinateView(
    pcms::test::StandardEvalCoords2D(), pcms::csys::CylindricalRZ::Create());
  REQUIRE_THROWS_WITH(space->GetGeometry()->Bind(coords.coordinate_view),
                      ContainsSubstring("is not the geometry's coordinate"));
}

TEST_CASE("Sources without derivatives reject the gradient")
{
  auto lib = Omega_h::Library{};
  auto mesh = pcms::test::BuildSquareMesh(lib, 8);
  auto space = pcms::test::MakeP1Space(mesh);
  auto bound =
    pcms::test::Bind(*space->GetGeometry(), pcms::test::StandardEvalCoords2D());
  REQUIRE_FALSE(space->Supports(pcms::DifferentialOperator::Gradient));
  REQUIRE_THROWS_WITH(space->CreatePointEvaluator<Real>(
                        bound, pcms::DifferentialOperator::Gradient),
                      ContainsSubstring("Gradient is not supported"));
}

TEST_CASE("Binding splits found and missing points and applies the policy")
{
  auto lib = Omega_h::Library{};
  auto mesh = pcms::test::BuildSquareMesh(lib, 8);
  auto space = pcms::test::MakeP1Space(mesh);
  auto field = space->CreateFunction<Real>("f");
  pcms::test::SetField(
    field,
    OMEGA_H_LAMBDA(Real x, Real y) { return pcms::test::linear_f(x, y); });
  // the second point is outside the unit square
  const std::vector<Real> pts = {0.5, 0.5, 2.0, 2.0, 0.25, 0.75};

  REQUIRE_THROWS_WITH(pcms::test::Bind(*space->GetGeometry(), pts),
                      ContainsSubstring("outside the source's domain"));
  REQUIRE_THROWS_WITH(
    pcms::test::Bind(*space->GetGeometry(), pts,
                     OutOfBoundsPolicy{OutOfBoundsMode::NEAREST_BOUNDARY}),
    ContainsSubstring("NEAREST_BOUNDARY"));

  auto bound = pcms::test::Bind(*space->GetGeometry(), pts,
                                OutOfBoundsPolicy{OutOfBoundsMode::FILL, -1.0});
  const auto found = Kokkos::create_mirror_view_and_copy(HostMemorySpace(),
                                                         bound->FoundPoints());
  const auto missing = Kokkos::create_mirror_view_and_copy(
    HostMemorySpace(), bound->MissingPoints());
  REQUIRE(found.extent(0) == 2);
  REQUIRE(found(0) == 0);
  REQUIRE(found(1) == 2);
  REQUIRE(missing.extent(0) == 1);
  REQUIRE(missing(0) == 1);

  pcms::test::CheckEvaluationWithFill(
    *space->CreatePointEvaluator<Real>(bound), field, pts, {true, false, true},
    OMEGA_H_LAMBDA(Real x, Real y) { return pcms::test::linear_f(x, y); },
    -1.0);
}
