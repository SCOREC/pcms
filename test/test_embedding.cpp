#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <pcms/utility/coordinate_systems/cylindrical.hpp>
#include "embedding_test_utils.hpp"

using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;
using pcms::DeviceMemorySpace;
using pcms::Real;
using pcms::test::Annulus;

namespace
{

// the analytic Jacobian J(j, i) = d x_j / d xi_i of the annulus
std::array<std::array<Real, 2>, 2> AnnulusJacobian(Real xi1, Real xi2)
{
  const Real rho = Annulus::Rho(xi1);
  const Real theta = Annulus::Theta(xi2);
  return {{{0.5 * std::cos(theta), -2.0 * rho * std::sin(theta)},
           {0.5 * std::sin(theta), 2.0 * rho * std::cos(theta)}}};
}

const std::vector<Real> Params = {0.31, 0.42, 0.77, 0.13, 0.5, 0.95, 0.05, 0.6};
// the position is a spline of 32 x 32 cells: a smooth map is reproduced to
// roughly its fourth derivative times h^4
constexpr Real SplineTol = 2e-6;

struct Fixture
{
  std::shared_ptr<const pcms::Function<Real>> position =
    pcms::test::MakeSplinePosition(pcms::test::MakeUnitSquareGrid(32),
                                   Annulus::Position);
  std::shared_ptr<const pcms::CoordinateMap> map =
    pcms::CreateCoordinateMap(position);
  pcms::test::DeviceCoordinates xi = pcms::test::CreateDeviceCoordinateView(
    Params, map->GetSourceCoordinateSystem(), 2);
};

} // namespace

TEST_CASE("The coordinate map of a position maps parameter points and gives "
          "its Jacobian",
          "[embedding]")
{
  const Fixture f;
  CHECK(f.map->GetTargetCoordinateSystem()->Kind() == "cylindrical_rz");
  CHECK(f.map->GetSourceCoordinateSystem()->Kind() == "xi");
  const auto bound = f.map->Bind(f.xi.coordinate_view);
  const auto mapped = pcms::test::CopyCoordinatesToHost(bound->Mapped());
  const auto jacobian_view = bound->Jacobian();
  Kokkos::View<Real***, pcms::HostMemorySpace> jacobian("j", 4, 2, 2);
  Kokkos::View<Real***, DeviceMemorySpace> device("j", 4, 2, 2);
  Kokkos::parallel_for(
    "copy", Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, 4),
    KOKKOS_LAMBDA(const pcms::LO p) {
      for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
          device(p, i, j) = jacobian_view(p, i, j);
        }
      }
    });
  Kokkos::deep_copy(jacobian, device);
  for (int p = 0; p < 4; ++p) {
    const Real xi1 = Params[std::size_t{2} * p];
    const Real xi2 = Params[2 * p + 1];
    const auto x = Annulus::Position(xi1, xi2);
    const auto j = AnnulusJacobian(xi1, xi2);
    CAPTURE(p);
    CHECK_THAT(mapped(p, 0), WithinAbs(x[0], SplineTol));
    CHECK_THAT(mapped(p, 1), WithinAbs(x[1], SplineTol));
    for (int a = 0; a < 2; ++a) {
      for (int b = 0; b < 2; ++b) {
        CHECK_THAT(jacobian(p, a, b), WithinAbs(j[a][b], 100 * SplineTol));
      }
    }
  }
}

TEST_CASE("Coordinates pass through value transforms", "[embedding]")
{
  const Fixture f;
  const auto bound = f.map->Bind(f.xi.coordinate_view);
  const auto coordinates =
    pcms::ValueBasis::Coordinates(pcms::csys::CreateCartesian(3));
  CHECK(bound->ValueTransforms(coordinates).empty());
  const pcms::ValueBasis vector{pcms::csys::CreateCylindricalRZ(),
                                pcms::ComponentScaling::Physical,
                                pcms::values::Vector};
  for (const auto* transform : bound->ValueTransforms(vector)) {
    CHECK(transform->OutputBasis(coordinates) == coordinates);
  }
  // coordinates are not scalars, and coordinates of different coordinate
  // systems differ
  CHECK_FALSE(pcms::detail::SameValueBasis(coordinates, pcms::ValueBasis{}));
  CHECK_FALSE(pcms::detail::SameValueBasis(
    coordinates,
    pcms::ValueBasis::Coordinates(pcms::csys::CreateCylindricalRThetaZ())));
  CHECK(pcms::detail::SameValueBasis(
    f.position->GetData().GetValueBasis(),
    pcms::ValueBasis::Coordinates(pcms::csys::CreateCylindricalRZ())));
}

TEST_CASE("A function of several components defaults to independent "
          "components",
          "[embedding]")
{
  const auto space =
    pcms::CreateLagrangeSpace(pcms::test::MakeUnitSquareGrid(2), 1, 2);
  CHECK(space->CreateFunction<Real>().GetData().GetValueBasis().kind ==
        pcms::ValueKind::Array);
  CHECK_THROWS_WITH(space->CreateFunction<Real>("x", pcms::ValueBasis{}),
                    ContainsSubstring("a scalar has one component"));
}

TEST_CASE("The coordinate map of a position refuses what it cannot do",
          "[embedding]")
{
  const auto geometry = pcms::test::MakeUnitSquareGrid(4);
  CHECK_THROWS_WITH(pcms::CreateCoordinateMap(nullptr),
                    ContainsSubstring("must not be null"));
  // values that are not coordinates
  const auto space = pcms::CreateSplineSpace(geometry, 2);
  const auto scalars = std::make_shared<pcms::Function<Real>>(
    space->CreateFunction<Real>("x", {.kind = pcms::ValueKind::Array}));
  CHECK_THROWS_WITH(pcms::CreateCoordinateMap(scalars),
                    ContainsSubstring("must be coordinates"));
  // two components are not the coordinates of a 3D coordinate system
  CHECK_THROWS_WITH(
    space->CreateFunction<Real>(
      "x", pcms::ValueBasis::Coordinates(pcms::csys::CreateCartesian(3))),
    ContainsSubstring("require 3 components"));
}
