#ifndef PCMS_TEST_EMBEDDING_TEST_UTILS_HPP
#define PCMS_TEST_EMBEDDING_TEST_UTILS_HPP

#include <pcms/backends/uniform_grid/geometry.hpp>
#include <pcms/field/embedding.hpp>
#include <pcms/field/function_spaces.hpp>
#include <pcms/utility/coordinate_systems/cartesian.hpp>
#include <pcms/utility/coordinate_systems/cylindrical.hpp>
#include <pcms/utility/coordinate_systems/parameter.hpp>
#include "field_test_utils.hpp"
#include <array>
#include <cmath>

namespace pcms::test
{

// A sector of an annulus around (Rc, 0) in the (R, Z) plane, the image of the
// unit square: rho = 0.5 + 0.5 xi1, theta = 0.2 + 2 xi2,
// (R, Z) = (Rc + rho cos theta, rho sin theta). Its parameter metric is
// diag(0.25, 4 rho^2) and its area is 0.75.
struct Annulus
{
  static constexpr Real Rc = 2.0;
  static Real Rho(Real xi1) { return 0.5 + 0.5 * xi1; }
  static Real Theta(Real xi2) { return 0.2 + 2.0 * xi2; }
  static std::array<Real, 2> Position(Real xi1, Real xi2)
  {
    const Real rho = Rho(xi1);
    return {Rc + rho * std::cos(Theta(xi2)), rho * std::sin(Theta(xi2))};
  }
  // the exact integral of R over the sector, the volume of revolution per
  // radian
  static Real MomentOfR()
  {
    return Rc * 0.75 +
           (1.0 - 0.125) / 3.0 * (std::sin(Theta(1.0)) - std::sin(Theta(0.0)));
  }
};

inline std::shared_ptr<const Geometry> MakeUnitSquareGrid(int divisions)
{
  UniformGrid<2> grid;
  grid.bot_left = {0.0, 0.0};
  grid.edge_length = {1.0, 1.0};
  grid.divisions = {divisions, divisions};
  return uniform_grid::CreateGeometry(grid,
                                      csys::CreateParameterSystem("xi", 2));
}

/// The spline position Function of a 2D map xi -> position(xi) on `geometry`,
/// whose values are points of `system`.
template <typename Map>
std::shared_ptr<const Function<Real>> MakeSplinePosition(
  const std::shared_ptr<const Geometry>& geometry, Map position,
  std::shared_ptr<const CoordinateSystem> system = csys::CreateCylindricalRZ())
{
  const auto space = CreateSplineSpace(geometry, 2);
  auto function = std::make_shared<Function<Real>>(space->CreateFunction<Real>(
    "position", ValueBasis::Coordinates(std::move(system))));
  const auto xi = CopyCoordinatesToHost(
    space->GetLayout()->GetDOFHolderCoordinates().GetValues());
  Kokkos::View<Real**, HostMemorySpace> values("position", xi.extent(0), 2);
  for (size_t i = 0; i < xi.extent(0); ++i) {
    const auto x = position(xi(i, 0), xi(i, 1));
    values(i, 0) = x[0];
    values(i, 1) = x[1];
  }
  function->SetDOFHolderDataHost(ValueView<const Real, HostMemorySpace>(
    function->GetData().GetValueBasis(), MakeConstRank2View(values)));
  return function;
}

} // namespace pcms::test

#endif // PCMS_TEST_EMBEDDING_TEST_UTILS_HPP
