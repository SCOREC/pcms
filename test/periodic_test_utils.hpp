#ifndef PCMS_TEST_PERIODIC_TEST_UTILS_HPP
#define PCMS_TEST_PERIODIC_TEST_UTILS_HPP

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <Kokkos_Core.hpp>
#include <pcms/field/function_spaces.hpp>
#include <pcms/utility/coordinate_systems/cartesian.hpp>
#include <pcms/utility/arrays.hpp>
#include <pcms/utility/detail/arrays.hpp>
#include <cmath>
#include <memory>
#include <numbers>
#include <vector>

// Helpers of the tests of the periodic and product spaces: evaluation at
// points given row by row, the adjoint identity and the evaluation matrix.
namespace pcms::test
{

inline Kokkos::View<Real**, DeviceMemorySpace> DevicePoints(
  const std::vector<std::vector<Real>>& points)
{
  Kokkos::View<Real**, HostMemorySpace> host("points", points.size(),
                                             points.at(0).size());
  for (std::size_t p = 0; p < points.size(); ++p) {
    for (std::size_t d = 0; d < points[p].size(); ++d) {
      host(p, d) = points[p][d];
    }
  }
  return pcms::detail::ToDeviceView(host);
}

inline void SetDofs(Field<Real>& field, const std::vector<Real>& c)
{
  Kokkos::View<Real**, HostMemorySpace> values("c", c.size(), 1);
  for (std::size_t i = 0; i < c.size(); ++i) {
    values(i, 0) = c[i];
  }
  field.SetDOFHolderDataHost(ValueView<const Real, HostMemorySpace>(
    field.GetData().GetValueBasis(), MakeConstRank2View(values)));
}

inline std::vector<Real> GetDofs(const Field<Real>& field)
{
  const auto values = field.GetDOFHolderDataHost().GetValues();
  std::vector<Real> out(values.extent(0));
  for (std::size_t i = 0; i < out.size(); ++i) {
    out[i] = values(i, 0);
  }
  return out;
}

/// The values (one row per point, columns of the operator's output) of the
/// field at the points, through a bound evaluator.
inline std::vector<std::vector<Real>> EvaluateAt(
  const FunctionSpace& space, const Field<Real>& field,
  const std::vector<std::vector<Real>>& points,
  DifferentialOperator op = DifferentialOperator::Value)
{
  const auto device = DevicePoints(points);
  const auto bound =
    space.GetGeometry()->Bind(CoordinateView<DeviceMemorySpace>(
      space.GetGeometry()->GetCoordinateSystem(), MakeConstRank2View(device)));
  const auto evaluator = space.CreatePointEvaluator<Real>(bound, op);
  const auto basis = evaluator->OutputBasis(field.GetData().GetValueBasis());
  const auto columns =
    op == DifferentialOperator::Value
      ? 1
      : space.GetGeometry()->GetCoordinateSystem()->Dimension();
  Kokkos::View<Real**, DeviceMemorySpace> out("out", points.size(), columns);
  evaluator->Evaluate(
    field, ValueView<Real, DeviceMemorySpace>(basis, MakeRank2View(out)));
  const auto host = Kokkos::create_mirror_view_and_copy(HostMemorySpace{}, out);
  std::vector<std::vector<Real>> result(points.size());
  for (std::size_t p = 0; p < points.size(); ++p) {
    for (int c = 0; c < columns; ++c) {
      result[p].push_back(host(p, c));
    }
  }
  return result;
}

/// Deterministic pseudo-random numbers in (-1, 1)
inline std::vector<Real> Pseudo(std::size_t n, int seed)
{
  std::vector<Real> out(n);
  for (std::size_t i = 0; i < n; ++i) {
    out[i] = std::sin(12.9898 * (i + 1) + 78.233 * seed) * 0.99;
  }
  return out;
}

/// Checks <E x, y> = <x, E^T y> for random x and y, and that the matrix of
/// GetMatrix is E.
inline void CheckAdjointAndMatrix(
  const std::shared_ptr<const FunctionSpace>& space,
  const std::vector<std::vector<Real>>& points)
{
  const auto num_dofs =
    static_cast<std::size_t>(space->GetLayout()->GetNumDofHolder());
  auto u = space->CreateFunction<Real>();
  const auto x = Pseudo(num_dofs, 1);
  SetDofs(u, x);
  const auto Ex = EvaluateAt(*space, u, points);
  const auto y = Pseudo(points.size(), 2);

  const auto device = DevicePoints(points);
  const auto bound =
    space->GetGeometry()->Bind(CoordinateView<DeviceMemorySpace>(
      space->GetGeometry()->GetCoordinateSystem(), MakeConstRank2View(device)));
  const auto evaluator = space->CreatePointEvaluator<Real>(bound);
  REQUIRE(evaluator->SupportsTranspose());
  Kokkos::View<Real**, HostMemorySpace> host_y("y", points.size(), 1);
  for (std::size_t p = 0; p < points.size(); ++p) {
    host_y(p, 0) = y[p];
  }
  const auto device_y = pcms::detail::ToDeviceView(host_y);
  auto functional = space->CreateFunction<Real>();
  evaluator->EvaluateTranspose(
    ValueView<const Real, DeviceMemorySpace>(
      functional.GetData().GetValueBasis(), MakeConstRank2View(device_y)),
    functional);
  const auto ETy = GetDofs(functional);
  Real lhs = 0;
  Real rhs = 0;
  for (std::size_t p = 0; p < points.size(); ++p) {
    lhs += Ex[p][0] * y[p];
  }
  for (std::size_t i = 0; i < num_dofs; ++i) {
    rhs += x[i] * ETy[i];
  }
  CHECK(lhs == Catch::Approx(rhs).epsilon(1e-12).margin(1e-12));

  // accumulates: a second transpose doubles the functional
  evaluator->EvaluateTranspose(
    ValueView<const Real, DeviceMemorySpace>(
      functional.GetData().GetValueBasis(), MakeConstRank2View(device_y)),
    functional);
  const auto twice = GetDofs(functional);
  for (std::size_t i = 0; i < num_dofs; ++i) {
    CHECK(twice[i] == Catch::Approx(2 * ETy[i]).margin(1e-13));
  }

  // E from the matrix: columns are gids, which are the holders in these
  // layouts
  REQUIRE(evaluator->SupportsMatrix());
  const auto matrix = evaluator->GetMatrix();
  const auto offsets =
    Kokkos::create_mirror_view_and_copy(HostMemorySpace{}, matrix.row_offsets);
  const auto columns =
    Kokkos::create_mirror_view_and_copy(HostMemorySpace{}, matrix.columns);
  const auto values =
    Kokkos::create_mirror_view_and_copy(HostMemorySpace{}, matrix.values);
  const auto layout_gids = space->GetLayout()->GetGidsHost();
  std::vector<Real> by_gid(num_dofs, 0.0);
  for (std::size_t i = 0; i < num_dofs; ++i) {
    by_gid[layout_gids(i)] = x[i];
  }
  for (std::size_t p = 0; p < points.size(); ++p) {
    Real sum = 0;
    for (LO j = offsets(p); j < offsets(p + 1); ++j) {
      sum += values(j) * by_gid[columns(j)];
    }
    CHECK(sum == Catch::Approx(Ex[p][0]).epsilon(1e-12).margin(1e-12));
  }
}

} // namespace pcms::test

#endif // PCMS_TEST_PERIODIC_TEST_UTILS_HPP
