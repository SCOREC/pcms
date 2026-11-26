#ifndef PCMS_TEST_EVALUATOR_TEST_UTILS_HPP
#define PCMS_TEST_EVALUATOR_TEST_UTILS_HPP

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
#include "pcms/field/point_rows.hpp"
#include "field_test_utils.hpp"

// Helpers of the evaluator tests: DOFs and points given row by row,
// evaluation at points, the adjoint identity and the evaluation rows.
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
  return pcms::test::ToDevice(host);
}

/// Sets the DOFs of `field` from `c`, holder by holder, components fastest.
inline void SetDofs(Field<Real>& field, const std::vector<Real>& c)
{
  const auto nc =
    static_cast<std::size_t>(field.GetLayout().GetNumComponents());
  Kokkos::View<Real**, HostMemorySpace> values("c", c.size() / nc, nc);
  for (std::size_t i = 0; i < c.size(); ++i) {
    values(i / nc, i % nc) = c[i];
  }
  field.SetDOFHolderDataHost(ValueView<const Real, HostMemorySpace>(
    field.GetData().GetValueBasis(), MakeConstRank2View(values)));
}

/// @return the DOFs of `field`, holder by holder, components fastest
inline std::vector<Real> GetDofs(const Field<Real>& field)
{
  const auto values = field.GetDOFHolderDataHost().GetValues();
  std::vector<Real> out;
  for (std::size_t i = 0; i < values.extent(0); ++i) {
    for (std::size_t c = 0; c < values.extent(1); ++c) {
      out.push_back(values(i, c));
    }
  }
  return out;
}

/// The values (one row per point, columns of the operator's output) of the
/// field at the points, through a bound evaluator.
inline std::vector<std::vector<Real>> EvaluateAt(
  const FunctionSpace& space, const Field<Real>& field,
  const std::vector<std::vector<Real>>& points,
  DifferentialOperator op = DifferentialOperator::Value,
  OutOfBoundsPolicy policy = {})
{
  const auto device = DevicePoints(points);
  const auto bound = space.GetGeometry()->Bind(
    CoordinateView<DeviceMemorySpace>(
      space.GetGeometry()->GetCoordinateSystem(), MakeConstRank2View(device)),
    policy);
  const auto evaluator = space.CreatePointEvaluator<Real>(bound, op);
  const auto basis = evaluator->OutputBasis(field.GetData().GetValueBasis());
  const auto columns =
    space.GetLayout()->GetNumComponents() *
    (op == DifferentialOperator::Value
       ? 1
       : space.GetGeometry()->GetCoordinateSystem()->Dimension());
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

/// Checks <E x, y> = <x, E^T y> for random x and y, and that the evaluation
/// rows are E.
inline void CheckAdjointAndRows(
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
  const auto device_y = pcms::test::ToDevice(host_y);
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

  // overwrites: a second transpose gives the same functional
  evaluator->EvaluateTranspose(
    ValueView<const Real, DeviceMemorySpace>(
      functional.GetData().GetValueBasis(), MakeConstRank2View(device_y)),
    functional);
  const auto twice = GetDofs(functional);
  for (std::size_t i = 0; i < num_dofs; ++i) {
    CHECK(twice[i] == Catch::Approx(ETy[i]).margin(1e-13));
  }

  // E from the evaluation rows
  const auto* source = pcms::detail::AsRowsSource(*evaluator);
  REQUIRE(source != nullptr);
  const auto& rows = source->Rows();
  const auto offsets =
    Kokkos::create_mirror_view_and_copy(HostMemorySpace{}, rows.offsets);
  const auto holders =
    Kokkos::create_mirror_view_and_copy(HostMemorySpace{}, rows.holders);
  const auto values =
    Kokkos::create_mirror_view_and_copy(HostMemorySpace{}, rows.values);
  for (std::size_t p = 0; p < points.size(); ++p) {
    Real sum = 0;
    for (LO j = offsets(p); j < offsets(p + 1); ++j) {
      sum += values(j) * x[holders(j)];
    }
    CHECK(sum == Catch::Approx(Ex[p][0]).epsilon(1e-12).margin(1e-12));
  }
}

} // namespace pcms::test

#endif // PCMS_TEST_EVALUATOR_TEST_UTILS_HPP
