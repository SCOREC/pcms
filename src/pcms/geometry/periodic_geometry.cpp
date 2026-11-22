#include "pcms/geometry/periodic_geometry.hpp"
#include "pcms/geometry/detail/geometry.hpp"
#include "pcms/utility/arrays.hpp"
#include "pcms/utility/detail/arrays.hpp"
#include "pcms/utility/assert.hpp"
#include <cmath>
#include <variant>
#include <vector>

#include "periodic_binding.hpp"

namespace pcms
{

void detail::GaussLegendre(int order, std::vector<Real>& nodes,
                           std::vector<Real>& weights)
{
  const int n = order / 2 + 1;
  nodes.assign(n, 0);
  weights.assign(n, 0);
  for (int i = 0; i < n; ++i) {
    Real x = std::cos(std::numbers::pi_v<Real> * (i + 0.75) / (n + 0.5));
    Real dp = 1;
    for (int it = 0; it < 100; ++it) {
      Real p0 = 1;
      Real p1 = x;
      for (int k = 2; k <= n; ++k) {
        const Real p2 = ((2 * k - 1) * x * p1 - (k - 1) * p0) / k;
        p0 = p1;
        p1 = p2;
      }
      dp = n * (x * p1 - p0) / (x * x - 1);
      const Real dx = p1 / dp;
      x -= dx;
      if (std::abs(dx) < 1e-16) {
        break;
      }
    }
    nodes[n - 1 - i] = x;
    weights[n - 1 - i] = 2 / ((1 - x * x) * dp * dp);
  }
}

PeriodicGeometry::PeriodicGeometry(
  std::shared_ptr<const CoordinateSystem> system, LO num_cells, Real period,
  Real origin)
  : system_(pcms::detail::ResolveCoordinateSystem(std::move(system), 1)),
    num_cells_(num_cells),
    period_(period),
    origin_(origin)
{
  if (system_->Dimension() != 1 || num_cells_ < 1 || !(period_ > 0)) {
    throw pcms_error("CreatePeriodicGeometry: needs a one-dimensional system, "
                     "a positive number of cells and a positive period");
  }
}

std::shared_ptr<BoundPoints> PeriodicGeometry::BindImpl(
  const BindRequest& request) const
{
  return std::make_shared<detail::PeriodicBinding>(
    shared_from_this(), request, num_cells_, origin_, period_);
}

Quadrature PeriodicGeometry::CreateQuadrature(
  const pcms::detail::QuadratureRule& rule) const
{
  const auto* trapezoid = std::get_if<detail::TrapezoidRule>(&rule);
  const auto* gauss_rule = std::get_if<detail::GaussRule>(&rule);
  if (trapezoid == nullptr && gauss_rule == nullptr) {
    return Geometry::CreateQuadrature(rule);
  }
  // per cell, the trapezoid rule's one node or the Gauss nodes, on [-1, 1]
  std::vector<Real> nodes = {-1};
  std::vector<Real> gauss = {2};
  const LO num_cells =
    trapezoid != nullptr ? trapezoid->num_points : num_cells_;
  if (gauss_rule != nullptr) {
    detail::GaussLegendre(gauss_rule->order, nodes, gauss);
  }
  const auto n = static_cast<int>(nodes.size());
  const Real h = period_ / num_cells;
  Kokkos::View<Real**, HostMemorySpace> points(
    "periodic_quadrature_points", static_cast<std::size_t>(num_cells) * n, 1);
  Kokkos::View<Real*, HostMemorySpace> weights(
    "periodic_quadrature_weights", static_cast<std::size_t>(num_cells) * n);
  for (LO c = 0; c < num_cells; ++c) {
    for (int i = 0; i < n; ++i) {
      points(c * n + i, 0) = origin_ + (c + 0.5 * (1 + nodes[i])) * h;
      weights(c * n + i) = 0.5 * h * gauss[i];
    }
  }
  const auto device_points = pcms::detail::ToDeviceView(points);
  return detail::BindQuadrature(
    *this,
    CoordinateView<DeviceMemorySpace>(system_,
                                      MakeConstRank2View(device_points)),
    Kokkos::create_mirror_view_and_copy(DeviceMemorySpace{}, weights));
}

std::shared_ptr<const Geometry> CreatePeriodicGeometry(
  std::shared_ptr<const CoordinateSystem> system, int num_cells,
  std::optional<Real> period, Real origin)
{
  if (system != nullptr && system->Dimension() == 1) {
    const auto declared = system->Period(0);
    if (declared && period && *declared != *period) {
      throw pcms_error("CreatePeriodicGeometry: the period differs from the "
                       "period the coordinate system declares for its axis");
    }
    if (declared) {
      period = declared;
    }
  }
  return std::make_shared<PeriodicGeometry>(
    std::move(system), num_cells, period.value_or(2 * std::numbers::pi_v<Real>),
    origin);
}

Quadrature CreateTrapezoidQuadrature(
  const std::shared_ptr<const Geometry>& geometry, int num_points)
{
  if (num_points < 1) {
    throw pcms_error("CreateTrapezoidQuadrature: needs a positive number of "
                     "points");
  }
  return detail::CreateQuadrature(*geometry, detail::TrapezoidRule{num_points});
}

} // namespace pcms
