#include "pcms/geometry/periodic_geometry.hpp"
#include "pcms/utility/arrays.hpp"
#include "pcms/utility/detail/arrays.hpp"
#include "pcms/utility/assert.hpp"
#include <cmath>
#include <variant>
#include <vector>

#include "periodic_binding.hpp"

namespace pcms
{

namespace
{

/// Gauss-Legendre nodes (ascending) and weights on [-1, 1].
void GaussLegendre(int n, std::vector<Real>& nodes, std::vector<Real>& weights)
{
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

} // namespace

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

Quadrature PeriodicGeometry::MakeQuadrature(
  Kokkos::View<Real**, HostMemorySpace> points,
  Kokkos::View<Real*, HostMemorySpace> weights) const
{
  const auto device_points = pcms::detail::ToDeviceView(points);
  const auto device_weights =
    Kokkos::create_mirror_view_and_copy(DeviceMemorySpace{}, weights);
  return Quadrature{Bind(CoordinateView<DeviceMemorySpace>(
                      system_, MakeConstRank2View(device_points))),
                    device_weights,
                    {}};
}

Quadrature PeriodicGeometry::CreateQuadrature(
  const pcms::detail::QuadratureRule& rule) const
{
  if (const auto* trapezoid = std::get_if<detail::TrapezoidRule>(&rule)) {
    const int num_points = trapezoid->num_points;
    Kokkos::View<Real**, HostMemorySpace> points("trapezoid_points", num_points,
                                                 1);
    Kokkos::View<Real*, HostMemorySpace> weights("trapezoid_weights",
                                                 num_points);
    for (int i = 0; i < num_points; ++i) {
      points(i, 0) = origin_ + period_ * i / num_points;
      weights(i) = period_ / num_points;
    }
    return MakeQuadrature(points, weights);
  }
  const auto* gauss_rule = std::get_if<detail::GaussRule>(&rule);
  if (gauss_rule == nullptr) {
    return Geometry::CreateQuadrature(rule);
  }
  const int n = std::max(1, gauss_rule->order / 2 + 1);
  std::vector<Real> nodes;
  std::vector<Real> gauss;
  GaussLegendre(n, nodes, gauss);
  const Real h = period_ / num_cells_;
  Kokkos::View<Real**, HostMemorySpace> points(
    "gauss_points", static_cast<std::size_t>(num_cells_) * n, 1);
  Kokkos::View<Real*, HostMemorySpace> weights(
    "gauss_weights", static_cast<std::size_t>(num_cells_) * n);
  for (LO c = 0; c < num_cells_; ++c) {
    for (int i = 0; i < n; ++i) {
      points(c * n + i, 0) = origin_ + (c + 0.5 * (1 + nodes[i])) * h;
      weights(c * n + i) = 0.5 * h * gauss[i];
    }
  }
  return MakeQuadrature(points, weights);
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
