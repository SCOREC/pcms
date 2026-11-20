#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <Kokkos_Core.hpp>
#include <pcms/utility/coordinate_systems/cartesian.hpp>
#include <pcms/utility/coordinate_systems/cylindrical.hpp>
#include <pcms/utility/coordinate_systems/parameter.hpp>
#include "field_test_utils.hpp"
#include <numbers>
#include <optional>

using Catch::Matchers::WithinAbs;
using pcms::DeviceMemorySpace;
using pcms::LO;
using pcms::Real;

namespace
{

// The metric and volume factor of `system` at the points, on the host.
struct BoundOnHost
{
  std::vector<Real> volume;
  std::vector<Real> metric; // [point][i][j]
  bool metric_broadcast;
};

BoundOnHost BindOnHost(const pcms::CoordinateSystem& system,
                       const std::vector<Real>& points)
{
  const int dim = system.Dimension();
  const auto n = static_cast<LO>(points.size() / dim);
  const auto coords = pcms::test::CreateDeviceRank2View(points, dim);
  const auto bound = system.Bind(pcms::MakeConstRank2View(coords));
  Kokkos::View<Real*, DeviceMemorySpace> volume("volume", n);
  Kokkos::View<Real***, DeviceMemorySpace> metric("metric", n, dim, dim);
  const auto v = bound->VolumeFactor();
  const auto g = bound->Metric();
  Kokkos::parallel_for(
    "copy", Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, n),
    KOKKOS_LAMBDA(const LO p) {
      volume(p) = v(p);
      for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
          metric(p, i, j) = g(p, i, j);
        }
      }
    });
  const auto v_host =
    Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace(), volume);
  const auto g_host =
    Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace(), metric);
  BoundOnHost out{{}, {}, g.stride(0) == 0};
  for (LO p = 0; p < n; ++p) {
    out.volume.push_back(v_host(p));
    for (int i = 0; i < dim; ++i) {
      for (int j = 0; j < dim; ++j) {
        out.metric.push_back(g_host(p, i, j));
      }
    }
  }
  return out;
}

} // namespace

TEST_CASE("Bound coordinate systems give the metric and the volume factor",
          "[coordinate_system][metric]")
{
  SECTION("Cartesian: identity, allocated once")
  {
    const auto bound =
      BindOnHost(*pcms::csys::CreateCartesian(2), {0.1, 0.2, 3.0, -4.0});
    CHECK(bound.metric_broadcast);
    CHECK(bound.volume == std::vector<Real>{1.0, 1.0});
    CHECK(bound.metric ==
          std::vector<Real>{1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0});
  }
  SECTION("cylindrical (R, theta, Z): diag(1, R^2, 1), volume factor R")
  {
    const auto system = pcms::csys::CreateCylindricalRThetaZ();
    const auto bound = BindOnHost(*system, {2.0, 0.3, 1.0, 3.0, -1.0, 0.0});
    CHECK_FALSE(bound.metric_broadcast);
    CHECK(bound.volume == std::vector<Real>{2.0, 3.0});
    CHECK(bound.metric[4] == 4.0);
    CHECK(bound.metric[9 + 4] == 9.0);
    CHECK(bound.metric[1] == 0.0);
    CHECK(system->Period(1) == std::optional<Real>(2.0 * std::numbers::pi));
    CHECK_FALSE(system->Period(0).has_value());
  }
  SECTION("cylindrical (R, Z): identity metric, volume factor R")
  {
    const auto bound =
      BindOnHost(*pcms::csys::CreateCylindricalRZ(), {2.0, 0.3, 3.0, 1.0});
    CHECK(bound.metric_broadcast);
    CHECK(bound.volume == std::vector<Real>{2.0, 3.0});
    CHECK(bound.metric ==
          std::vector<Real>{1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0});
  }
}

TEST_CASE("Parameter systems are labels without a metric",
          "[coordinate_system][metric]")
{
  const auto system = pcms::csys::CreateParameterSystem(
    "r_theta_zeta", 3, {std::nullopt, 2.0 * std::numbers::pi, 1.5});
  CHECK(system->Kind() == "r_theta_zeta");
  CHECK(system->Dimension() == 3);
  CHECK_FALSE(system->HasMetric());
  CHECK_FALSE(system->IsOrthogonal());
  CHECK_FALSE(system->Period(0).has_value());
  CHECK(system->Period(1) == std::optional<Real>(2.0 * std::numbers::pi));
  CHECK(system->Period(2) == std::optional<Real>(1.5));
  CHECK(pcms::csys::CreateCartesian(2)->HasMetric());
  const auto coords = pcms::test::CreateDeviceRank2View({0.1, 0.2, 0.3}, 3);
  CHECK_THROWS_AS(system->Bind(pcms::MakeConstRank2View(coords)),
                  pcms::pcms_error);
  // equal by identity
  CHECK(*system == *system);
  CHECK_FALSE(*system == *pcms::csys::CreateParameterSystem(
                           "r_theta_zeta", 3,
                           {std::nullopt, 2.0 * std::numbers::pi, 1.5}));
  CHECK_FALSE(*pcms::csys::CreateParameterSystem("x", 2) ==
              *pcms::csys::CreateCartesian(2));
  CHECK_THROWS(pcms::csys::CreateParameterSystem("x", 0));
  CHECK_THROWS(pcms::csys::CreateParameterSystem("x", 2, {1.0}));
  CHECK_THROWS(pcms::csys::CreateParameterSystem("x", 1, {-1.0}));
}
