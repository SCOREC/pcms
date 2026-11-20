#include "pcms/utility/coordinate_system.hpp"
#include "pcms/utility/coordinate_systems/cartesian.hpp"
#include "pcms/utility/coordinate_systems/cylindrical.hpp"
#include "pcms/utility/coordinate_systems/parameter.hpp"
#include "pcms/utility/assert.hpp"
#include "pcms/utility/detail/memory_spaces.hpp"
#include <numbers>
#include <string>
#include <utility>
#include <vector>

namespace pcms
{

BoundCoordinateSystem::MetricView BoundCoordinateSystem::Metric() const noexcept
{
  using Extents = Kokkos::dextents<LO, 3>;
  const LO stride0 =
    metric_.extent(0) == 1 ? 0 : static_cast<LO>(metric_.stride(0));
  return MetricView(
    metric_.data(),
    Kokkos::layout_stride::mapping<Extents>(
      Extents(num_points_, metric_.extent(1), metric_.extent(2)),
      std::array<LO, 3>{stride0, static_cast<LO>(metric_.stride(1)),
                        static_cast<LO>(metric_.stride(2))}));
}

BoundCoordinateSystem::FactorView BoundCoordinateSystem::VolumeFactor()
  const noexcept
{
  using Extents = Kokkos::dextents<LO, 1>;
  return FactorView(
    volume_factor_.data(),
    Kokkos::layout_stride::mapping<Extents>(
      Extents(num_points_),
      std::array<LO, 1>{volume_factor_.extent(0) == 1
                          ? 0
                          : static_cast<LO>(volume_factor_.stride(0))}));
}

namespace
{

// The same metric (the identity) and volume factor 1 at every point.
std::unique_ptr<BoundCoordinateSystem> BindIdentity(int dim, LO num_points)
{
  Kokkos::View<Real***, DeviceMemorySpace> metric("identity_metric", 1, dim,
                                                  dim);
  Kokkos::View<Real*, DeviceMemorySpace> volume("unit_volume_factor", 1);
  Kokkos::deep_copy(volume, 1.0);
  Kokkos::parallel_for(
    "identity_metric",
    Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, dim),
    KOKKOS_LAMBDA(const int i) { metric(0, i, i) = 1.0; });
  return std::make_unique<BoundCoordinateSystem>(num_points, metric, volume);
}

// The identity metric of the (R, Z) plane or the diagonal metric
// (1, R^2, 1) of (R, theta, Z), and the volume factor R (per radian of
// theta for the plane) at points whose first coordinate is R.
std::unique_ptr<BoundCoordinateSystem> BindCylindrical(
  Rank2View<const Real, DeviceMemorySpace> points, bool with_theta)
{
  const LO n = points.extent(0);
  Kokkos::View<Real*, DeviceMemorySpace> volume("cylindrical_volume_factor", n);
  if (!with_theta) {
    detail::ForEach(
      "cylindrical_volume_factor", n,
      KOKKOS_LAMBDA(const LO p) { volume(p) = points(p, 0); });
    Kokkos::View<Real***, DeviceMemorySpace> metric("identity_metric", 1, 2, 2);
    Kokkos::parallel_for(
      "identity_metric",
      Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, 2),
      KOKKOS_LAMBDA(const int i) { metric(0, i, i) = 1.0; });
    return std::make_unique<BoundCoordinateSystem>(n, metric, volume);
  }
  Kokkos::View<Real***, DeviceMemorySpace> metric("cylindrical_metric", n, 3,
                                                  3);
  detail::ForEach(
    "cylindrical_metric", n, KOKKOS_LAMBDA(const LO p) {
      const Real r = points(p, 0);
      metric(p, 0, 0) = 1.0;
      metric(p, 1, 1) = r * r;
      metric(p, 2, 2) = 1.0;
      volume(p) = r;
    });
  return std::make_unique<BoundCoordinateSystem>(n, metric, volume);
}

class Cartesian final : public CoordinateSystem
{
public:
  explicit Cartesian(int dim) noexcept : dim_(dim) {}

  [[nodiscard]] std::string_view Kind() const noexcept override
  {
    return dim_ == 0 ? "cartesian (dimension deferred)" : "cartesian";
  }
  [[nodiscard]] int Dimension() const noexcept override { return dim_; }
  [[nodiscard]] bool IsOrthogonal() const noexcept override { return true; }
  [[nodiscard]] bool IsNormalized() const noexcept override { return true; }
  [[nodiscard]] std::unique_ptr<BoundCoordinateSystem> Bind(
    Rank2View<const Real, DeviceMemorySpace> points) const override
  {
    return BindIdentity(dim_, points.extent(0));
  }
  [[nodiscard]] bool operator==(
    const CoordinateSystem& other) const noexcept override
  {
    const auto* p = dynamic_cast<const Cartesian*>(&other);
    return p != nullptr && p->dim_ == dim_;
  }

private:
  int dim_;
};

class CylindricalRThetaZ final : public CoordinateSystem
{
public:
  [[nodiscard]] std::string_view Kind() const noexcept override
  {
    return "cylindrical_rthetaz";
  }
  [[nodiscard]] int Dimension() const noexcept override { return 3; }
  [[nodiscard]] bool IsOrthogonal() const noexcept override { return true; }
  [[nodiscard]] bool IsNormalized() const noexcept override { return false; }
  [[nodiscard]] std::optional<Real> Period(int axis) const noexcept override
  {
    return axis == 1 ? std::optional<Real>(2.0 * std::numbers::pi)
                     : std::nullopt;
  }
  [[nodiscard]] std::unique_ptr<BoundCoordinateSystem> Bind(
    Rank2View<const Real, DeviceMemorySpace> points) const override
  {
    return BindCylindrical(points, true);
  }
  [[nodiscard]] bool operator==(
    const CoordinateSystem& other) const noexcept override
  {
    return dynamic_cast<const CylindricalRThetaZ*>(&other) != nullptr;
  }
};

// The poloidal plane of an axisymmetric system: its metric is the identity,
// and its volume factor is R, the volume per radian of the toroidal angle, so
// that quadrature weights integrate volumes of revolution (times 2 pi).
class CylindricalRZ final : public CoordinateSystem
{
public:
  [[nodiscard]] std::string_view Kind() const noexcept override
  {
    return "cylindrical_rz";
  }
  [[nodiscard]] int Dimension() const noexcept override { return 2; }
  [[nodiscard]] bool IsOrthogonal() const noexcept override { return true; }
  [[nodiscard]] bool IsNormalized() const noexcept override { return true; }
  [[nodiscard]] std::unique_ptr<BoundCoordinateSystem> Bind(
    Rank2View<const Real, DeviceMemorySpace> points) const override
  {
    return BindCylindrical(points, false);
  }
  [[nodiscard]] bool operator==(
    const CoordinateSystem& other) const noexcept override
  {
    return dynamic_cast<const CylindricalRZ*>(&other) != nullptr;
  }
};

// Coordinates that are labels only: no metric.
class ParameterSystem final : public CoordinateSystem
{
public:
  ParameterSystem(std::string name, int dimension,
                  std::vector<std::optional<Real>> periods)
    : name_(std::move(name)),
      dimension_(dimension),
      periods_(std::move(periods))
  {
  }
  [[nodiscard]] std::string_view Kind() const noexcept override
  {
    return name_;
  }
  [[nodiscard]] int Dimension() const noexcept override { return dimension_; }
  [[nodiscard]] bool IsOrthogonal() const noexcept override { return false; }
  [[nodiscard]] bool IsNormalized() const noexcept override { return false; }
  [[nodiscard]] bool HasMetric() const noexcept override { return false; }
  [[nodiscard]] std::optional<Real> Period(int axis) const noexcept override
  {
    if (axis < 0 || static_cast<std::size_t>(axis) >= periods_.size()) {
      return std::nullopt;
    }
    return periods_[static_cast<std::size_t>(axis)];
  }
  [[nodiscard]] std::unique_ptr<BoundCoordinateSystem> Bind(
    Rank2View<const Real, DeviceMemorySpace>) const override
  {
    throw pcms_error("CoordinateSystem::Bind: the parameter system '" + name_ +
                     "' has no metric; the metric of a geometry over it is "
                     "the geometry's");
  }
  [[nodiscard]] bool operator==(
    const CoordinateSystem& other) const noexcept override
  {
    return this == &other;
  }

private:
  std::string name_;
  int dimension_;
  std::vector<std::optional<Real>> periods_;
};

} // namespace

namespace csys
{

std::shared_ptr<const CoordinateSystem> CreateParameterSystem(
  std::string name, int dimension, std::vector<std::optional<Real>> periods)
{
  if (dimension < 1 || dimension > 6) {
    throw pcms_error(
      "csys::CreateParameterSystem: dimension must be 1..6, got " +
      std::to_string(dimension));
  }
  if (!periods.empty() &&
      periods.size() != static_cast<std::size_t>(dimension)) {
    throw pcms_error("csys::CreateParameterSystem: give no periods or one "
                     "entry per axis");
  }
  for (const auto& period : periods) {
    if (period.has_value() && !(*period > 0)) {
      throw pcms_error("csys::CreateParameterSystem: periods must be positive");
    }
  }
  return std::make_shared<const ParameterSystem>(std::move(name), dimension,
                                                 std::move(periods));
}

std::shared_ptr<const CoordinateSystem> CreateCartesian(int dim)
{
  if (dim < 1 || dim > 3) {
    throw pcms_error("csys::CreateCartesian: dimension must be 1..3, got " +
                     std::to_string(dim));
  }
  return std::make_shared<const Cartesian>(dim);
}

std::shared_ptr<const CoordinateSystem> CreateCylindricalRThetaZ()
{
  return std::make_shared<const CylindricalRThetaZ>();
}

std::shared_ptr<const CoordinateSystem> CreateCylindricalRZ()
{
  return std::make_shared<const CylindricalRZ>();
}

} // namespace csys

namespace detail
{

std::shared_ptr<const CoordinateSystem> DeferredCartesian()
{
  return std::make_shared<const Cartesian>(0);
}

bool IsCartesian(const CoordinateSystem& system) noexcept
{
  return dynamic_cast<const Cartesian*>(&system) != nullptr;
}

std::shared_ptr<const CoordinateSystem> ResolveCoordinateSystem(
  std::shared_ptr<const CoordinateSystem> system, int dim)
{
  if (system == nullptr) {
    throw pcms_error(
      "ResolveCoordinateSystem: coordinate system must not be null");
  }
  // To avoid requiring the user specifying the dimension twice in
  // FunctionSpace constructors, we have a default
  // Cartesian constructor that constructs with dimension 0.
  if (IsCartesian(*system) && system->Dimension() == 0) {
    return csys::CreateCartesian(dim);
  }
  if (system->Dimension() != dim) {
    throw pcms_error("ResolveCoordinateSystem: coordinate system '" +
                     std::string(system->Kind()) + "' has " +
                     std::to_string(system->Dimension()) +
                     " coordinate columns but the data has " +
                     std::to_string(dim));
  }
  return system;
}

} // namespace detail
} // namespace pcms
