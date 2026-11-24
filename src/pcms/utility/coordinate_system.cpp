#include "pcms/utility/coordinate_system.hpp"
#include "pcms/utility/coordinate_systems/cartesian.hpp"
#include "pcms/utility/coordinate_systems/cylindrical.hpp"
#include "pcms/utility/coordinate_systems/parameter.hpp"
#include "pcms/utility/assert.hpp"
#include "pcms/utility/detail/memory_spaces.hpp"
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

/// @throws pcms_error if a periodic axis has an axis outside `dimension` or
///         a length that is not positive
std::vector<PeriodicAxis> CheckPeriodicAxes(
  std::vector<PeriodicAxis> periodic_axes, int dimension,
  const std::string& caller)
{
  for (const auto& id : periodic_axes) {
    if (id.axis < 0 || id.axis >= dimension || !(id.length > 0)) {
      throw pcms_error(caller + ": a periodic axis needs an axis of the "
                                "system and a positive length");
    }
  }
  return periodic_axes;
}

class Cartesian final : public CoordinateSystem
{
public:
  explicit Cartesian(int dim, std::vector<PeriodicAxis> periodic_axes = {})
    : dim_(dim), periodic_axes_(std::move(periodic_axes))
  {
  }

  [[nodiscard]] std::string_view Kind() const noexcept override
  {
    return dim_ == 0 ? "cartesian (dimension deferred)" : "cartesian";
  }
  [[nodiscard]] int Dimension() const noexcept override { return dim_; }
  [[nodiscard]] bool IsOrthonormal() const noexcept override { return true; }
  [[nodiscard]] std::span<const PeriodicAxis> PeriodicAxes()
    const noexcept override
  {
    return periodic_axes_;
  }
  [[nodiscard]] std::unique_ptr<BoundCoordinateSystem> Bind(
    Rank2View<const Real, DeviceMemorySpace> points) const override
  {
    return BindIdentity(dim_, points.extent(0));
  }
  [[nodiscard]] bool operator==(
    const CoordinateSystem& other) const noexcept override
  {
    const auto* p = dynamic_cast<const Cartesian*>(&other);
    return p != nullptr && p->dim_ == dim_ &&
           p->periodic_axes_ == periodic_axes_;
  }

private:
  int dim_;
  std::vector<PeriodicAxis> periodic_axes_;
};

class CylindricalRThetaZ final : public CoordinateSystem
{
public:
  explicit CylindricalRThetaZ(Real phi_period)
    : periodic_axes_{PeriodicAxis{1, 0, phi_period, nullptr}}
  {
  }
  [[nodiscard]] std::string_view Kind() const noexcept override
  {
    return "cylindrical_rthetaz";
  }
  [[nodiscard]] int Dimension() const noexcept override { return 3; }
  [[nodiscard]] bool IsOrthonormal() const noexcept override { return false; }
  [[nodiscard]] std::span<const PeriodicAxis> PeriodicAxes()
    const noexcept override
  {
    return periodic_axes_;
  }
  [[nodiscard]] std::unique_ptr<BoundCoordinateSystem> Bind(
    Rank2View<const Real, DeviceMemorySpace> points) const override
  {
    return BindCylindrical(points, true);
  }
  [[nodiscard]] bool operator==(
    const CoordinateSystem& other) const noexcept override
  {
    const auto* p = dynamic_cast<const CylindricalRThetaZ*>(&other);
    return p != nullptr && p->periodic_axes_ == periodic_axes_;
  }

private:
  std::vector<PeriodicAxis> periodic_axes_;
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
  [[nodiscard]] bool IsOrthonormal() const noexcept override { return true; }
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
                  std::vector<PeriodicAxis> periodic_axes)
    : name_(std::move(name)),
      dimension_(dimension),
      periodic_axes_(std::move(periodic_axes))
  {
  }
  [[nodiscard]] std::string_view Kind() const noexcept override
  {
    return name_;
  }
  [[nodiscard]] int Dimension() const noexcept override { return dimension_; }
  [[nodiscard]] bool IsOrthonormal() const noexcept override { return false; }
  [[nodiscard]] bool HasMetric() const noexcept override { return false; }
  [[nodiscard]] std::span<const PeriodicAxis> PeriodicAxes()
    const noexcept override
  {
    return periodic_axes_;
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
  std::vector<PeriodicAxis> periodic_axes_;
};

} // namespace

namespace csys
{

std::shared_ptr<const CoordinateSystem> CreateParameterSystem(
  std::string name, int dimension, std::vector<PeriodicAxis> periodic_axes)
{
  if (dimension < 1 || dimension > 6) {
    throw pcms_error(
      "csys::CreateParameterSystem: dimension must be 1..6, got " +
      std::to_string(dimension));
  }
  return std::make_shared<const ParameterSystem>(
    std::move(name), dimension,
    CheckPeriodicAxes(std::move(periodic_axes), dimension,
                      "csys::CreateParameterSystem"));
}

std::shared_ptr<const CoordinateSystem> CreateCartesian(
  int dim, std::vector<PeriodicAxis> periodic_axes)
{
  if (dim < 1 || dim > 3) {
    throw pcms_error("csys::CreateCartesian: dimension must be 1..3, got " +
                     std::to_string(dim));
  }
  return std::make_shared<const Cartesian>(
    dim,
    CheckPeriodicAxes(std::move(periodic_axes), dim, "csys::CreateCartesian"));
}

std::shared_ptr<const CoordinateSystem> CreateCylindricalRThetaZ(
  Real phi_period)
{
  if (!(phi_period > 0)) {
    throw pcms_error(
      "csys::CreateCylindricalRThetaZ: the phi period must be positive");
  }
  return std::make_shared<const CylindricalRThetaZ>(phi_period);
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
