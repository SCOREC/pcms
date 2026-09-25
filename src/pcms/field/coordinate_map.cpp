#include "pcms/field/coordinate_map.hpp"
#include "pcms/field/coordinate_systems/cartesian.hpp"
#include "pcms/field/coordinate_systems/cylindrical.hpp"
#include "pcms/utility/assert.h"
#include "pcms/utility/profile.h"
#include <Kokkos_Core.hpp>
#include <string>

namespace pcms
{

namespace
{

using ExecutionSpace = DeviceMemorySpace::execution_space;

const ValueBasis& CartesianPhysicalVector()
{
  static const ValueBasis basis{csys::Cartesian::Create(3),
                                ComponentScaling::Physical,
                                VarianceSignature{Variance::Contravariant}};
  return basis;
}

const ValueBasis& CylindricalPhysicalVector()
{
  static const ValueBasis basis{csys::CylindricalRThetaZ::Create(),
                                ComponentScaling::Physical,
                                VarianceSignature{Variance::Contravariant}};
  return basis;
}

Kokkos::View<Real**, DeviceMemorySpace> AllocateCoordinates(size_t num_points,
                                                            int dim)
{
  return Kokkos::View<Real**, DeviceMemorySpace>(
    Kokkos::view_alloc(Kokkos::WithoutInitializing, "mapped_points"),
    num_points, static_cast<size_t>(dim));
}

Kokkos::View<Real*, DeviceMemorySpace> AllocateTrig(const std::string& label,
                                                    size_t num_points)
{
  return Kokkos::View<Real*, DeviceMemorySpace>(
    Kokkos::view_alloc(Kokkos::WithoutInitializing, label), num_points);
}

// Kernels as free functions: CUDA extended lambdas may not be defined inside
// member functions of classes with private members.
void BindCartesianToCylindrical(
  Rank2View<const Real, DeviceMemorySpace> in,
  Kokkos::View<Real**, DeviceMemorySpace> out,
  Kokkos::View<Real*, DeviceMemorySpace> cos_theta,
  Kokkos::View<Real*, DeviceMemorySpace> sin_theta)
{
  Kokkos::parallel_for(
    "cartesian_to_cylindrical_bind",
    Kokkos::RangePolicy<ExecutionSpace>(0, static_cast<LO>(in.extent(0))),
    KOKKOS_LAMBDA(const LO i) {
      const Real x = in(i, 0);
      const Real y = in(i, 1);
      const Real theta = Kokkos::atan2(y, x);
      out(i, 0) = Kokkos::sqrt(x * x + y * y);
      out(i, 1) = theta;
      out(i, 2) = in(i, 2);
      cos_theta(i) = Kokkos::cos(theta);
      sin_theta(i) = Kokkos::sin(theta);
    });
}

void BindCylindricalToCartesian(
  Rank2View<const Real, DeviceMemorySpace> in,
  Kokkos::View<Real**, DeviceMemorySpace> out,
  Kokkos::View<Real*, DeviceMemorySpace> cos_theta,
  Kokkos::View<Real*, DeviceMemorySpace> sin_theta)
{
  Kokkos::parallel_for(
    "cylindrical_to_cartesian_bind",
    Kokkos::RangePolicy<ExecutionSpace>(0, static_cast<LO>(in.extent(0))),
    KOKKOS_LAMBDA(const LO i) {
      const Real r = in(i, 0);
      const Real theta = in(i, 1);
      const Real c = Kokkos::cos(theta);
      const Real s = Kokkos::sin(theta);
      out(i, 0) = r * c;
      out(i, 1) = r * s;
      out(i, 2) = in(i, 2);
      cos_theta(i) = c;
      sin_theta(i) = s;
    });
}

void RotateCylindricalToCartesianValues(
  const Kokkos::View<Real*, DeviceMemorySpace>& c,
  const Kokkos::View<Real*, DeviceMemorySpace>& s,
  Rank2View<const Real, DeviceMemorySpace> in,
  Rank2View<Real, DeviceMemorySpace> out, LO n)
{
  Kokkos::parallel_for(
    "cylindrical_to_cartesian_values",
    Kokkos::RangePolicy<ExecutionSpace>(0, n), KOKKOS_LAMBDA(const LO i) {
      const Real vr = in(i, 0);
      const Real vt = in(i, 1);
      out(i, 0) = vr * c(i) - vt * s(i);
      out(i, 1) = vr * s(i) + vt * c(i);
      out(i, 2) = in(i, 2);
    });
}

void RotateCartesianToCylindricalValues(
  const Kokkos::View<Real*, DeviceMemorySpace>& c,
  const Kokkos::View<Real*, DeviceMemorySpace>& s,
  Rank2View<const Real, DeviceMemorySpace> in,
  Rank2View<Real, DeviceMemorySpace> out, LO n)
{
  Kokkos::parallel_for(
    "cartesian_to_cylindrical_values",
    Kokkos::RangePolicy<ExecutionSpace>(0, n), KOKKOS_LAMBDA(const LO i) {
      const Real vx = in(i, 0);
      const Real vy = in(i, 1);
      out(i, 0) = vx * c(i) + vy * s(i);
      out(i, 1) = -vx * s(i) + vy * c(i);
      out(i, 2) = in(i, 2);
    });
}

class BoundCylindricalRotation : public BoundCoordinateMap
{
public:
  BoundCylindricalRotation(std::shared_ptr<const CoordinateSystem> system,
                           Kokkos::View<Real**, DeviceMemorySpace> mapped,
                           Kokkos::View<Real*, DeviceMemorySpace> cos_theta,
                           Kokkos::View<Real*, DeviceMemorySpace> sin_theta,
                           const ValueBasis& stored, const ValueBasis& output)
    : BoundCoordinateMap(std::move(system), std::move(mapped)),
      cos_theta_(std::move(cos_theta)),
      sin_theta_(std::move(sin_theta)),
      stored_(stored),
      output_(output)
  {
  }

  [[nodiscard]] ValueBasis OutputBasis(const ValueBasis& stored) const override
  {
    if (!SameValueBasis(stored, stored_)) {
      throw pcms_error(
        "BoundCoordinateMap::OutputBasis: this map re-expresses values stored "
        "in '" +
        std::string(stored_.system->Kind()) +
        "' physical contravariant components; it cannot re-express the basis "
        "the field declares");
    }
    return output_;
  }

protected:
  Kokkos::View<Real*, DeviceMemorySpace> cos_theta_;
  Kokkos::View<Real*, DeviceMemorySpace> sin_theta_;

private:
  ValueBasis stored_;
  ValueBasis output_;
};

class BoundCartesianToCylindrical final : public BoundCylindricalRotation
{
public:
  using BoundCylindricalRotation::BoundCylindricalRotation;

  void TransformValues(ValueView<const Real, DeviceMemorySpace> in,
                       ValueView<Real, DeviceMemorySpace> out) const override
  {
    PCMS_FUNCTION_TIMER;
    ValidateTransformViews("CartesianToCylindrical::TransformValues", in, out);
    RotateCylindricalToCartesianValues(cos_theta_, sin_theta_, in.GetValues(),
                                       out.GetValues(), NumPoints());
  }
};

class BoundCylindricalToCartesian final : public BoundCylindricalRotation
{
public:
  using BoundCylindricalRotation::BoundCylindricalRotation;

  void TransformValues(ValueView<const Real, DeviceMemorySpace> in,
                       ValueView<Real, DeviceMemorySpace> out) const override
  {
    PCMS_FUNCTION_TIMER;
    ValidateTransformViews("CylindricalToCartesian::TransformValues", in, out);
    RotateCartesianToCylindricalValues(cos_theta_, sin_theta_, in.GetValues(),
                                       out.GetValues(), NumPoints());
  }
};

} // namespace

BoundCoordinateMap::BoundCoordinateMap(
  std::shared_ptr<const CoordinateSystem> system,
  Kokkos::View<Real**, DeviceMemorySpace> mapped, PointStatusView status)
  : system_(std::move(system)), mapped_(std::move(mapped)), status_(status)
{
}

CoordinateView<DeviceMemorySpace> BoundCoordinateMap::MappedPoints() const
{
  return CoordinateView<DeviceMemorySpace>(system_,
                                           MakeConstRank2View(mapped_));
}

ValueBasis BoundCoordinateMap::OutputBasis(const ValueBasis&) const
{
  throw pcms_error("BoundCoordinateMap::OutputBasis: this coordinate map maps "
                   "points only and cannot re-express component values");
}

void BoundCoordinateMap::TransformValues(
  ValueView<const Real, DeviceMemorySpace>,
  ValueView<Real, DeviceMemorySpace>) const
{
  throw pcms_error("BoundCoordinateMap::TransformValues: this coordinate map "
                   "maps points only and cannot re-express component values");
}

void BoundCoordinateMap::ValidateTransformViews(
  const char* who, const ValueView<const Real, DeviceMemorySpace>& in,
  const ValueView<Real, DeviceMemorySpace>& out) const
{
  detail::CheckWrittenValueBasis(who, out.GetBasis(),
                                 OutputBasis(in.GetBasis()));
  if (static_cast<LO>(in.extent(0)) != NumPoints() ||
      static_cast<LO>(out.extent(0)) != NumPoints()) {
    throw pcms_error(std::string(who) +
                     ": value counts do not match the bound point set");
  }
  if (in.extent(1) != out.extent(1)) {
    throw pcms_error(std::string(who) + ": input and output widths differ");
  }
}

void CoordinateMap::ValidateSourcePoints(
  const CoordinateView<DeviceMemorySpace>& points, const char* who) const
{
  if (!SameCoordinateSystem(points.GetCoordinateSystem(),
                            GetSourceCoordinateSystem())) {
    throw pcms_error(std::string(who) +
                     ": points are tagged with coordinate system '" +
                     std::string(points.GetCoordinateSystem()->Kind()) +
                     "' but this map's source coordinate system is '" +
                     std::string(GetSourceCoordinateSystem()->Kind()) +
                     "' (coordinate systems compare by object identity)");
  }
}

std::shared_ptr<const CoordinateSystem>
CartesianToCylindrical::GetSourceCoordinateSystem() const noexcept
{
  return csys::Cartesian::Create(3);
}

std::shared_ptr<const CoordinateSystem>
CartesianToCylindrical::GetTargetCoordinateSystem() const noexcept
{
  return csys::CylindricalRThetaZ::Create();
}

std::unique_ptr<BoundCoordinateMap> CartesianToCylindrical::Bind(
  const CoordinateView<DeviceMemorySpace>& query_points) const
{
  PCMS_FUNCTION_TIMER;
  ValidateSourcePoints(query_points, "CartesianToCylindrical");
  const auto in = query_points.GetValues();
  const auto n = in.extent(0);
  auto mapped = AllocateCoordinates(n, 3);
  auto cos_theta = AllocateTrig("bind_cos_theta", n);
  auto sin_theta = AllocateTrig("bind_sin_theta", n);
  BindCartesianToCylindrical(in, mapped, cos_theta, sin_theta);
  // The kernel reads the caller's point buffer; Bind's contract is that the
  // caller may release it as soon as Bind returns.
  ExecutionSpace().fence("CartesianToCylindrical::Bind: complete binding");
  return std::make_unique<BoundCartesianToCylindrical>(
    GetTargetCoordinateSystem(), std::move(mapped), std::move(cos_theta),
    std::move(sin_theta), CylindricalPhysicalVector(),
    CartesianPhysicalVector());
}

std::shared_ptr<const CoordinateSystem>
CylindricalToCartesian::GetSourceCoordinateSystem() const noexcept
{
  return csys::CylindricalRThetaZ::Create();
}

std::shared_ptr<const CoordinateSystem>
CylindricalToCartesian::GetTargetCoordinateSystem() const noexcept
{
  return csys::Cartesian::Create(3);
}

std::unique_ptr<BoundCoordinateMap> CylindricalToCartesian::Bind(
  const CoordinateView<DeviceMemorySpace>& query_points) const
{
  PCMS_FUNCTION_TIMER;
  ValidateSourcePoints(query_points, "CylindricalToCartesian");
  const auto in = query_points.GetValues();
  const auto n = in.extent(0);
  auto mapped = AllocateCoordinates(n, 3);
  auto cos_theta = AllocateTrig("bind_cos_theta", n);
  auto sin_theta = AllocateTrig("bind_sin_theta", n);
  BindCylindricalToCartesian(in, mapped, cos_theta, sin_theta);
  ExecutionSpace().fence("CylindricalToCartesian::Bind: complete binding");
  return std::make_unique<BoundCylindricalToCartesian>(
    GetTargetCoordinateSystem(), std::move(mapped), std::move(cos_theta),
    std::move(sin_theta), CartesianPhysicalVector(),
    CylindricalPhysicalVector());
}

} // namespace pcms
