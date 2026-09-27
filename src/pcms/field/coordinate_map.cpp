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

/// Direction the query points travel through a cylindrical map.
enum class PointDirection
{
  /// (x, y, z) -> (r, theta, z), with theta from atan2(y, x) in (-pi, pi]
  ToCylindrical,
  /// (r, theta, z) -> (x, y, z)
  ToCartesian
};

// Kernels as free functions: CUDA extended lambdas may not be defined inside
// member functions of classes with private members.

/// @brief Maps query points between Cartesian and cylindrical coordinates and
///        caches the angle each point sits at.
///
/// @param in query points, one per row, in the map's source coordinate system:
///        (x, y, z) for ToCylindrical, (r, theta, z) for ToCartesian
/// @param out receives the mapped points, one per row of `in`, in the map's
///        target coordinate system: (r, theta, z) for ToCylindrical,
///        (x, y, z) for ToCartesian. Must be preallocated with the extents of
///        `in`
/// @param cos_theta receives cos(theta) at each point; one entry per row of
///        `in`
/// @param sin_theta receives sin(theta) at each point; one entry per row of
///        `in`
/// @param direction which way the points are mapped
///
/// @details This is the per-point work of Bind for both cylindrical maps. It
/// produces the mapped coordinates the source is evaluated at and, in the same
/// pass, the cos/sin of each point's angle. The angle is what the bound map
/// later needs to rotate vector components between the (e_r, e_theta, e_z) and
/// (e_x, e_y, e_z) bases, so caching it here lets every later TransformValues
/// call, one per evaluated field, reuse it instead of recomputing trig. z
/// passes through unchanged. The kernel is launched asynchronously and reads
/// `in`; CoordinateMap::Bind fences before the caller may release it.
void BindCylindricalPoints(Rank2View<const Real, DeviceMemorySpace> in,
                           Kokkos::View<Real**, DeviceMemorySpace> out,
                           Kokkos::View<Real*, DeviceMemorySpace> cos_theta,
                           Kokkos::View<Real*, DeviceMemorySpace> sin_theta,
                           PointDirection direction)
{
  Kokkos::parallel_for(
    "cylindrical_map_bind",
    Kokkos::RangePolicy<ExecutionSpace>(0, static_cast<LO>(in.extent(0))),
    KOKKOS_LAMBDA(const LO i) {
      Real c, s;
      if (direction == PointDirection::ToCylindrical) {
        const Real x = in(i, 0);
        const Real y = in(i, 1);
        const Real theta = Kokkos::atan2(y, x);
        c = Kokkos::cos(theta);
        s = Kokkos::sin(theta);
        out(i, 0) = Kokkos::sqrt(x * x + y * y);
        out(i, 1) = theta;
      } else {
        const Real r = in(i, 0);
        c = Kokkos::cos(in(i, 1));
        s = Kokkos::sin(in(i, 1));
        out(i, 0) = r * c;
        out(i, 1) = r * s;
      }
      out(i, 2) = in(i, 2);
      cos_theta(i) = c;
      sin_theta(i) = s;
    });
}

void RotateValues(const Kokkos::View<Real*, DeviceMemorySpace>& c,
                  const Kokkos::View<Real*, DeviceMemorySpace>& s, Real sign,
                  Rank2View<Real, DeviceMemorySpace> values, LO n)
{
  Kokkos::parallel_for(
    "cylindrical_map_rotate_values", Kokkos::RangePolicy<ExecutionSpace>(0, n),
    KOKKOS_LAMBDA(const LO i) {
      const Real v0 = values(i, 0);
      const Real v1 = values(i, 1);
      const Real si = sign * s(i);
      values(i, 0) = v0 * c(i) - v1 * si;
      values(i, 1) = v0 * si + v1 * c(i);
    });
}

/// Rotates physical vector components about z by +theta (cylindrical to
/// Cartesian) or -theta (Cartesian to cylindrical) at each bound point.
class BoundCylindricalRotation final : public BoundCoordinateMap
{
public:
  BoundCylindricalRotation(std::shared_ptr<const CoordinateSystem> system,
                           Kokkos::View<Real**, DeviceMemorySpace> mapped,
                           Kokkos::View<Real*, DeviceMemorySpace> cos_theta,
                           Kokkos::View<Real*, DeviceMemorySpace> sin_theta,
                           const ValueBasis& stored, const ValueBasis& output,
                           Real sign)
    : BoundCoordinateMap(std::move(system), std::move(mapped)),
      cos_theta_(std::move(cos_theta)),
      sin_theta_(std::move(sin_theta)),
      stored_(stored),
      output_(output),
      sign_(sign)
  {
  }

  /// Values already in the output basis pass through unchanged: components
  /// with respect to the output basis do not depend on where the point is
  /// expressed, so only the points needed mapping.
  [[nodiscard]] ValueBasis OutputBasis(const ValueBasis& stored) const override
  {
    if (!SameValueBasis(stored, stored_) && !SameValueBasis(stored, output_)) {
      throw pcms_error(
        "BoundCoordinateMap::OutputBasis: this map re-expresses values stored "
        "in '" +
        std::string(stored_.system->Kind()) + "' or '" +
        std::string(output_.system->Kind()) +
        "' physical contravariant components; it cannot re-express the basis "
        "the field declares");
    }
    return output_;
  }

  ValueView<Real, DeviceMemorySpace> TransformValues(
    ValueView<Real, DeviceMemorySpace> values) const override
  {
    PCMS_FUNCTION_TIMER;
    const ValueBasis output = ValidateTransformValues(
      "BoundCylindricalRotation::TransformValues", values);
    if (!SameValueBasis(values.GetBasis(), output_)) {
      RotateValues(cos_theta_, sin_theta_, sign_, values.GetValues(),
                   NumPoints());
    }
    return ValueView<Real, DeviceMemorySpace>(output, values.GetValues());
  }

private:
  Kokkos::View<Real*, DeviceMemorySpace> cos_theta_;
  Kokkos::View<Real*, DeviceMemorySpace> sin_theta_;
  ValueBasis stored_;
  ValueBasis output_;
  Real sign_;
};

/// Values run opposite to points: a map into cylindrical coordinates rotates
/// cylindrical components into Cartesian ones, and vice versa.
std::unique_ptr<BoundCoordinateMap> BindCylindrical(
  const CoordinateView<DeviceMemorySpace>& query_points,
  std::shared_ptr<const CoordinateSystem> target, PointDirection direction)
{
  const auto in = query_points.GetValues();
  const auto n = in.extent(0);
  auto mapped = AllocateCoordinates(n, 3);
  auto cos_theta = AllocateTrig("bind_cos_theta", n);
  auto sin_theta = AllocateTrig("bind_sin_theta", n);
  BindCylindricalPoints(in, mapped, cos_theta, sin_theta, direction);
  const bool to_cylindrical = direction == PointDirection::ToCylindrical;
  const auto& stored =
    to_cylindrical ? CylindricalPhysicalVector() : CartesianPhysicalVector();
  const auto& output =
    to_cylindrical ? CartesianPhysicalVector() : CylindricalPhysicalVector();
  return std::make_unique<BoundCylindricalRotation>(
    std::move(target), std::move(mapped), std::move(cos_theta),
    std::move(sin_theta), stored, output, to_cylindrical ? 1.0 : -1.0);
}

} // namespace

BoundCoordinateMap::BoundCoordinateMap(
  std::shared_ptr<const CoordinateSystem> system,
  Kokkos::View<Real**, DeviceMemorySpace> mapped)
  : system_(std::move(system)), mapped_(std::move(mapped))
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

ValueView<Real, DeviceMemorySpace> BoundCoordinateMap::TransformValues(
  ValueView<Real, DeviceMemorySpace>) const
{
  throw pcms_error("BoundCoordinateMap::TransformValues: this coordinate map "
                   "maps points only and cannot re-express component values");
}

ValueBasis BoundCoordinateMap::ValidateTransformValues(
  const char* who, const ValueView<Real, DeviceMemorySpace>& values) const
{
  ValueBasis output = OutputBasis(values.GetBasis());
  if (static_cast<LO>(values.extent(0)) != NumPoints()) {
    throw pcms_error(std::string(who) +
                     ": value count does not match the bound point set");
  }
  return output;
}

CoordinateMap::CoordinateMap(std::shared_ptr<const CoordinateSystem> source,
                             std::shared_ptr<const CoordinateSystem> target)
  : source_(std::move(source)), target_(std::move(target))
{
  if (source_ == nullptr || target_ == nullptr) {
    throw pcms_error("CoordinateMap: coordinate systems must not be null");
  }
}

std::unique_ptr<BoundCoordinateMap> CoordinateMap::Bind(
  const CoordinateView<DeviceMemorySpace>& query_points) const
{
  PCMS_FUNCTION_TIMER;
  if (!SameCoordinateSystem(query_points.GetCoordinateSystem(), source_)) {
    throw pcms_error("CoordinateMap::Bind: points are tagged with coordinate "
                     "system '" +
                     std::string(query_points.GetCoordinateSystem()->Kind()) +
                     "' but this map's source coordinate system is '" +
                     std::string(source_->Kind()) + "'");
  }
  auto bound = BindImpl(query_points);
  // BindImpl kernels read the caller's point buffer, which the caller may
  // release as soon as Bind returns.
  ExecutionSpace().fence("CoordinateMap::Bind: complete binding");
  if (bound == nullptr ||
      bound->NumPoints() !=
        static_cast<LO>(query_points.GetValues().extent(0)) ||
      !SameCoordinateSystem(bound->MappedPoints().GetCoordinateSystem(),
                            target_)) {
    throw pcms_error("CoordinateMap::Bind: BindImpl must return one mapped "
                     "point per query point, in the target coordinate system");
  }
  return bound;
}

CartesianToCylindrical::CartesianToCylindrical()
  : CoordinateMap(csys::Cartesian::Create(3),
                  csys::CylindricalRThetaZ::Create())
{
}

std::unique_ptr<BoundCoordinateMap> CartesianToCylindrical::BindImpl(
  const CoordinateView<DeviceMemorySpace>& query_points) const
{
  return BindCylindrical(query_points, GetTargetCoordinateSystem(),
                         PointDirection::ToCylindrical);
}

CylindricalToCartesian::CylindricalToCartesian()
  : CoordinateMap(csys::CylindricalRThetaZ::Create(),
                  csys::Cartesian::Create(3))
{
}

std::unique_ptr<BoundCoordinateMap> CylindricalToCartesian::BindImpl(
  const CoordinateView<DeviceMemorySpace>& query_points) const
{
  return BindCylindrical(query_points, GetTargetCoordinateSystem(),
                         PointDirection::ToCartesian);
}

} // namespace pcms
