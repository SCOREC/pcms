#include "pcms/field/coordinate_map.hpp"
#include "pcms/utility/coordinate_systems/cartesian.hpp"
#include "pcms/utility/coordinate_systems/cylindrical.hpp"
#include "pcms/utility/assert.hpp"
#include "pcms/utility/detail/profile.hpp"
#include <Kokkos_Core.hpp>
#include <string>

namespace pcms
{

namespace
{

using ExecutionSpace = DeviceMemorySpace::execution_space;

// The kernels of the maps between coordinate systems: the image of a point
// and the Jacobian jacobian[3 * j + i] = d out_j / d x_i there. On the axis
// R = 0 the angle is the one atan2 gives and the Jacobian is written with
// R = 1 in its angular entries, the convention of a zero scale factor in
// ValueTransform: physical components then rotate by that angle.

// (x, y, z) -> (R, theta, Z), theta = atan2(y, x) in (-pi, pi]
struct CartesianToCylindricalKernel
{
  KOKKOS_FUNCTION void operator()(const Kokkos::Array<Real, 3>& x,
                                  Kokkos::Array<Real, 3>& out,
                                  Kokkos::Array<Real, 9>& jacobian) const
  {
    const Real r = Kokkos::sqrt(x[0] * x[0] + x[1] * x[1]);
    const Real theta = Kokkos::atan2(x[1], x[0]);
    const Real c = Kokkos::cos(theta);
    const Real s = Kokkos::sin(theta);
    const Real r_safe = r > 0 ? r : Real(1);
    out[0] = r;
    out[1] = theta;
    out[2] = x[2];
    jacobian = {c, s, 0, -s / r_safe, c / r_safe, 0, 0, 0, 1};
  }
};

// (R, theta, Z) -> (x, y, z)
struct CylindricalToCartesianKernel
{
  KOKKOS_FUNCTION void operator()(const Kokkos::Array<Real, 3>& x,
                                  Kokkos::Array<Real, 3>& out,
                                  Kokkos::Array<Real, 9>& jacobian) const
  {
    const Real c = Kokkos::cos(x[1]);
    const Real s = Kokkos::sin(x[1]);
    const Real r_safe = x[0] != 0 ? x[0] : Real(1);
    out[0] = x[0] * c;
    out[1] = x[0] * s;
    out[2] = x[2];
    jacobian = {c, -r_safe * s, 0, s, r_safe * c, 0, 0, 0, 1};
  }
};

// S coordinates -> S + 1, coordinate `axis` of the image fixed at `value`
template <int S>
struct InsertionKernel
{
  int axis;
  Real value;

  KOKKOS_FUNCTION void operator()(
    const Kokkos::Array<Real, S>& x, Kokkos::Array<Real, S + 1>& out,
    Kokkos::Array<Real, static_cast<std::size_t>(S + 1) * S>& jacobian) const
  {
    for (int j = 0, i = 0; j <= S; ++j) {
      out[j] = j == axis ? value : x[i];
      for (int k = 0; k < S; ++k) {
        jacobian[S * j + k] = j != axis && k == i ? 1 : 0;
      }
      i += j == axis ? 0 : 1;
    }
  }
};

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

const ValueTransform& BoundCoordinateMap::ValueTransform() const
{
  if (!value_transform_) {
    // values flow against the map: from the target system, where the points
    // were located, to the source system of the query points; the scale
    // factors of both come from the target's metric at the mapped points
    const auto metric = system_->HasMetric()
                          ? system_->Bind(MakeConstRank2View(mapped_))
                          : nullptr;
    value_transform_.emplace(detail::Pullback(
      Jacobian(), system_,
      metric != nullptr ? metric->Metric() : detail::MetricView{}, source_));
  }
  return *value_transform_;
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
  if (!pcms::detail::SameCoordinateSystem(query_points.GetCoordinateSystem(),
                                          source_)) {
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
      !pcms::detail::SameCoordinateSystem(
        bound->MappedPoints().GetCoordinateSystem(), target_)) {
    throw pcms_error("CoordinateMap::Bind: BindImpl must return one mapped "
                     "point per query point, in the target coordinate system");
  }
  bound->source_ = source_;
  return bound;
}

std::shared_ptr<const CoordinateMap> CreateCartesianToCylindrical()
{
  return CreateCoordinateMap(csys::CreateCartesian(3),
                             csys::CreateCylindricalRThetaZ(),
                             CartesianToCylindricalKernel{});
}

std::shared_ptr<const CoordinateMap> CreateCylindricalToCartesian()
{
  return CreateCoordinateMap(csys::CreateCylindricalRThetaZ(),
                             csys::CreateCartesian(3),
                             CylindricalToCartesianKernel{});
}

std::shared_ptr<const CoordinateMap> CreateCoordinateEmbedding(
  std::shared_ptr<const CoordinateSystem> source,
  std::shared_ptr<const CoordinateSystem> target, int axis, Real value)
{
  if (source == nullptr || target == nullptr) {
    throw pcms_error(
      "CreateCoordinateEmbedding: coordinate systems must not be null");
  }
  const int dim = source->Dimension();
  if (target->Dimension() != dim + 1 || dim > 2 || axis < 0 || axis > dim) {
    throw pcms_error("CreateCoordinateEmbedding: the target must have one "
                     "coordinate more than the source (at most 3), and the "
                     "axis must be one of its coordinates");
  }
  if (dim == 1) {
    return CreateCoordinateMap(std::move(source), std::move(target),
                               InsertionKernel<1>{axis, value});
  }
  return CreateCoordinateMap(std::move(source), std::move(target),
                             InsertionKernel<2>{axis, value});
}

} // namespace pcms
