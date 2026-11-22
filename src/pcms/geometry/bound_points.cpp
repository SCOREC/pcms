#include "pcms/geometry/bound_points.hpp"
#include "pcms/geometry/geometry.hpp"
#include "pcms/field/coordinate_map.hpp"
#include "pcms/utility/arrays.hpp"
#include "pcms/utility/assert.hpp"
#include <algorithm>
#include <numeric>
#include <ranges>
#include <string>

namespace pcms
{
namespace
{

LO CountFound(Kokkos::View<const bool*, DeviceMemorySpace> found)
{
  const auto n = static_cast<LO>(found.extent(0));
  if (pcms::detail::IsSmallHostBatch(found.extent(0))) {
    return static_cast<LO>(std::ranges::count_if(
      std::views::iota(LO{0}, n), [&](LO i) { return found(i); }));
  }
  LO num_found = 0;
  Kokkos::parallel_reduce(
    "BoundPoints::CountFound",
    Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, n),
    KOKKOS_LAMBDA(const LO i, LO& count) { count += found(i) ? 1 : 0; },
    num_found);
  return num_found;
}

void SplitPoints(Kokkos::View<const bool*, DeviceMemorySpace> found,
                 Kokkos::View<LO*, DeviceMemorySpace> found_points,
                 Kokkos::View<LO*, DeviceMemorySpace> missing_points)
{
  if (pcms::detail::IsSmallHostBatch(found.extent(0))) {
    // the found and missing indices, each in increasing order
    std::ranges::partition_copy(
      std::views::iota(LO{0}, static_cast<LO>(found.extent(0))),
      found_points.data(), missing_points.data(),
      [&](LO i) { return found(i); });
    return;
  }
  Kokkos::parallel_scan(
    "BoundPoints::SplitPoints",
    Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(
      0, static_cast<LO>(found.extent(0))),
    KOKKOS_LAMBDA(const LO i, LO& num_found_before, const bool final) {
      if (final) {
        if (found(i)) {
          found_points(num_found_before) = i;
        } else {
          missing_points(i - num_found_before) = i;
        }
      }
      num_found_before += found(i) ? 1 : 0;
    });
}

/// Copies the query coordinates into the BoundPoints's points, whatever the
/// layouts. (A free function: CUDA does not allow device lambdas in private
/// member functions.)
template <typename PointsView>
void CopyPoints(PointsView points,
                Rank2View<const Real, DeviceMemorySpace> coords)
{
  const auto dim = static_cast<int>(coords.extent(1));
  pcms::detail::ForEach(
    "BoundPoints::CopyPoints", static_cast<LO>(coords.extent(0)),
    KOKKOS_LAMBDA(const LO i) {
      for (int d = 0; d < dim; ++d) {
        points(i, d) = coords(i, d);
      }
    });
}

} // namespace

BoundPoints::BoundPoints(std::shared_ptr<const Geometry> geometry,
                         const BindRequest& request)
  : geometry_(std::move(geometry))
{
  Reset(request);
}

void BoundPoints::Reset(const BindRequest& request)
{
  policy_ = request.policy;
  factors_.reset();
  map_ = request.map;
  coordinate_map_ = request.coordinate_map;
  system_ = request.query_coords.GetCoordinateSystem();
  const auto coords = request.coords.GetValues();
  num_points_ = static_cast<LO>(request.query_coords.GetValues().extent(0));
  PCMS_ALWAYS_ASSERT(static_cast<LO>(coords.extent(0)) == num_points_);
  detail::ReserveRows(points_, "bound_points", num_points_, coords.extent(1));
  auto points = Kokkos::subview(points_, Kokkos::make_pair(LO{0}, num_points_),
                                Kokkos::ALL);
  CopyPoints(points, coords);
  // every point is found until SetFound
  detail::ReserveRows(point_indices_, "bound_point_indices", num_points_);
  num_found_ = num_points_;
  const auto all =
    Kokkos::subview(point_indices_, Kokkos::make_pair(LO{0}, num_points_));
  if (pcms::detail::IsSmallHostBatch(num_points_)) {
    std::iota(all.data(), all.data() + num_points_, LO{0});
  } else {
    pcms::detail::iota_view<LO>(all);
  }
}

CoordinateView<DeviceMemorySpace> BoundPoints::LocatedPoints() const
{
  return CoordinateView<DeviceMemorySpace>(
    map_ != nullptr ? map_->MappedPoints().GetCoordinateSystem() : system_,
    MakeConstRank2View(LocatedPointValues()));
}

Kokkos::View<const Real**, DeviceMemorySpace> BoundPoints::LocatedPointValues()
  const noexcept
{
  return Kokkos::subview(points_, Kokkos::make_pair(LO{0}, num_points_),
                         Kokkos::ALL);
}

void BoundPoints::ReplaceLocatedPoints(
  Kokkos::View<const Real**, DeviceMemorySpace> points)
{
  PCMS_ALWAYS_ASSERT(static_cast<LO>(points.extent(0)) == num_points_);
  CopyPoints(Kokkos::subview(points_, Kokkos::make_pair(LO{0}, num_points_),
                             Kokkos::ALL),
             MakeConstRank2View(points));
}

void BoundPoints::Rebind(RebindKey, const BindRequest& request)
{
  Reset(request);
  RebindImpl();
}

bool BoundPoints::Unshared() const noexcept
{
  return points_.use_count() <= 1 && point_indices_.use_count() <= 1 &&
         UnsharedImpl();
}

void BoundPoints::SetFound(Kokkos::View<const bool*, DeviceMemorySpace> found)
{
  PCMS_ALWAYS_ASSERT(static_cast<LO>(found.extent(0)) == NumPoints());
  const LO num_found = CountFound(found);
  const LO num_missing = NumPoints() - num_found;
  if (num_missing > 0 && policy_.mode == OutOfBoundsMode::ERROR) {
    throw pcms_error("Bind: " + std::to_string(num_missing) +
                     " query point(s) are outside the source's domain");
  }
  // the found points first, then the missing ones, in the storage Reset sized
  SplitPoints(
    found, Kokkos::subview(point_indices_, Kokkos::make_pair(0, num_found)),
    Kokkos::subview(point_indices_, Kokkos::make_pair(num_found, num_points_)));
  num_found_ = num_found;
}

detail::GeometricFactors BoundPoints::ComputeFactors() const
{
  return detail::SystemFactors(*geometry_->GetCoordinateSystem(),
                               LocatedPointValues());
}

const detail::GeometricFactors& BoundPoints::Factors() const
{
  if (!factors_) {
    factors_.emplace(ComputeFactors());
  }
  return *factors_;
}

View<3, const Real, DeviceMemorySpace, Kokkos::layout_stride>
BoundPoints::Jacobian() const
{
  return Factors().Jacobian();
}

BoundCoordinateSystem::MetricView BoundPoints::Metric() const
{
  return Factors().Metric();
}

BoundCoordinateSystem::FactorView BoundPoints::VolumeFactor() const
{
  return Factors().VolumeFactor();
}

detail::GeometricFactors detail::SystemFactors(
  const CoordinateSystem& system,
  Kokkos::View<const Real**, DeviceMemorySpace> points)
{
  if (!system.HasMetric()) {
    throw pcms_error("BoundPoints: the geometric factors need a metric; the "
                     "coordinate system '" +
                     std::string(system.Kind()) + "' has none");
  }
  const int dim = system.Dimension();
  Kokkos::View<Real***, DeviceMemorySpace> identity("identity_jacobian", 1, dim,
                                                    dim);
  Kokkos::parallel_for(
    "identity_jacobian",
    Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, dim),
    KOKKOS_LAMBDA(const int i) { identity(0, i, i) = 1.0; });
  return GeometricFactors(static_cast<LO>(points.extent(0)), identity,
                          system.Bind(MakeConstRank2View(points)));
}

} // namespace pcms
