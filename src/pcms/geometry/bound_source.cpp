#include "pcms/geometry/bound_source.hpp"
#include "pcms/geometry/geometry.hpp"
#include "pcms/field/coordinate_map.hpp"
#include "pcms/utility/arrays.h"
#include "pcms/utility/assert.h"
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
  if (detail::IsSmallHostBatch(found.extent(0))) {
    return static_cast<LO>(std::ranges::count_if(
      std::views::iota(LO{0}, n), [&](LO i) { return found(i); }));
  }
  LO num_found = 0;
  Kokkos::parallel_reduce(
    "BoundSource::CountFound",
    Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, n),
    KOKKOS_LAMBDA(const LO i, LO& count) { count += found(i) ? 1 : 0; },
    num_found);
  return num_found;
}

void SplitPoints(Kokkos::View<const bool*, DeviceMemorySpace> found,
                 Kokkos::View<LO*, DeviceMemorySpace> found_points,
                 Kokkos::View<LO*, DeviceMemorySpace> missing_points)
{
  if (detail::IsSmallHostBatch(found.extent(0))) {
    // the found and missing indices, each in increasing order
    std::ranges::partition_copy(
      std::views::iota(LO{0}, static_cast<LO>(found.extent(0))),
      found_points.data(), missing_points.data(),
      [&](LO i) { return found(i); });
    return;
  }
  Kokkos::parallel_scan(
    "BoundSource::SplitPoints",
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

/// Copies the query coordinates into the bound source's points, whatever the
/// layouts. (A free function: CUDA does not allow device lambdas in private
/// member functions.)
template <typename PointsView>
void CopyPoints(PointsView points,
                Rank2View<const Real, DeviceMemorySpace> coords)
{
  const auto dim = static_cast<int>(coords.extent(1));
  detail::ForEach(
    "BoundSource::CopyPoints", static_cast<LO>(coords.extent(0)),
    KOKKOS_LAMBDA(const LO i) {
      for (int d = 0; d < dim; ++d) {
        points(i, d) = coords(i, d);
      }
    });
}

} // namespace

BoundSource::BoundSource(std::shared_ptr<const Geometry> geometry,
                         const BindRequest& request)
  : geometry_(std::move(geometry))
{
  Reset(request);
}

void BoundSource::Reset(const BindRequest& request)
{
  policy_ = request.policy;
  map_ = request.map;
  coordinate_map_ = request.coordinate_map;
  system_ = request.query_coords.GetCoordinateSystem();
  const auto coords = request.coords.GetValues();
  num_points_ = static_cast<LO>(request.query_coords.GetValues().extent(0));
  PCMS_ALWAYS_ASSERT(static_cast<LO>(coords.extent(0)) == num_points_);
  if (map_ == nullptr) {
    detail::ReserveRows(points_, "bound_points", num_points_, coords.extent(1));
    auto points = Kokkos::subview(
      points_, Kokkos::make_pair(LO{0}, num_points_), Kokkos::ALL);
    CopyPoints(points, coords);
  }
  // every point is found until SetFound
  detail::ReserveRows(point_indices_, "bound_point_indices", num_points_);
  num_found_ = num_points_;
  const auto all =
    Kokkos::subview(point_indices_, Kokkos::make_pair(LO{0}, num_points_));
  if (detail::IsSmallHostBatch(num_points_)) {
    std::iota(all.data(), all.data() + num_points_, LO{0});
  } else {
    iota_view<LO>(all);
  }
}

Kokkos::View<const Real**, DeviceMemorySpace> BoundSource::LocatedPointValues()
  const noexcept
{
  if (map_ != nullptr) {
    return map_->Mapped();
  }
  return Kokkos::subview(points_, Kokkos::make_pair(LO{0}, num_points_),
                         Kokkos::ALL);
}

void BoundSource::Rebind(RebindKey, const BindRequest& request)
{
  Reset(request);
  RebindImpl();
}

bool BoundSource::Unshared() const noexcept
{
  return points_.use_count() <= 1 && point_indices_.use_count() <= 1 &&
         UnsharedImpl();
}

void BoundSource::SetFound(Kokkos::View<const bool*, DeviceMemorySpace> found)
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

} // namespace pcms
