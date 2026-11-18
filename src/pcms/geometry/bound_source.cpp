#include "pcms/geometry/bound_source.hpp"
#include "pcms/geometry/geometry.hpp"
#include "pcms/field/coordinate_map.hpp"
#include "pcms/utility/arrays.h"
#include "pcms/utility/assert.h"
#include <string>

namespace pcms
{
namespace
{

Kokkos::View<Real**, DeviceMemorySpace> CopyPoints(
  Rank2View<const Real, DeviceMemorySpace> coords)
{
  Kokkos::View<Real**, DeviceMemorySpace> points(
    Kokkos::view_alloc(Kokkos::WithoutInitializing, "bound_points"),
    coords.extent(0), coords.extent(1));
  ConvertMismatchLayoutView2D(points, coords);
  return points;
}

Kokkos::View<LO*, DeviceMemorySpace> AllPoints(LO n)
{
  Kokkos::View<LO*, DeviceMemorySpace> all(
    Kokkos::view_alloc(Kokkos::WithoutInitializing, "bound_found_points"), n);
  iota_view<LO>(all);
  return all;
}

LO CountFound(Kokkos::View<const bool*, DeviceMemorySpace> found)
{
  LO num_found = 0;
  Kokkos::parallel_reduce(
    "BoundSource::CountFound",
    Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(
      0, static_cast<LO>(found.extent(0))),
    KOKKOS_LAMBDA(const LO i, LO& count) { count += found(i) ? 1 : 0; },
    num_found);
  return num_found;
}

void SplitPoints(Kokkos::View<const bool*, DeviceMemorySpace> found,
                 Kokkos::View<LO*, DeviceMemorySpace> found_points,
                 Kokkos::View<LO*, DeviceMemorySpace> missing_points)
{
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

} // namespace

BoundSource::BoundSource(std::shared_ptr<const Geometry> geometry,
                         const BindRequest& request)
  : geometry_(std::move(geometry)),
    policy_(request.policy),
    map_(request.map),
    coordinate_map_(request.coordinate_map),
    system_(request.query_coords.GetCoordinateSystem()),
    num_points_(static_cast<LO>(request.query_coords.GetValues().extent(0))),
    located_(map_ != nullptr ? map_->Mapped()
                             : Kokkos::View<const Real**, DeviceMemorySpace>(
                                 CopyPoints(request.coords.GetValues()))),
    found_(AllPoints(num_points_))
{
  PCMS_ALWAYS_ASSERT(static_cast<LO>(request.coords.GetValues().extent(0)) ==
                     num_points_);
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
  Kokkos::View<LO*, DeviceMemorySpace> found_points(
    Kokkos::view_alloc(Kokkos::WithoutInitializing, "bound_found_points"),
    num_found);
  Kokkos::View<LO*, DeviceMemorySpace> missing_points(
    Kokkos::view_alloc(Kokkos::WithoutInitializing, "bound_missing_points"),
    num_missing);
  SplitPoints(found, found_points, missing_points);
  found_ = found_points;
  missing_ = missing_points;
}

} // namespace pcms
