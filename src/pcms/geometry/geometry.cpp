#include "pcms/geometry/geometry.hpp"
#include "pcms/field/coordinate_map.hpp"
#include "pcms/utility/assert.h"
#include <string>

namespace pcms
{

std::shared_ptr<const BoundSource> Geometry::Bind(
  CoordinateView<DeviceMemorySpace> coords, OutOfBoundsPolicy policy,
  std::shared_ptr<const CoordinateMap> map) const
{
  return Bind(std::move(coords), OwnerHints{}, policy, std::move(map));
}

std::shared_ptr<const BoundSource> Geometry::Bind(
  CoordinateView<DeviceMemorySpace> coords, const OwnerHints& owners,
  OutOfBoundsPolicy policy, std::shared_ptr<const CoordinateMap> map) const
{
  const auto& query_system = coords.GetCoordinateSystem();
  if (map != nullptr) {
    if (!SameCoordinateSystem(query_system, map->GetSourceCoordinateSystem()) ||
        !SameCoordinateSystem(map->GetTargetCoordinateSystem(),
                              GetCoordinateSystem())) {
      throw pcms_error(
        "Geometry::Bind: the query points are in '" +
        std::string(query_system->Kind()) + "' and the map takes '" +
        std::string(map->GetSourceCoordinateSystem()->Kind()) + "' to '" +
        std::string(map->GetTargetCoordinateSystem()->Kind()) +
        "'; it must take the query coordinate system to the geometry's '" +
        std::string(GetCoordinateSystem()->Kind()) + "'");
    }
  } else if (!SameCoordinateSystem(query_system, GetCoordinateSystem())) {
    throw pcms_error("Geometry::Bind: query coordinate system '" +
                     std::string(query_system->Kind()) +
                     "' is not the geometry's coordinate system '" +
                     std::string(GetCoordinateSystem()->Kind()) +
                     "'; pass a CoordinateMap from one to the other as "
                     "Bind's map argument");
  }
  if (policy.mode == OutOfBoundsMode::NEAREST_BOUNDARY &&
      !SupportsNearestBoundary()) {
    throw pcms_error("Geometry::Bind: this geometry does not support "
                     "OutOfBoundsMode::NEAREST_BOUNDARY");
  }
  // a hint naming another geometry, or not one entity per point, is not one;
  // nor is a hint given with a map, as it describes the unmapped points
  const OwnerHints hint =
    map == nullptr &&
        TakesHints(owners, static_cast<LO>(coords.GetValues().extent(0)))
      ? owners
      : OwnerHints{};
  std::shared_ptr<BoundSource> bound;
  if (map != nullptr) {
    std::shared_ptr<const BoundCoordinateMap> bound_map = map->Bind(coords);
    const CoordinateView<DeviceMemorySpace> mapped = bound_map->MappedPoints();
    bound = BindImpl(
      BindRequest{mapped, policy, hint, std::move(bound_map), coords, map});
  } else {
    bound = BindImpl(BindRequest{coords, policy, hint, nullptr, coords});
  }
  PCMS_ALWAYS_ASSERT(bound != nullptr);
  return bound;
}

SupportResults Geometry::Neighbourhood(CoordinateView<DeviceMemorySpace>, int,
                                       const NeighbourhoodOptions&,
                                       const OwnerHints&) const
{
  throw pcms_error("Geometry::Neighbourhood: this geometry has no "
                   "neighbourhood search");
}

DifferentialView Geometry::Differential(const BoundSource&) const
{
  return {};
}

void detail::CheckNeighbourhoodArguments(
  const Geometry& geometry, const CoordinateView<DeviceMemorySpace>& coords,
  const NeighbourhoodOptions& options, const std::string& caller)
{
  if (!SameCoordinateSystem(coords.GetCoordinateSystem(),
                            geometry.GetCoordinateSystem())) {
    throw pcms_error(caller + ": query coordinate system '" +
                     std::string(coords.GetCoordinateSystem()->Kind()) +
                     "' is not the geometry's coordinate system '" +
                     std::string(geometry.GetCoordinateSystem()->Kind()) + "'");
  }
  // the searches measure Euclidean distances
  if (!HasIdentityMetric(*geometry.GetCoordinateSystem())) {
    throw pcms_error(
      caller +
      ": requires a coordinate system whose metric is the "
      "identity (orthogonal basis and unit scale factors); got '" +
      std::string(geometry.GetCoordinateSystem()->Kind()) + "'");
  }
  if (options.max_supports != 0 &&
      options.max_supports != 3 * options.min_supports) {
    throw pcms_error(caller + ": max_supports must be 0 or 3 * min_supports (" +
                     std::to_string(3 * options.min_supports) + "), got " +
                     std::to_string(options.max_supports) +
                     "; other bounds are not honoured yet");
  }
}

void detail::CheckEntityIds(Kokkos::View<const LO*, DeviceMemorySpace> ids,
                            LO num_entities, const std::string& caller,
                            const std::string& entity)
{
  LO min_id = 0;
  LO max_id = 0;
  Kokkos::parallel_reduce(
    "CheckEntityIds",
    Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(
      0, static_cast<LO>(ids.extent(0))),
    KOKKOS_LAMBDA(const LO i, LO& lo, LO& hi) {
      lo = ids(i) < lo ? ids(i) : lo;
      hi = ids(i) > hi ? ids(i) : hi;
    },
    Kokkos::Min<LO>(min_id), Kokkos::Max<LO>(max_id));
  if (ids.extent(0) > 0 && (min_id < 0 || max_id >= num_entities)) {
    throw pcms_error(caller + ": hinted " + entity + " ids span [" +
                     std::to_string(min_id) + ", " + std::to_string(max_id) +
                     "], outside [0, " + std::to_string(num_entities) + ")");
  }
}

} // namespace pcms
