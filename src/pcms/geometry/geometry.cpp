#include "pcms/geometry/geometry.hpp"
#include "pcms/geometry/detail/geometry.hpp"
#include "pcms/field/coordinate_map.hpp"
#include "pcms/utility/assert.hpp"
#include "pcms/utility/detail/memory_spaces.hpp"
#include "pcms/utility/detail/periodic.hpp"
#include <cmath>
#include <numbers>
#include <string>
#include <vector>

namespace pcms
{
namespace
{

/// The checks Bind makes of the query coordinate system and the policy.
void CheckBindArguments(const Geometry& geometry,
                        const CoordinateView<DeviceMemorySpace>& coords,
                        const OutOfBoundsPolicy& policy,
                        const CoordinateMap* map)
{
  const auto& query_system = coords.GetCoordinateSystem();
  const auto& system = geometry.GetCoordinateSystem();
  if (map != nullptr) {
    if (!pcms::detail::SameCoordinateSystem(query_system,
                                            map->GetSourceCoordinateSystem()) ||
        !pcms::detail::SameCoordinateSystem(map->GetTargetCoordinateSystem(),
                                            system)) {
      throw pcms_error(
        "Geometry::Bind: the query points are in '" +
        std::string(query_system->Kind()) + "' and the map takes '" +
        std::string(map->GetSourceCoordinateSystem()->Kind()) + "' to '" +
        std::string(map->GetTargetCoordinateSystem()->Kind()) +
        "'; it must take the query coordinate system to the geometry's '" +
        std::string(system->Kind()) + "'");
    }
  } else if (!pcms::detail::SameCoordinateSystem(query_system, system) &&
             !pcms::detail::SameCoordinateSystem(
               query_system, geometry.GetParameterSystem())) {
    throw pcms_error("Geometry::Bind: query coordinate system '" +
                     std::string(query_system->Kind()) +
                     "' is not the geometry's coordinate system '" +
                     std::string(system->Kind()) +
                     "'; pass a CoordinateMap from one to the other as "
                     "Bind's map argument, or points in the geometry's "
                     "parameter system");
  }
  if (policy.mode == OutOfBoundsMode::NEAREST_BOUNDARY &&
      !geometry.SupportsNearestBoundary()) {
    throw pcms_error("Geometry::Bind: this geometry does not support "
                     "OutOfBoundsMode::NEAREST_BOUNDARY");
  }
}

/// Moves coordinate `axis` into [lower, lower + length) and counts the lengths
/// each point moved down.
void ReduceAxis(Kokkos::View<Real**, DeviceMemorySpace> points,
                Kokkos::View<Real*, DeviceMemorySpace> wraps, int axis,
                Real lower, Real length)
{
  pcms::detail::ForEach(
    "ReduceAxis", static_cast<LO>(points.extent(0)), KOKKOS_LAMBDA(const LO p) {
      Real n = 0;
      points(p, axis) =
        pcms::detail::WrapToPeriod(points(p, axis), lower, length, n);
      wraps(p) = n;
    });
}

} // namespace

Kokkos::View<Real**, DeviceMemorySpace> detail::ReduceToPeriod(
  const CoordinateView<DeviceMemorySpace>& coords)
{
  const auto identifications = coords.GetCoordinateSystem()->Identifications();
  if (identifications.empty()) {
    return {};
  }
  const auto in = coords.GetValues();
  const auto n = static_cast<LO>(in.extent(0));
  const auto dim = static_cast<int>(in.extent(1));
  Kokkos::View<Real**, DeviceMemorySpace> points(
    Kokkos::view_alloc(Kokkos::WithoutInitializing, "reduced_points"), n, dim);
  pcms::detail::ForEach(
    "CopyReduced", n, KOKKOS_LAMBDA(const LO p) {
      for (int d = 0; d < dim; ++d) {
        points(p, d) = in(p, d);
      }
    });
  Kokkos::View<Real*, DeviceMemorySpace> wraps(
    Kokkos::view_alloc(Kokkos::WithoutInitializing, "identification_wraps"), n);
  bool mapped = false;
  for (const auto& id : identifications) {
    ReduceAxis(points, wraps, id.axis, id.lower, id.length);
    if (id.map != nullptr) {
      id.map->Apply(MakeRank2View(points), MakeRank1View(wraps));
      mapped = true;
    }
  }
  // a map can move an axis that an earlier translation reduced
  if (mapped) {
    for (const auto& id : identifications) {
      if (id.map == nullptr) {
        ReduceAxis(points, wraps, id.axis, id.lower, id.length);
      }
    }
  }
  return points;
}

namespace
{

/// @return `coords`, or its points `reduced` if not empty
CoordinateView<DeviceMemorySpace> Located(
  const CoordinateView<DeviceMemorySpace>& coords,
  const Kokkos::View<Real**, DeviceMemorySpace>& reduced)
{
  if (reduced.data() == nullptr) {
    return coords;
  }
  return CoordinateView<DeviceMemorySpace>(coords.GetCoordinateSystem(),
                                           MakeConstRank2View(reduced));
}

} // namespace

std::shared_ptr<const BoundPoints> Geometry::Bind(
  CoordinateView<DeviceMemorySpace> coords, OutOfBoundsPolicy policy,
  std::shared_ptr<const CoordinateMap> map) const
{
  return Bind(std::move(coords), OwnerHints{}, policy, std::move(map));
}

std::shared_ptr<const BoundPoints> Geometry::Bind(
  CoordinateView<DeviceMemorySpace> coords, const OwnerHints& owners,
  OutOfBoundsPolicy policy, std::shared_ptr<const CoordinateMap> map) const
{
  CheckBindArguments(*this, coords, policy, map.get());
  // a hint naming another geometry, or not one entity per point, is not one;
  // nor is a hint given with a map, as it describes the unmapped points
  const OwnerHints hint =
    map == nullptr &&
        TakesHints(owners, static_cast<LO>(coords.GetValues().extent(0)))
      ? owners
      : OwnerHints{};
  std::shared_ptr<BoundPoints> bound;
  if (map != nullptr) {
    std::shared_ptr<const BoundCoordinateMap> bound_map = map->Bind(coords);
    const CoordinateView<DeviceMemorySpace> mapped = bound_map->MappedPoints();
    const auto reduced = detail::ReduceToPeriod(mapped);
    bound = BindImpl(BindRequest{Located(mapped, reduced), policy, hint,
                                 std::move(bound_map), coords, map});
  } else {
    const auto reduced = detail::ReduceToPeriod(coords);
    bound = BindImpl(
      BindRequest{Located(coords, reduced), policy, hint, nullptr, coords});
  }
  PCMS_ALWAYS_ASSERT(bound != nullptr);
  return bound;
}

std::shared_ptr<const BoundPoints> Geometry::Rebind(
  std::shared_ptr<const BoundPoints> previous,
  CoordinateView<DeviceMemorySpace> coords) const
{
  if (previous == nullptr) {
    throw pcms_error("Geometry::Rebind: the previous BoundPoints must not be "
                     "null");
  }
  if (previous->GetGeometry().get() != this) {
    throw pcms_error("Geometry::Rebind: the previous BoundPoints belongs to "
                     "another geometry");
  }
  if (previous->GetCoordinateMap() != nullptr || previous.use_count() != 1 ||
      !previous->Unshared()) {
    const OutOfBoundsPolicy policy = previous->Policy();
    std::shared_ptr<const CoordinateMap> map = previous->GetCoordinateMap();
    previous = Bind(std::move(coords), policy, std::move(map));
    return previous;
  }
  CheckBindArguments(*this, coords, previous->Policy(), nullptr);
  // nothing else refers to `previous`, so writing it is not observable
  auto reused = std::const_pointer_cast<BoundPoints>(std::move(previous));
  const auto reduced = detail::ReduceToPeriod(coords);
  reused->Rebind(
    RebindKey{},
    BindRequest{
      Located(coords, reduced), reused->Policy(), {}, nullptr, coords});
  return reused;
}

SupportResults Geometry::Neighbourhood(CoordinateView<DeviceMemorySpace>, int,
                                       const NeighbourhoodOptions&,
                                       const OwnerHints&) const
{
  throw pcms_error("Geometry::Neighbourhood: this geometry has no "
                   "neighbourhood search");
}

Quadrature Geometry::CreateQuadrature(const detail::QuadratureRule&) const
{
  throw pcms_error("CreateQuadrature: this geometry has no quadrature of "
                   "that rule");
}

Quadrature detail::CreateQuadrature(const Geometry& geometry,
                                    const QuadratureRule& rule)
{
  return geometry.CreateQuadrature(rule);
}

Quadrature CreateGaussQuadrature(
  const std::shared_ptr<const Geometry>& geometry, int order)
{
  if (order < 0) {
    throw pcms_error("CreateGaussQuadrature: the order must not be negative");
  }
  return detail::CreateQuadrature(*geometry, detail::GaussRule{order});
}

Quadrature CreateTrapezoidQuadrature(
  const std::shared_ptr<const Geometry>& geometry, int num_points)
{
  if (num_points < 1) {
    throw pcms_error("CreateTrapezoidQuadrature: needs a positive number of "
                     "points");
  }
  return detail::CreateQuadrature(*geometry, detail::TrapezoidRule{num_points});
}

Quadrature CreateMonteCarloQuadrature(
  const std::shared_ptr<const Geometry>& geometry, int samples_per_cell,
  std::uint64_t seed)
{
  if (samples_per_cell < 1) {
    throw pcms_error("CreateMonteCarloQuadrature: the samples per cell must "
                     "be positive");
  }
  return detail::CreateQuadrature(
    *geometry, detail::MonteCarloRule{samples_per_cell, seed});
}

Quadrature MapQuadrature(const std::shared_ptr<const Geometry>& geometry,
                         const Quadrature& quadrature)
{
  const auto& system = quadrature.points->GetGeometry()->GetCoordinateSystem();
  const auto parameter = geometry->GetParameterSystem();
  if (parameter == nullptr ||
      !detail::SameCoordinateSystem(system, parameter)) {
    throw pcms_error("MapQuadrature: the quadrature is on a geometry in '" +
                     std::string(system->Kind()) +
                     "', not in the geometry's parameter system");
  }
  const auto points = quadrature.points->LocatedPoints();
  Kokkos::View<Real*, DeviceMemorySpace> weights(
    Kokkos::view_alloc(Kokkos::WithoutInitializing, "mapped_weights"),
    quadrature.weights.extent(0));
  Kokkos::deep_copy(weights, quadrature.weights);
  return detail::WeightedQuadrature(*geometry, points, weights);
}

void detail::GaussLegendre(int order, std::vector<Real>& nodes,
                           std::vector<Real>& weights)
{
  const int n = order / 2 + 1;
  nodes.assign(n, 0);
  weights.assign(n, 0);
  for (int i = 0; i < n; ++i) {
    Real x = std::cos(std::numbers::pi_v<Real> * (i + 0.75) / (n + 0.5));
    Real dp = 1;
    for (int it = 0; it < 100; ++it) {
      Real p0 = 1;
      Real p1 = x;
      for (int k = 2; k <= n; ++k) {
        const Real p2 = ((2 * k - 1) * x * p1 - (k - 1) * p0) / k;
        p0 = p1;
        p1 = p2;
      }
      dp = n * (x * p1 - p0) / (x * x - 1);
      const Real dx = p1 / dp;
      x -= dx;
      if (std::abs(dx) < 1e-16) {
        break;
      }
    }
    nodes[n - 1 - i] = x;
    weights[n - 1 - i] = 2 / ((1 - x * x) * dp * dp);
  }
}

void detail::CheckNeighbourhoodArguments(
  const Geometry& geometry, const CoordinateView<DeviceMemorySpace>& coords,
  const NeighbourhoodOptions& options, const std::string& caller)
{
  if (!pcms::detail::SameCoordinateSystem(coords.GetCoordinateSystem(),
                                          geometry.GetCoordinateSystem())) {
    throw pcms_error(caller + ": query coordinate system '" +
                     std::string(coords.GetCoordinateSystem()->Kind()) +
                     "' is not the geometry's coordinate system '" +
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

void detail::ScaleByVolumeFactor(Kokkos::View<Real*, DeviceMemorySpace> weights,
                                 const BoundPoints& bound)
{
  const auto volume = bound.VolumeFactor();
  if (volume.extent(0) == 0) {
    return;
  }
  // a constant factor 1 leaves the weights as they are
  if (volume.stride(0) == 0) {
    Real factor = 0;
    Kokkos::deep_copy(factor,
                      Kokkos::View<const Real, DeviceMemorySpace,
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
                        volume.data_handle()));
    if (factor == 1) {
      return;
    }
  }
  pcms::detail::ForEach(
    "ScaleByVolumeFactor", static_cast<LO>(weights.extent(0)),
    KOKKOS_LAMBDA(const LO i) { weights(i) *= volume(i); });
  Kokkos::fence("ScaleByVolumeFactor");
}

std::shared_ptr<const BoundPoints> detail::BindOn(const Quadrature& quadrature,
                                                  const Geometry& geometry,
                                                  OutOfBoundsPolicy policy)
{
  if (quadrature.points->GetGeometry().get() == &geometry) {
    return quadrature.points;
  }
  for (const auto& bound : quadrature.also) {
    if (bound->GetGeometry().get() == &geometry) {
      return bound;
    }
  }
  return geometry.Bind(quadrature.points->LocatedPoints(), policy);
}

Quadrature detail::WeightedQuadrature(
  const Geometry& geometry, const CoordinateView<DeviceMemorySpace>& points,
  Kokkos::View<Real*, DeviceMemorySpace> weights, const OwnerHints& owners)
{
  auto bound = geometry.Bind(points, owners);
  ScaleByVolumeFactor(weights, *bound);
  return Quadrature{std::move(bound), weights, {}};
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
