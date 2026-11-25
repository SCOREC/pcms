#ifndef PCMS_GEOMETRY_DISTRIBUTED_ROUTING_PLAN_HPP
#define PCMS_GEOMETRY_DISTRIBUTED_ROUTING_PLAN_HPP

#include "pcms/geometry/bound_points.hpp"
#include "pcms/utility/memory_spaces.hpp"
#include "pcms/utility/types.hpp"

#include <Kokkos_Core.hpp>
#include <mpi.h>

#include <cstddef>
#include <vector>

namespace pcms
{

/**
 * \brief Which rank answers for each of this rank's query points, and which
 * points this rank answers for the others.
 *
 * Built by a DistributedBinding from the owners its geometry located; valid
 * for the query points and the partition it was built from.
 */
struct DistributedRoutingPlan
{
  /**
   * Communicator of the partition the plan was built on. Every exchange that
   * uses the plan runs on it.
   */
  MPI_Comm comm = MPI_COMM_NULL;
  /** Spatial dimension of the query points. */
  int dim = 0;
  /** Number of query points on this rank. */
  LO n_points = 0;

  /** Query points this rank owns, as indices into the query set, ascending. */
  Kokkos::View<LO*, DeviceMemorySpace> own_points;

  /**
   * Query points another rank owns, grouped by owner rank in ascending order.
   */
  Kokkos::View<LO*, DeviceMemorySpace> remote_points;
  /** How many of `remote_points` go to each rank. */
  std::vector<int> remote_counts;

  /**
   * Coordinates other ranks sent this rank to evaluate, `[n][dim]`, grouped
   * by requesting rank.
   */
  Kokkos::View<Real**, DeviceMemorySpace> incoming_coords;
  /** How many of `incoming_coords` came from each rank. */
  std::vector<int> incoming_counts;

  /** Query points no rank could locate, ascending. */
  Kokkos::View<LO*, DeviceMemorySpace> unlocated_points;

  /**
   * Claimed owners of `own_points` and of the incoming points, in their order,
   * as entities of the local geometry: what the owners already found, for
   * the local binds to verify instead of searching. Empty when an owner's
   * entity is unknown.
   */
  OwnerHints own_hints;
  OwnerHints incoming_hints;

  /**
   * True when this rank owns every one of its query points: nothing is sent
   * away and nothing is unlocated. Other ranks may still send points here.
   */
  [[nodiscard]] bool AllOwnedHere() const noexcept
  {
    return static_cast<LO>(own_points.extent(0)) == n_points;
  }

  /** How many query points get their value from another rank. */
  [[nodiscard]] std::size_t RoutedPointCount() const noexcept
  {
    return remote_points.extent(0);
  }

  /** How many points other ranks asked this rank to evaluate. */
  [[nodiscard]] std::size_t IncomingPointCount() const noexcept
  {
    return incoming_coords.extent(0);
  }
};

namespace detail
{

// `counts` times `factor`, for sending `factor` scalars per point.
inline std::vector<int> Scaled(const std::vector<int>& counts, int factor)
{
  std::vector<int> out(counts.size());
  for (size_t r = 0; r < counts.size(); ++r) {
    out[r] = counts[r] * factor;
  }
  return out;
}

} // namespace detail

} // namespace pcms

#endif // PCMS_GEOMETRY_DISTRIBUTED_ROUTING_PLAN_HPP
