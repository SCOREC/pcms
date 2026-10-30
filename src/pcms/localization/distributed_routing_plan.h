#ifndef PCMS_LOCALIZATION_DISTRIBUTED_ROUTING_PLAN_H
#define PCMS_LOCALIZATION_DISTRIBUTED_ROUTING_PLAN_H

#include "pcms/localization/distributed_point_search.h"
#include "pcms/utility/arrays.h"
#include "pcms/utility/memory_spaces.h"
#include "pcms/utility/mpi_utilities.h"
#include "pcms/utility/types.h"

#include <Kokkos_Core.hpp>
#include <Omega_h_comm.hpp>

#include <cstddef>
#include <vector>

namespace pcms
{

/**
 * \brief Which rank answers for each of this rank's query points, and which
 * points this rank answers for the others.
 *
 * This is the distributed counterpart of a serial evaluator's localization
 * hint. A serial hint records where each point is inside the local mesh, in a
 * form only its own backend understands. A routing plan records only which
 * rank owns each point. The owner then localizes the point with whatever
 * backend it uses, so the plan does not depend on the field or the function
 * space. One plan can be shared by evaluators of several function spaces on
 * the same mesh partition.
 *
 * Built by `BuildRoutingPlan`, which is collective. The plan is valid for the
 * query points and the mesh partition it was built from.
 */
struct DistributedRoutingPlan
{
  /// Communicator of the mesh the plan was built on. Every exchange that uses
  /// the plan runs on it.
  Omega_h::CommPtr comm;
  /// Spatial dimension of the query points.
  int dim = 0;
  /// Number of query points on this rank.
  LO n_points = 0;

  /// Query points this rank owns, as indices into the query set, ascending.
  Kokkos::View<LO*, DeviceMemorySpace> own_points;

  /// Query points another rank owns, grouped by owner rank in ascending order.
  Kokkos::View<LO*, DeviceMemorySpace> remote_points;
  /// How many of `remote_points` go to each rank.
  std::vector<int> remote_counts;

  /// Coordinates other ranks sent this rank to evaluate, `[n][dim]`, grouped
  /// by requesting rank.
  Kokkos::View<Real**, DeviceMemorySpace> incoming_coords;
  /// How many of `incoming_coords` came from each rank.
  std::vector<int> incoming_counts;

  /// Query points no rank could locate, ascending.
  Kokkos::View<LO*, DeviceMemorySpace> unlocated_points;

  /// True when this rank owns every one of its query points: nothing is sent
  /// away and nothing is unlocated. Other ranks may still send points here.
  [[nodiscard]] bool AllOwnedHere() const noexcept
  {
    return static_cast<LO>(own_points.extent(0)) == n_points;
  }

  /// How many query points get their value from another rank.
  [[nodiscard]] std::size_t RoutedPointCount() const noexcept
  {
    return remote_points.extent(0);
  }

  /// How many points other ranks asked this rank to evaluate.
  [[nodiscard]] std::size_t IncomingPointCount() const noexcept
  {
    return incoming_coords.extent(0);
  }
};

namespace detail
{

/// `counts` times `factor`, for sending `factor` scalars per point.
inline std::vector<int> Scaled(const std::vector<int>& counts, int factor)
{
  std::vector<int> out(counts.size());
  for (size_t r = 0; r < counts.size(); ++r) {
    out[r] = counts[r] * factor;
  }
  return out;
}

} // namespace detail

/**
 * \brief Locate `points` over the whole decomposition and send each point
 * another rank owns to that rank. Collective.
 *
 * Every rank must call it, with its own points (possibly none).
 */
template <int Dim>
[[nodiscard]] DistributedRoutingPlan BuildRoutingPlan(
  DistributedPointSearch<Dim>& search,
  Kokkos::View<const Real* [Dim], DeviceMemorySpace> points)
{
  using Search = DistributedPointSearch<Dim>;
  Omega_h::Mesh& mesh = search.GetMesh();

  DistributedRoutingPlan plan;
  plan.comm = mesh.comm();
  plan.dim = Dim;
  plan.n_points = static_cast<LO>(points.extent(0));

  const auto located = search.Locate(points);

  // Split the query indices by owner. The search result is host data, so this
  // is host work on indices only; coordinates stay on the device.
  const int rank = search.Rank();
  const int nranks = plan.comm->size();
  std::vector<std::vector<LO>> by_owner(static_cast<size_t>(nranks));
  std::vector<LO> own;
  std::vector<LO> unlocated;
  for (LO i = 0; i < plan.n_points; ++i) {
    const auto& r = located[static_cast<size_t>(i)];
    if (r.status != Search::ResultStatus::LOCATED) {
      unlocated.push_back(i);
    } else if (r.owner_rank == rank) {
      own.push_back(i);
    } else {
      by_owner[static_cast<size_t>(r.owner_rank)].push_back(i);
    }
  }
  std::vector<LO> remote;
  plan.remote_counts.assign(static_cast<size_t>(nranks), 0);
  for (int r = 0; r < nranks; ++r) {
    const auto& points_of_r = by_owner[static_cast<size_t>(r)];
    remote.insert(remote.end(), points_of_r.begin(), points_of_r.end());
    plan.remote_counts[static_cast<size_t>(r)] =
      static_cast<int>(points_of_r.size());
  }
  plan.own_points = CopyToDevice(own);
  plan.remote_points = CopyToDevice(remote);
  plan.unlocated_points = CopyToDevice(unlocated);

  // Send each remote point to its owner, once. Only the remote points'
  // coordinates leave the device, because MPI takes host buffers.
  const auto send = CopyToHostRowMajor(GatherRows(points, plan.remote_points));
  std::vector<int> incoming_scalar_counts;
  const auto incoming = ExchangeRecords(plan.comm->get_impl(),
                                        detail::Scaled(plan.remote_counts, Dim),
                                        send, incoming_scalar_counts);
  plan.incoming_counts.resize(incoming_scalar_counts.size());
  for (size_t r = 0; r < plan.incoming_counts.size(); ++r) {
    plan.incoming_counts[r] = incoming_scalar_counts[r] / Dim;
  }
  plan.incoming_coords = CopyRowMajorToDevice(incoming, Dim);
  return plan;
}

} // namespace pcms

#endif // PCMS_LOCALIZATION_DISTRIBUTED_ROUTING_PLAN_H
