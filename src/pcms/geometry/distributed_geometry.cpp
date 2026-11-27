#include "pcms/geometry/distributed_geometry.hpp"
#include "pcms/geometry/distributed_binding.hpp"
#include "pcms/utility/arrays.hpp"
#include "pcms/utility/assert.hpp"
#include "pcms/utility/detail/arrays.hpp"
#include "pcms/utility/mpi_type.hpp"
#include "pcms/utility/mpi_utilities.hpp"

#include <mpi.h>

#include <exception>
#include <string>
#include <vector>

namespace pcms
{

namespace
{

// Throws on every rank if `error` is non-empty on any rank, so a failure on
// one rank never leaves the others blocked in the next collective.
void ThrowIfAnyRankFailed(MPI_Comm comm, const std::string& context,
                          const std::string& error)
{
  const int local = error.empty() ? 0 : 1;
  int failed = 0;
  MPI_Allreduce(&local, &failed, 1, MPI_INT, MPI_MAX, comm);
  if (failed != 0) {
    throw pcms_error(
      context + ": " +
      (error.empty() ? std::string("failed on another rank") : error));
  }
}

// Per point of `bound`, 1 if it is found, else 0.
std::vector<char> FoundFlags(const BoundPoints& bound)
{
  std::vector<char> found(static_cast<size_t>(bound.NumPoints()), 1);
  const auto missing = Kokkos::create_mirror_view_and_copy(
    HostMemorySpace{}, bound.MissingPoints());
  for (size_t k = 0; k < missing.extent(0); ++k) {
    found[static_cast<size_t>(missing(k))] = 0;
  }
  return found;
}

std::vector<LO> ToHostVector(Kokkos::View<const LO*, DeviceMemorySpace> view)
{
  const auto h = Kokkos::create_mirror_view_and_copy(HostMemorySpace{}, view);
  return std::vector<LO>(h.data(), h.data() + h.extent(0));
}

// Hints naming `local`'s entities of dimension `entity_dim`, or none if an
// entity is unknown.
OwnerHints HintsOf(const std::shared_ptr<const Geometry>& local, int entity_dim,
                   const std::vector<LO>& entities)
{
  for (const LO e : entities) {
    if (e < 0) {
      return {};
    }
  }
  OwnerHints hints;
  hints.geometry = local;
  hints.entity_dim = entity_dim;
  hints.entities = pcms::detail::CopyToDevice(entities);
  return hints;
}

// Sends each of `points` another rank owns to that rank. Collective.
DistributedRoutingPlan BuildRoutingPlan(
  const DistributedGeometry& geometry,
  Kokkos::View<const Real**, DeviceMemorySpace> points,
  const std::vector<PointOwner>& located)
{
  DistributedRoutingPlan plan;
  plan.comm = geometry.GetComm();
  plan.dim = static_cast<int>(points.extent(1));
  plan.n_points = static_cast<LO>(points.extent(0));
  if (located.size() != static_cast<size_t>(plan.n_points)) {
    throw pcms_error("DistributedGeometry::LocateOwners: returned " +
                     std::to_string(located.size()) + " owners for " +
                     std::to_string(plan.n_points) + " points");
  }

  // Split the query indices by owner.
  int rank = 0;
  int nranks = 1;
  MPI_Comm_rank(plan.comm, &rank);
  MPI_Comm_size(plan.comm, &nranks);
  std::vector<std::vector<LO>> by_owner(static_cast<size_t>(nranks));
  std::vector<LO> own;
  std::vector<LO> own_entities;
  std::vector<LO> unlocated;
  for (LO i = 0; i < plan.n_points; ++i) {
    const auto& o = located[static_cast<size_t>(i)];
    if (o.rank < 0) {
      unlocated.push_back(i);
    } else if (o.rank == rank) {
      own.push_back(i);
      own_entities.push_back(o.entity);
    } else {
      by_owner[static_cast<size_t>(o.rank)].push_back(i);
    }
  }
  std::vector<LO> remote;
  // the entity each remote point's owner located it in, in its numbering
  std::vector<LO> remote_entities;
  plan.remote_counts.assign(static_cast<size_t>(nranks), 0);
  for (int r = 0; r < nranks; ++r) {
    const auto& points_of_r = by_owner[static_cast<size_t>(r)];
    remote.insert(remote.end(), points_of_r.begin(), points_of_r.end());
    for (const LO i : points_of_r) {
      remote_entities.push_back(located[static_cast<size_t>(i)].entity);
    }
    plan.remote_counts[static_cast<size_t>(r)] =
      static_cast<int>(points_of_r.size());
  }
  plan.own_points = pcms::detail::CopyToDevice(own);
  plan.remote_points = pcms::detail::CopyToDevice(remote);
  plan.unlocated_points = pcms::detail::CopyToDevice(unlocated);

  // Send each remote point to its owner, once.
  const auto send = pcms::detail::CopyToHostRowMajor(
    pcms::detail::GatherRows(points, plan.remote_points));
  std::vector<int> incoming_scalar_counts;
  const auto incoming =
    ExchangeRecords(plan.comm, detail::Scaled(plan.remote_counts, plan.dim),
                    send, incoming_scalar_counts);
  plan.incoming_counts.resize(incoming_scalar_counts.size());
  for (size_t r = 0; r < plan.incoming_counts.size(); ++r) {
    plan.incoming_counts[r] =
      plan.dim > 0 ? incoming_scalar_counts[r] / plan.dim : 0;
  }
  plan.incoming_coords = pcms::detail::CopyRowMajorToDevice(incoming, plan.dim);

  // Each owner gets back the entity it located each of its incoming points
  // in, so its local bind verifies that entity instead of searching.
  const auto incoming_entities = ExchangeRecordsWithKnownCounts(
    plan.comm, plan.remote_counts, remote_entities, plan.incoming_counts);
  const auto& local = geometry.GetLocalGeometry();
  const int entity_dim = geometry.OwnerEntityDim();
  plan.own_hints = HintsOf(local, entity_dim, own_entities);
  plan.incoming_hints = HintsOf(local, entity_dim, incoming_entities);
  return plan;
}

} // namespace

DistributedGeometry::DistributedGeometry(std::shared_ptr<const Geometry> local,
                                         MPI_Comm comm)
  : local_(std::move(local)), comm_(comm)
{
  if (local_ == nullptr) {
    throw pcms_error(
      "DistributedGeometry: the local geometry must not be null");
  }
}

std::shared_ptr<BoundPoints> DistributedGeometry::BindImpl(
  const BindRequest& request) const
{
  return std::make_shared<DistributedBinding>(
    std::static_pointer_cast<const DistributedGeometry>(shared_from_this()),
    request);
}

DistributedBinding::DistributedBinding(
  std::shared_ptr<const DistributedGeometry> geometry,
  const BindRequest& request)
  : BoundPoints(geometry, request), distributed_(*geometry)
{
  Locate(request.owners);
}

void DistributedBinding::Locate(const OwnerHints& owners)
{
  const auto points = LocatedPointValues();
  plan_ = BuildRoutingPlan(distributed_, points,
                           distributed_.LocateOwners(points, owners));

  // The owner of a point found it with the same local geometry, so the local
  // binds find it again; they fill rather than throw, and the policy is
  // applied below, on every rank together.
  OutOfBoundsPolicy fill;
  fill.mode = OutOfBoundsMode::FILL;
  fill.fill_value = Policy().fill_value;
  const auto& local = *distributed_.GetLocalGeometry();
  const auto& system = local.GetCoordinateSystem();
  std::string error;
  try {
    const auto own_coords = pcms::detail::GatherRows(points, plan_.own_points);
    // the elements the routing found, verified rather than searched for
    own_bound_ = local.Bind(CoordinateView<DeviceMemorySpace>(
                              system, pcms::MakeConstRank2View(own_coords)),
                            plan_.own_hints, fill);
    incoming_bound_ =
      local.Bind(CoordinateView<DeviceMemorySpace>(
                   system, pcms::MakeConstRank2View(plan_.incoming_coords)),
                 plan_.incoming_hints, fill);
  } catch (const std::exception& e) {
    error = e.what();
  }
  ThrowIfAnyRankFailed(plan_.comm, "DistributedGeometry::Bind", error);

  // Each owner tells the requester which of its points it found.
  const auto remote_found = ExchangeRecordsWithKnownCounts(
    plan_.comm, plan_.incoming_counts, FoundFlags(*incoming_bound_),
    plan_.remote_counts);
  const auto own_found = FoundFlags(*own_bound_);

  const auto n = static_cast<size_t>(NumPoints());
  Kokkos::View<bool*, HostMemorySpace> found_h("found", n);
  Kokkos::deep_copy(found_h, false);
  const auto own = ToHostVector(plan_.own_points);
  const auto remote = ToHostVector(plan_.remote_points);
  for (size_t k = 0; k < own.size(); ++k) {
    found_h(static_cast<size_t>(own[k])) = own_found[k] != 0;
  }
  for (size_t k = 0; k < remote.size(); ++k) {
    found_h(static_cast<size_t>(remote[k])) = remote_found[k] != 0;
  }
  LO missing = 0;
  for (size_t i = 0; i < n; ++i) {
    missing += found_h(i) ? 0 : 1;
  }

  // Under ERROR every rank throws if a point is missing on any rank.
  LO any_missing = 0;
  MPI_Allreduce(&missing, &any_missing, 1, utility::GetMPIType(missing),
                MPI_MAX, plan_.comm);
  if (Policy().mode == OutOfBoundsMode::ERROR && any_missing > 0) {
    throw pcms_error(missing > 0
                       ? "DistributedGeometry::Bind: a query point lies "
                         "outside the domain and the out-of-bounds policy is "
                         "ERROR"
                       : "DistributedGeometry::Bind: a query point of "
                         "another rank lies outside the domain and the "
                         "out-of-bounds policy is ERROR");
  }
  const auto found =
    Kokkos::create_mirror_view_and_copy(DeviceMemorySpace{}, found_h);
  SetFound(found);
}

} // namespace pcms
