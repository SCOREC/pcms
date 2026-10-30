#include "pcms/localization/distributed_point_search.h"
#include "pcms/utility/mpi_utilities.h"

#include <Omega_h_array.hpp>

#include <mpi.h>
#include <stdexcept>
#include <vector>

namespace pcms
{

namespace
{

/// A query point sent to a candidate rank.
template <int Dim>
struct LocateRequest
{
  LO point; // index on the requesting rank
  Real coords[Dim];
};

/// A candidate rank's answer for a point it owns.
template <int Dim>
struct LocateReply
{
  LO point; // index on the requesting rank
  LO element_id;
  GO global_element_id;
  Real barycentric[Dim + 1];
};

/// `divisions`, after checking it is the same on every rank. The check is
/// collective, so every rank reaches the same verdict and throws together.
LO AgreedGridDivisions(Omega_h::Mesh& mesh, LO divisions)
{
  const LO lowest = mesh.comm()->allreduce(divisions, Omega_h_Op::OMEGA_H_MIN);
  const LO highest = mesh.comm()->allreduce(divisions, Omega_h_Op::OMEGA_H_MAX);
  if (lowest != highest) {
    throw std::invalid_argument(
      "DistributedPointSearch: grid_divisions must be the same on every rank");
  }
  if (lowest < 0) {
    throw std::invalid_argument(
      "DistributedPointSearch: grid_divisions must not be negative");
  }
  return lowest;
}

} // namespace

template <int Dim>
DistributedPointSearch<Dim>::DistributedPointSearch(
  Omega_h::Mesh& mesh, const LocalSearch& local_search, Real tolerance_reach,
  LO grid_divisions)
  : mesh_(mesh),
    local_search_(local_search),
    routing_(mesh,
             CandidateRankTable<Dim>::RoutingGrid(
               mesh, AgreedGridDivisions(mesh, grid_divisions)),
             tolerance_reach),
    tolerance_reach_(tolerance_reach)
{
  if (mesh.dim() != Dim) {
    throw std::invalid_argument(
      "DistributedPointSearch: mesh dimension does not match the search");
  }
  rank_ = mesh.comm()->rank();

  // Ownership is consulted once per point, so read it once. `nghost_layers() >
  // 0` is also the guard for using `Mesh::owned` at all: on a mesh that has not
  // been partitioned Omega_h refuses to answer an ownership query on more than
  // one rank, and without ghosts every local element is owned anyway.
  has_ghost_layer_ = mesh.nghost_layers() > 0;
  if (has_ghost_layer_) {
    const auto owned = Omega_h::HostRead<Omega_h::I8>(mesh.owned(mesh.dim()));
    owned_elements_.assign(owned.begin(), owned.end());
  } else {
    owned_elements_.assign(static_cast<size_t>(mesh.nelems()), 1);
  }

  // Global ids identify an element across ranks. They are what lets the answer
  // to "which element" be checked across ranks; a mesh without
  // them still works, the id is just unavailable.
  if (mesh.has_tag(mesh.dim(), "global")) {
    const auto globals = Omega_h::HostRead<GO>(mesh.globals(mesh.dim()));
    global_ids_.assign(globals.begin(), globals.end());
  }
}

template <int Dim>
std::vector<typename DistributedPointSearch<Dim>::LocalClassification>
DistributedPointSearch<Dim>::ClassifyLocally(
  Kokkos::View<const Real* [Dim]> points) const
{
  const auto n_points = static_cast<LO>(points.extent(0));
  std::vector<LocalClassification> classification(
    static_cast<size_t>(n_points));

  if (n_points == 0) {
    return classification;
  }

  const auto results = local_search_(points);
  const auto owners = local_search_.GetOwningElementIds(results);

  const auto results_h =
    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, results);
  const auto owners_h =
    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, owners);

  for (LO i = 0; i < n_points; ++i) {
    auto& out = classification[static_cast<size_t>(i)];
    const LO element_id = owners_h(i);
    // `IsOutside` is the search's own notion of "no element contains this
    // point": it covers both an empty candidate cell and a point the search
    // could only resolve to a nearest entity.
    if (results_h(i).IsOutside() || element_id < 0) {
      out.status = LocalStatus::NOT_FOUND_HERE;
      out.element_id = -1;
      continue;
    }
    out.element_id = element_id;
    out.barycentric = results_h(i).parametric_coords;
    out.status = OwnsLocalElement(element_id) ? LocalStatus::LOCATED_HERE
                                              : LocalStatus::LOCATED_ELSEWHERE;
  }
  return classification;
}

template <int Dim>
std::vector<typename DistributedPointSearch<Dim>::Result>
DistributedPointSearch<Dim>::Locate(Kokkos::View<const Real* [Dim]> points)
{
  const auto n_points = static_cast<LO>(points.extent(0));
  std::vector<Result> results(static_cast<size_t>(n_points));

  const auto classified = ClassifyLocally(points);
  const auto points_h =
    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, points);
  const auto& grid = routing_.GetGrid();
  const int nranks = mesh_.comm()->size();

  // ---- 1. local search, 2. candidates from the global grid ----
  std::vector<std::vector<LocateRequest<Dim>>> requests_by_rank(
    static_cast<size_t>(nranks));
  last_routed_ = 0;
  for (LO i = 0; i < n_points; ++i) {
    auto& r = results[static_cast<size_t>(i)];
    const auto& c = classified[static_cast<size_t>(i)];
    if (c.status == LocalStatus::LOCATED_HERE) {
      r.status = ResultStatus::LOCATED;
      r.owner_rank = rank_;
      r.element_id = c.element_id;
      r.global_element_id = GlobalElementId(c.element_id);
      r.barycentric = c.barycentric;
      continue;
    }
    LocateRequest<Dim> request{};
    request.point = i;
    Omega_h::Vector<Dim> point;
    for (int d = 0; d < Dim; ++d) {
      point[d] = points_h(i, d);
      request.coords[d] = points_h(i, d);
    }
    // Outside the global bounding box no element can contain the point, so it
    // stays NOT_FOUND and no rank is asked about it.
    if (!grid.IsPointInBounds(point, tolerance_reach_)) {
      continue;
    }
    ++last_routed_;
    // This rank already answered for anything it owns, so it never asks
    // itself. That also stops a ghost hit from bouncing straight back here.
    for (const int candidate : routing_.CandidatesForPoint(point)) {
      if (candidate != rank_) {
        requests_by_rank[static_cast<size_t>(candidate)].push_back(request);
      }
    }
  }

  // The exchange is collective, so whether to run it is decided on the global
  // total: a rank with nothing to send must still join if any other rank sends.
  LO send_total = 0;
  for (const auto& requests : requests_by_rank) {
    send_total += static_cast<LO>(requests.size());
  }
  if (mesh_.comm()->allreduce(send_total, Omega_h_Op::OMEGA_H_SUM) == 0) {
    return results;
  }

  // ---- 3. send to candidates, which search locally ----
  auto comm = mesh_.comm()->get_impl();
  std::vector<int> send_counts(static_cast<size_t>(nranks));
  std::vector<LocateRequest<Dim>> send;
  for (int r = 0; r < nranks; ++r) {
    const auto& requests = requests_by_rank[static_cast<size_t>(r)];
    send_counts[static_cast<size_t>(r)] = static_cast<int>(requests.size());
    send.insert(send.end(), requests.begin(), requests.end());
  }
  std::vector<int> recv_counts;
  const auto received = ExchangeRecords(comm, send_counts, send, recv_counts);

  Kokkos::View<Real* [Dim]> remote_points("remote query points",
                                          received.size());
  auto remote_h = Kokkos::create_mirror_view(remote_points);
  for (size_t t = 0; t < received.size(); ++t) {
    for (int d = 0; d < Dim; ++d) {
      remote_h(t, d) = received[t].coords[d];
    }
  }
  Kokkos::deep_copy(remote_points, remote_h);
  const auto remote_classified = ClassifyLocally(remote_points);

  // Reply only for points this rank owns. A candidate that cannot place a
  // point has nothing to say about it.
  std::vector<int> reply_counts(static_cast<size_t>(nranks), 0);
  std::vector<LocateReply<Dim>> replies;
  size_t t = 0;
  for (int r = 0; r < nranks; ++r) {
    for (int k = 0; k < recv_counts[static_cast<size_t>(r)]; ++k, ++t) {
      const auto& c = remote_classified[t];
      if (c.status != LocalStatus::LOCATED_HERE) {
        continue;
      }
      LocateReply<Dim> reply{};
      reply.point = received[t].point;
      reply.element_id = c.element_id;
      reply.global_element_id = GlobalElementId(c.element_id);
      for (int d = 0; d <= Dim; ++d) {
        reply.barycentric[d] = c.barycentric[d];
      }
      replies.push_back(reply);
      ++reply_counts[static_cast<size_t>(r)];
    }
  }

  // ---- 4. receive the answers; the lowest owning rank wins ----
  std::vector<int> reply_recv_counts;
  const auto answers =
    ExchangeRecords(comm, reply_counts, replies, reply_recv_counts);
  // Answers arrive grouped by source rank in ascending order, so the first
  // answer seen for a point is the lowest rank's.
  t = 0;
  for (int r = 0; r < nranks; ++r) {
    for (int k = 0; k < reply_recv_counts[static_cast<size_t>(r)]; ++k, ++t) {
      const auto& a = answers[t];
      if (a.point < 0 || a.point >= n_points) {
        continue;
      }
      auto& result = results[static_cast<size_t>(a.point)];
      if (result.status == ResultStatus::LOCATED) {
        continue;
      }
      result.status = ResultStatus::LOCATED;
      result.owner_rank = r;
      result.element_id = a.element_id;
      result.global_element_id = a.global_element_id;
      for (int d = 0; d <= Dim; ++d) {
        result.barycentric[d] = a.barycentric[d];
      }
    }
  }
  return results;
}

template class DistributedPointSearch<2>;
template class DistributedPointSearch<3>;

} // namespace pcms