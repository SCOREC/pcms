#include <mpi.h>
#include <Kokkos_Core.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_mesh.hpp>

#include "pcms/coupler/field_exchange_planner.h"
#include "pcms/coupler/field_serializer.h"
#include "pcms/coupler/overlap_mask.h"
#include "pcms/field/data/simple.h"
#include "pcms/field/function_space/lagrange.h"
#include "pcms/utility/arrays.h"
#include "pcms/utility/types.h"

#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <set>
#include <vector>

using pcms::GO;
using pcms::LO;
using pcms::Real;

namespace
{
void require(bool cond, const char* msg)
{
  if (!cond) {
    int rank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::fprintf(stderr, "[rank %d] distributed coupling check failed: %s\n",
                 rank, msg);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
}

// Builds a one-block GID message advertising `gids`, i.e. the structure the
// exchange planner produces for a single destination.
std::vector<pcms::GO> MakeGidMessage(const std::vector<pcms::GO>& gids)
{
  require(!gids.empty(), "cannot build an empty GID message");
  std::vector<pcms::GO> msg(
    static_cast<size_t>(pcms::ent_offsets_len) + gids.size(), 0);
  // Every DOF holder in these tests sits on a mesh vertex (entity dim 0), so
  // the per-entity-dimension header counts equal the total GID count from slot
  // one onwards.
  for (int e = 1; e < pcms::ent_offsets_len; ++e) {
    msg[static_cast<size_t>(e)] = static_cast<pcms::GO>(gids.size());
  }
  std::copy(gids.begin(), gids.end(), msg.begin() + pcms::ent_offsets_len);
  return msg;
}

// In-message layout for a one-block message (one sender group) received by
// `rank`. The planner reads only the entries describing this rank's segment of
// the receive buffer: its start and its end.
redev::InMessageLayout MakeInMessageLayout(int rank, int nproc, size_t msg_size)
{
  redev::InMessageLayout layout;
  layout.srcRanks.assign(static_cast<size_t>(nproc), 0);
  layout.offset.assign(static_cast<size_t>(nproc) + 1, 0);
  layout.offset[static_cast<size_t>(rank)] = 0;
  layout.offset[static_cast<size_t>(rank) + 1] =
    static_cast<redev::LO>(msg_size);
  return layout;
}

// ---------------------------------------------------------------------------
// Checks
// ---------------------------------------------------------------------------

// Shared state for the distributed coupling checks: the (ghosted) mesh, its
// order-1 Lagrange space and the rank counts. Everything else the checks need
// (counts, GIDs, masks) is derived from the layout.
struct DistributedCouplingContext
{
  Omega_h::Mesh& mesh;
  std::shared_ptr<pcms::LagrangeFunctionSpace> space;
  int rank;
  int nproc;

  const pcms::FieldLayout& Layout() const { return *space->GetLayout(); }
  LO NumOwned() const { return Layout().GetNumOwnedDofHolder(); }
  LO NumLocal() const { return Layout().GetNumLocalDofHolder(); }
  int NumComponents() const { return Layout().GetNumComponents(); }
};

// Exercises the owned-only serialization contract used by distributed coupling:
// only owned (rank-exclusive) DOF holders are put on the wire, and received
// values are scattered back into the local (owned + ghost) field.
void check_serializer_round_trip(const DistributedCouplingContext& ctx)
{
  const pcms::FieldLayout& layout = ctx.Layout();
  const LO n_owned = ctx.NumOwned();
  const LO n_local = ctx.NumLocal();
  const LO n_comp = ctx.NumComponents();
  const auto owned_to_local = layout.GetOwnedToLocalHost();
  require(n_owned > 0, "rank owns at least one DOF holder");
  require(n_owned >= 2, "rank owns at least two DOF holders");
  require(n_local > n_owned, "ghosted mesh has local > owned");

  // Set local field data: value = local index.
  auto field = ctx.space->CreateFunction<Real>();
  std::vector<Real> local_data(static_cast<size_t>(n_local));
  for (LO i = 0; i < n_local; ++i) {
    local_data[static_cast<size_t>(i)] = static_cast<Real>(i);
  }
  field.SetDOFHolderDataHost(pcms::Rank2View<const Real, pcms::HostMemorySpace>(
    local_data.data(), n_local, n_comp));

  // Owned-indexed permutation: the last owned holder is "outside the overlap
  // region" (perm = -1); the rest get compact buffer slots 0..n_participating.
  const LO n_participating = n_owned - 1;
  std::vector<LO> permutation(static_cast<size_t>(n_owned), -1);
  for (LO o = 0; o < n_participating; ++o) {
    permutation[static_cast<size_t>(o)] = o;
  }

  // Serialize: only participating owned data is written to the buffer.
  pcms::FieldSerializer<Real> serializer;
  std::vector<Real> buffer(static_cast<size_t>(n_participating) * n_comp);
  const int sent =
    serializer.Serialize(field.GetData(), layout, pcms::make_array_view(buffer),
                         pcms::make_const_array_view(permutation));
  require(sent == static_cast<int>(n_participating) * n_comp,
          "serialized size matches participating holders");

  auto owned_field =
    pcms::FlattenToRank1View(field.GetOwnedDOFHolderDataHost());
  require(owned_field.size() == static_cast<size_t>(n_owned),
          "owned field size");
  for (LO o = 0; o < n_participating; ++o) {
    require(buffer[static_cast<size_t>(o)] ==
              owned_field[static_cast<size_t>(o)],
            "serialized value == owned value");
  }

  // Deserialize into a fresh field: participating owned values land at their
  // local positions; ghost and non-participating owned slots stay zero.
  auto target = ctx.space->CreateFunction<Real>();
  serializer.Deserialize(target.GetData(), layout,
                         pcms::make_const_array_view(buffer),
                         pcms::make_const_array_view(permutation));
  auto target_data = pcms::FlattenToRank1View(target.GetDOFHolderDataHost());
  require(target_data.size() == static_cast<size_t>(n_local),
          "target local size");

  for (LO o = 0; o < n_participating; ++o) {
    const LO local = owned_to_local(o);
    require(target_data[static_cast<size_t>(local)] ==
              local_data[static_cast<size_t>(local)],
            "participating owned value restored");
  }

  const LO non_participating_local = owned_to_local(n_owned - 1);
  require(target_data[static_cast<size_t>(non_participating_local)] == Real(0),
          "non-participating owned holder is zero");

  auto owned_mask = layout.GetOwnedHost();
  for (LO i = 0; i < n_local; ++i) {
    if (!owned_mask[static_cast<size_t>(i)]) {
      require(target_data[static_cast<size_t>(i)] == Real(0),
              "ghost holder is zero");
    }
  }
}

// Exercises the distributed exchange end to end: an ExchangePlan packed into a
// GID message, GID matching on the receive side, and the payload round trip.
// The redev transport is replaced by MPI collectives (Alltoall/Alltoallv), and
// the exchange partners are the ranks that actually own this rank's ghost DOF
// holders, so the check holds for any partitioning and rank count.
void check_distributed_exchange(const DistributedCouplingContext& ctx)
{
  Omega_h::Mesh& mesh = ctx.mesh;
  const pcms::FieldLayout& layout = ctx.Layout();
  const int rank = ctx.rank;
  const int nproc = ctx.nproc;
  const LO n_owned = ctx.NumOwned();
  const LO n_local = ctx.NumLocal();
  const LO n_comp = ctx.NumComponents();
  const auto owned_to_local = layout.GetOwnedToLocalHost();
  const auto owned_mask = layout.GetOwnedHost();
  auto owned_gids = layout.GetOwnedGidsHost();
  auto local_gids = layout.GetGidsHost();
  auto owner_remotes = mesh.ask_owners(0);
  auto owner_ranks = Omega_h::HostRead<Omega_h::I32>(owner_remotes.ranks);

  // The real exchange partners: this rank's ghost DOF holders, grouped by the
  // rank that owns them. A ghost owned by rank n is a DOF this rank would have
  // to obtain from n in a coupled exchange.
  std::vector<int> neighbor_set;
  std::vector<std::vector<GO>> ghosts_from(static_cast<size_t>(nproc));
  {
    std::set<int> seen;
    for (LO i = 0; i < n_local; ++i) {
      if (!owned_mask[static_cast<size_t>(i)]) {
        const int owner = static_cast<int>(owner_ranks[static_cast<size_t>(i)]);
        ghosts_from[static_cast<size_t>(owner)].push_back(
          static_cast<GO>(local_gids(i)));
        seen.insert(owner);
      }
    }
    neighbor_set.assign(seen.begin(), seen.end());
  }
  require(!neighbor_set.empty(),
          "the ghosted mesh should have at least one neighbour");

  // Send-side plan: the planner is driven by a synthetic peer partition that
  // routes every owned holder to one actual neighbour. The routing is synthetic
  // (a single mesh has no peer partition to query), but it drives the same
  // partition-based plan construction a coupled application uses.
  const int destination = neighbor_set.front();
  std::vector<int> ptn_ranks = {destination, rank};
  std::vector<double> ptn_cuts = {0.0, 2.0}; // heap layout: cuts[1] is the cut
  redev::RCBPtn peer_ptn(2, ptn_ranks, ptn_cuts);
  redev::Partition exchange_partition{peer_ptn};

  // Every owned holder participates; the unmapped (-1) path is exercised on the
  // receive side, where the advertised list is a strict subset of this rank's
  // owned holders.
  pcms::OverlapMask overlap_mask(static_cast<size_t>(n_owned));

  pcms::GenericFieldExchangePlanner planner;
  pcms::ExchangePlan send_plan =
    planner.BuildExchangePlan(layout, exchange_partition, &overlap_mask);

  // Plan structure: a single destination (an actual neighbour) described by a
  // valid CSR, and a permutation that compacts the owned holders in owned
  // order.
  require(send_plan.dest_ranks.size() == 1,
          "send plan should route to a single destination");
  require(send_plan.dest_ranks[0] == destination,
          "send plan destination should be an actual neighbour");
  require(send_plan.offsets.size() == 2, "send plan CSR size");
  require(send_plan.offsets[0] == 0, "send plan CSR start");
  require(send_plan.offsets[1] == static_cast<redev::LO>(send_plan.msg_size),
          "send plan CSR end != message size");
  require(send_plan.msg_size == static_cast<size_t>(n_owned),
          "every owned holder should participate");
  require(send_plan.permutation.size() == static_cast<size_t>(n_owned),
          "send plan permutation is owned-indexed");
  LO next_slot = 0;
  for (LO o = 0; o < n_owned; ++o) {
    require(send_plan.permutation[static_cast<size_t>(o)] == next_slot,
            "send plan permutation must compact owned holders in owned order");
    ++next_slot;
  }

  // GID message: one block whose header counts the exchanged holders per mesh
  // entity dimension and whose payload lists the owned GIDs.
  const size_t gid_header =
    send_plan.dest_ranks.size() * static_cast<size_t>(pcms::ent_offsets_len);
  std::vector<pcms::GO> gid_msg(send_plan.msg_size + gid_header, 0);
  planner.FillGidMessage(layout, send_plan,
                         pcms::Rank1View<pcms::GO, pcms::HostMemorySpace>(
                           gid_msg.data(), gid_msg.size()));
  require(gid_msg[0] == 0, "GID message header must start at zero");
  for (int e = 1; e < pcms::ent_offsets_len; ++e) {
    require(gid_msg[static_cast<size_t>(e)] ==
              static_cast<pcms::GO>(send_plan.msg_size),
            "GID message header must count the vertex holders");
  }
  for (LO o = 0; o < n_owned; ++o) {
    const auto slot =
      static_cast<size_t>(send_plan.permutation[static_cast<size_t>(o)]);
    require(gid_msg[gid_header + slot] == static_cast<pcms::GO>(owned_gids(o)),
            "GID message payload != owned GIDs");
  }

  // Exchange each ghost GID with its owner. Two-phase all-to-all: first the
  // per-neighbour counts, then the GID lists themselves. Afterwards ghost_recv
  // holds, grouped by sender rank, every GID this rank owns that some neighbour
  // admits holding -- i.e. the owned DOFs this rank must send values for.
  std::vector<int> send_counts(static_cast<size_t>(nproc), 0);
  std::vector<int> recv_counts(static_cast<size_t>(nproc), 0);
  std::vector<int> sdispls(static_cast<size_t>(nproc), 0);
  std::vector<int> rdispls(static_cast<size_t>(nproc), 0);
  for (int n = 0; n < nproc; ++n) {
    send_counts[static_cast<size_t>(n)] =
      static_cast<int>(ghosts_from[static_cast<size_t>(n)].size());
  }
  MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT,
               MPI_COMM_WORLD);
  std::vector<GO> ghost_send;
  for (int n = 0; n < nproc; ++n) {
    sdispls[static_cast<size_t>(n)] = static_cast<int>(ghost_send.size());
    ghost_send.insert(ghost_send.end(),
                      ghosts_from[static_cast<size_t>(n)].begin(),
                      ghosts_from[static_cast<size_t>(n)].end());
  }
  LO total_recv = 0;
  for (int n = 0; n < nproc; ++n) {
    rdispls[static_cast<size_t>(n)] = static_cast<int>(total_recv);
    total_recv += recv_counts[static_cast<size_t>(n)];
  }
  std::vector<GO> ghost_recv(static_cast<size_t>(total_recv), 0);
  MPI_Alltoallv(ghost_send.data(), send_counts.data(), sdispls.data(),
                MPI_INT64_T, ghost_recv.data(), recv_counts.data(),
                rdispls.data(), MPI_INT64_T, MPI_COMM_WORLD);
  require(total_recv > 0,
          "no neighbour admits holding any of this rank's owned DOFs");

  // The received GID list is what the receive plan is built from: it gathers
  // into a single block (one sender's message) every GID some neighbour admits
  // holding, then matches those against the holders this rank owns.
  const std::vector<pcms::GO> recv_gid_message = MakeGidMessage(ghost_recv);
  pcms::ExchangePlan recv_plan = planner.BuildReceivePlan(
    layout,
    pcms::GlobalIDView<pcms::HostMemorySpace>(recv_gid_message.data(),
                                              recv_gid_message.size()),
    rank, nproc, MakeInMessageLayout(rank, nproc, recv_gid_message.size()));

  require(recv_plan.permutation.size() == static_cast<size_t>(n_owned),
          "receive plan permutation is owned-indexed");
  require(recv_plan.msg_size == static_cast<size_t>(total_recv),
          "receive plan payload length != advertised GID count");
  require(recv_plan.dest_ranks.size() == 1 && recv_plan.offsets.size() == 2,
          "receive plan should describe a single sender");
  require(recv_plan.offsets[0] == 0 &&
            recv_plan.offsets[1] == static_cast<redev::LO>(recv_plan.msg_size),
          "receive plan CSR");
  // A rank's owned DOF can be a ghost on several neighbours (e.g. a corner
  // vertex shared by four ranks), so the same GID may legitimately be
  // advertised more than once. Matching only cares about the distinct set.
  const std::set<GO> advertised_to_me(ghost_recv.begin(), ghost_recv.end());
  LO n_matched = 0;
  for (LO o = 0; o < n_owned; ++o) {
    const bool advertised =
      advertised_to_me.count(static_cast<GO>(owned_gids(o))) > 0;
    const LO slot = recv_plan.permutation[static_cast<size_t>(o)];
    require((slot >= 0) == advertised,
            "receive plan matched the wrong owned holders");
    if (advertised) {
      ++n_matched;
      require(slot < static_cast<LO>(recv_plan.msg_size),
              "receive plan buffer index out of range");
    }
  }
  require(n_matched == static_cast<LO>(advertised_to_me.size()),
          "receive plan must match every advertised GID against an owned "
          "holder");
  require(n_matched < n_owned,
          "every owned holder was advertised: the unmatched-permutation path "
          "was not exercised");

  // Local packing: Serialize this rank's owned values (the GIDs themselves)
  // with the receive plan permutation and round-trip them through Deserialize
  // to confirm each matched holder maps onto its slot (checked just below). The
  // cross-rank broadcast after it sends the same values back to every
  // advertiser via the transpose of the GID exchange.
  auto exchange_field = ctx.space->CreateFunction<Real>();
  std::vector<Real> exchange_vals(static_cast<size_t>(n_local), Real(0));
  for (LO o = 0; o < n_owned; ++o) {
    exchange_vals[static_cast<size_t>(owned_to_local(o))] =
      static_cast<Real>(owned_gids(o));
  }
  exchange_field.SetDOFHolderDataHost(
    pcms::Rank2View<const Real, pcms::HostMemorySpace>(exchange_vals.data(),
                                                       n_local, n_comp));

  pcms::FieldSerializer<Real> exchange_serializer;
  std::vector<Real> payload(recv_plan.msg_size, Real(0));
  const int packed = exchange_serializer.Serialize(
    exchange_field.GetData(), layout, pcms::make_array_view(payload),
    pcms::make_const_array_view(recv_plan.permutation));
  require(packed == static_cast<int>(recv_plan.msg_size),
          "packed payload size != receive plan message size");

  // Local check of the packing: a matched holder's owned value lands at its
  // plan slot, and holders the plan did not match stay untouched.
  auto echo_field = ctx.space->CreateFunction<Real>();
  exchange_serializer.Deserialize(
    echo_field.GetData(), layout, pcms::make_const_array_view(payload),
    pcms::make_const_array_view(recv_plan.permutation));
  auto echo_data = pcms::FlattenToRank1View(echo_field.GetDOFHolderDataHost());
  for (LO o = 0; o < n_owned; ++o) {
    const LO local = owned_to_local(o);
    if (recv_plan.permutation[static_cast<size_t>(o)] >= 0) {
      require(echo_data[static_cast<size_t>(local)] ==
                static_cast<Real>(owned_gids(o)),
              "packed value for a matched holder != its owned value");
    } else {
      require(echo_data[static_cast<size_t>(local)] == Real(0),
              "unmatched holder must stay untouched by deserialization");
    }
  }
  for (LO i = 0; i < n_local; ++i) {
    if (!owned_mask[static_cast<size_t>(i)]) {
      require(echo_data[static_cast<size_t>(i)] == Real(0),
              "deserialization must not write ghost holders");
    }
  }

  // Cross-rank broadcast: the owner sends its value for every advertised GID
  // back to each rank that advertised it, in that rank's own order. This is the
  // same broadcast Omega_h's sync_array performs for ghost values -- it must
  // serve every advertiser, because a corner DOF may be requested by several
  // neighbours.
  std::vector<Real> return_buf(ghost_recv.size(), Real(0));
  for (size_t k = 0; k < ghost_recv.size(); ++k) {
    return_buf[k] = static_cast<Real>(ghost_recv[k]);
  }
  std::vector<Real> payload_recv(ghost_send.size(), Real(0));
  MPI_Alltoallv(return_buf.data(), recv_counts.data(), rdispls.data(),
                MPI_DOUBLE, payload_recv.data(), send_counts.data(),
                sdispls.data(), MPI_DOUBLE, MPI_COMM_WORLD);
  for (size_t k = 0; k < ghost_send.size(); ++k) {
    require(payload_recv[k] == static_cast<Real>(ghost_send[k]),
            "values returned for ghost holders != the owning rank's values");
  }
}
} // namespace

// Exercises the distributed coupling data path on a ghosted mesh: the
// owned-only serialization contract (check_serializer_round_trip) and the
// exchange plan, GID matching and payload round trip a coupler drives
// (check_distributed_exchange).
int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  int result = 0;
  {
    Omega_h::Library lib(&argc, &argv);
    auto world = lib.world();
    const int rank = world->rank();
    const int nproc = world->size();
    require(nproc >= 2, "test requires at least 2 MPI ranks");

    const int n = 4 * nproc;
    Omega_h::Mesh mesh = Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1.0, 1.0,
                                            0.0, /*nx=*/n, /*ny=*/n, /*nz=*/0,
                                            /*symmetric=*/false);
    mesh.set_parting(OMEGA_H_GHOSTED, 1, false);

    auto space = pcms::LagrangeFunctionSpace::FromMesh(
      mesh, /*order=*/1, /*num_components=*/1,
      pcms::CoordinateSystem::Cartesian, "global",
      pcms::LagrangeFunctionSpace::Backend::OmegaH);
    const DistributedCouplingContext ctx{mesh, space, rank, nproc};

    check_serializer_round_trip(ctx);

    check_distributed_exchange(ctx);

    if (rank == 0) {
      std::printf("distributed coupling test passed (nproc=%d)\n", nproc);
    }
  }
  MPI_Finalize();
  return result;
}
