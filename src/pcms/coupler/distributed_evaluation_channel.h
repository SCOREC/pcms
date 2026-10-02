#ifndef PCMS_DISTRIBUTED_EVALUATION_CHANNEL_H
#define PCMS_DISTRIBUTED_EVALUATION_CHANNEL_H
#include <algorithm>
#include <array>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <string>
#include <vector>
#include <Kokkos_Core.hpp>
#include <redev.h>
#include "pcms/coupler/distributed_evaluation_partition.h"
#include "pcms/field/distributed_partition.h"
#include "pcms/field/field_layout.h"
#include "pcms/utility/arrays.h"
#include "pcms/utility/assert.h"
#include "pcms/utility/memory_spaces.h"
#include "pcms/utility/types.h"

namespace pcms
{

/**
 * @brief Resolves distributed-field query points across ranks via a standalone
 * evaluation server.
 *
 * DistributedEvaluationChannel<T> resolves query points that a rank's local
 * evaluator could not localize by routing them through a standalone,
 * single-rank evaluation server over a redev::Channel, instead of through
 * direct peer-to-peer MPI. The server holds no field data: it only knows the
 * routing partition (EvaluationPartitionQuery) and relays points to whichever
 * field rank(s) that partition names as candidate owners, then relays their
 * answers back to the asking rank. See EvaluationPartitionQuery /
 * BoundingBoxEvaluationPartition for the routing hint itself.
 *
 * The routing partition is derived from the field's own mesh, not from a
 * separate partition file: before the first round each field rank publishes
 * its own owned bounding box (PublishOwnedBounds, driven automatically by
 * DistributedPointEvaluator) and the server assembles the partition from those
 * boxes (ReceiveRoutingPartition). Because the partition comes from exactly
 * the distribution the field ranks evaluate against, the two sides cannot
 * disagree, and no mesh or partition file is needed on the server.
 * EvaluationPartitionQuery remains the extension point: a different routing
 * source (e.g. a redev::Partition) can be supplied to RunServerRound later
 * without changing the protocol.
 *
 * One instance is created per (field, value type) pair, on both the field
 * job's ranks (role = client) and the standalone server process (role =
 * server; must be launched with exactly one rank -- see
 * pcms/tools/distributed_field_eval_server.cpp). Both sides must construct
 * their instance with the same `name` and call the round methods
 * collectively and in the same relative order.
 */
template <typename T>
class DistributedEvaluationChannel
{
public:
  // points: query coordinates local to this rank, [n][dim].
  // values: per-point results, [n][num_components]; entries for points that
  //   remain unresolved are left unmodified.
  // resolved: [n]; on input, true marks points already resolved locally (and
  //   thus skipped) and false marks points needing server-mediated
  //   resolution.
  using LocalEvaluate =
    std::function<void(Rank2View<const Real, HostMemorySpace> query_points,
                       Rank2View<T, HostMemorySpace> out_values,
                       Rank1View<bool, HostMemorySpace> out_resolved)>;

  /**
   * @brief Constructs one endpoint of the distributed-evaluation channel.
   *
   * @param expected_rounds declares how many resolution rounds this side will
   * run, and exists only as a diagnostic guard against the otherwise-silent
   * failure mode where the two sides disagree about the round count (see
   * CheckDeclaredRoundCount, which the setup exchange drives). On the client it
   * is the total number of Evaluate() calls the application intends to make; on
   * the server it is the number of RunServerRound() calls the driver will make.
   * -1 (the default) means "not declared", which skips the check. Note this
   * count is unrelated to the number of field ranks: one round fans out across
   * every rank.
   */
  DistributedEvaluationChannel(redev::Redev& redev, redev::Channel& channel,
                               MPI_Comm mpi_comm, const std::string& name,
                               int dim, int num_components,
                               int expected_rounds = -1)
    : redev_(redev),
      channel_(channel),
      mpi_comm_(mpi_comm),
      name_(name),
      dim_(dim),
      num_components_(num_components),
      expected_rounds_(expected_rounds),
      query_coords_comm_(
        channel.CreateComm<Real>(name + "_query_coords", mpi_comm)),
      query_meta_comm_(
        channel.CreateComm<LO>(name + "_query_meta", mpi_comm)),
      answer_values_comm_(
        channel.CreateComm<T>(name + "_answer_values", mpi_comm)),
      answer_meta_comm_(
        channel.CreateComm<LO>(name + "_answer_meta", mpi_comm)),
      answer_resolved_comm_(
        channel.CreateComm<int8_t>(name + "_answer_resolved", mpi_comm)),
      partition_coords_comm_(
        channel.CreateComm<Real>(name + "_partition_coords", mpi_comm)),
      partition_flags_comm_(
        channel.CreateComm<LO>(name + "_partition_flags", mpi_comm))
  {
    PCMS_ALWAYS_ASSERT(redev_.GetProcessType() == redev::ProcessType::Client ||
                       redev_.GetProcessType() == redev::ProcessType::Server);
  }

  /**
   * @brief Client-side half of the one-time channel setup exchange.
   *
   * Called once before the first resolution round (DistributedPointEvaluator
   * does it automatically). Every field rank reports its own owned bounding
   * box, derived from its own FieldLayout, and the server assembles the routing
   * partition from those. The same exchange also carries each side's declared
   * round count, so a mismatch is reported immediately rather than hanging
   * mid-run (see CheckDeclaredRoundCount). The boxes are what make the routing
   * partition consistent with the field's real mesh distribution by
   * construction: derived from the mesh, shipped over the channel, never read
   * from a separate file. Each rank contributes exactly one box, so no message
   * is ever empty. See DistributedBoundingBoxPartition::LocalOwnedBounds.
   */
  void PublishOwnedBounds(const FieldLayout& layout)
  {
    PCMS_ALWAYS_ASSERT(redev_.GetProcessType() == redev::ProcessType::Client);
    if (setup_exchanged_) {
      return;
    }
    setup_exchanged_ = true;

    const int dim = layout.GetDimension();
    PCMS_ALWAYS_ASSERT(dim > 0 && dim <= 3);
    std::array<Real, 3> box_min{};
    std::array<Real, 3> box_max{};
    bool owns_nothing = false;
    DistributedBoundingBoxPartition::LocalOwnedBounds(layout, dim, box_min,
                                                      box_max, owns_nothing);

    // One segment to the single server rank: coords are
    // [min[0..dim-1], max[0..dim-1]] and flags are
    // [owns_nothing, dim, declared_round_count].
    redev::LOs dest{0};
    redev::LOs coord_offsets{0, 2 * dim};
    redev::LOs flag_offsets{0, 3};
    partition_coords_comm_.SetOutMessageLayout(dest, coord_offsets);
    partition_flags_comm_.SetOutMessageLayout(dest, flag_offsets);

    std::vector<Real> coords(static_cast<size_t>(2 * dim));
    for (int d = 0; d < dim; ++d) {
      coords[static_cast<size_t>(d)] = box_min[static_cast<size_t>(d)];
      coords[static_cast<size_t>(dim + d)] = box_max[static_cast<size_t>(d)];
    }
    redev::LOs flags{owns_nothing ? 1 : 0, static_cast<redev::LO>(dim),
                     static_cast<redev::LO>(expected_rounds_)};

    // c2s: this rank's box and declared round count.
    channel_scope_.BeginSendCommunicationPhase(*this);
    partition_coords_comm_.Send(coords.data(), redev::Mode::Synchronous);
    partition_flags_comm_.Send(flags.data(), redev::Mode::Synchronous);
    channel_scope_.EndSendCommunicationPhase(*this);

    // s2c: the server's declared round count. Sent to client rank 0 only, so
    // the other client ranks see an empty message and skip the comparison.
    channel_scope_.BeginReceiveCommunicationPhase(*this);
    const auto server_rounds =
      partition_flags_comm_.Recv(redev::Mode::Synchronous);
    channel_scope_.EndReceiveCommunicationPhase(*this);
    CheckDeclaredRoundCount(
      expected_rounds_,
      server_rounds.empty() ? -1 : static_cast<int>(server_rounds[0]));
  }

  /**
   * @brief Server-side half of the one-time channel setup exchange.
   *
   * Called once before the round loop (the server driver does it explicitly).
   * Receives every field rank's box, returns the routing partition built from
   * them, and participates in the round-count handshake (see
   * CheckDeclaredRoundCount). Segments arrive in sender-rank order (redev
   * concatenates by rank), which is exactly the rank-major layout
   * DistributedBoundingBoxPartition::FromGatheredBounds expects.
   */
  [[nodiscard]] std::unique_ptr<EvaluationPartitionQuery>
  ReceiveRoutingPartition(Real padding = 1.0e-8)
  {
    PCMS_ALWAYS_ASSERT(redev_.GetProcessType() == redev::ProcessType::Server);
    if (setup_exchanged_) {
      return nullptr;
    }
    setup_exchanged_ = true;

    // c2s: every field rank's box and declared round count.
    channel_scope_.BeginReceiveCommunicationPhase(*this);
    auto coords = partition_coords_comm_.Recv(redev::Mode::Synchronous);
    auto flags = partition_flags_comm_.Recv(redev::Mode::Synchronous);
    channel_scope_.EndReceiveCommunicationPhase(*this);

    const size_t num_ranks = flags.size() / 3;
    PCMS_ALWAYS_ASSERT(num_ranks > 0);
    const int dim = static_cast<int>(flags[1]);
    PCMS_ALWAYS_ASSERT(dim > 0 && dim <= 3);
    // The channel was constructed with the dimension the driver expects; the
    // shipped boxes must agree with it or the protocol would misread messages.
    PCMS_ALWAYS_ASSERT(dim == dim_);
    PCMS_ALWAYS_ASSERT(coords.size() ==
                       num_ranks * static_cast<size_t>(2 * dim));

    std::vector<Real> mins(num_ranks * static_cast<size_t>(dim));
    std::vector<Real> maxs(num_ranks * static_cast<size_t>(dim));
    std::vector<char> owns_nothing(num_ranks, 0);
    const int client_rounds = static_cast<int>(flags[2]);
    for (size_t s = 0; s < num_ranks; ++s) {
      owns_nothing[s] = static_cast<char>(flags[s * 3] != 0);
      PCMS_ALWAYS_ASSERT(static_cast<int>(flags[s * 3 + 1]) == dim);
      // Every field rank must agree on the declared round count.
      PCMS_ALWAYS_ASSERT(static_cast<int>(flags[s * 3 + 2]) == client_rounds);
      for (int d = 0; d < dim; ++d) {
        mins[s * static_cast<size_t>(dim) + static_cast<size_t>(d)] =
          coords[s * static_cast<size_t>(2 * dim) + static_cast<size_t>(d)];
        maxs[s * static_cast<size_t>(dim) + static_cast<size_t>(d)] =
          coords[s * static_cast<size_t>(2 * dim) + static_cast<size_t>(dim) +
                 static_cast<size_t>(d)];
      }
    }

    // s2c: declare our own round count back (to client rank 0 only).
    redev::LOs dest{0};
    redev::LOs offsets{0, 1};
    partition_flags_comm_.SetOutMessageLayout(dest, offsets);
    redev::LO mine = static_cast<redev::LO>(expected_rounds_);
    channel_scope_.BeginSendCommunicationPhase(*this);
    partition_flags_comm_.Send(&mine, redev::Mode::Synchronous);
    channel_scope_.EndSendCommunicationPhase(*this);

    CheckDeclaredRoundCount(expected_rounds_, client_rounds);

    return std::make_unique<BoundingBoxEvaluationPartition>(
      DistributedBoundingBoxPartition::FromGatheredBounds(
        dim, padding, std::move(mins), std::move(maxs),
        std::move(owns_nothing)));
  }

  // Client-only. Every field rank must call this collectively once per
  // resolution round: it both asks the server about its own unresolved
  // points and, in the same round, answers any points the server redirects to
  // it as a candidate owner. Updates values/resolved in place for any point
  // the server relays back with resolved == true.
  void ResolveUnresolved(Rank2View<const Real, HostMemorySpace> points,
                        Rank2View<T, HostMemorySpace> values,
                        Rank1View<bool, HostMemorySpace> resolved,
                        const LocalEvaluate& local_evaluate)
  {
    PCMS_ALWAYS_ASSERT(redev_.GetProcessType() == redev::ProcessType::Client);
    // The one-time setup exchange (PublishOwnedBounds) performs the
    // round-count handshake; it must have run before any round.
    PCMS_ALWAYS_ASSERT(setup_exchanged_);
    int my_rank = 0;
    MPI_Comm_rank(mpi_comm_, &my_rank);
    const LO n = static_cast<LO>(points.extent(0));

    std::vector<LO> ask_local_index;
    std::vector<Real> ask_coords;
    for (LO i = 0; i < n; ++i) {
      if (resolved(i)) {
        continue;
      }
      ask_local_index.push_back(i);
      for (int d = 0; d < dim_; ++d) {
        ask_coords.push_back(points(i, d));
      }
    }
    const LO num_ask = static_cast<LO>(ask_local_index.size());
    std::vector<LO> ask_meta(static_cast<size_t>(num_ask) * 2);
    for (LO k = 0; k < num_ask; ++k) {
      ask_meta[static_cast<size_t>(k) * 2] = my_rank;
      ask_meta[static_cast<size_t>(k) * 2 + 1] = ask_local_index[static_cast<size_t>(k)];
    }

    // Phase 1 (ask): send every still-unresolved point to the single server
    // rank, as one contiguous segment.
    redev::LOs ask_dest{0};
    redev::LOs ask_coord_offsets{0, num_ask * dim_};
    redev::LOs ask_meta_offsets{0, num_ask * 2};
    query_coords_comm_.SetOutMessageLayout(ask_dest, ask_coord_offsets);
    query_meta_comm_.SetOutMessageLayout(ask_dest, ask_meta_offsets);
    channel_scope_.BeginSendCommunicationPhase(*this);
    query_coords_comm_.Send(ask_coords.data(), redev::Mode::Synchronous);
    query_meta_comm_.Send(ask_meta.data(), redev::Mode::Synchronous);
    channel_scope_.EndSendCommunicationPhase(*this);

    // Phase 2 (redirect): receive whatever points the server thinks this rank
    // might own (its own asks included, if it is its own candidate), evaluate
    // them locally, and send the answers back.
    channel_scope_.BeginReceiveCommunicationPhase(*this);
    auto redirected_coords = query_coords_comm_.Recv(redev::Mode::Synchronous);
    auto redirected_meta = query_meta_comm_.Recv(redev::Mode::Synchronous);
    channel_scope_.EndReceiveCommunicationPhase(*this);

    const LO num_redirected =
      static_cast<LO>(redirected_meta.size()) / 2;
    Rank2View<const Real, HostMemorySpace> redirected_points_view(
      redirected_coords.data(), num_redirected, dim_);
    std::vector<T> answer_values(static_cast<size_t>(num_redirected) *
                                static_cast<size_t>(num_components_));
    std::vector<bool> answer_resolved_bits(static_cast<size_t>(num_redirected));
    if (num_redirected > 0) {
      Rank2View<T, HostMemorySpace> answer_values_view(
        answer_values.data(), num_redirected, num_components_);
      // Rank1View<bool,...> needs real contiguous bool storage; std::vector
      // <bool> is bit-packed and not addressable, so stage through a
      // Kokkos::View (a plain array of bool, unlike std::vector<bool>).
      Kokkos::View<bool*, HostMemorySpace> resolved_storage(
        "resolved_storage", static_cast<size_t>(num_redirected));
      Rank1View<bool, HostMemorySpace> answer_resolved_view(
        resolved_storage.data(), num_redirected);
      local_evaluate(redirected_points_view, answer_values_view,
                    answer_resolved_view);
      for (LO k = 0; k < num_redirected; ++k) {
        answer_resolved_bits[static_cast<size_t>(k)] = resolved_storage(k);
      }
    }

    // Phase 3 (answer): send the (asker, local_index, value, resolved) tuples
    // back to the server, one per redirected point.
    std::vector<LO> answer_meta(static_cast<size_t>(num_redirected) * 2);
    std::vector<int8_t> answer_resolved_flags(static_cast<size_t>(num_redirected));
    for (LO k = 0; k < num_redirected; ++k) {
      answer_meta[static_cast<size_t>(k) * 2] =
        redirected_meta[static_cast<size_t>(k) * 2];
      answer_meta[static_cast<size_t>(k) * 2 + 1] =
        redirected_meta[static_cast<size_t>(k) * 2 + 1];
      answer_resolved_flags[static_cast<size_t>(k)] =
        answer_resolved_bits[static_cast<size_t>(k)] ? 1 : 0;
    }
    // Single destination (server rank 0), one contiguous segment.
    redev::LOs answer_dest{0};
    redev::LOs answer_val_offsets{0, num_redirected * num_components_};
    redev::LOs answer_meta_offsets{0, num_redirected * 2};
    redev::LOs answer_res_offsets{0, num_redirected};
    answer_values_comm_.SetOutMessageLayout(answer_dest, answer_val_offsets);
    answer_meta_comm_.SetOutMessageLayout(answer_dest, answer_meta_offsets);
    answer_resolved_comm_.SetOutMessageLayout(answer_dest, answer_res_offsets);
    channel_scope_.BeginSendCommunicationPhase(*this);
    answer_values_comm_.Send(answer_values.data(), redev::Mode::Synchronous);
    answer_meta_comm_.Send(answer_meta.data(), redev::Mode::Synchronous);
    answer_resolved_comm_.Send(answer_resolved_flags.data(),
                              redev::Mode::Synchronous);
    channel_scope_.EndSendCommunicationPhase(*this);

    // Phase 4 (relay): receive the final answers to this rank's own asks and
    // merge them in.
    channel_scope_.BeginReceiveCommunicationPhase(*this);
    auto relay_values = answer_values_comm_.Recv(redev::Mode::Synchronous);
    auto relay_index = answer_meta_comm_.Recv(redev::Mode::Synchronous);
    auto relay_resolved = answer_resolved_comm_.Recv(redev::Mode::Synchronous);
    channel_scope_.EndReceiveCommunicationPhase(*this);

    const LO num_relay = static_cast<LO>(relay_resolved.size());
    for (LO k = 0; k < num_relay; ++k) {
      if (relay_resolved[static_cast<size_t>(k)] == 0) {
        continue;
      }
      const LO local_index = relay_index[static_cast<size_t>(k)];
      if (resolved(local_index)) {
        continue;
      }
      for (int c = 0; c < num_components_; ++c) {
        values(local_index, c) =
          relay_values[static_cast<size_t>(k) *
                      static_cast<size_t>(num_components_) +
                      static_cast<size_t>(c)];
      }
      resolved(local_index) = true;
    }
  }

  // Server-only. Must be called collectively once per resolution round, after
  // every field rank has called ResolveUnresolved for that same round.
  void RunServerRound(const EvaluationPartitionQuery& partition)
  {
    PCMS_ALWAYS_ASSERT(redev_.GetProcessType() == redev::ProcessType::Server);
    // The one-time setup exchange (ReceiveRoutingPartition) performs the
    // round-count handshake; it must have run before any round.
    PCMS_ALWAYS_ASSERT(setup_exchanged_);

    // Phase 1 (ask): collect every field rank's unresolved points.
    channel_scope_.BeginReceiveCommunicationPhase(*this);
    auto ask_coords = query_coords_comm_.Recv(redev::Mode::Synchronous);
    auto ask_meta = query_meta_comm_.Recv(redev::Mode::Synchronous);
    channel_scope_.EndReceiveCommunicationPhase(*this);

    const LO num_ask = static_cast<LO>(ask_meta.size()) / 2;
    std::vector<int> owner_dest;
    std::vector<Real> redirect_coords;
    std::vector<LO> redirect_meta;
    std::vector<int> candidates;
    for (LO k = 0; k < num_ask; ++k) {
      candidates.clear();
      std::array<Real, 3> pt{0.0, 0.0, 0.0};
      for (int d = 0; d < dim_; ++d) {
        pt[static_cast<size_t>(d)] =
          ask_coords[static_cast<size_t>(k) * static_cast<size_t>(dim_) +
                    static_cast<size_t>(d)];
      }
      partition.GetCandidateRanks(pt, candidates);
      for (int owner : candidates) {
        owner_dest.push_back(owner);
        for (int d = 0; d < dim_; ++d) {
          redirect_coords.push_back(pt[static_cast<size_t>(d)]);
        }
        redirect_meta.push_back(ask_meta[static_cast<size_t>(k) * 2]);
        redirect_meta.push_back(ask_meta[static_cast<size_t>(k) * 2 + 1]);
      }
    }

    // Sort by destination rank so dest/offsets describe contiguous segments,
    // as SetOutMessageLayout requires.
    const LO num_redirect = static_cast<LO>(owner_dest.size());
    std::vector<LO> order(static_cast<size_t>(num_redirect));
    for (LO k = 0; k < num_redirect; ++k) {
      order[static_cast<size_t>(k)] = k;
    }
    std::stable_sort(order.begin(), order.end(), [&](LO a, LO b) {
      return owner_dest[static_cast<size_t>(a)] <
            owner_dest[static_cast<size_t>(b)];
    });

    std::vector<Real> sorted_coords(static_cast<size_t>(num_redirect) *
                                   static_cast<size_t>(dim_));
    std::vector<LO> sorted_meta(static_cast<size_t>(num_redirect) * 2);
    redev::LOs dest;
    std::vector<LO> coord_counts;
    std::vector<LO> meta_counts;
    int prev_dest = -1;
    for (LO out_k = 0; out_k < num_redirect; ++out_k) {
      const LO in_k = order[static_cast<size_t>(out_k)];
      for (int d = 0; d < dim_; ++d) {
        sorted_coords[static_cast<size_t>(out_k) * static_cast<size_t>(dim_) +
                      static_cast<size_t>(d)] =
          redirect_coords[static_cast<size_t>(in_k) * static_cast<size_t>(dim_) +
                          static_cast<size_t>(d)];
      }
      sorted_meta[static_cast<size_t>(out_k) * 2] =
        redirect_meta[static_cast<size_t>(in_k) * 2];
      sorted_meta[static_cast<size_t>(out_k) * 2 + 1] =
        redirect_meta[static_cast<size_t>(in_k) * 2 + 1];
      const int d = owner_dest[static_cast<size_t>(in_k)];
      if (d != prev_dest) {
        dest.push_back(d);
        coord_counts.push_back(0);
        meta_counts.push_back(0);
        prev_dest = d;
      }
      coord_counts.back() += dim_;
      meta_counts.back() += 2;
    }
    // dest[i] owns [offsets[i], offsets[i+1]); building the offsets from the
    // per-destination counts keeps offsets.size() == dest.size() + 1 exactly
    // (an empty dest still yields the single leading zero SetOutMessageLayout
    // needs). Seeding the offsets with a literal 0 *and* pushing the first
    // block's start would leave a duplicated leading zero, which makes redev
    // read a zero-length segment for the first destination and send nothing.
    redev::LOs coord_offsets = MakeCsrOffsets(coord_counts);
    redev::LOs meta_offsets = MakeCsrOffsets(meta_counts);

    // Phase 2 (redirect): send each point to its candidate owner rank(s).
    query_coords_comm_.SetOutMessageLayout(dest, coord_offsets);
    query_meta_comm_.SetOutMessageLayout(dest, meta_offsets);
    channel_scope_.BeginSendCommunicationPhase(*this);
    query_coords_comm_.Send(sorted_coords.data(), redev::Mode::Synchronous);
    query_meta_comm_.Send(sorted_meta.data(), redev::Mode::Synchronous);
    channel_scope_.EndSendCommunicationPhase(*this);

    // Phase 3 (answer): collect every owner rank's computed values.
    channel_scope_.BeginReceiveCommunicationPhase(*this);
    auto in_values = answer_values_comm_.Recv(redev::Mode::Synchronous);
    auto in_meta = answer_meta_comm_.Recv(redev::Mode::Synchronous);
    auto in_resolved = answer_resolved_comm_.Recv(redev::Mode::Synchronous);
    channel_scope_.EndReceiveCommunicationPhase(*this);

    const LO num_answers = static_cast<LO>(in_resolved.size());
    std::vector<LO> order2(static_cast<size_t>(num_answers));
    for (LO k = 0; k < num_answers; ++k) {
      order2[static_cast<size_t>(k)] = k;
    }
    std::stable_sort(order2.begin(), order2.end(), [&](LO a, LO b) {
      return in_meta[static_cast<size_t>(a) * 2] <
            in_meta[static_cast<size_t>(b) * 2];
    });
    std::vector<T> relay_values(static_cast<size_t>(num_answers) *
                               static_cast<size_t>(num_components_));
    std::vector<LO> relay_index(static_cast<size_t>(num_answers));
    std::vector<int8_t> relay_resolved(static_cast<size_t>(num_answers));
    redev::LOs relay_dest;
    std::vector<LO> val_counts;
    std::vector<LO> idx_counts;
    prev_dest = -1;
    for (LO out_k = 0; out_k < num_answers; ++out_k) {
      const LO in_k = order2[static_cast<size_t>(out_k)];
      for (int c = 0; c < num_components_; ++c) {
        relay_values[static_cast<size_t>(out_k) *
                    static_cast<size_t>(num_components_) +
                    static_cast<size_t>(c)] =
          in_values[static_cast<size_t>(in_k) *
                   static_cast<size_t>(num_components_) +
                   static_cast<size_t>(c)];
      }
      relay_index[static_cast<size_t>(out_k)] =
        in_meta[static_cast<size_t>(in_k) * 2 + 1];
      relay_resolved[static_cast<size_t>(out_k)] =
        in_resolved[static_cast<size_t>(in_k)];
      const int asker = in_meta[static_cast<size_t>(in_k) * 2];
      if (asker != prev_dest) {
        relay_dest.push_back(asker);
        val_counts.push_back(0);
        idx_counts.push_back(0);
        prev_dest = asker;
      }
      val_counts.back() += num_components_;
      idx_counts.back() += 1;
    }
    // See the note in the redirect phase: the offsets must be exactly
    // dest.size() + 1 long, with one entry per answer per destination.
    redev::LOs relay_val_offsets = MakeCsrOffsets(val_counts);
    redev::LOs relay_idx_offsets = MakeCsrOffsets(idx_counts);
    redev::LOs relay_res_offsets = relay_idx_offsets;

    // Phase 4 (relay): send final answers back to the original asking ranks.
    answer_values_comm_.SetOutMessageLayout(relay_dest, relay_val_offsets);
    answer_meta_comm_.SetOutMessageLayout(relay_dest, relay_idx_offsets);
    answer_resolved_comm_.SetOutMessageLayout(relay_dest, relay_res_offsets);
    channel_scope_.BeginSendCommunicationPhase(*this);
    answer_values_comm_.Send(relay_values.data(), redev::Mode::Synchronous);
    answer_meta_comm_.Send(relay_index.data(), redev::Mode::Synchronous);
    answer_resolved_comm_.Send(relay_resolved.data(), redev::Mode::Synchronous);
    channel_scope_.EndSendCommunicationPhase(*this);
  }

private:
  /**
   * @brief Aborts with a clear message if the two sides' declared round counts
   * disagree.
   *
   * The single setup exchange carries each side's count across. A mismatch
   * would eventually hang (the server blocked in BeginStep waiting for a round
   * that never comes, or a client blocked on a server that has already exited),
   * so it aborts now instead. -1 on either side means "not declared" and skips
   * the check. Called by both directions within PublishOwnedBounds /
   * ReceiveRoutingPartition, on the rank that observes both numbers; the
   * message is printed only there.
   */
  void CheckDeclaredRoundCount(int mine, int theirs)
  {
    if (mine < 0 || theirs < 0 || mine == theirs) {
      return;
    }
    int rank = 0;
    MPI_Comm_rank(mpi_comm_, &rank);
    const char* side =
      redev_.GetProcessType() == redev::ProcessType::Server ? "server"
                                                           : "client";
    std::fprintf(
      stderr,
      "\n[pcms DistributedEvaluationChannel '%s'] declared round counts "
      "do not match: this %s will run %d round(s), the other side "
      "declared %d.\n"
      "The round count is the total number of Evaluate()/batch calls "
      "serviced, and is independent of the number of field ranks (one "
      "round fans out across every rank). Refusing to continue, because "
      "a mismatch would hang or abort mid-run.\n",
      name_.c_str(), side, mine, theirs);
    std::fflush(stderr);
    MPI_Abort(mpi_comm_, EXIT_FAILURE);
  }

  /**
   * @brief Builds CSR offsets [0, c0, c0+c1, ...] from per-destination item
   * counts, so offsets.size() == counts.size() + 1 as SetOutMessageLayout
   * requires.
   */
  static redev::LOs MakeCsrOffsets(const std::vector<LO>& counts)
  {
    redev::LOs offsets(counts.size() + 1, 0);
    for (std::size_t i = 0; i < counts.size(); ++i) {
      offsets[i + 1] = offsets[i] + counts[i];
    }
    return offsets;
  }

  /**
   * @brief Forwards the channel's Begin/EndCommunicationPhase calls.
   *
   * Kept as a nested helper only so the many phase transitions above read as
   * symmetric pairs.
   */
  struct ChannelScope
  {
    void BeginSendCommunicationPhase(DistributedEvaluationChannel& self)
    {
      self.channel_.BeginSendCommunicationPhase();
    }
    void EndSendCommunicationPhase(DistributedEvaluationChannel& self)
    {
      self.channel_.EndSendCommunicationPhase();
    }
    void BeginReceiveCommunicationPhase(DistributedEvaluationChannel& self)
    {
      self.channel_.BeginReceiveCommunicationPhase();
    }
    void EndReceiveCommunicationPhase(DistributedEvaluationChannel& self)
    {
      self.channel_.EndReceiveCommunicationPhase();
    }
  };

  redev::Redev& redev_;
  redev::Channel& channel_;
  MPI_Comm mpi_comm_;
  std::string name_;
  int dim_;
  int num_components_;
  int expected_rounds_ = -1;
  bool setup_exchanged_ = false;
  ChannelScope channel_scope_;
  redev::BidirectionalComm<Real> query_coords_comm_;
  redev::BidirectionalComm<LO> query_meta_comm_;
  redev::BidirectionalComm<T> answer_values_comm_;
  redev::BidirectionalComm<LO> answer_meta_comm_;
  redev::BidirectionalComm<int8_t> answer_resolved_comm_;
  redev::BidirectionalComm<Real> partition_coords_comm_;
  redev::BidirectionalComm<LO> partition_flags_comm_;
};

} // namespace pcms

#endif // PCMS_DISTRIBUTED_EVALUATION_CHANNEL_H
