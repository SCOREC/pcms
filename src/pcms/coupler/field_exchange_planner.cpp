#include "pcms/coupler/field_exchange_planner.h"
#include "partition.h"
#include "pcms/field/gid_permutation.hpp"
#include "pcms/utility/assert.h"
#include "pcms/utility/inclusive_scan.h"
#include "pcms/utility/profile.h"
#include <algorithm>
#include <map>

namespace pcms
{

namespace
{

struct PartitionMapping
{
  std::vector<LO> indices;
  EntOffsetsArray ent_offsets;

  PartitionMapping() { ent_offsets.fill(0); }
};

using ReversePartitionMap2 = std::map<pcms::LO, PartitionMapping>;

struct OutMsg
{
  redev::LOs dest;
  redev::LOs offset;
};

static size_t GetMessageBlockIndex(
  LO permutation_entry, Rank1View<const redev::LO, HostMemorySpace> offsets)
{
  for (size_t i = 0; i < offsets.size() - 1; ++i) {
    if (permutation_entry >= offsets(i) && permutation_entry < offsets(i + 1)) {
      return i;
    }
  }
  PCMS_ALWAYS_ASSERT(false);
  return 0;
}

static OutMsg ConstructOutMessage(const ReversePartitionMap2& reverse_partition)
{
  PCMS_FUNCTION_TIMER;
  OutMsg out;
  redev::LOs counts;
  counts.reserve(reverse_partition.size());
  out.dest.reserve(reverse_partition.size());
  for (const auto& rank : reverse_partition) {
    out.dest.push_back(rank.first);
    counts.push_back(rank.second.indices.size());
  }
  out.offset.resize(counts.size() + 1);
  out.offset[0] = 0;
  pcms::inclusive_scan(counts.begin(), counts.end(),
                       std::next(out.offset.begin(), 1));
  return out;
}

static redev::LOs ConstructPermutation(
  const ReversePartitionMap2& reverse_partition, size_t num_entries,
  int& length)
{
  PCMS_FUNCTION_TIMER;
  // Holders that do not participate in this exchange (non-owned, or owned but
  // outside the overlap region) stay unmapped; serialization skips them.
  redev::LOs permutation = MakeUnmappedPermutation(num_entries);
  LO entry = 0;
  for (const auto& rank : reverse_partition) {
    for (unsigned e = 0; e < rank.second.ent_offsets.size() - 1; ++e) {
      const int start = rank.second.ent_offsets[e];
      const int end = rank.second.ent_offsets[e + 1];

      for (int i = start; i < end; ++i) {
        auto index = rank.second.indices[i];
        PCMS_ALWAYS_ASSERT(static_cast<size_t>(index) < permutation.size());
        permutation[index] = entry++;
      }
    }
  }
  length = entry;
  return permutation;
}

// Parses a received GID message -- a concatenation of per-sender
// [entity-offset header | GIDs] blocks -- into the permutation mapping each
// local holder to the compact (header-excluded) buffer index of its GID. This
// owns the on-the-wire message format; the GID-to-holder matching is delegated
// to pcms::AppendGidPermutation. Holders absent from the message are left
// unmapped. Returns the total payload length (GID count excluding headers) via
// payload_length.
static redev::LOs ConstructPermutation(
  GlobalIDView<HostMemorySpace> local_gids,
  GlobalIDView<HostMemorySpace> received_msg,
  const EntOffsetsArray& ent_offsets, int& payload_length)
{
  PCMS_FUNCTION_TIMER;
  GidToIndexMaps gid_to_buffer_index;
  size_t offset = 0;
  LO payload_offset = 0;
  while (true) {
    // Guard the header read before dereferencing its last entry (the payload
    // length), so a truncated message aborts rather than reading out of bounds.
    PCMS_ALWAYS_ASSERT(offset + ent_offsets_len <= received_msg.size());
    int message_length = received_msg(offset + ent_offsets_len - 1);

    PCMS_ALWAYS_ASSERT(offset + ent_offsets_len + message_length - 1 <
                       received_msg.size());

    for (size_t e = 0; e < ent_offsets_len - 1; ++e) {
      size_t start = received_msg(offset + e);
      size_t end = received_msg(offset + e + 1);

      for (size_t i = start; i < end; ++i) {
        const GO gid = received_msg(offset + ent_offsets_len + i);
        gid_to_buffer_index[e][gid] = payload_offset + i;
      }
    }

    payload_offset += message_length;
    offset += message_length + ent_offsets_len;
    if (offset >= received_msg.size())
      break;
  }

  redev::LOs permutation;
  permutation.reserve(local_gids.size());
  AppendGidPermutation(local_gids, ent_offsets, gid_to_buffer_index,
                       permutation,
                       /*allow_missing=*/true);

  REDEV_ALWAYS_ASSERT(permutation.size() == local_gids.size());
  payload_length = payload_offset;
  return permutation;
}

static OutMsg ConstructOutMessage(int rank, int nproc,
                                  const redev::InMessageLayout& in)
{
  PCMS_FUNCTION_TIMER;
  REDEV_ALWAYS_ASSERT(!in.srcRanks.empty());
  auto nAppProcs = in.srcRanks.size() / static_cast<size_t>(nproc);
  redev::LOs senderDeg(nAppProcs);
  for (size_t i = 0; i < nAppProcs - 1; ++i) {
    senderDeg[i] =
      in.srcRanks[(i + 1) * nproc + rank] - in.srcRanks[i * nproc + rank];
  }
  const auto totInMsgs = in.offset[rank + 1] - in.offset[rank];
  senderDeg[nAppProcs - 1] =
    totInMsgs - in.srcRanks[(nAppProcs - 1) * nproc + rank];
  OutMsg out;
  for (size_t i = 0; i < nAppProcs; ++i) {
    if (senderDeg[i] > 0) {
      REDEV_ALWAYS_ASSERT(senderDeg[i] > ent_offsets_len);
      senderDeg[i] -= ent_offsets_len;
      out.dest.push_back(i);
    }
  }
  redev::GO sum = 0;
  for (auto deg : senderDeg) {
    if (deg > 0) {
      out.offset.push_back(sum);
      sum += deg;
    }
  }
  out.offset.push_back(sum);
  return out;
}

static ReversePartitionMap2 BuildReversePartitionMap(
  const FieldLayout& layout, const redev::Partition& partition,
  const OverlapMask& overlap_mask)
{
  PCMS_FUNCTION_TIMER;
  auto owned = layout.GetOwnedHost();
  auto owned_to_local = layout.GetOwnedToLocalHost();
  auto class_dims = layout.GetDOFHolderClassificationDimensionsHost();
  auto class_ids = layout.GetDOFHolderClassificationIdsHost();
  auto coords = layout.GetOwnedDOFHolderCoordinates().GetValues();
  auto ent_offsets = layout.GetEntOffsets();
  int mesh_dim = static_cast<int>(coords.extent(1));

  // move coords to host
  Kokkos::View<Real**, DeviceMemorySpace> coords_device(
    "coords_device", coords.extent(0), coords.extent(1));
  Kokkos::parallel_for(
    "copy_coords", Kokkos::RangePolicy<>(0, coords.extent(0)),
    KOKKOS_LAMBDA(int i) {
      for (unsigned d = 0; d < coords.extent(1); ++d) {
        coords_device(i, d) = coords(i, d);
      }
    });
  auto coords_host =
    Kokkos::create_mirror_view_and_copy(HostMemorySpace(), coords_device);

  ReversePartitionMap2 reverse_partition;
  const LO n_owned = static_cast<LO>(coords.extent(0));
  std::array<Real, 3> coord{};
  auto overlap_mask_view = overlap_mask.GetMask(layout);

  for (LO owned_index = 0; owned_index < n_owned; ++owned_index) {
    // Classification dims/ids, entity offsets, and the ownership mask are still
    // local-indexed, so map the owned index back to its local index.
    const LO local_index =
      owned_to_local.size() == 0 ? owned_index : owned_to_local(owned_index);

    // Skip holders not owned by this rank. This is required for xgc
    if (!owned(local_index))
      continue;

    if (!overlap_mask_view[owned_index])
      continue;

    for (int d = 0; d < mesh_dim; ++d)
      coord[d] = coords_host(owned_index, d);
    for (int d = mesh_dim; d < 3; ++d)
      coord[d] = 0.0;

    int mesh_ent_dim = GetMeshEntityDim(local_index, ent_offsets);
    LO class_dim = class_dims[local_index];
    LO class_id = class_ids[local_index];

    auto dr = std::visit(GetRank{class_id, class_dim, coord}, partition);
    reverse_partition[dr].indices.emplace_back(owned_index);

    for (size_t e = static_cast<size_t>(mesh_ent_dim) + 1; e < ent_offsets_len;
         ++e) {
      reverse_partition[dr].ent_offsets[e] += 1;
    }
  }
  return reverse_partition;
}

} // namespace

ExchangePlan GenericFieldExchangePlanner::BuildExchangePlan(
  const FieldLayout& layout, const redev::Partition& partition,
  const OverlapMask* overlap_mask) const
{
  PCMS_FUNCTION_TIMER;
  PCMS_ALWAYS_ASSERT(overlap_mask != nullptr);
  auto gids = layout.GetOwnedGidsHost();

  const ReversePartitionMap2 reverse_partition =
    BuildReversePartitionMap(layout, partition, *overlap_mask);

  ExchangePlan plan;

  auto out_msg = ConstructOutMessage(reverse_partition);
  plan.dest_ranks = std::move(out_msg.dest);
  plan.offsets = std::move(out_msg.offset);

  int length = 0;
  plan.permutation =
    ConstructPermutation(reverse_partition, gids.size(), length);
  plan.msg_size = static_cast<size_t>(length);

  return plan;
}

ExchangePlan GenericFieldExchangePlanner::BuildReceivePlan(
  const FieldLayout& layout, GlobalIDView<HostMemorySpace> received_gids,
  int rank, int nproc, const redev::InMessageLayout& in_message_layout) const
{
  PCMS_FUNCTION_TIMER;
  auto gids = layout.GetOwnedGidsHost();
  auto ent_offsets = layout.GetOwnedEntOffsets();

  ExchangePlan plan;
  auto out_msg = ConstructOutMessage(rank, nproc, in_message_layout);
  plan.dest_ranks = std::move(out_msg.dest);
  plan.offsets = std::move(out_msg.offset);
  int length = 0;
  plan.permutation =
    ConstructPermutation(gids, received_gids, ent_offsets, length);
  plan.msg_size = static_cast<size_t>(length);
  return plan;
}

void GenericFieldExchangePlanner::FillGidMessage(
  const FieldLayout& layout, const ExchangePlan& plan,
  Rank1View<GO, HostMemorySpace> gid_message) const
{
  PCMS_FUNCTION_TIMER;
  const size_t header_size = plan.dest_ranks.size() * ent_offsets_len;
  PCMS_ALWAYS_ASSERT(static_cast<size_t>(gid_message.size()) ==
                     plan.msg_size + header_size);

  auto gids = layout.GetOwnedGidsHost();
  auto owned_to_local = layout.GetOwnedToLocalHost();
  auto ent_offsets = layout.GetEntOffsets();
  auto offsets = Rank1View<const redev::LO, HostMemorySpace>(
    plan.offsets.data(), plan.offsets.size());

  std::vector<EntOffsetsArray> per_rank_offsets(plan.dest_ranks.size());

  for (LO owned_index = 0; owned_index < static_cast<LO>(gids.size());
       ++owned_index) {
    LO perm_index = plan.permutation[owned_index];
    // Owned holders outside the overlap region carry the sentinel and have no
    // slot in the message.
    if (perm_index < 0)
      continue;
    auto block_index = GetMessageBlockIndex(perm_index, offsets);
    const auto gid_index =
      perm_index + static_cast<LO>((block_index + 1) * ent_offsets_len);
    gid_message(gid_index) = gids(owned_index);

    const LO local_index =
      owned_to_local.size() == 0 ? owned_index : owned_to_local(owned_index);
    int mesh_ent_dim = GetMeshEntityDim(local_index, ent_offsets);
    for (size_t e = static_cast<size_t>(mesh_ent_dim) + 1; e < ent_offsets_len;
         ++e) {
      per_rank_offsets[block_index][e] += 1;
    }
  }

  for (size_t block_index = 0; block_index < per_rank_offsets.size();
       ++block_index) {
    auto header_offset = static_cast<size_t>(plan.offsets[block_index]) +
                         block_index * ent_offsets_len;
    for (size_t e = 0; e < ent_offsets_len; ++e) {
      gid_message(header_offset + e) =
        static_cast<GO>(per_rank_offsets[block_index][e]);
    }
  }
}

} // namespace pcms
