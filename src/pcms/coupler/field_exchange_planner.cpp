#include "pcms/coupler/field_exchange_planner.h"
#include "partition.h"
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

// Returns the mesh entity dimension for DOF local_index based on the
// entity offsets array.  ent_offsets[d]..ent_offsets[d+1] is the range of
// DOF indices belonging to mesh entity dimension d.
static int GetMeshEntityDim(LO local_index, const EntOffsetsArray& ent_offsets)
{
  for (int d = 0; d < ent_offsets_len - 1; ++d) {
    if (local_index >= static_cast<LO>(ent_offsets[d]) &&
        local_index < static_cast<LO>(ent_offsets[d + 1])) {
      return d;
    }
  }
  return ent_offsets_len - 1;
}

static size_t GetMessageBlockIndex(
  LO permutation_entry, Rank1View<const redev::LO, HostMemorySpace> offsets)
{
  auto begin = offsets.data_handle();
  auto end = begin + offsets.size();
  auto it = std::upper_bound(begin, end, permutation_entry);
  PCMS_ALWAYS_ASSERT(it != begin);
  return static_cast<size_t>(std::distance(begin, it - 1));
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
    counts.push_back(rank.second.indices.size() +
                     rank.second.ent_offsets.size());
  }
  out.offset.resize(counts.size() + 1);
  out.offset[0] = 0;
  pcms::inclusive_scan(counts.begin(), counts.end(),
                       std::next(out.offset.begin(), 1));
  return out;
}

static redev::LOs ConstructPermutation(
  const ReversePartitionMap2& reverse_partition, size_t num_entries,
  int* length)
{
  PCMS_FUNCTION_TIMER;
  redev::LOs permutation(num_entries);
  LO entry = 0;
  for (const auto& rank : reverse_partition) {
    entry += ent_offsets_len;

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
  *length = entry;
  return permutation;
}

static redev::LOs ConstructPermutation(
  GlobalIDView<HostMemorySpace> local_gids,
  GlobalIDView<HostMemorySpace> received_msg,
  const EntOffsetsArray& ent_offsets)
{
  PCMS_FUNCTION_TIMER;
  std::array<std::map<pcms::GO, pcms::LO>, 4> gid_to_buffer_index;
  size_t offset = 0;
  while (true) {
    GlobalIDView<HostMemorySpace> received_offsets(
      received_msg.data_handle() + offset, ent_offsets_len);
    int length = received_offsets[received_offsets.size() - 1];
    GlobalIDView<HostMemorySpace> received_gids(
      received_msg.data_handle() + offset + ent_offsets_len, length);

    PCMS_ALWAYS_ASSERT(offset + ent_offsets_len + length - 1 <
                       received_msg.size());

    for (size_t e = 0; e < received_offsets.size() - 1; ++e) {
      size_t start = received_offsets[e];
      size_t end = received_offsets[e + 1];

      for (size_t i = start; i < end; ++i) {
        gid_to_buffer_index[e][received_gids[i]] = offset + ent_offsets_len + i;
      }
    }

    offset += length + ent_offsets_len;
    if (offset >= received_msg.size())
      break;
  }

  redev::LOs permutation;
  permutation.reserve(local_gids.size());
  for (size_t e = 0; e < ent_offsets.size() - 1; ++e) {
    const auto start = ent_offsets[e];
    const auto end = ent_offsets[e + 1];

    for (size_t i = start; i < end; ++i)
      permutation.push_back(gid_to_buffer_index[e][local_gids[i]]);
  }

  REDEV_ALWAYS_ASSERT(permutation.size() == local_gids.size());
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
    if (senderDeg[i] > 0)
      out.dest.push_back(i);
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
  auto class_dims = layout.GetDOFHolderClassificationDimensionsHost();
  auto class_ids = layout.GetDOFHolderClassificationIdsHost();
  auto coords = layout.GetDOFHolderCoordinates().GetValues();
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
  LO n = static_cast<LO>(owned.extent(0));
  std::array<Real, 3> coord{};
  auto overlap_mask_view = overlap_mask.GetMask(layout);

  for (LO local_index = 0; local_index < n; ++local_index) {
    if (!owned[local_index])
      continue;

    if (!overlap_mask_view[local_index])
      continue;

    for (int d = 0; d < mesh_dim; ++d)
      coord[d] = coords_host(local_index, d);
    for (int d = mesh_dim; d < 3; ++d)
      coord[d] = 0.0;

    int mesh_ent_dim = GetMeshEntityDim(local_index, ent_offsets);
    LO class_dim = class_dims[local_index];
    LO class_id = class_ids[local_index];

    auto dr = std::visit(GetRank{class_id, class_dim, coord}, partition);
    reverse_partition[dr].indices.emplace_back(local_index);

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
  auto gids = layout.GetGidsHost();

  const ReversePartitionMap2 reverse_partition =
    BuildReversePartitionMap(layout, partition, *overlap_mask);

  ExchangePlan plan;

  auto out_msg = ConstructOutMessage(reverse_partition);
  plan.dest_ranks = std::move(out_msg.dest);
  plan.offsets = std::move(out_msg.offset);

  int length = 0;
  plan.permutation =
    ConstructPermutation(reverse_partition, gids.size(), &length);
  plan.msg_size = static_cast<size_t>(length);

  return plan;
}

ExchangePlan GenericFieldExchangePlanner::BuildReceivePlan(
  const FieldLayout& layout, GlobalIDView<HostMemorySpace> received_gids,
  int rank, int nproc, const redev::InMessageLayout& in_message_layout) const
{
  PCMS_FUNCTION_TIMER;
  auto gids = layout.GetGidsHost();
  auto ent_offsets = layout.GetEntOffsets();

  ExchangePlan plan;
  auto out_msg = ConstructOutMessage(rank, nproc, in_message_layout);
  plan.dest_ranks = std::move(out_msg.dest);
  plan.offsets = std::move(out_msg.offset);
  plan.permutation = ConstructPermutation(gids, received_gids, ent_offsets);
  plan.msg_size = received_gids.size();
  return plan;
}

void GenericFieldExchangePlanner::FillGidMessage(
  const FieldLayout& layout, const ExchangePlan& plan,
  Rank1View<GO, HostMemorySpace> gid_message) const
{
  PCMS_FUNCTION_TIMER;
  PCMS_ALWAYS_ASSERT(static_cast<size_t>(gid_message.size()) == plan.msg_size);

  auto gids = layout.GetGidsHost();
  auto owned = layout.GetOwnedHost();
  auto ent_offsets = layout.GetEntOffsets();
  auto offsets = Rank1View<const redev::LO, HostMemorySpace>(
    plan.offsets.data(), plan.offsets.size());

  std::vector<EntOffsetsArray> per_rank_offsets(plan.dest_ranks.size());

  for (LO local_index = 0; local_index < static_cast<LO>(gids.size());
       ++local_index) {
    if (!owned[local_index])
      continue;

    LO perm_index = plan.permutation[local_index];
    gid_message[perm_index] = gids[local_index];

    auto block_index = GetMessageBlockIndex(perm_index, offsets);
    int mesh_ent_dim = GetMeshEntityDim(local_index, ent_offsets);
    for (size_t e = static_cast<size_t>(mesh_ent_dim) + 1; e < ent_offsets_len;
         ++e) {
      per_rank_offsets[block_index][e] += 1;
    }
  }

  for (size_t block_index = 0; block_index < per_rank_offsets.size();
       ++block_index) {
    auto header_offset = static_cast<size_t>(plan.offsets[block_index]);
    for (size_t e = 0; e < ent_offsets_len; ++e) {
      gid_message[header_offset + e] =
        static_cast<GO>(per_rank_offsets[block_index][e]);
    }
  }
}

} // namespace pcms
