#include "field_layout_communicator.h"

namespace pcms
{

namespace field_layout_communicator
{
// reverse partition is a map that has the partition rank as a key
// and the values are an vector where each entry is the index into
// the array of data to send
OutMsg ConstructOutMessage(const ReversePartitionMap2& reverse_partition)
{
  PCMS_FUNCTION_TIMER;
  OutMsg out;
  redev::LOs counts;
  counts.reserve(reverse_partition.size());
  out.dest.clear();
  out.dest.reserve(reverse_partition.size());
  // number of entries for each rank
  for (auto& rank : reverse_partition) {
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

size_t count_entries(const ReversePartitionMap2& reverse_partition)
{
  PCMS_FUNCTION_TIMER;
  size_t num_entries = 0;
  for (const auto& v : reverse_partition) {
    num_entries += v.second.indices.size();
  }
  return num_entries;
}

// note this function can be parallelized by making use of the offsets
redev::LOs ConstructPermutation(const ReversePartitionMap2& reverse_partition,
                                size_t num_entries, int* length)
{
  PCMS_FUNCTION_TIMER;
  redev::LOs permutation(num_entries);
  LO entry = 0;
  for (auto& rank : reverse_partition) {
    entry += ent_offsets_len;

    for (int e = 0; e < rank.second.ent_offsets.size() - 1; ++e) {
      int start = rank.second.ent_offsets[e];
      int end = rank.second.ent_offsets[e + 1];

      for (int i = start; i < end; ++i) {
        LO index = rank.second.indices[i];
        PCMS_ALWAYS_ASSERT(index < permutation.size());
        permutation[index] = entry++;
      }
    }
  }
  *length = entry;
  return permutation;
}

/**
 *
 * @param local_gids local gids are the mesh GIDs in local mesh iteration order
 * @param received_gids received GIDs are the GIDS in the order of the incomming
 * message1
 * @return permutation array such that GIDS(Permutation[i]) = msgs
 */
redev::LOs ConstructPermutation(GlobalIDView<HostMemorySpace> local_gids,
                                GlobalIDView<HostMemorySpace> received_msg,
                                EntOffsetsArray ent_offsets)
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

    for (int e = 0; e < received_offsets.size() - 1; ++e) {
      size_t start = received_offsets[e];
      size_t end = received_offsets[e + 1];

      for (int i = start; i < end; ++i) {
        gid_to_buffer_index[e][received_gids[i]] = offset + ent_offsets_len + i;
      }
    }

    offset += length + ent_offsets_len;
    if (offset >= received_msg.size())
      break;
  }

  redev::LOs permutation;
  permutation.reserve(local_gids.size());
  for (int e = 0; e < ent_offsets.size() - 1; ++e) {
    size_t start = ent_offsets[e];
    size_t end = ent_offsets[e + 1];

    for (int i = start; i < end; ++i) {
      permutation.push_back(gid_to_buffer_index[e][local_gids[i]]);
    }
  }

  REDEV_ALWAYS_ASSERT(permutation.size() == local_gids.size());
  return permutation;
}

OutMsg ConstructOutMessage(int rank, int nproc,
                           const redev::InMessageLayout& in)
{
  PCMS_FUNCTION_TIMER;
  REDEV_ALWAYS_ASSERT(!in.srcRanks.empty());
  // auto nAppProcs =
  // Omega_h::divide_no_remainder(in.srcRanks.size(),static_cast<size_t>(nproc));
  auto nAppProcs = in.srcRanks.size() / static_cast<size_t>(nproc);
  // build dest and offsets arrays from incoming message metadata
  redev::LOs senderDeg(nAppProcs);
  for (size_t i = 0; i < nAppProcs - 1; i++) {
    senderDeg[i] =
      in.srcRanks[(i + 1) * nproc + rank] - in.srcRanks[i * nproc + rank];
  }
  const auto totInMsgs = in.offset[rank + 1] - in.offset[rank];
  senderDeg[nAppProcs - 1] =
    totInMsgs - in.srcRanks[(nAppProcs - 1) * nproc + rank];
  OutMsg out;
  for (size_t i = 0; i < nAppProcs; i++) {
    if (senderDeg[i] > 0) {
      out.dest.push_back(i);
    }
  }
  redev::GO sum = 0;
  for (auto deg : senderDeg) { // exscan over values > 0
    if (deg > 0) {
      out.offset.push_back(sum);
      sum += deg;
    }
  }
  out.offset.push_back(sum);
  return out;
}

} // namespace field_layout_communicator
} // namespace pcms