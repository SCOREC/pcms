#include "pcms/field/gid_permutation.hpp"
#include "pcms/utility/assert.h"
#include <algorithm>

namespace pcms
{

namespace
{

// Sentinel permutation entry for a local holder that does not participate in an
// exchange (a non-owned/ghost holder, or an owned holder outside the coupled
// overlap region whose GID is never transmitted). It is deliberately negative
// so serialization/deserialization can skip it, and any accidental use as a
// buffer index is an obvious out-of-bounds access rather than a silent read or
// write of holder 0.
constexpr LO UnmappedPermutationIndex = -1;

} // namespace

std::vector<LO> MakeUnmappedPermutation(std::size_t size)
{
  return std::vector<LO>(size, UnmappedPermutationIndex);
}

void AppendGidPermutation(GlobalIDView<HostMemorySpace> local_gids,
                          const EntOffsetsArray& ent_offsets,
                          const GidToIndexMaps& gid_to_index,
                          std::vector<LO>& permutation, bool allow_missing)
{
  for (size_t e = 0; e < ent_offsets.size() - 1; ++e) {
    for (size_t i = ent_offsets[e]; i < ent_offsets[e + 1]; ++i) {
      const auto entry = gid_to_index[e].find(local_gids(i));
      if (entry != gid_to_index[e].end()) {
        permutation.push_back(entry->second);
        continue;
      }
      PCMS_ALWAYS_ASSERT(allow_missing);
      permutation.push_back(UnmappedPermutationIndex);
    }
  }
}

std::vector<LO> BuildSortedGidPermutation(GlobalIDView<HostMemorySpace> gids,
                                          const EntOffsetsArray& ent_offsets)
{
  GidToIndexMaps gid_to_index;
  LO next_index = 0;
  for (size_t e = 0; e < ent_offsets.size() - 1; ++e) {
    std::vector<GO> sorted_gids;
    sorted_gids.reserve(ent_offsets[e + 1] - ent_offsets[e]);
    for (size_t i = ent_offsets[e]; i < ent_offsets[e + 1]; ++i) {
      sorted_gids.push_back(gids(i));
    }
    std::sort(sorted_gids.begin(), sorted_gids.end());
    PCMS_ALWAYS_ASSERT(
      std::adjacent_find(sorted_gids.begin(), sorted_gids.end()) ==
      sorted_gids.end());
    for (const GO gid : sorted_gids) {
      gid_to_index[e].emplace(gid, next_index++);
    }
  }
  std::vector<LO> permutation;
  permutation.reserve(gids.size());
  AppendGidPermutation(gids, ent_offsets, gid_to_index, permutation,
                       /*allow_missing=*/false);
  return permutation;
}

} // namespace pcms
