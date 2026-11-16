#ifndef PCMS_FIELD_GID_PERMUTATION_HPP
#define PCMS_FIELD_GID_PERMUTATION_HPP

#include "pcms/field/field_layout.h"
#include <array>
#include <cstddef>
#include <map>
#include <vector>

namespace pcms
{

// Maps a GID to its position in an exchange buffer, grouped by entity-dimension
// block: index e holds the GIDs classified to mesh entity dimension e.
using GidToIndexMaps = std::array<std::map<GO, LO>, ent_offsets_len - 1>;

// Returns a permutation of `size` entries with every local holder initially
// unmapped. Unmapped entries hold a negative sentinel so that serialization
// skips them and any accidental use as a buffer index is an obvious
// out-of-bounds access rather than a silent read or write of holder 0. Callers
// overwrite the entries for holders that participate in the exchange.
std::vector<LO> MakeUnmappedPermutation(std::size_t size);

// Appends, in entity-offset (local holder) order, the position each local GID
// maps to according to gid_to_index. When allow_missing is true, a GID absent
// from gid_to_index is left unmapped (negative sentinel): a local layout may
// contain ghost holders and holders outside the overlap region, which are never
// exchanged and which (de)serialization skips. When allow_missing is false --
// e.g. the sorted permutation below, whose maps are built from these very GIDs
// -- a missing GID is a hard error.
void AppendGidPermutation(GlobalIDView<HostMemorySpace> local_gids,
                          const EntOffsetsArray& ent_offsets,
                          const GidToIndexMaps& gid_to_index,
                          std::vector<LO>& permutation,
                          bool allow_missing = false);

// Builds a permutation that maps each local holder to a contiguous index
// ordered by GID within each entity block. Every GID must be unique within its
// block.
std::vector<LO> BuildSortedGidPermutation(GlobalIDView<HostMemorySpace> gids,
                                          const EntOffsetsArray& ent_offsets);

} // namespace pcms

#endif // PCMS_FIELD_GID_PERMUTATION_HPP
