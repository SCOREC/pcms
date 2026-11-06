#include "pcms/field/field_layout.h"
#include "pcms/field/gid_permutation.hpp"
#include "pcms/utility/assert.h"

namespace pcms
{

void FieldLayout::BuildGlobalToLocalPermutation(
  Kokkos::View<const GO*, DeviceMemorySpace> local_gids,
  const EntOffsetsArray& local_ent_offsets)
{
  // The sort runs on the host; the GIDs are copied there only for it.
  const auto gids_host =
    Kokkos::create_mirror_view_and_copy(HostMemorySpace(), local_gids);
  const auto permutation = BuildSortedGidPermutation(
    GlobalIDView<HostMemorySpace>(gids_host.data(), gids_host.size()),
    local_ent_offsets);
  global_to_local_host_ = Kokkos::View<LO*, HostMemorySpace>(
    "global_to_local_host", permutation.size());
  for (size_t i = 0; i < permutation.size(); ++i)
    global_to_local_host_(i) = permutation[i];

  global_to_local_ =
    Kokkos::View<LO*, DeviceMemorySpace>("global_to_local", permutation.size());
  Kokkos::deep_copy(global_to_local_, global_to_local_host_);
}

Kokkos::View<const LO*, HostMemorySpace>
FieldLayout::GetGlobalToLocalPermutationHost() const
{
  // A non-empty view means BuildGlobalToLocalPermutation() ran in the derived
  // layout's constructor; an empty one means the layout never built it.
  PCMS_ALWAYS_ASSERT(global_to_local_host_.size() > 0);
  return global_to_local_host_;
}

Kokkos::View<const LO*, DeviceMemorySpace>
FieldLayout::GetGlobalToLocalPermutation() const
{
  PCMS_ALWAYS_ASSERT(global_to_local_.size() > 0);
  return global_to_local_;
}

} // namespace pcms
