#include "pcms/field/field_layout.hpp"
#include "pcms/field/gid_permutation.hpp"
#include "pcms/utility/assert.hpp"

namespace pcms
{

void FieldLayout::BuildGlobalToLocalPermutation()
{
  const auto permutation =
    BuildSortedGidPermutation(GetGidsHost(), GetEntOffsets());
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
