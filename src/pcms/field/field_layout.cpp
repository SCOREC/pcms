#include "pcms/field/field_layout.hpp"
#include "pcms/geometry/geometry.hpp"
#include "pcms/field/gid_permutation.hpp"
#include "pcms/utility/assert.hpp"

namespace pcms
{
using field::BuildSortedGidPermutation;
} // namespace pcms

namespace pcms
{

OwnerHints DOFHolderOwners(const FieldLayout& layout)
{
  const int entity_dim = layout.GetDOFHolderEntityDim();
  if (layout.GetGeometry() == nullptr || entity_dim < 0 || layout.HasGhosts()) {
    return {};
  }
  // without ghosts, DOF holder i is the i-th DOF holder of the layout's one
  // populated dimension, whose block starts at that dimension's offset
  const LO n = layout.GetNumDofHolder();
  const auto first = static_cast<LO>(layout.GetEntOffsets()[entity_dim]);
  Kokkos::View<LO*, DeviceMemorySpace> entities(
    Kokkos::view_alloc(Kokkos::WithoutInitializing, "dof_holder_owners"), n);
  Kokkos::parallel_for(
    "DOFHolderOwners",
    Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, n),
    KOKKOS_LAMBDA(const LO i) { entities(i) = i - first; });
  return OwnerHints{layout.GetGeometry(), entity_dim, entities};
}

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
