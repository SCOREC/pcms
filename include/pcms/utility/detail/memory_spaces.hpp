#ifndef PCMS_UTILITY_DETAIL_MEMORY_SPACES_HPP
#define PCMS_UTILITY_DETAIL_MEMORY_SPACES_HPP

#include "pcms/utility/memory_spaces.hpp"

namespace pcms::detail
{

/// Kernels over at most this many items run as plain loops when the device
/// memory is host memory: opening a parallel region (about a microsecond with
/// OpenMP) costs more than such a loop, and dominates a bind of a few points.
inline constexpr std::size_t SmallHostBatchSize = 256;

/// @return whether a kernel over `n` items on DeviceMemorySpace runs as a
///         plain loop on the calling thread: the memory is host-accessible
///         and `n` is at most SmallHostBatchSize
[[nodiscard]] constexpr bool IsSmallHostBatch(std::size_t n) noexcept
{
  return Kokkos::SpaceAccessibility<Kokkos::HostSpace,
                                    DeviceMemorySpace>::accessible &&
         n <= detail::SmallHostBatchSize;
}

/// Runs `functor(i)` for every i in [0, n): as a plain loop on the calling
/// thread for a small host batch (IsSmallHostBatch), else as a kernel on
/// DeviceMemorySpace's execution space. For the kernels of a bind or an
/// evaluation, whose work is proportional to the number of query points.
template <typename Index, typename Functor>
void ForEach(const char* label, Index n, const Functor& functor)
{
  if (detail::IsSmallHostBatch(static_cast<std::size_t>(n))) {
    for (Index i = 0; i < n; ++i) {
      functor(i);
    }
    return;
  }
  Kokkos::parallel_for(
    label, Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, n),
    functor);
}

} // namespace pcms::detail

#endif // PCMS_UTILITY_DETAIL_MEMORY_SPACES_HPP
