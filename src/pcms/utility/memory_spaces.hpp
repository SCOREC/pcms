#ifndef PCMS_COUPLING_MEMORY_SPACES_H
#define PCMS_COUPLING_MEMORY_SPACES_H
#include <Kokkos_Core.hpp>
#include <cstddef>

namespace pcms
{
// struct CudaMemorySpace{};
// struct HostMemorySpace{};
// using CudaMemorySpace = Kokkos::Cuda;
using HostMemorySpace = Kokkos::HostSpace;
using DefaultExecutionSpace = Kokkos::DefaultExecutionSpace;

// DeviceMemorySpace is the memory space of the default GPU execution space when
// Kokkos is built with a GPU backend. When no GPU backend is present it aliases
// HostMemorySpace so that generic code using DeviceMemorySpace compiles without
// change. PCMS_HAS_DISTINCT_DEVICE_MEMORY_SPACE is defined only when the device
// memory space actually differs from the host memory space; use it to guard
// device-only overloads.
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP) ||               \
  defined(KOKKOS_ENABLE_SYCL)
#define PCMS_HAS_DISTINCT_DEVICE_MEMORY_SPACE
using DeviceMemorySpace = Kokkos::DefaultExecutionSpace::memory_space;
#else
using DeviceMemorySpace = HostMemorySpace;
#endif

namespace detail
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
         n <= SmallHostBatchSize;
}

/// Runs `functor(i)` for every i in [0, n): as a plain loop on the calling
/// thread for a small host batch (IsSmallHostBatch), else as a kernel on
/// DeviceMemorySpace's execution space. For the kernels of a bind or an
/// evaluation, whose work is proportional to the number of query points.
template <typename Index, typename Functor>
void ForEach(const char* label, Index n, const Functor& functor)
{
  if (IsSmallHostBatch(static_cast<std::size_t>(n))) {
    for (Index i = 0; i < n; ++i) {
      functor(i);
    }
    return;
  }
  Kokkos::parallel_for(
    label, Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, n),
    functor);
}

} // namespace detail

} // namespace pcms

#endif // PCMS_COUPLING_MEMORY_SPACES_H
