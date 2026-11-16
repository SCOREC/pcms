#ifndef PCMS_COUPLING_MEMORY_SPACES_H
#define PCMS_COUPLING_MEMORY_SPACES_H
#include <Kokkos_Core.hpp>

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

} // namespace pcms

#endif // PCMS_COUPLING_MEMORY_SPACES_H
