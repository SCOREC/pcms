#ifndef PCMS_COUPLING_MEMORY_SPACES_H
#define PCMS_COUPLING_MEMORY_SPACES_H
#include <Kokkos_Core.hpp>

namespace pcms {
//struct CudaMemorySpace{};
//struct HostMemorySpace{};
//using CudaMemorySpace = Kokkos::Cuda;
using  HostMemorySpace = Kokkos::HostSpace;
using DefaultExecutionSpace = Kokkos::DefaultExecutionSpace;

}

#endif // PCMS_COUPLING_MEMORY_SPACES_H
