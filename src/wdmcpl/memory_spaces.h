#ifndef WDM_COUPLING_MEMORY_SPACES_H
#define WDM_COUPLING_MEMORY_SPACES_H
#include <Kokkos_Core.hpp>

namespace wdmcpl {
//struct CudaMemorySpace{};
//struct HostMemorySpace{};
//using CudaMemorySpace = Kokkos::Cuda;
using  HostMemorySpace = Kokkos::HostSpace;
using DefaultExecutionSpace = Kokkos::DefaultExecutionSpace;

}

#endif // WDM_COUPLING_MEMORY_SPACES_H
