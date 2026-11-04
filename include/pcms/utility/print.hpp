#ifndef PCMS_UTILITY_PRINT_HPP
#define PCMS_UTILITY_PRINT_HPP
#include "pcms/configuration.hpp"
#include <cstdio>

#include <Kokkos_Core.hpp>

namespace pcms
{

#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__) ||               \
  defined(__SYCL_DEVICE_ONLY__)
#define PCMS_ACTIVE_GPU_EXECUTION
#endif

FILE* getStdout();
FILE* getStderr();

void setStdout(FILE* out);
void setStderr(FILE* err);

template <typename... Args>
void printError(const char* fmt, const Args&... args)
{
#if defined(PCMS_ENABLE_PRINT)
  fprintf(getStdout(), fmt, args...);
#endif
}

template <typename... Args>
KOKKOS_INLINE_FUNCTION void printInfo(const char* fmt, const Args&... args)
{
#if defined(PCMS_ENABLE_PRINT) && !defined(PCMS_ACTIVE_GPU_EXECUTION)
  fprintf(getStdout(), fmt, args...);
#endif
}

template <typename... Args>
KOKKOS_INLINE_FUNCTION void printDebugInfo(const char* fmt, const Args&... args)
{
#if !defined(NDEBUG) && defined(PCMS_ENABLE_PRINT)
#if !defined(PCMS_ACTIVE_GPU_EXECUTION)
  fprintf(getStdout(), fmt, args...);
#else // For GPU execution
  printf(fmt, args...);
#endif
#else
  (void)fmt;
  ((void)args, ...);
#endif
}

} // namespace pcms

#endif // PCMS_UTILITY_PRINT_HPP
