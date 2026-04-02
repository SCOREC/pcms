#ifndef PCMS_PRINT_H
#define PCMS_PRINT_H
#include "pcms/configuration.h"
#include <cstdio>

#include <Kokkos_Core.hpp>

namespace pcms
{

#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__) ||               \
  defined(__SYCL_DEVICE_ONLY__)
#define ACTIVE_GPU_EXECUTION
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
#if defined(PCMS_ENABLE_PRINT) && !defined(ACTIVE_GPU_EXECUTION)
  fprintf(getStdout(), fmt, args...);
#endif
}

template <typename... Args>
KOKKOS_INLINE_FUNCTION void printDebugInfo(const char* fmt, const Args&... args)
{
#if !defined(NDEBUG) && defined(PCMS_ENABLE_PRINT)
#if !defined(ACTIVE_GPU_EXECUTION)
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

#endif // PCMS_PRINT_H
