#ifndef PCMS_UTILITY_CHECK_HPP
#define PCMS_UTILITY_CHECK_HPP

// Debug checks usable in device code, with the semantics of Omega_h's
// OMEGA_H_CHECK / OMEGA_H_CHECK_PRINTF in a Kokkos build: with NDEBUG the
// condition is evaluated and ignored; otherwise a failed check prints the
// message (PCMS_CHECK_PRINTF) and asserts.

#include <Kokkos_Core.hpp>

#include <cassert>

#ifdef NDEBUG
#define PCMS_CHECK(cond) static_cast<void>(cond)
#define PCMS_CHECK_PRINTF(cond, format, ...) static_cast<void>(cond)
#else
#define PCMS_CHECK(cond) assert(cond)
#define PCMS_CHECK_PRINTF(cond, format, ...)                                   \
  do {                                                                         \
    if (!(cond)) {                                                             \
      Kokkos::printf("ERROR: " format "\n", __VA_ARGS__);                      \
    }                                                                          \
    assert(cond);                                                              \
  } while (0)
#endif

#endif // PCMS_UTILITY_CHECK_HPP
