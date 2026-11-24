#ifndef PCMS_UTILITY_DETAIL_RANDOM_HPP
#define PCMS_UTILITY_DETAIL_RANDOM_HPP

#include "pcms/utility/memory_spaces.hpp"
#include "pcms/utility/types.hpp"
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <cstdint>

namespace pcms::detail
{

/// @brief A generator whose stream depends only on `seed` and the entity's
///        global id: the same samples for an entity on any number of ranks,
///        partition, thread count or device. The state is splitmix64 of
///        (seed, gid + 1), which decorrelates neighbouring ids; never zero
///        (XorShift's absorbing state).
KOKKOS_INLINE_FUNCTION Kokkos::Random_XorShift64<DefaultExecutionSpace>
EntityGenerator(std::uint64_t seed, GO gid)
{
  std::uint64_t z = seed + 0x9E3779B97F4A7C15ULL * (gid + 1);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  z ^= z >> 31;
  return Kokkos::Random_XorShift64<DefaultExecutionSpace>(z != 0 ? z : 1);
}

} // namespace pcms::detail

#endif // PCMS_UTILITY_DETAIL_RANDOM_HPP
