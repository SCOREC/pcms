#ifndef PCMS_UTILITY_ASSERT_HPP
#define PCMS_UTILITY_ASSERT_HPP
#include <stdexcept>
#include "pcms/utility/detail/assert.hpp"

// https://stackoverflow.com/questions/16683146/can-macros-be-overloaded-by-number-of-arguments
#define PCMS_CAT(A, B) A##B
#define PCMS_SELECT(NAME, NUM) PCMS_CAT(NAME##_, NUM)

#define PCMS_GET_COUNT(_1, _2, _3, _4, _5, _6 /* ad nauseam */, COUNT, ...)    \
  COUNT
#define PCMS_VA_SIZE(...) PCMS_GET_COUNT(__VA_ARGS__, 6, 5, 4, 3, 2, 1)

#define PCMS_VA_SELECT(NAME, ...)                                              \
  PCMS_SELECT(NAME, PCMS_VA_SIZE(__VA_ARGS__))(__VA_ARGS__)

#define PCMS_ALWAYS_ASSERT(...) PCMS_VA_SELECT(PCMS_ALWAYS_ASSERT, __VA_ARGS__)
// PCMS_ALWAYS_ASSERT(cond) reports the rank in MPI_COMM_WORLD when MPI is
// initialized. The forms that take a communicator, PCMS_ALWAYS_ASSERT(cond,
// comm) and PCMS_ALWAYS_ASSERT(cond, comm, msg), are for code that uses MPI:
// they expand to MPI_Comm_rank, so the caller includes <mpi.h> (this header
// does not).
#define PCMS_ALWAYS_ASSERT_1(pcms_macro_cond)                                  \
  /* NOLINTBEGIN(cppcoreguidelines-avoid-do-while) */                          \
  do {                                                                         \
    if (!(pcms_macro_cond)) {                                                  \
      pcms::detail::AssertFail(#pcms_macro_cond, __FILE__, __LINE__, nullptr); \
    }                                                                          \
  } while (0) /* NOLINTEND(cppcoreguidelines-avoid-do-while) */
#define PCMS_ALWAYS_ASSERT_2(pcms_macro_cond, pcms_macro_comm)                 \
  PCMS_ALWAYS_ASSERT_3(pcms_macro_cond, pcms_macro_comm, nullptr)
#define PCMS_ALWAYS_ASSERT_3(pcms_macro_cond, pcms_macro_comm, pcms_macro_msg) \
  /* NOLINTBEGIN(cppcoreguidelines-avoid-do-while) */                          \
  do {                                                                         \
    if (!(pcms_macro_cond)) {                                                  \
      int pcms_macro_rank = -1;                                                \
      MPI_Comm_rank(pcms_macro_comm, &pcms_macro_rank);                        \
      pcms::detail::AssertFailOnRank(#pcms_macro_cond, __FILE__, __LINE__,     \
                                     pcms_macro_rank, pcms_macro_msg);         \
    }                                                                          \
  } while (0) /* NOLINTEND(cppcoreguidelines-avoid-do-while) */

namespace pcms
{

struct pcms_error : std::runtime_error
{
  using std::runtime_error::runtime_error;
};

} // namespace pcms

#endif // PCMS_UTILITY_ASSERT_HPP
