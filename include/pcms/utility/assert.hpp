#ifndef PCMS_COUPLING_ASSERT_H
#define PCMS_COUPLING_ASSERT_H
#include <stdexcept>
#include <mpi.h>

// https://stackoverflow.com/questions/16683146/can-macros-be-overloaded-by-number-of-arguments
#define PCMS_CAT(A, B) A##B
#define PCMS_SELECT(NAME, NUM) PCMS_CAT(NAME##_, NUM)

#define PCMS_GET_COUNT(_1, _2, _3, _4, _5, _6 /* ad nauseam */, COUNT, ...)    \
  COUNT
#define PCMS_VA_SIZE(...) PCMS_GET_COUNT(__VA_ARGS__, 6, 5, 4, 3, 2, 1)

#define PCMS_VA_SELECT(NAME, ...)                                              \
  PCMS_SELECT(NAME, PCMS_VA_SIZE(__VA_ARGS__))(__VA_ARGS__)

#define PCMS_ALWAYS_ASSERT(...) PCMS_VA_SELECT(PCMS_ALWAYS_ASSERT, __VA_ARGS__)
#define PCMS_ALWAYS_ASSERT_1(cond) PCMS_ALWAYS_ASSERT_2(cond, MPI_COMM_WORLD)
#define PCMS_ALWAYS_ASSERT_2(pcms_macro_cond, pcms_macro_comm)                                                 \
  /* NOLINTBEGIN(modernize-avoid-c-arrays,cppcoreguidelines-avoid-do-while,cppcoreguidelines-avoid-c-arrays)*/ \
  do {                                                                                                         \
    if (!(pcms_macro_cond)) {                                                                                  \
      int pcms_macro_rank = -1;                                                                                \
      MPI_Comm_rank(pcms_macro_comm, &pcms_macro_rank);                                                        \
      char pcms_macro_omsg[2048];                                                                              \
      snprintf(pcms_macro_omsg, 2048, "%s failed at %s + %d on rank %d\n",                                     \
               #pcms_macro_cond, __FILE__, __LINE__, pcms_macro_rank);                                         \
      pcms::Pcms_Assert_Fail(pcms_macro_omsg);                                                                 \
    }                                                                                                          \
  } while (                                                                                                    \
    0) /* NOLINTEND(modernize-avoid-c-arrays,cppcoreguidelines-avoid-do-while,cppcoreguidelines-avoid-c-arrays)*/
#define PCMS_ALWAYS_ASSERT_3(pcms_macro_cond, pcms_macro_comm, pcms_macro_msg)                                 \
  /* NOLINTBEGIN(modernize-avoid-c-arrays,cppcoreguidelines-avoid-do-while,cppcoreguidelines-avoid-c-arrays)*/ \
  do {                                                                                                         \
    if (!(pcms_macro_cond)) {                                                                                  \
      int pcms_macro_rank = -1;                                                                                \
      MPI_Comm_rank(pcms_macro_comm, &pcms_macro_rank);                                                        \
      char pcms_macro_omsg[2048];                                                                              \
      snprintf(pcms_macro_omsg, 2048, "%s failed at %s + %d on rank %d: %s\n",                                 \
               #pcms_macro_cond, __FILE__, __LINE__, pcms_macro_rank,                                          \
               pcms_macro_msg);                                                                                \
      pcms::Pcms_Assert_Fail(pcms_macro_omsg);                                                                 \
    }                                                                                                          \
  } while (                                                                                                    \
    0) /* NOLINTEND(modernize-avoid-c-arrays,cppcoreguidelines-avoid-do-while,cppcoreguidelines-avoid-c-arrays)*/

namespace pcms
{

struct pcms_error : std::runtime_error
{
  using std::runtime_error::runtime_error;
};

// from scorec/core/pcu_fail.h
void Pcms_Assert_Fail(const char* msg) __attribute__((noreturn));
} // namespace pcms

#endif // PCMS_COUPLING_ASSERT_H
