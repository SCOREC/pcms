#ifndef PCMS_COUPLING_ASSERT_H
#define PCMS_COUPLING_ASSERT_H
#include <cassert>

#define PCMS_ALWAYS_ASSERT(cond) do {          \
  if (! (cond)) {                               \
    char omsg[2048];                            \
    snprintf(omsg, 2048, "%s failed at %s + %d \n",    \
      #cond, __FILE__, __LINE__);               \
    pcms::Pcms_Assert_Fail(omsg);             \
  }                                             \
} while (0)


namespace pcms {
//from scorec/core/pcu_fail.h
void Pcms_Assert_Fail(const char* msg) __attribute__ ((noreturn));
}


#endif // PCMS_COUPLING_ASSERT_H
