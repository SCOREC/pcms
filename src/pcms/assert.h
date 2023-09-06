#ifndef WDM_COUPLING_ASSERT_H
#define WDM_COUPLING_ASSERT_H
#include <cassert>

#define WDMCPL_ALWAYS_ASSERT(cond) do {          \
  if (! (cond)) {                               \
    char omsg[2048];                            \
    snprintf(omsg, 2048, "%s failed at %s + %d \n",    \
      #cond, __FILE__, __LINE__);               \
    pcms::Wdmcpl_Assert_Fail(omsg);             \
  }                                             \
} while (0)


namespace pcms {
//from scorec/core/pcu_fail.h
void Wdmcpl_Assert_Fail(const char* msg) __attribute__ ((noreturn));
}


#endif // WDM_COUPLING_ASSERT_H
