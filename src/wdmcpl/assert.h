#ifndef WDM_COUPLING_ASSERT_H
#define WDM_COUPLING_ASSERT_H
#include <cassert>

#define WDMCPL_ALWAYS_ASSERT(cond) do {          \
  if (! (cond)) {                               \
    char omsg[2048];                            \
    sprintf(omsg, "%s failed at %s + %d \n",    \
      #cond, __FILE__, __LINE__);               \
    wdmcpl::Wdmcpl_Assert_Fail(omsg);             \
  }                                             \
} while (0)


namespace wdmcpl {
//from scorec/core/pcu_fail.h
void Wdmcpl_Assert_Fail(const char* msg) __attribute__ ((noreturn));
}


#endif // WDM_COUPLING_ASSERT_H
