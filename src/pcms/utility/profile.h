#ifndef PCMS_SRC_PCMS_PROFILE_H
#define PCMS_SRC_PCMS_PROFILE_H
#include <perfstubs_api/timer.h>

// perfstubs currently has a buffer overflow
// I am disabling this, as it could cause surprising issues until I have time to
// understand and resolve the problem see
// https://github.com/SCOREC/pcms/issues/284
// #define PCMS_FUNCTION_TIMER PERFSTUBS_SCOPED_TIMER_FUNC()
#define PCMS_FUNCTION_TIMER

#endif // PCMS_SRC_PCMS_PROFILE_H
