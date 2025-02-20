//
// Created by hasanm4 on 2/17/25.
//

#ifndef PCMS_INTERPOLATOR_CAPI_H
#define PCMS_INTERPOLATOR_CAPI_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

struct PcmsInterpolatorHandle { void* pointer; };
typedef struct PcmsInterpolatorHandle PcmsInterpolatorHandle;

/*
enum for interpolation type
*/

PcmsInterpolatorHandle pcms_create_interpolator(void* oh_mesh, double radius);
void pcms_destroy_interpolator(PcmsInterpolatorHandle interpolator);



#ifdef __cplusplus
}
#endif

#endif // PCMS_INTERPOLATOR_CAPI_H
