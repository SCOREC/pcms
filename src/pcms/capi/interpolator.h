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

struct PcmsInterpolatorOHMeshHandle {
  void* mesh_handle;
  void* lib_handle;
};
typedef struct PcmsInterpolatorOHMeshHandle PcmsInterpolatorOHMeshHandle;

/*
enum for interpolation type
*/

PcmsInterpolatorHandle pcms_create_interpolator(PcmsInterpolatorOHMeshHandle oh_mesh, double radius);
void pcms_destroy_interpolator(PcmsInterpolatorHandle interpolator);

PcmsInterpolatorOHMeshHandle read_oh_mesh(const char* filename);
void release_oh_mesh(PcmsInterpolatorOHMeshHandle oh_mesh_handle);


void pcms_interpolate(PcmsInterpolatorHandle interpolator, double* input, int input_size, double* output, int output_size);



#ifdef __cplusplus
}
#endif

#endif // PCMS_INTERPOLATOR_CAPI_H
