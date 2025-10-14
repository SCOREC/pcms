//
// Created by hasanm4 on 2/17/25.
//

#ifndef PCMS_INTERPOLATOR_CAPI_H
#define PCMS_INTERPOLATOR_CAPI_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

struct PcmsInterpolatorHandle
{
  void* pointer;
};
typedef struct PcmsInterpolatorHandle PcmsInterpolatorHandle;

struct PcmsPointBasedInterpolatorHandle
{
  void* pointer;
};
typedef struct PcmsPointBasedInterpolatorHandle
  PcmsPointBasedInterpolatorHandle;

struct PcmsInterpolatorOHMeshHandle
{
  void* mesh_handle;
  void* lib_handle;
};
typedef struct PcmsInterpolatorOHMeshHandle PcmsInterpolatorOHMeshHandle;

/*
enum for interpolation type
*/

PcmsInterpolatorHandle pcms_create_interpolator(
  PcmsInterpolatorOHMeshHandle oh_mesh, double radius);
PcmsPointBasedInterpolatorHandle pcms_create_point_based_interpolator(
  void* source_points, int source_points_size, void* target_points,
  int target_points_size, double radius, int degree, int min_req_supports,
  double lambda, double decay_factor);
PcmsPointBasedInterpolatorHandle pcms_create_degas2xgc_interpolator(
  const char* xgc_mesh_filename, const char* dg2_mesh_filename, double radius,
  int degree, int min_req_supports, double lambda, double decay_factor);
PcmsPointBasedInterpolatorHandle pcms_create_degas2xgcnode_interpolator(
  void* target_points, int target_points_size, const char* dg2_mesh_filename,
  double radius, void* dg2_elem_count, int degree, int min_req_supports,
  double lambda, double decay_factor);
PcmsPointBasedInterpolatorHandle pcms_create_xgcnodedegas2_interpolator(
  const char* dg2_mesh_filename, void* source_points, int source_points_size,
  double radius, void* dg2_elem_count, int degree, int min_req_supports,
  double lambda, double decay_factor);
void pcms_destroy_interpolator(PcmsInterpolatorHandle interpolator);
void pcms_destroy_point_based_interpolator(
  PcmsPointBasedInterpolatorHandle interpolator);

PcmsInterpolatorOHMeshHandle read_oh_mesh(const char* filename);
void release_oh_mesh(PcmsInterpolatorOHMeshHandle oh_mesh_handle);

void pcms_interpolate(PcmsInterpolatorHandle interpolator, void* input,
                      int input_size, void* output, int output_size);
void pcms_interpolate_point_based(PcmsPointBasedInterpolatorHandle interpolator,
                                  void* input, int input_size, void* output,
                                  int output_size);

#ifdef __cplusplus
}
#endif

#endif // PCMS_INTERPOLATOR_CAPI_H
