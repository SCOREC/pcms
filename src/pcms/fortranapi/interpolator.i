%module pcms_interpolator
%{
#include "pcms/capi/interpolator.h"
#include "pcms/capi/kokkos.h"
%}
%include <../external/flibhpc/include/mpi.i>
%include <stdint.i>
%include <typemaps.i>

struct PcmsInterpolatorOHMeshHandle
{
  void* mesh_handle;
  void* lib_handle;
};
typedef struct PcmsInterpolatorOHMeshHandle PcmsInterpolatorOHMeshHandle;

struct PcmsInterpolatorHandle {
  void* pointer;
};
typedef struct PcmsInterpolatorHandle PcmsInterpolatorHandle;



PcmsInterpolatorHandle pcms_create_point_based_interpolator(void* source_points, int source_points_size,
                                                                      void* target_points, int target_points_size, double radius, int degree, int min_req_supports, double lambda, double decay_factor);
PcmsInterpolatorHandle pcms_create_degas2xgc_interpolator(const char* xgc_mesh_filename, const char* dg2_mesh_filename, double radius, int degree, int min_req_supports, double lambda, double decay_factor);
PcmsInterpolatorHandle pcms_create_degas2xgcnode_interpolator(void* target_points, int target_points_size,
                                                                const char* dg2_mesh_filename, double radius, void* dg2_elem_count, int degree, int min_req_supports, double lambda, double decay_factor);
PcmsInterpolatorHandle pcms_create_xgcnodedegas2_interpolator(const char* dg2_mesh_filename, void* source_points, int source_points_size,
                                                                double radius, void* dg2_elem_count, int degree, int min_req_supports, double lambda, double decay_factor);
PcmsInterpolatorHandle pcms_create_interpolator(PcmsInterpolatorOHMeshHandle oh_mesh, double radius);
void pcms_destroy_interpolator(PcmsInterpolatorHandle interpolator);

void pcms_kokkos_initialize_without_args();
void pcms_kokkos_finalize();

PcmsInterpolatorOHMeshHandle read_oh_mesh(const char* filename);
void release_oh_mesh(PcmsInterpolatorOHMeshHandle oh_mesh);

void pcms_interpolate(PcmsInterpolatorHandle interpolator, void* input, int input_size, void* output, int output_size);
