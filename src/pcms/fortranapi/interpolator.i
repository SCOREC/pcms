%module pcms_interpolator
%{
#include "pcms/capi/interpolator.h"
#include "pcms/capi/kokkos.h"
%}
%include <../external/flibhpc/include/mpi.i>
%include <stdint.i>
%include <typemaps.i>
%import "mesh.i"


struct PcmsInterpolatorHandle {
  void* pointer;
};
typedef struct PcmsInterpolatorHandle PcmsInterpolatorHandle;



PcmsInterpolatorHandle pcms_create_point_based_interpolator(void* source_points, int source_points_size,
                                                                      void* target_points, int target_points_size, double radius, int degree, int min_req_supports, double lambda, double decay_factor);
PcmsInterpolatorHandle pcms_create_degas2xgcnode_interpolator(void* target_points, int target_points_size,
                                                                const char* dg2_mesh_filename, double radius, void* dg2_elem_count, int degree, int min_req_supports, double lambda, double decay_factor);
PcmsInterpolatorHandle pcms_create_xgcnodedegas2_interpolator(const char* dg2_mesh_filename, void* source_points, int source_points_size,
                                                                double radius, void* dg2_elem_count, int degree, int min_req_supports, double lambda, double decay_factor);
PcmsInterpolatorHandle pcms_create_interpolator(PcmsOmegaHMeshHandle oh_mesh, double radius);
void pcms_destroy_interpolator(PcmsInterpolatorHandle interpolator);

void pcms_kokkos_initialize_without_args();
void pcms_kokkos_finalize();

void pcms_interpolate(PcmsInterpolatorHandle interpolator, void* input, int input_size, void* output, int output_size);


// --- Conservative Projection (mesh-intersection-based) ---

struct PcmsConservativeProjectionHandle {
  void* pointer;
};
typedef struct PcmsConservativeProjectionHandle PcmsConservativeProjectionHandle;

PcmsConservativeProjectionHandle pcms_create_conservative_projection(
    const char* source_mesh_name, int source_order,
    const char* target_mesh_name, int target_order);
int pcms_conservative_projection_get_source_size(
    PcmsConservativeProjectionHandle projection);
int pcms_conservative_projection_get_target_size(
    PcmsConservativeProjectionHandle projection);
void pcms_conservative_projection_apply(
    PcmsConservativeProjectionHandle projection,
    void* source_data, int source_size,
    void* target_data, int target_size);
void pcms_destroy_conservative_projection(
    PcmsConservativeProjectionHandle projection);
