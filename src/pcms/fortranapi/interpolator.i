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



PcmsInterpolatorHandle pcms_create_interpolator(PcmsInterpolatorOHMeshHandle oh_mesh, double radius);
void pcms_destroy_interpolator(PcmsInterpolatorHandle interpolator);

void pcms_kokkos_initialize_without_args();
void pcms_kokkos_finalize();

PcmsInterpolatorOHMeshHandle read_oh_mesh(const char* filename);
void release_oh_mesh(PcmsInterpolatorOHMeshHandle oh_mesh);

void pcms_interpolate(PcmsInterpolatorHandle interpolator, double* input, int input_size, double* output, int output_size);