%module pcms_mesh
%{
#include "pcms/capi/mesh.h"
%}
%include <../external/flibhpc/include/mpi.i>
%include <stdint.i>
%include <typemaps.i>

struct PcmsOmegaHLibraryHandle
{
  void* lib_handle;
};
typedef struct PcmsOmegaHLibraryHandle PcmsOmegaHLibraryHandle;

struct PcmsOmegaHMeshHandle
{
  void* mesh_handle;
};
typedef struct PcmsOmegaHMeshHandle PcmsOmegaHMeshHandle;


PcmsOmegaHLibraryHandle pcms_create_omega_h_library();
PcmsOmegaHMeshHandle pcms_create_omega_h_mesh(const char* filename, PcmsOmegaHLibraryHandle oh_lib_handle);
void pcms_destroy_omega_h_mesh(PcmsOmegaHMeshHandle oh_mesh);
void pcms_destroy_omega_h_library(PcmsOmegaHLibraryHandle oh_lib_handle);
