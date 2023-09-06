%module pcms
%{
#include "pcms/capi/client.h"
#include "pcms/capi/kokkos.h"
%}
%include <../external/flibhpc/include/mpi.i>
%include <stdint.i>
%include <typemaps.i>

%fortrancallback("%s") in_overlap_func;
extern "C" {
int8_t in_overlap_func(int dimension, int id);
}


enum WdmCplAdapterType
{
  WDMCPL_ADAPTER_XGC,
  WDMCPL_ADAPTER_OMEGAH,
  WDMCPL_ADAPTER_GENE,
  WDMCPL_ADAPTER_GEM
};
enum WdmCplType
{
  WDMCPL_FLOAT,
  WDMCPL_DOUBLE,
  WDMCPL_INT,
  WDMCPL_LONG_INT
};


WdmCplClientHandle* pcms_create_client(const char* name, MPI_Comm comm);
void pcms_destroy_client(WdmCplClientHandle*);

WdmCplReverseClassificationHandle* pcms_load_reverse_classification(
  const char* file, MPI_Comm comm);
void pcms_destroy_reverse_classification(WdmCplReverseClassificationHandle*);

int pcms_reverse_classification_count_verts(
  WdmCplReverseClassificationHandle*);

typedef int8_t (*in_overlap_function)(int, int);
WdmCplFieldAdapterHandle* pcms_create_xgc_field_adapter(
  const char* name, MPI_Comm plane_comm, void* data, int size, WdmCplType data_type,
  const WdmCplReverseClassificationHandle* rc, in_overlap_function in_overlap);

WdmCplFieldAdapterHandle* pcms_create_dummy_field_adapter();

void pcms_destroy_field_adapter(WdmCplFieldAdapterHandle*);

WdmCplFieldHandle* pcms_add_field(WdmCplClientHandle* client_handle,
                                    const char* name,
                                    WdmCplFieldAdapterHandle* adapter_handle,
                                    int participates);
void pcms_send_field_name(WdmCplClientHandle*, const char* name);
void pcms_receive_field_name(WdmCplClientHandle*, const char* name);

void pcms_send_field(WdmCplFieldHandle*);
void pcms_receive_field(WdmCplFieldHandle*);

void pcms_begin_send_phase(WdmCplClientHandle*);
void pcms_end_send_phase(WdmCplClientHandle*);
void pcms_begin_receive_phase(WdmCplClientHandle*);
void pcms_end_receive_phase(WdmCplClientHandle*);

void pcms_kokkos_initialize_without_args();
void pcms_kokkos_finalize();
