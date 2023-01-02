%module wdmcpl
%{
#include "wdmcpl/capi/client.h"
%}
%include <../external/flibhpc/include/mpi.i>
%include <stdint.i>
%include <typemaps.i>

%fortrancallback("%s") in_overlap_func;
extern "C" {
int8_t in_overlap_func(int dimension, int id);
}

// TODO try with bodies to improve names once other stuff works..
struct WdmCplClient;
struct WdmCplXgcAdapter;
struct WdmCplOmegaHMeshHandle;
struct WdmCplReverseClassificationHandle;
struct WdmCplFieldHandle;

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

struct WdmCplFieldAdapter
{
  WdmCplAdapterType type;
  WdmCplType data_type;
  union
  {
    WdmCplXgcAdapter* xgc_;
    // OmegaHAdapter* omega_h;
    // GENEAdapter* gene_;
    // GEMAdapter* gem_;
  };
};

WdmCplClient* wdmcpl_create_client(const char* name, MPI_Comm comm);
void wdmcpl_destroy_client(WdmCplClient*);

WdmCplReverseClassificationHandle* wdmcpl_load_reverse_classification(
  const char* file, MPI_Comm comm);
void wdmcpl_destroy_reverse_classification(WdmCplReverseClassificationHandle*);

int wdmcpl_reverse_classification_count_verts(
  WdmCplReverseClassificationHandle*);

typedef int8_t (*in_overlap_function)(int, int);
WdmCplFieldAdapter* wdmcpl_create_xgc_field_adapter(
  const char* name, void* data, int size, WdmCplType data_type,
  const WdmCplReverseClassificationHandle* rc, in_overlap_function in_overlap);

void wdmcpl_destroy_field_adapter(WdmCplFieldAdapter*);

WdmCplFieldHandle* wdmcpl_add_field(WdmCplClient* client_handle, const char* name,
                      WdmCplFieldAdapter* adapter_handle);
void wdmcpl_send_field_name(WdmCplClient*, const char* name);
void wdmcpl_receive_field_name(WdmCplClient*, const char* name);

void wdmcpl_send_field(WdmCplFieldHandle*);
void wdmcpl_receive_field(WdmCplFieldHandle*);


