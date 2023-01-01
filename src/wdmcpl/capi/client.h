#ifndef WDM_COUPLING_CAPI_CLIENT_H
#define WDM_COUPLING_CAPI_CLIENT_H
#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif
struct WdmCplClient;
typedef struct WdmCplClient WdmCplClient;
// TODO move into XGC specific header
struct WdmCplXgcAdapter;
typedef struct WdmCplXgcAdapter WdmCplXgcAdapter;
struct WdmCplOmegaHMeshHandle;
typedef struct WdmCplOmegaHMeshHandle WdmCplOmegaHMeshHandle;
struct WdmCplReverseClassificationHandle;
typedef struct WdmCplReverseClassificationHandle WdmCplReverseClassificationHandle;

struct WdmCplFieldHandle;
typedef struct WdmCplFieldHandle WdmCplFieldHandle;

enum WdmCplAdapterType
{
  WDMCPL_ADAPTER_XGC,
//#ifdef WDMCPL_HAS_OMEGA_H
  WDMCPL_ADAPTER_OMEGAH,
//#endif
  WDMCPL_ADAPTER_GENE,
  WDMCPL_ADAPTER_GEM
};
typedef enum WdmCplAdapterType WdmCplAdapterType;
enum WdmCplType
{
  WDMCPL_FLOAT,
  WDMCPL_DOUBLE,
  WDMCPL_INT,
  WDMCPL_LONG_INT
};
typedef enum WdmCplType WdmCplType;

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
typedef struct WdmCplFieldAdapter WdmCplFieldAdapter;

WdmCplClient* wdmcpl_create_client(const char* name, MPI_Comm comm);
void wdmcpl_destroy_client(WdmCplClient*);

// returns a pointer to a handle to a reverse classification object
WdmCplReverseClassificationHandle* wdmcpl_load_reverse_classification(
  const char* file, MPI_Comm comm);
void wdmcpl_destroy_reverse_classification(WdmCplReverseClassificationHandle*);

// this function is helpful for test cases so we can compute the total number of
// vertexes in the mesh without reading additional files. This function is not
// likely to be needed for production cases
int wdmcpl_reverse_classification_count_verts(
  WdmCplReverseClassificationHandle*);

// takes in overlap function takes a geometric dimension and a geometric id
// C doesn't have a builtin bool type, so we use int for compatability with C++
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
#ifdef __cplusplus
}
#endif

#endif // WDM_COUPLING_CLIENT_H
