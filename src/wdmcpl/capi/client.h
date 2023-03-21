#ifndef WDM_COUPLING_CAPI_CLIENT_H
#define WDM_COUPLING_CAPI_CLIENT_H
#include <mpi.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif
struct WdmCplClientHandle;
typedef struct WdmCplClientHandle WdmCplClientHandle;
struct WdmCplOmegaHMeshHandle;
typedef struct WdmCplOmegaHMeshHandle WdmCplOmegaHMeshHandle;
struct WdmCplReverseClassificationHandle;
typedef struct WdmCplReverseClassificationHandle WdmCplReverseClassificationHandle;
struct WdmCplFieldAdapterHandle;
typedef struct WdmCplFieldAdapterHandle WdmCplFieldAdapterHandle;
struct WdmCplFieldHandle;
typedef struct WdmCplFieldHandle WdmCplFieldHandle;

enum WdmCplAdapterType
{
  WDMCPL_ADAPTER_XGC,
  WDMCPL_ADAPTER_OMEGAH,
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

WdmCplClientHandle* wdmcpl_create_client(const char* name, MPI_Comm comm);
void wdmcpl_destroy_client(WdmCplClientHandle*);

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
WdmCplFieldAdapterHandle* wdmcpl_create_xgc_field_adapter(
  const char* name, MPI_Comm plane_comm, void* data, int size, WdmCplType data_type,
  const WdmCplReverseClassificationHandle* rc, in_overlap_function in_overlap);

WdmCplFieldAdapterHandle* wdmcpl_create_dummy_field_adapter();

void wdmcpl_destroy_field_adapter(WdmCplFieldAdapterHandle*);

WdmCplFieldHandle* wdmcpl_add_field(WdmCplClientHandle* client_handle,
                                    const char* name,
                                    WdmCplFieldAdapterHandle* adapter_handle,
                                    int participates);
void wdmcpl_send_field_name(WdmCplClientHandle*, const char* name);
void wdmcpl_receive_field_name(WdmCplClientHandle*, const char* name);

void wdmcpl_send_field(WdmCplFieldHandle*);
void wdmcpl_receive_field(WdmCplFieldHandle*);

void wdmcpl_begin_send_phase(WdmCplClientHandle*);
void wdmcpl_end_send_phase(WdmCplClientHandle*);
void wdmcpl_begin_receive_phase(WdmCplClientHandle*);
void wdmcpl_end_receive_phase(WdmCplClientHandle*);
#ifdef __cplusplus
}
#endif

#endif // WDM_COUPLING_CLIENT_H
