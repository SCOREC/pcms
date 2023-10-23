#ifndef PCMS_COUPLING_CAPI_CLIENT_H
#define PCMS_COUPLING_CAPI_CLIENT_H
#include <mpi.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

struct PcmsClientHandle { void* pointer; };
typedef struct PcmsClientHandle PcmsClientHandle;
struct PcmsOmegaHMeshHandle { void* pointer; };
typedef struct PcmsOmegaHMeshHandle PcmsOmegaHMeshHandle;
struct PcmsReverseClassificationHandle { void* pointer; };
typedef struct PcmsReverseClassificationHandle PcmsReverseClassificationHandle;
struct PcmsFieldAdapterHandle { void* pointer; };
typedef struct PcmsFieldAdapterHandle PcmsFieldAdapterHandle;
struct PcmsFieldHandle { void* pointer; };
typedef struct PcmsFieldHandle PcmsFieldHandle;

enum PcmsAdapterType
{
  PCMS_ADAPTER_XGC,
  PCMS_ADAPTER_OMEGAH,
  PCMS_ADAPTER_GENE,
  PCMS_ADAPTER_GEM
};
typedef enum PcmsAdapterType PcmsAdapterType;
enum PcmsType
{
  PCMS_FLOAT,
  PCMS_DOUBLE,
  PCMS_INT,
  PCMS_LONG_INT
};
typedef enum PcmsType PcmsType;

//change to a struct holding a pointer
PcmsClientHandle pcms_create_client(const char* name, MPI_Comm comm);
void pcms_destroy_client(PcmsClientHandle);

// returns a pointer to a handle to a reverse classification object
PcmsReverseClassificationHandle pcms_load_reverse_classification(
  const char* file, MPI_Comm comm);
void pcms_destroy_reverse_classification(PcmsReverseClassificationHandle);

// this function is helpful for test cases so we can compute the total number of
// vertexes in the mesh without reading additional files. This function is not
// likely to be needed for production cases
int pcms_reverse_classification_count_verts(
  PcmsReverseClassificationHandle);

// takes in overlap function takes a geometric dimension and a geometric id
// C doesn't have a builtin bool type, so we use int for compatability with C++
typedef int8_t (*in_overlap_function)(int, int);
PcmsFieldAdapterHandle pcms_create_xgc_field_adapter(
  const char* name, MPI_Comm plane_comm, void* data, int size, PcmsType data_type,
  const PcmsReverseClassificationHandle rc, in_overlap_function in_overlap);

PcmsFieldAdapterHandle pcms_create_dummy_field_adapter();

void pcms_destroy_field_adapter(PcmsFieldAdapterHandle);

PcmsFieldHandle pcms_add_field(PcmsClientHandle client_handle,
                                    const char* name,
                                    PcmsFieldAdapterHandle adapter_handle,
                                    int participates);
void pcms_send_field_name(PcmsClientHandle, const char* name);
void pcms_receive_field_name(PcmsClientHandle, const char* name);

void pcms_send_field(PcmsFieldHandle);
void pcms_receive_field(PcmsFieldHandle);

void pcms_begin_send_phase(PcmsClientHandle);
void pcms_end_send_phase(PcmsClientHandle);
void pcms_begin_receive_phase(PcmsClientHandle);
void pcms_end_receive_phase(PcmsClientHandle);
#ifdef __cplusplus
}
#endif

#endif // PCMS_COUPLING_CLIENT_H
