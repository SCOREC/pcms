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

enum WdmCplAdapterType
{
  WDMCPL_ADAPTER_XGC,
#ifdef WDMCPL_HAS_OMEGA_H
  WDMCPL_ADAPTER_OMEGAH,
#endif
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

// FIXME communicate partition type in underlying library
WdmCplClient* wdmcpl_create_client(const char* name, MPI_Comm comm);
void wdmcpl_destroy_client(WdmCplClient*);

WdmCplFieldAdapter* wdmcpl_create_xgc_field_adapter(const char* name,
                                                    void* data, int size,
                                                    WdmCplType data_type,
                                                    char* classification_file);

void wdmcpl_destroy_field_adapter(WdmCplFieldAdapter*);

void wdmcpl_add_field(WdmCplClient* client_handle, const char* name,
                      WdmCplFieldAdapter* adapter_handle);
void wdmcpl_send_field(WdmCplClient*, const char* name);
void wdmcpl_receive_field(WdmCplClient*, const char* name);
#ifdef __cplusplus
}
#endif

#endif // WDM_COUPLING_CLIENT_H
