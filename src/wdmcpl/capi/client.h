#ifndef WDM_COUPLING_CAPI_CLIENT_H
#define WDM_COUPLING_CAPI_CLIENT_H
#include <mpi.h>

extern "C" {
struct WdmCplClient;
typedef struct WdmCplClient WdmCplClient;
// TODO move into XGC specific header
struct WdmCplXgcAdapter;
typedef struct WdmCplXgcAdapter WdmCplXgcAdapter;

struct WdmCplOmegaHMeshHandle;
typedef struct WdmCplOmegaHMeshHandle WdmCplOmegaHMeshHandle;

struct WdmCplAdapterType_s
{
  enum Type
  {
    XGC,
#ifdef WDMCPL_HAS_OMEGA_H
    OmegaH,
#endif
    GENE,
    GEM
  };
  // empty struct in C++ is 1 byte and 0 byte in C
  char dummy;
};
typedef WdmCplAdapterType_s::Type WdmCplAdapterType;
struct WdmCplType_s {
  enum Type {
    FLOAT,
    DOUBLE,
    INT,
    LONG_INT
  };
  // empty struct in C++ is 1 byte and 0 byte in C
  char dummy;
};
typedef WdmCplType_s::Type WdmCplType;

struct WdmCplFieldAdapter
{
  WdmCplAdapterType type;
  WdmCplType data_type;
  /*
  enum Type
  {
    XGC,
    OmegaH,
    GENE,
    GEM
  } type;
  enum DataType
  {
    FLOAT,
    DOUBLE,
    INT,
    LONG_INT
  } data_type;
   */
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
[[nodiscard]] WdmCplClient* wdmcpl_create_client(const char* name,
                                                 MPI_Comm comm);
void wdmcpl_destroy_client(WdmCplClient*);

WdmCplFieldAdapter* wdmcpl_create_xgc_field_adapter(
  const char* name, void* data, int size,
  WdmCplType data_type, char* classification_file);

void wdmcpl_destroy_field_adapter(WdmCplFieldAdapter*);

void wdmcpl_add_field(WdmCplClient* client_handle, const char* name,
                      WdmCplFieldAdapter* adapter_handle);
void wdmcpl_send_field(WdmCplClient*, const char* name);
void wdmcpl_receive_field(WdmCplClient*, const char* name);
};

#endif // WDM_COUPLING_CLIENT_H
