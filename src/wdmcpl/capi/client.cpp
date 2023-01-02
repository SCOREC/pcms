#include "client.h"
#include "wdmcpl.h"
#include "wdmcpl/xgc_field_adapter.h"
// #ifdef WDMCPL_HAS_OMEGA_H
//   #include "wdmcpl/omega_h_field.h"
// #endif
#include <fstream>
#include "wdmcpl/xgc_reverse_classification.h"
extern "C" {
// TODO move into XGC specific header
struct WdmCplXgcAdapter;
typedef struct WdmCplXgcAdapter WdmCplXgcAdapter;
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
};
namespace wdmcpl
{

template <typename T>
static WdmCplFieldHandle* wdmcpl_add_field_t(WdmCplClientHandle* client_handle, const char* name,
                               WdmCplFieldAdapter* adapter_handle)
{
  auto* client = reinterpret_cast<wdmcpl::CouplerClient*>(client_handle);
  CoupledField* field = nullptr;
  switch (adapter_handle->type) {
    case WDMCPL_ADAPTER_XGC: {
      auto* adapter =
        reinterpret_cast<wdmcpl::XGCFieldAdapter<T>*>(adapter_handle->xgc_);
      field = client->AddField(name, *adapter);
      break;
    }
#ifdef WDMCPL_HAS_OMEGA_H
    case WDMCPL_ADAPTER_OMEGAH: {
      if constexpr (std::is_same_v<T, double> || std::is_same_v<T, int>) {
        auto* adapter = reinterpret_cast<wdmcpl::OmegaHFieldAdapter<T>*>(
          adapter_handle->xgc_);
        field = client->AddField(name, *adapter);
      } else {
        // error trying to use omega_h adapter with unsupported data type
        Wdmcpl_Assert_Fail("Omega_h can only use double or int fields.");
      }
      break;
    }
#endif
    case WDMCPL_ADAPTER_GENE:
    case WDMCPL_ADAPTER_GEM:
      Wdmcpl_Assert_Fail("Unhandled field adapter type in C interface\n");
  }
  return reinterpret_cast<WdmCplFieldHandle*>(field);
}

template <typename F>
auto ApplyToTypes(WdmCplType type, const F&& f)
{
  switch (type) {
    case WDMCPL_DOUBLE: {
      return f(double{});
    }
    case WDMCPL_FLOAT: {
      return f(float{});
    }
    case WDMCPL_INT: {
      return f(int{});
    }
    case WDMCPL_LONG_INT: {
      return f(long{});
    }
  }
  // quash compiler warning
  using dummy_type = typename std::result_of_t<decltype(f)(int)>;
  return static_cast<dummy_type>(nullptr);
}

} // namespace wdmcpl

[[nodiscard]] WdmCplClientHandle* wdmcpl_create_client(const char* name,
                                                 MPI_Comm comm)
{
  auto* client = new wdmcpl::CouplerClient(name, comm);
  return reinterpret_cast<WdmCplClientHandle*>(client);
}
void wdmcpl_destroy_client(WdmCplClientHandle* client)
{
  if (client != nullptr)
    delete reinterpret_cast<wdmcpl::CouplerClient*>(client);
}
WdmCplReverseClassificationHandle* wdmcpl_load_reverse_classification(
  const char* file, MPI_Comm comm)
{
  std::filesystem::path filepath{file};
  auto* rc = new wdmcpl::ReverseClassificationVertex{
    wdmcpl::ReadReverseClassificationVertex(filepath, comm)};
  return reinterpret_cast<WdmCplReverseClassificationHandle*>(rc);
}
void wdmcpl_destroy_reverse_classification(WdmCplReverseClassificationHandle* rc)
{
  if(rc != nullptr)
    delete reinterpret_cast<wdmcpl::ReverseClassificationVertex*>(rc);
}

WdmCplFieldHandle* wdmcpl_add_field(WdmCplClientHandle* client_handle, const char* name,
                      WdmCplFieldAdapterHandle* adapter_handle)
{

  WdmCplFieldHandle* field = nullptr;
  auto* adapter = reinterpret_cast<WdmCplFieldAdapter*>(adapter_handle);
  wdmcpl::ApplyToTypes(adapter->data_type, [&](auto d) {
    using T = decltype(d);
    field = wdmcpl::wdmcpl_add_field_t<T>(client_handle, name, adapter);
  });
  return field;
}
void wdmcpl_send_field_name(WdmCplClientHandle* client_handle, const char* name)
{
  auto* client = reinterpret_cast<wdmcpl::CouplerClient*>(client_handle);
  client->SendField(name);
}
void wdmcpl_receive_field_name(WdmCplClientHandle* client_handle, const char* name)
{
  auto* client = reinterpret_cast<wdmcpl::CouplerClient*>(client_handle);
  client->ReceiveField(name);
}
void wdmcpl_send_field(WdmCplFieldHandle* field_handle) {
  auto * field = reinterpret_cast<wdmcpl::CoupledField*>(field_handle);
  field->Send();
}
void wdmcpl_receive_field(WdmCplFieldHandle* field_handle) {
  auto * field = reinterpret_cast<wdmcpl::CoupledField*>(field_handle);
  field->Receive();
}
WdmCplFieldAdapterHandle* wdmcpl_create_xgc_field_adapter(
  const char* name, void* data, int size, WdmCplType data_type,
  const WdmCplReverseClassificationHandle* rc, in_overlap_function in_overlap)
{
  // auto classification = wdmcpl::ReadXGCNodeClassification(infile);
  auto* field_adapter = new WdmCplFieldAdapter;
  field_adapter->type = WDMCPL_ADAPTER_XGC;
  WDMCPL_ALWAYS_ASSERT(rc != nullptr);
  auto* reverse_classification =
    reinterpret_cast<const wdmcpl::ReverseClassificationVertex*>(rc);
  switch (data_type) {
    case WDMCPL_DOUBLE: field_adapter->data_type = WDMCPL_DOUBLE; break;
    case WDMCPL_FLOAT: field_adapter->data_type = WDMCPL_FLOAT; break;
    case WDMCPL_INT: field_adapter->data_type = WDMCPL_INT; break;
    case WDMCPL_LONG_INT: field_adapter->data_type = WDMCPL_LONG_INT; break;
  }
  WDMCPL_ALWAYS_ASSERT(reverse_classification != nullptr);
  wdmcpl::ApplyToTypes(data_type, [&](auto d) {
    using T = decltype(d);
    wdmcpl::ScalarArrayView<T, wdmcpl::HostMemorySpace> data_view(
      reinterpret_cast<T*>(data), size);
    auto* adapter =
      new wdmcpl::XGCFieldAdapter<T>(name, data_view, *reverse_classification, in_overlap);
    // auto* adapter = new wdmcpl::XGCFieldAdapter<T>(
    //   name, std::move(classification.dimension),
    //   std::move(classification.geometric_id), data_view);
    field_adapter->xgc_ = reinterpret_cast<WdmCplXgcAdapter*>(adapter);
  });
  return reinterpret_cast<WdmCplFieldAdapterHandle*>(field_adapter);
}

static void wdmcpl_destroy_xgc_adapter(WdmCplXgcAdapter* xgc_adapter,
                                       WdmCplType data_type)
{
  if (xgc_adapter == nullptr) {
    return;
  }
  wdmcpl::ApplyToTypes(data_type, [&xgc_adapter](auto d) {
    using T = decltype(d);
    delete reinterpret_cast<wdmcpl::XGCFieldAdapter<T>*>(xgc_adapter);
    xgc_adapter = nullptr;
  });
}
void wdmcpl_destroy_field_adapter(WdmCplFieldAdapterHandle* adapter_handle)
{
  auto* adapter = reinterpret_cast<WdmCplFieldAdapter*>(adapter_handle);
  switch (adapter->type) {
    case (WDMCPL_ADAPTER_XGC):
      wdmcpl_destroy_xgc_adapter(adapter->xgc_, adapter->data_type);
      break;
#ifdef WDMCPL_HAS_OMEGA_H
    case (WDMCPL_ADAPTER_OMEGAH): break;
#endif
    case (WDMCPL_ADAPTER_GENE): break;
    case (WDMCPL_ADAPTER_GEM): break;
  }

  if (adapter != nullptr) {
    delete adapter;
    adapter = nullptr;
  }
}
int wdmcpl_reverse_classification_count_verts(
  WdmCplReverseClassificationHandle* rc)
{
  auto* reverse_classification =
    reinterpret_cast<const wdmcpl::ReverseClassificationVertex*>(rc);
  return std::accumulate(reverse_classification->begin(),
                         reverse_classification->end(), 0,
                         [](auto current, const auto& verts) {
                           return current + verts.second.size();
                         });
}
