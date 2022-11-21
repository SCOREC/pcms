#include "client.h"
#include "wdmcpl.h"
#include "wdmcpl/xgc_field_adapter.h"
#ifdef WDMCPL_HAS_OMEGA_H
  #include "wdmcpl/omega_h_field.h"
#endif
#include <fstream>

namespace wdmcpl
{

template <typename T>
static void wdmcpl_add_field_t(WdmCplClient* client_handle, const char* name,
                               WdmCplFieldAdapter* adapter_handle)
{
  auto* client = reinterpret_cast<wdmcpl::CouplerClient*>(client_handle);
  switch (adapter_handle->type) {
    case WDMCPL_ADAPTER_XGC: {
      auto* adapter =
        reinterpret_cast<wdmcpl::XGCFieldAdapter<T>*>(adapter_handle->xgc_);
      client->AddField(name, *adapter);
      break;
    }
#ifdef WDMCPL_HAS_OMEGA_H
    case WdmCplAdapterType::OmegaH: {
      auto* adapter =
        reinterpret_cast<wdmcpl::OmegaHFieldAdapter<T>*>(adapter_handle->xgc_);
      client->AddField(name, *adapter);
      break;
    }
#endif
    case WDMCPL_ADAPTER_GENE:
    case WDMCPL_ADAPTER_GEM:
      std::cerr << "Unhandled field adapter type in C interface\n";
      std::abort();
  }
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

extern "C" {

[[nodiscard]] WdmCplClient* wdmcpl_create_client(const char* name,
                                                 MPI_Comm comm)
{
  auto* client = new wdmcpl::CouplerClient(name, comm);
  return reinterpret_cast<WdmCplClient*>(client);
}
void wdmcpl_destroy_client(WdmCplClient* client)
{
  if (client != nullptr)
    delete reinterpret_cast<wdmcpl::CouplerClient*>(client);
}

void wdmcpl_add_field(WdmCplClient* client_handle, const char* name,
                      WdmCplFieldAdapter* adapter_handle)
{
  wdmcpl::ApplyToTypes(adapter_handle->data_type, [&](auto d) {
    using T = decltype(d);
    wdmcpl::wdmcpl_add_field_t<T>(client_handle, name, adapter_handle);
  });
}
void wdmcpl_send_field(WdmCplClient* client_handle, const char* name)
{
  auto* client = reinterpret_cast<wdmcpl::CouplerClient*>(client_handle);
  client->SendField(name);
}
void wdmcpl_receive_field(WdmCplClient* client_handle, const char* name)
{
  auto* client = reinterpret_cast<wdmcpl::CouplerClient*>(client_handle);
  client->ReceiveField(name);
}
WdmCplFieldAdapter* wdmcpl_create_xgc_field_adapter(const char* name,
                                                    void* data, int size,
                                                    WdmCplType data_type,
                                                    char* classification_file)
{
  std::ifstream infile(classification_file);
  auto classification = wdmcpl::ReadXGCNodeClassification(infile);
  WdmCplFieldAdapter* field_adapter = new WdmCplFieldAdapter;
  field_adapter->type = WDMCPL_ADAPTER_XGC;
  switch (data_type) {
    case WDMCPL_DOUBLE: field_adapter->data_type = WDMCPL_DOUBLE; break;
    case WDMCPL_FLOAT: field_adapter->data_type = WDMCPL_FLOAT; break;
    case WDMCPL_INT: field_adapter->data_type = WDMCPL_INT; break;
    case WDMCPL_LONG_INT: field_adapter->data_type = WDMCPL_LONG_INT; break;
  }
  wdmcpl::ApplyToTypes(data_type, [&](auto d) {
    using T = decltype(d);
    wdmcpl::ScalarArrayView<T, wdmcpl::HostMemorySpace> data_view(
      reinterpret_cast<T*>(data), size);
    auto* adapter = new wdmcpl::XGCFieldAdapter<T>(
      name, std::move(classification.dimension),
      std::move(classification.geometric_id), data_view);
    field_adapter->xgc_ = reinterpret_cast<WdmCplXgcAdapter*>(adapter);
  });
  return field_adapter;
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
void wdmcpl_destroy_field_adapter(WdmCplFieldAdapter* adapter)
{
  switch (adapter->type) {
    case (WDMCPL_ADAPTER_XGC):
      wdmcpl_destroy_xgc_adapter(adapter->xgc_, adapter->data_type);
      break;
#ifdef WDMCPL_HAS_OMEGA_H
    case (WdmCplAdapterType::OmegaH): break;
#endif
    case (WDMCPL_ADAPTER_GENE): break;
    case (WDMCPL_ADAPTER_GEM): break;
  }

  if (adapter != nullptr) {
    delete adapter;
    adapter = nullptr;
  }
}
}