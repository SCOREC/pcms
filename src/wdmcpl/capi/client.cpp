#include "client.h"
#include "wdmcpl.h"
#include "wdmcpl/xgc_field_adapter.h"
#include <variant>
#include <redev_variant_tools.h>
// #ifdef WDMCPL_HAS_OMEGA_H
//   #include "wdmcpl/omega_h_field.h"
// #endif
#include <fstream>
#include "wdmcpl/xgc_reverse_classification.h"
#include "wdmcpl/dummy_field_adapter.h"
namespace wdmcpl
{
// Note that we have a closed set of types that can be used in the C interface
using FieldAdapterVariant =
  std::variant<std::monostate, wdmcpl::XGCFieldAdapter<double>,
               wdmcpl::XGCFieldAdapter<float>, wdmcpl::XGCFieldAdapter<int>,
               wdmcpl::XGCFieldAdapter<long>,
               wdmcpl::DummyFieldAdapter
//#ifdef WDMCPL_HAS_OMEGA_H
//               ,
//               wdmcpl::OmegaHFieldAdapter<double>,
//               wdmcpl::OmegaHFieldAdapter<int>
//#endif
               >;

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
  //std::filesystem::path filepath{file};
  auto* rc = new wdmcpl::ReverseClassificationVertex{
    wdmcpl::ReadReverseClassificationVertex(file, comm)};
  return reinterpret_cast<WdmCplReverseClassificationHandle*>(rc);
}
void wdmcpl_destroy_reverse_classification(
  WdmCplReverseClassificationHandle* rc)
{
  if (rc != nullptr)
    delete reinterpret_cast<wdmcpl::ReverseClassificationVertex*>(rc);
}

struct AddFieldVariantOperators {
  AddFieldVariantOperators(const char* name, wdmcpl::CouplerClient* client, int participates)
  : name_(name), client_(client), participates_(participates)
  {
  }

  [[nodiscard]]
  wdmcpl::CoupledField* operator()(const std::monostate&) const noexcept { return nullptr; }
  template <typename FieldAdapter>
  [[nodiscard]]
  wdmcpl::CoupledField* operator()(const FieldAdapter& field_adapter) const noexcept {
        return client_->AddField(name_, field_adapter, participates_);
  }

  const char* name_;
  wdmcpl::CouplerClient* client_;
  bool participates_;
};

WdmCplFieldHandle* wdmcpl_add_field(WdmCplClientHandle* client_handle,
                                    const char* name,
                                    WdmCplFieldAdapterHandle* adapter_handle,
                                    int participates)
{

  auto* adapter =
    reinterpret_cast<wdmcpl::FieldAdapterVariant*>(adapter_handle);
  auto* client = reinterpret_cast<wdmcpl::CouplerClient*>(client_handle);
  WDMCPL_ALWAYS_ASSERT(client != nullptr);
  WDMCPL_ALWAYS_ASSERT(adapter != nullptr);
  //wdmcpl::CoupledField* field = std::visit(
  //  redev::overloaded{
  //    [](const std::monostate&) -> wdmcpl::CoupledField* { return nullptr; },
  //    [&name, &client, participates](const auto& field_adapter) {
  //      return client->AddField(name, field_adapter, participates);
  //    }},
  //  *adapter);
  wdmcpl::CoupledField* field = std::visit(AddFieldVariantOperators{name, client, participates},*adapter);
  return reinterpret_cast<WdmCplFieldHandle*>(field);
}
void wdmcpl_send_field_name(WdmCplClientHandle* client_handle, const char* name)
{
  auto* client = reinterpret_cast<wdmcpl::CouplerClient*>(client_handle);
  WDMCPL_ALWAYS_ASSERT(client != nullptr);
  client->SendField(name);
}
void wdmcpl_receive_field_name(WdmCplClientHandle* client_handle,
                               const char* name)
{
  auto* client = reinterpret_cast<wdmcpl::CouplerClient*>(client_handle);
  WDMCPL_ALWAYS_ASSERT(client != nullptr);
  client->ReceiveField(name);
}
void wdmcpl_send_field(WdmCplFieldHandle* field_handle)
{
  auto* field = reinterpret_cast<wdmcpl::CoupledField*>(field_handle);
  WDMCPL_ALWAYS_ASSERT(field != nullptr);
  field->Send();
}
void wdmcpl_receive_field(WdmCplFieldHandle* field_handle)
{
  auto* field = reinterpret_cast<wdmcpl::CoupledField*>(field_handle);
  WDMCPL_ALWAYS_ASSERT(field != nullptr);
  field->Receive();
}
template <typename T>
void wdmcpl_create_xgc_field_adapter_t(
  const char* name, MPI_Comm comm, void* data, int size,
  const wdmcpl::ReverseClassificationVertex& reverse_classification,
  in_overlap_function in_overlap, wdmcpl::FieldAdapterVariant& field_adapter)
{
  WDMCPL_ALWAYS_ASSERT((size >0) ? (data!=nullptr) : true);
  wdmcpl::ScalarArrayView<T, wdmcpl::HostMemorySpace> data_view(
    reinterpret_cast<T*>(data), size);
  field_adapter.emplace<wdmcpl::XGCFieldAdapter<T>>(
    name, comm, data_view, reverse_classification, in_overlap);
}
WdmCplFieldAdapterHandle* wdmcpl_create_xgc_field_adapter(
  const char* name, MPI_Comm comm, void* data, int size, WdmCplType data_type,
  const WdmCplReverseClassificationHandle* rc, in_overlap_function in_overlap)
{
  auto* field_adapter = new wdmcpl::FieldAdapterVariant{};
  WDMCPL_ALWAYS_ASSERT(rc != nullptr);
  auto* reverse_classification =
    reinterpret_cast<const wdmcpl::ReverseClassificationVertex*>(rc);
  WDMCPL_ALWAYS_ASSERT(reverse_classification != nullptr);
  switch (data_type) {
    case WDMCPL_DOUBLE:
      wdmcpl_create_xgc_field_adapter_t<double>(
        name, comm, data, size, *reverse_classification, in_overlap, *field_adapter);
      break;
    case WDMCPL_FLOAT:
      wdmcpl_create_xgc_field_adapter_t<float>(
        name, comm, data, size, *reverse_classification, in_overlap, *field_adapter);
      break;
    case WDMCPL_INT:
      wdmcpl_create_xgc_field_adapter_t<int>(
        name, comm, data, size, *reverse_classification, in_overlap, *field_adapter);
      break;
    case WDMCPL_LONG_INT:
      wdmcpl_create_xgc_field_adapter_t<long int>(
        name, comm, data, size, *reverse_classification, in_overlap, *field_adapter);
      break;
    default:
      printf("tyring to create XGC adapter with invalid type! %d", data_type);
      std::abort();
  }
  return reinterpret_cast<WdmCplFieldAdapterHandle*>(field_adapter);
}
WdmCplFieldAdapterHandle* wdmcpl_create_dummy_field_adapter() {
  auto* field_adapter = new wdmcpl::FieldAdapterVariant{wdmcpl::DummyFieldAdapter{}};
  return reinterpret_cast<WdmCplFieldAdapterHandle*>(field_adapter);
}

void wdmcpl_destroy_field_adapter(WdmCplFieldAdapterHandle* adapter_handle)
{
  auto* adapter =
    reinterpret_cast<wdmcpl::FieldAdapterVariant*>(adapter_handle);
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
  WDMCPL_ALWAYS_ASSERT(reverse_classification != nullptr);
  return std::accumulate(reverse_classification->begin(),
                         reverse_classification->end(), 0,
                         [](auto current, const auto& verts) {
                           return current + verts.second.size();
                         });
}
void wdmcpl_begin_send_phase(WdmCplClientHandle* h) {
  auto* client = reinterpret_cast<wdmcpl::CouplerClient*>(h);
  WDMCPL_ALWAYS_ASSERT(client != nullptr);
  client ->BeginSendPhase();
}
void wdmcpl_end_send_phase(WdmCplClientHandle* h) {
  auto* client = reinterpret_cast<wdmcpl::CouplerClient*>(h);
  WDMCPL_ALWAYS_ASSERT(client != nullptr);
  client ->EndSendPhase();
}
void wdmcpl_begin_receive_phase(WdmCplClientHandle* h ) {
  auto* client = reinterpret_cast<wdmcpl::CouplerClient*>(h);
  WDMCPL_ALWAYS_ASSERT(client != nullptr);
  client ->BeginReceivePhase();
}
void wdmcpl_end_receive_phase(WdmCplClientHandle* h) {
  auto* client = reinterpret_cast<wdmcpl::CouplerClient*>(h);
  WDMCPL_ALWAYS_ASSERT(client != nullptr);
  client ->EndReceivePhase();
}
