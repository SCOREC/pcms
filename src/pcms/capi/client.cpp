#include "client.h"
#include "pcms.h"
#include "pcms/xgc_field_adapter.h"
#include <variant>
#include <redev_variant_tools.h>
// #ifdef WDMCPL_HAS_OMEGA_H
//   #include "pcms/omega_h_field.h"
// #endif
#include <fstream>
#include "pcms/xgc_reverse_classification.h"
#include "pcms/dummy_field_adapter.h"
namespace pcms
{
// Note that we have a closed set of types that can be used in the C interface
using FieldAdapterVariant =
  std::variant<std::monostate, pcms::XGCFieldAdapter<double>,
               pcms::XGCFieldAdapter<float>, pcms::XGCFieldAdapter<int>,
               pcms::XGCFieldAdapter<long>,
               pcms::DummyFieldAdapter
//#ifdef WDMCPL_HAS_OMEGA_H
//               ,
//               pcms::OmegaHFieldAdapter<double>,
//               pcms::OmegaHFieldAdapter<int>
//#endif
               >;

} // namespace pcms

[[nodiscard]] WdmCplClientHandle* pcms_create_client(const char* name,
                                                       MPI_Comm comm)
{
  auto* client = new pcms::CouplerClient(name, comm);
  return reinterpret_cast<WdmCplClientHandle*>(client);
}
void pcms_destroy_client(WdmCplClientHandle* client)
{
  if (client != nullptr)
    delete reinterpret_cast<pcms::CouplerClient*>(client);
}
WdmCplReverseClassificationHandle* pcms_load_reverse_classification(
  const char* file, MPI_Comm comm)
{
  //std::filesystem::path filepath{file};
  auto* rc = new pcms::ReverseClassificationVertex{
    pcms::ReadReverseClassificationVertex(file, comm)};
  return reinterpret_cast<WdmCplReverseClassificationHandle*>(rc);
}
void pcms_destroy_reverse_classification(
  WdmCplReverseClassificationHandle* rc)
{
  if (rc != nullptr)
    delete reinterpret_cast<pcms::ReverseClassificationVertex*>(rc);
}

struct AddFieldVariantOperators {
  AddFieldVariantOperators(const char* name, pcms::CouplerClient* client, int participates)
  : name_(name), client_(client), participates_(participates)
  {
  }

  [[nodiscard]]
  pcms::CoupledField* operator()(const std::monostate&) const noexcept { return nullptr; }
  template <typename FieldAdapter>
  [[nodiscard]]
  pcms::CoupledField* operator()(const FieldAdapter& field_adapter) const noexcept {
        return client_->AddField(name_, field_adapter, participates_);
  }

  const char* name_;
  pcms::CouplerClient* client_;
  bool participates_;
};

WdmCplFieldHandle* pcms_add_field(WdmCplClientHandle* client_handle,
                                    const char* name,
                                    WdmCplFieldAdapterHandle* adapter_handle,
                                    int participates)
{

  auto* adapter =
    reinterpret_cast<pcms::FieldAdapterVariant*>(adapter_handle);
  auto* client = reinterpret_cast<pcms::CouplerClient*>(client_handle);
  WDMCPL_ALWAYS_ASSERT(client != nullptr);
  WDMCPL_ALWAYS_ASSERT(adapter != nullptr);
  //pcms::CoupledField* field = std::visit(
  //  redev::overloaded{
  //    [](const std::monostate&) -> pcms::CoupledField* { return nullptr; },
  //    [&name, &client, participates](const auto& field_adapter) {
  //      return client->AddField(name, field_adapter, participates);
  //    }},
  //  *adapter);
  pcms::CoupledField* field = std::visit(AddFieldVariantOperators{name, client, participates},*adapter);
  return reinterpret_cast<WdmCplFieldHandle*>(field);
}
void pcms_send_field_name(WdmCplClientHandle* client_handle, const char* name)
{
  auto* client = reinterpret_cast<pcms::CouplerClient*>(client_handle);
  WDMCPL_ALWAYS_ASSERT(client != nullptr);
  client->SendField(name);
}
void pcms_receive_field_name(WdmCplClientHandle* client_handle,
                               const char* name)
{
  auto* client = reinterpret_cast<pcms::CouplerClient*>(client_handle);
  WDMCPL_ALWAYS_ASSERT(client != nullptr);
  client->ReceiveField(name);
}
void pcms_send_field(WdmCplFieldHandle* field_handle)
{
  auto* field = reinterpret_cast<pcms::CoupledField*>(field_handle);
  WDMCPL_ALWAYS_ASSERT(field != nullptr);
  field->Send();
}
void pcms_receive_field(WdmCplFieldHandle* field_handle)
{
  auto* field = reinterpret_cast<pcms::CoupledField*>(field_handle);
  WDMCPL_ALWAYS_ASSERT(field != nullptr);
  field->Receive();
}
template <typename T>
void pcms_create_xgc_field_adapter_t(
  const char* name, MPI_Comm comm, void* data, int size,
  const pcms::ReverseClassificationVertex& reverse_classification,
  in_overlap_function in_overlap, pcms::FieldAdapterVariant& field_adapter)
{
  WDMCPL_ALWAYS_ASSERT((size >0) ? (data!=nullptr) : true);
  pcms::ScalarArrayView<T, pcms::HostMemorySpace> data_view(
    reinterpret_cast<T*>(data), size);
  field_adapter.emplace<pcms::XGCFieldAdapter<T>>(
    name, comm, data_view, reverse_classification, in_overlap);
}
WdmCplFieldAdapterHandle* pcms_create_xgc_field_adapter(
  const char* name, MPI_Comm comm, void* data, int size, WdmCplType data_type,
  const WdmCplReverseClassificationHandle* rc, in_overlap_function in_overlap)
{
  auto* field_adapter = new pcms::FieldAdapterVariant{};
  WDMCPL_ALWAYS_ASSERT(rc != nullptr);
  auto* reverse_classification =
    reinterpret_cast<const pcms::ReverseClassificationVertex*>(rc);
  WDMCPL_ALWAYS_ASSERT(reverse_classification != nullptr);
  switch (data_type) {
    case WDMCPL_DOUBLE:
      pcms_create_xgc_field_adapter_t<double>(
        name, comm, data, size, *reverse_classification, in_overlap, *field_adapter);
      break;
    case WDMCPL_FLOAT:
      pcms_create_xgc_field_adapter_t<float>(
        name, comm, data, size, *reverse_classification, in_overlap, *field_adapter);
      break;
    case WDMCPL_INT:
      pcms_create_xgc_field_adapter_t<int>(
        name, comm, data, size, *reverse_classification, in_overlap, *field_adapter);
      break;
    case WDMCPL_LONG_INT:
      pcms_create_xgc_field_adapter_t<long int>(
        name, comm, data, size, *reverse_classification, in_overlap, *field_adapter);
      break;
    default:
      printf("tyring to create XGC adapter with invalid type! %d", data_type);
      std::abort();
  }
  return reinterpret_cast<WdmCplFieldAdapterHandle*>(field_adapter);
}
WdmCplFieldAdapterHandle* pcms_create_dummy_field_adapter() {
  auto* field_adapter = new pcms::FieldAdapterVariant{pcms::DummyFieldAdapter{}};
  return reinterpret_cast<WdmCplFieldAdapterHandle*>(field_adapter);
}

void pcms_destroy_field_adapter(WdmCplFieldAdapterHandle* adapter_handle)
{
  auto* adapter =
    reinterpret_cast<pcms::FieldAdapterVariant*>(adapter_handle);
  if (adapter != nullptr) {
    delete adapter;
    adapter = nullptr;
  }
}
int pcms_reverse_classification_count_verts(
  WdmCplReverseClassificationHandle* rc)
{
  auto* reverse_classification =
    reinterpret_cast<const pcms::ReverseClassificationVertex*>(rc);
  WDMCPL_ALWAYS_ASSERT(reverse_classification != nullptr);
  return std::accumulate(reverse_classification->begin(),
                         reverse_classification->end(), 0,
                         [](auto current, const auto& verts) {
                           return current + verts.second.size();
                         });
}
void pcms_begin_send_phase(WdmCplClientHandle* h) {
  auto* client = reinterpret_cast<pcms::CouplerClient*>(h);
  WDMCPL_ALWAYS_ASSERT(client != nullptr);
  client ->BeginSendPhase();
}
void pcms_end_send_phase(WdmCplClientHandle* h) {
  auto* client = reinterpret_cast<pcms::CouplerClient*>(h);
  WDMCPL_ALWAYS_ASSERT(client != nullptr);
  client ->EndSendPhase();
}
void pcms_begin_receive_phase(WdmCplClientHandle* h ) {
  auto* client = reinterpret_cast<pcms::CouplerClient*>(h);
  WDMCPL_ALWAYS_ASSERT(client != nullptr);
  client ->BeginReceivePhase();
}
void pcms_end_receive_phase(WdmCplClientHandle* h) {
  auto* client = reinterpret_cast<pcms::CouplerClient*>(h);
  WDMCPL_ALWAYS_ASSERT(client != nullptr);
  client ->EndReceivePhase();
}
