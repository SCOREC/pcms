#include "client.h"
#include "pcms.h"
#include "pcms/xgc_field_adapter.h"
#include <variant>
#include <redev_variant_tools.h>
// #ifdef PCMS_HAS_OMEGA_H
//   #include "pcms/omega_h_field.h"
// #endif
#include <fstream>
#include "pcms/xgc_reverse_classification.h"
#include "pcms/dummy_field_adapter.h"
#include "pcms/print.h"
namespace pcms
{
// Note that we have a closed set of types that can be used in the C interface
using FieldAdapterVariant =
  std::variant<std::monostate, pcms::XGCFieldAdapter<double>,
               pcms::XGCFieldAdapter<float>, pcms::XGCFieldAdapter<int>,
               pcms::XGCFieldAdapter<long>,
               pcms::DummyFieldAdapter
//#ifdef PCMS_HAS_OMEGA_H
//               ,
//               pcms::OmegaHFieldAdapter<double>,
//               pcms::OmegaHFieldAdapter<int>
//#endif
               >;

} // namespace pcms

[[nodiscard]] PcmsClientHandle pcms_create_client(const char* name,
                                                       MPI_Comm comm)
{
  auto* client = new pcms::CouplerClient(name, comm);
  return {reinterpret_cast<void*>(client)};
}
void pcms_destroy_client(PcmsClientHandle client)
{
  if (client.pointer != nullptr)
    delete reinterpret_cast<pcms::CouplerClient*>(client.pointer);
}
PcmsReverseClassificationHandle pcms_load_reverse_classification(
  const char* file, MPI_Comm comm)
{
  //std::filesystem::path filepath{file};
  auto* rc = new pcms::ReverseClassificationVertex{
    pcms::ReadReverseClassificationVertex(file, comm)};
  return {reinterpret_cast<void*>(rc)};
}
void pcms_destroy_reverse_classification(
  PcmsReverseClassificationHandle rc)
{
  if (rc.pointer != nullptr)
    delete reinterpret_cast<pcms::ReverseClassificationVertex*>(rc.pointer);
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

PcmsFieldHandle pcms_add_field(PcmsClientHandle client_handle,
                                    const char* name,
                                    PcmsFieldAdapterHandle adapter_handle,
                                    int participates)
{

  auto* adapter =
    reinterpret_cast<pcms::FieldAdapterVariant*>(adapter_handle.pointer);
  auto* client = reinterpret_cast<pcms::CouplerClient*>(client_handle.pointer);
  PCMS_ALWAYS_ASSERT(client != nullptr);
  PCMS_ALWAYS_ASSERT(adapter != nullptr);
  //pcms::CoupledField* field = std::visit(
  //  redev::overloaded{
  //    [](const std::monostate&) -> pcms::CoupledField* { return nullptr; },
  //    [&name, &client, participates](const auto& field_adapter) {
  //      return client->AddField(name, field_adapter, participates);
  //    }},
  //  *adapter);
  pcms::CoupledField* field = std::visit(AddFieldVariantOperators{name, client, participates},*adapter);
  return {reinterpret_cast<void*>(field)};
}
void pcms_send_field_name(PcmsClientHandle client_handle, const char* name)
{
  auto* client = reinterpret_cast<pcms::CouplerClient*>(client_handle.pointer);
  PCMS_ALWAYS_ASSERT(client != nullptr);
  client->SendField(name);
}
void pcms_receive_field_name(PcmsClientHandle client_handle,
                               const char* name)
{
  auto* client = reinterpret_cast<pcms::CouplerClient*>(client_handle.pointer);
  PCMS_ALWAYS_ASSERT(client != nullptr);
  client->ReceiveField(name);
}
void pcms_send_field(PcmsFieldHandle field_handle)
{
  auto* field = reinterpret_cast<pcms::CoupledField*>(field_handle.pointer);
  PCMS_ALWAYS_ASSERT(field != nullptr);
  field->Send();
}
void pcms_receive_field(PcmsFieldHandle field_handle)
{
  auto* field = reinterpret_cast<pcms::CoupledField*>(field_handle.pointer);
  PCMS_ALWAYS_ASSERT(field != nullptr);
  field->Receive();
}
template <typename T>
void pcms_create_xgc_field_adapter_t(
  const char* name, MPI_Comm comm, void* data, int size,
  const pcms::ReverseClassificationVertex& reverse_classification,
  in_overlap_function in_overlap, pcms::FieldAdapterVariant& field_adapter)
{
  PCMS_ALWAYS_ASSERT((size >0) ? (data!=nullptr) : true);
  pcms::ScalarArrayView<T, pcms::HostMemorySpace> data_view(
    reinterpret_cast<T*>(data), size);
  field_adapter.emplace<pcms::XGCFieldAdapter<T>>(
    name, comm, data_view, reverse_classification, in_overlap);
}
PcmsFieldAdapterHandle pcms_create_xgc_field_adapter(
  const char* name, MPI_Comm comm, void* data, int size, PcmsType data_type,
  const PcmsReverseClassificationHandle rc, in_overlap_function in_overlap)
{
  auto* field_adapter = new pcms::FieldAdapterVariant{};
  PCMS_ALWAYS_ASSERT(rc.pointer != nullptr);
  auto* reverse_classification =
    reinterpret_cast<const pcms::ReverseClassificationVertex*>(rc.pointer);
  PCMS_ALWAYS_ASSERT(reverse_classification != nullptr);
  switch (data_type) {
    case PCMS_DOUBLE:
      pcms_create_xgc_field_adapter_t<double>(
        name, comm, data, size, *reverse_classification, in_overlap, *field_adapter);
      break;
    case PCMS_FLOAT:
      pcms_create_xgc_field_adapter_t<float>(
        name, comm, data, size, *reverse_classification, in_overlap, *field_adapter);
      break;
    case PCMS_INT:
      pcms_create_xgc_field_adapter_t<int>(
        name, comm, data, size, *reverse_classification, in_overlap, *field_adapter);
      break;
    case PCMS_LONG_INT:
      pcms_create_xgc_field_adapter_t<long int>(
        name, comm, data, size, *reverse_classification, in_overlap, *field_adapter);
      break;
    default:
      pcms::pPrintError("tyring to create XGC adapter with invalid type! %d", data_type);
      std::abort();
  }
  return {reinterpret_cast<void*>(field_adapter)};
}
PcmsFieldAdapterHandle pcms_create_dummy_field_adapter() {
  auto* field_adapter = new pcms::FieldAdapterVariant{pcms::DummyFieldAdapter{}};
  return {reinterpret_cast<void*>(field_adapter)};
}

void pcms_destroy_field_adapter(PcmsFieldAdapterHandle adapter_handle)
{
  auto* adapter =
    reinterpret_cast<pcms::FieldAdapterVariant*>(adapter_handle.pointer);
  if (adapter != nullptr) {
    delete adapter;
    adapter = nullptr;
  }
}
int pcms_reverse_classification_count_verts(
  PcmsReverseClassificationHandle rc)
{
  auto* reverse_classification =
    reinterpret_cast<const pcms::ReverseClassificationVertex*>(rc.pointer);
  PCMS_ALWAYS_ASSERT(reverse_classification != nullptr);
  return std::accumulate(reverse_classification->begin(),
                         reverse_classification->end(), 0,
                         [](auto current, const auto& verts) {
                           return current + verts.second.size();
                         });
}
void pcms_begin_send_phase(PcmsClientHandle h) {
  auto* client = reinterpret_cast<pcms::CouplerClient*>(h.pointer);
  PCMS_ALWAYS_ASSERT(client != nullptr);
  client ->BeginSendPhase();
}
void pcms_end_send_phase(PcmsClientHandle h) {
  auto* client = reinterpret_cast<pcms::CouplerClient*>(h.pointer);
  PCMS_ALWAYS_ASSERT(client != nullptr);
  client ->EndSendPhase();
}
void pcms_begin_receive_phase(PcmsClientHandle h ) {
  auto* client = reinterpret_cast<pcms::CouplerClient*>(h.pointer);
  PCMS_ALWAYS_ASSERT(client != nullptr);
  client ->BeginReceivePhase();
}
void pcms_end_receive_phase(PcmsClientHandle h) {
  auto* client = reinterpret_cast<pcms::CouplerClient*>(h.pointer);
  PCMS_ALWAYS_ASSERT(client != nullptr);
  client ->EndReceivePhase();
}