#include "client.h"
#include "pcms.h"
#include "pcms/field/field_factory.h"
#include "pcms/field/function_space/xgc.h"
#include "pcms/field/data/xgc.h"
#include "pcms/field/layout/xgc.h"
#include "pcms/coupler/serializer/xgc.h"
#include "pcms/discretization/discretization/xgc_reverse_classification.h"
#include "pcms/coupler/coupler.hpp"
#include "pcms/field/layout/empty.h"
#include "pcms/field/data/simple.h"
#include "pcms/utility/assert.h"
#include <map>
#include <memory>
#include <numeric>
#include <variant>
#include "pcms/field/coordinate_systems/cylindrical.hpp"

namespace pcms
{

namespace detail
{

template <typename T>
struct XGCFieldRegistration
{
  XGCFieldFactory function_space;
  MPI_Comm plane_comm;
  Rank1View<T, HostMemorySpace> data;

  XGCFieldRegistration(XGCFieldFactory fs, MPI_Comm comm,
                       Rank1View<T, HostMemorySpace> d)
    : function_space(std::move(fs)), plane_comm(comm), data(d)
  {
  }
};

struct DummyFieldRegistration
{};

class EmptyFieldFactory : public FieldFactory
{
public:
  EmptyFieldFactory(std::string layout_name,
                    std::shared_ptr<const CoordinateSystem> system)
  {
    auto layout = std::make_shared<EmptyFieldLayout>(std::move(system));
    layout->SetName(std::move(layout_name));
    layout_ = layout;
  }

  [[nodiscard]] std::shared_ptr<const FieldLayout> GetLayout()
    const noexcept override
  {
    return layout_;
  }

protected:
  [[nodiscard]] FieldVariant CreateFieldImpl(
    Type value_type, FieldMetadata metadata) const override
  {
    return apply_to_type(
      value_type, [this, metadata](auto tag) -> FieldVariant {
        using T = typename decltype(tag)::type;
        return WrapField<T>(
          layout_, std::make_unique<SimpleFieldData<T>>(layout_, metadata));
      });
  }

  [[nodiscard]] FieldVariant CreateFieldImpl(
    FieldDataVariant data) const override
  {
    return std::visit(
      [this](auto&& fd) -> FieldVariant {
        using FD = std::decay_t<decltype(fd)>;
        using T = typename FD::element_type::value_type;
        PCMS_ALWAYS_ASSERT(fd != nullptr);
        if (dynamic_cast<const SimpleFieldData<T>*>(fd.get()) == nullptr) {
          throw pcms_error(
            "EmptyFieldFactory::CreateField: requires SimpleFieldData");
        }
        if (fd->GetDOFHolderDataHost().size() !=
            detail::ExpectedFlatFieldDataSize(*layout_)) {
          throw pcms_error(
            "EmptyFieldFactory::CreateField: field data size does not match "
            "layout");
        }
        return WrapField<T>(layout_, std::forward<decltype(fd)>(fd));
      },
      std::move(data));
  }

private:
  std::shared_ptr<const EmptyFieldLayout> layout_;
};

} // namespace detail

struct ClientState
{
  std::unique_ptr<Coupler> coupler;
  Application* app = nullptr;
  using HandleVariant = std::variant<FieldHandle<double>, FieldHandle<float>,
                                     FieldHandle<int>, FieldHandle<pcms::GO>>;
  std::map<std::string, HandleVariant> field_handles;
};

using FieldAdapterVariant = std::variant<
  std::monostate, detail::XGCFieldRegistration<double>,
  detail::XGCFieldRegistration<float>, detail::XGCFieldRegistration<int>,
  detail::XGCFieldRegistration<pcms::GO>, detail::DummyFieldRegistration>;

template <typename T>
using ClientFieldHandle = FieldHandle<T>;

inline ClientState::HandleVariant RegisterField(Application& /*app*/,
                                                std::string name,
                                                const std::monostate&, bool)
{
  throw pcms_error("pcms_add_field: field adapter for field '" + name +
                   "' was never initialized");
}

template <typename T>
ClientState::HandleVariant RegisterField(
  Application& app, std::string name,
  const detail::XGCFieldRegistration<T>& registration, bool participates)
{
  // The layout's coupling identity is the field name
  // not the arbitrary adapter name given at creation.
  registration.function_space.SetLayoutName(name);
  auto field = registration.function_space.template CreateField<T>(
    std::move(name), std::make_unique<XGCFieldData<T>>(
                       registration.function_space.GetXGCLayout(),
                       FieldMetadata{}, registration.data));
  std::unique_ptr<FieldSerializer<T>> serializer =
    std::make_unique<XGCFieldSerializer<T>>(registration.plane_comm,
                                            participates);
  return ClientState::HandleVariant{
    app.AddField(std::move(field), std::move(serializer), participates)};
}

inline ClientState::HandleVariant RegisterField(
  Application& app, std::string name, const detail::DummyFieldRegistration&,
  bool participates)
{
  // FIXME this is essentially a hack for now since the only CAPI client
  // is XGC so we need to match the XGCLayout.
  // The CAPI needs to be completely reworked for PCMS2, but that
  // should happen after we finalize new coordinate tranformations
  // and higher level field control
  auto function_space =
    detail::EmptyFieldFactory{name, csys::CylindricalRZ::Create()};
  auto field =
    function_space.CreateField<int>(std::move(name), FieldMetadata{});
  std::unique_ptr<FieldSerializer<int>> serializer =
    std::make_unique<FieldSerializer<int>>();
  return ClientState::HandleVariant{
    app.AddField(std::move(field), std::move(serializer), participates)};
}

} // namespace pcms

[[nodiscard]] PcmsClientHandle pcms_create_client(const char* name,
                                                  MPI_Comm comm)
{
  auto* client = new pcms::ClientState{};
  client->coupler =
    std::make_unique<pcms::Coupler>(name, comm, false, redev::Partition{});
  client->app = client->coupler->AddApplication(name);
  return {reinterpret_cast<void*>(client),
          reinterpret_cast<void*>(client->app)};
}

void pcms_destroy_client(PcmsClientHandle client)
{
  delete reinterpret_cast<pcms::ClientState*>(client.couplerPointer);
}

PcmsReverseClassificationHandle pcms_load_reverse_classification(
  const char* file, MPI_Comm comm)
{
  auto* rc = new pcms::ReverseClassificationVertex{
    pcms::ReadReverseClassificationVertex(file, comm)};
  return {reinterpret_cast<void*>(rc)};
}

void pcms_destroy_reverse_classification(PcmsReverseClassificationHandle rc)
{
  delete reinterpret_cast<pcms::ReverseClassificationVertex*>(rc.pointer);
}

PcmsFieldHandle pcms_add_field(PcmsClientHandle client_handle, const char* name,
                               PcmsFieldAdapterHandle adapter_handle,
                               int participates)
{
  auto* client =
    reinterpret_cast<pcms::ClientState*>(client_handle.couplerPointer);
  auto* app = reinterpret_cast<pcms::Application*>(client_handle.appPointer);
  auto* adapter =
    reinterpret_cast<pcms::FieldAdapterVariant*>(adapter_handle.pointer);
  PCMS_ALWAYS_ASSERT(client != nullptr);
  PCMS_ALWAYS_ASSERT(app != nullptr);
  PCMS_ALWAYS_ASSERT(adapter != nullptr);

  auto handle = std::visit(
    [&](const auto& registration) {
      return pcms::RegisterField(*app, name, registration, participates);
    },
    *adapter);

  auto [it, inserted] =
    client->field_handles.try_emplace(name, std::move(handle));
  if (!inserted) {
    throw pcms::pcms_error("Field with this name already exists");
  }
  return {reinterpret_cast<void*>(&it->second)};
}

void pcms_send_field_name(PcmsClientHandle client_handle, const char* name)
{
  auto* app = reinterpret_cast<pcms::Application*>(client_handle.appPointer);
  PCMS_ALWAYS_ASSERT(app != nullptr);
  app->SendField(name);
}

void pcms_receive_field_name(PcmsClientHandle client_handle, const char* name)
{
  auto* app = reinterpret_cast<pcms::Application*>(client_handle.appPointer);
  PCMS_ALWAYS_ASSERT(app != nullptr);
  app->ReceiveField(name);
}

void pcms_send_field(PcmsFieldHandle field_handle)
{
  auto* field =
    reinterpret_cast<pcms::ClientState::HandleVariant*>(field_handle.pointer);
  PCMS_ALWAYS_ASSERT(field != nullptr);
  std::visit([](auto& typed_handle) { typed_handle.Send(); }, *field);
}

void pcms_receive_field(PcmsFieldHandle field_handle)
{
  auto* field =
    reinterpret_cast<pcms::ClientState::HandleVariant*>(field_handle.pointer);
  PCMS_ALWAYS_ASSERT(field != nullptr);
  std::visit([](auto& typed_handle) { typed_handle.Receive(); }, *field);
}

template <typename T>
void pcms_create_xgc_field_adapter_t(
  const char* /* name */, MPI_Comm comm, void* data, int size,
  const pcms::ReverseClassificationVertex& reverse_classification,
  in_overlap_function in_overlap, pcms::FieldAdapterVariant& field_adapter)
{
  PCMS_ALWAYS_ASSERT((size > 0) ? (data != nullptr) : true);
  auto function_space =
    pcms::XGCFieldFactory(reverse_classification, in_overlap, size);
  pcms::Rank1View<T, pcms::HostMemorySpace> data_view(
    reinterpret_cast<T*>(data), size);
  field_adapter.emplace<pcms::detail::XGCFieldRegistration<T>>(
    std::move(function_space), comm, data_view);
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
      pcms_create_xgc_field_adapter_t<double>(name, comm, data, size,
                                              *reverse_classification,
                                              in_overlap, *field_adapter);
      break;
    case PCMS_FLOAT:
      pcms_create_xgc_field_adapter_t<float>(name, comm, data, size,
                                             *reverse_classification,
                                             in_overlap, *field_adapter);
      break;
    case PCMS_INT:
      pcms_create_xgc_field_adapter_t<int>(name, comm, data, size,
                                           *reverse_classification, in_overlap,
                                           *field_adapter);
      break;
    case PCMS_LONG_INT:
      pcms_create_xgc_field_adapter_t<pcms::GO>(name, comm, data, size,
                                                *reverse_classification,
                                                in_overlap, *field_adapter);
      break;
    default:
      PCMS_ALWAYS_ASSERT(false, MPI_COMM_WORLD,
                         "tyring to create XGC adapter with invalid type!\n");
  }
  return {reinterpret_cast<void*>(field_adapter)};
}

PcmsFieldAdapterHandle pcms_create_dummy_field_adapter()
{
  auto* field_adapter =
    new pcms::FieldAdapterVariant{pcms::detail::DummyFieldRegistration{}};
  return {reinterpret_cast<void*>(field_adapter)};
}

void pcms_destroy_field_adapter(PcmsFieldAdapterHandle adapter_handle)
{
  delete reinterpret_cast<pcms::FieldAdapterVariant*>(adapter_handle.pointer);
}

int pcms_reverse_classification_count_verts(PcmsReverseClassificationHandle rc)
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

void pcms_begin_send_phase(PcmsClientHandle h)
{
  auto* app = reinterpret_cast<pcms::Application*>(h.appPointer);
  PCMS_ALWAYS_ASSERT(app != nullptr);
  app->BeginSendPhase();
}

void pcms_end_send_phase(PcmsClientHandle h)
{
  auto* app = reinterpret_cast<pcms::Application*>(h.appPointer);
  PCMS_ALWAYS_ASSERT(app != nullptr);
  app->EndSendPhase();
}

void pcms_begin_receive_phase(PcmsClientHandle h)
{
  auto* app = reinterpret_cast<pcms::Application*>(h.appPointer);
  PCMS_ALWAYS_ASSERT(app != nullptr);
  app->BeginReceivePhase();
}

void pcms_end_receive_phase(PcmsClientHandle h)
{
  auto* app = reinterpret_cast<pcms::Application*>(h.appPointer);
  PCMS_ALWAYS_ASSERT(app != nullptr);
  app->EndReceivePhase();
}
