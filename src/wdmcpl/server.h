#ifndef WDM_COUPLING_SERVER_H
#define WDM_COUPLING_SERVER_H
#include "wdmcpl/common.h"
#include "wdmcpl/field_communicator.h"
#include "wdmcpl/omega_h_field.h"
#include <map>
#include <typeinfo>

namespace wdmcpl
{
namespace detail
{
template <typename T, typename... Args>
auto& find_or_create_internal_field(
  const std::string& key, std::map<std::string, InternalField>& internal_fields,
  Args&&... args)
{
  auto [it, inserted] = internal_fields.try_emplace(
    key, std::in_place_type<OmegaHField<T, InternalCoordinateElement>>, key,
    std::forward<Args>(args)...);
  WDMCPL_ALWAYS_ASSERT(
    (std::holds_alternative<OmegaHField<T, InternalCoordinateElement>>(
      it->second)));
  return it->second;
}
} // namespace detail
using CombinerFunction = std::function<void(
  nonstd::span<const std::reference_wrapper<InternalField>>, InternalField&)>;

// TODO: come up with better name for this...Don't like CoupledFieldServer
// because it's necessarily tied to the Server of the xgc_coupler
class ConvertibleCoupledField
{
public:
  template <typename FieldAdapterT, typename CommT>
  ConvertibleCoupledField(const std::string& name, FieldAdapterT field_adapter,
                          FieldCommunicator<CommT> field_comm,
                          Omega_h::Mesh& internal_mesh,
                          TransferOptions native_to_internal,
                          TransferOptions internal_to_native,
                          Omega_h::Read<Omega_h::I8> internal_field_mask = {})
    : internal_field_{OmegaHField<typename FieldAdapterT::value_type,
                                  InternalCoordinateElement>(
        name + ".__internal__", internal_mesh, internal_field_mask)}
  {
    coupled_field_ = std::make_unique<CoupledFieldModel<FieldAdapterT, CommT>>(
      std::move(field_adapter), std::move(field_comm),
      std::move(native_to_internal), std::move(internal_to_native));
  }
  template <typename FieldAdapterT>
  ConvertibleCoupledField(const std::string& name,
                          FieldAdapterT field_adapter, MPI_Comm mpi_comm,
                          redev::Redev& redev,
                          redev::BidirectionalChannel& channel,
                          Omega_h::Mesh& internal_mesh,
                          TransferOptions native_to_internal,
                          TransferOptions internal_to_native,
                          Omega_h::Read<Omega_h::I8> internal_field_mask)
    : internal_field_{OmegaHField<typename FieldAdapterT::value_type,
                                  InternalCoordinateElement>(
        name + ".__internal__", internal_mesh, internal_field_mask)}
  {
    coupled_field_ =
      std::make_unique<CoupledFieldModel<FieldAdapterT, FieldAdapterT>>(
        name, std::move(field_adapter), mpi_comm, redev, channel,
        std::move(native_to_internal), std::move(internal_to_native));
  }

  void Send() { coupled_field_->Send(); }
  void Receive() { coupled_field_->Receive(); }
  void SyncNativeToInternal()
  {
    coupled_field_->SyncNativeToInternal(internal_field_);
  }
  void SyncInternalToNative()
  {
    coupled_field_->SyncInternalToNative(internal_field_);
  }
  [[nodiscard]] InternalField& GetInternalField() noexcept
  {
    return internal_field_;
  }
  [[nodiscard]] const InternalField& GetInternalField() const noexcept
  {
    return internal_field_;
  }
  template <typename T>
  [[nodiscard]]
  T* GetFieldAdapter() const {
    if(typeid(T) == coupled_field_->GetFieldAdapterType()) {
      auto* adapter = coupled_field_->GetFieldAdapter();
      return reinterpret_cast<T*>(adapter);
    }
    std::cerr<<"Requested type does not match field adapter type\n";
    std::abort();
  }
  struct CoupledFieldConcept
  {
    virtual void Send() = 0;
    virtual void Receive() = 0;
    virtual void SyncNativeToInternal(InternalField&) = 0;
    virtual void SyncInternalToNative(const InternalField&) = 0;
    [[nodiscard]]
    virtual const std::type_info& GetFieldAdapterType() const noexcept = 0;
    [[nodiscard]]
    virtual void* GetFieldAdapter() noexcept = 0;
    virtual ~CoupledFieldConcept() = default;
  };
  template <typename FieldAdapterT, typename CommT>
  struct CoupledFieldModel final : CoupledFieldConcept
  {
    using value_type = typename FieldAdapterT::value_type;

    CoupledFieldModel(FieldAdapterT&& field_adapter,
                      FieldCommunicator<CommT>&& comm,
                      TransferOptions&& native_to_internal,
                      TransferOptions&& internal_to_native)
      : field_adapter_(std::move(field_adapter)),
        comm_(std::move(comm)),
        native_to_internal_(std::move(native_to_internal)),
        internal_to_native_(std::move(internal_to_native)),
        type_info_(typeid(FieldAdapterT))
    {
    }
    CoupledFieldModel(const std::string& name, FieldAdapterT&& field_adapter,
                      MPI_Comm mpi_comm, redev::Redev& redev, redev::BidirectionalChannel& channel,
                      TransferOptions&& native_to_internal,
                      TransferOptions&& internal_to_native)
      : field_adapter_(std::move(field_adapter)),
        comm_(FieldCommunicator<FieldAdapterT>(name, mpi_comm, redev, channel,field_adapter_)),
        native_to_internal_(std::move(native_to_internal)),
        internal_to_native_(std::move(internal_to_native)),
        type_info_(typeid(FieldAdapterT))
    {
    }
    void Send() final { comm_.Send(); };
    void Receive() final { comm_.Receive(); };
    void SyncNativeToInternal(InternalField& internal_field) final
    {
      ConvertFieldAdapterToOmegaH(field_adapter_, internal_field,
                                  native_to_internal_.transfer_method,
                                  native_to_internal_.evaluation_method);
    };
    void SyncInternalToNative(const InternalField& internal_field) final
    {
      ConvertOmegaHToFieldAdapter(internal_field, field_adapter_,
                                  internal_to_native_.transfer_method,
                                  internal_to_native_.evaluation_method);
    };
    virtual const std::type_info& GetFieldAdapterType() const noexcept { return type_info_; }
    virtual void* GetFieldAdapter() noexcept { return reinterpret_cast<void*>(&field_adapter_); };

    FieldAdapterT field_adapter_;
    FieldCommunicator<CommT> comm_;
    TransferOptions native_to_internal_;
    TransferOptions internal_to_native_;
    const std::type_info& type_info_;
  };

private:
  std::unique_ptr<CoupledFieldConcept> coupled_field_;
  // even though we know the type of the internal field,
  // we store it as the InternalField variant since this avoids any copies
  // This comes at the cost of a slightly larger type with need to use the get<>
  // function
  InternalField internal_field_;
};
// TODO: strategy to merge Server/CLient Application and Fields
class Application
{
public:
  Application(std::string name, redev::Redev& rdv,
              MPI_Comm comm,
              redev::Redev& redev,
              Omega_h::Mesh& internal_mesh,
              adios2::Params params, redev::TransportType transport_type,
              std::string path)
    :
      mpi_comm_(comm),
      redev_(redev),
      channel_{rdv.CreateAdiosChannel(std::move(name), std::move(params),
                                      transport_type, std::move(path))},
      internal_mesh_{internal_mesh}
  {
  }
  // FIXME should take a file path for the parameters, not take adios2 params.
  // These fields are supposed to be agnostic to adios2...
  template <typename FieldAdapterT>
  ConvertibleCoupledField* AddField(
    std::string name, FieldAdapterT&& field_adapter,
    FieldTransferMethod to_field_transfer_method,
    FieldEvaluationMethod to_field_eval_method,
    FieldTransferMethod from_field_transfer_method,
    FieldEvaluationMethod from_field_eval_method,
    Omega_h::Read<Omega_h::I8> internal_field_mask = {})
  {
    auto [it, inserted] = fields_.template try_emplace(
      name, name, std::forward<FieldAdapterT>(field_adapter),
      mpi_comm_, redev_, channel_,
      internal_mesh_,
      TransferOptions{to_field_transfer_method, to_field_eval_method},
      TransferOptions{from_field_transfer_method, from_field_eval_method},
      internal_field_mask);
    if (!inserted) {
      std::cerr << "OHField with this name" << name << "already exists!\n";
      std::terminate();
    }
    return &(it->second);
  }
  void SendField(const std::string& name)
  {
    detail::find_or_error(name, fields_).Send();
  };
  void ReceiveField(const std::string& name)
  {
    detail::find_or_error(name, fields_).Receive();
  };

private:
  MPI_Comm mpi_comm_;
  redev::Redev& redev_;
  redev::BidirectionalChannel channel_;
  // map is used rather than unordered_map because we give pointers to the
  // internal data and rehash of unordered_map can cause pointer invalidation.
  // map is less cache friendly, but pointers are not invalidated.
  std::map<std::string, ConvertibleCoupledField> fields_;
  Omega_h::Mesh& internal_mesh_;
};
class GatherOperation
{
public:
  GatherOperation(std::vector<std::reference_wrapper<ConvertibleCoupledField>>
                    fields_to_gather,
                  InternalField& combined_field, CombinerFunction combiner)
    : coupled_fields_(std::move(fields_to_gather)),
      combined_field_(combined_field),
      combiner_(std::move(combiner))
  {
    internal_fields_.reserve(coupled_fields_.size());
    std::transform(coupled_fields_.begin(), coupled_fields_.end(),
                   std::back_inserter(internal_fields_),
                   [](ConvertibleCoupledField& fld) {
                     return std::ref(fld.GetInternalField());
                   });
  }
  void Run() const
  {
    for (auto& field : coupled_fields_) {
      field.get().Receive();
      field.get().SyncNativeToInternal();
    }
    combiner_(internal_fields_, combined_field_);
  };

private:
  std::vector<std::reference_wrapper<ConvertibleCoupledField>> coupled_fields_;
  std::vector<std::reference_wrapper<InternalField>> internal_fields_;
  InternalField& combined_field_;
  CombinerFunction combiner_;
};
class ScatterOperation
{
public:
  ScatterOperation(std::vector<std::reference_wrapper<ConvertibleCoupledField>>
                     fields_to_scatter,
                   InternalField& combined_field)
    : coupled_fields_(std::move(fields_to_scatter)),
      combined_field_{combined_field}
  {

    internal_fields_.reserve(coupled_fields_.size());
    std::transform(begin(coupled_fields_), end(coupled_fields_),
                   std::back_inserter(internal_fields_),
                   [](ConvertibleCoupledField& fld) {
                     return std::ref(fld.GetInternalField());
                   });
  }
  void Run() const
  {
    // possible we may need to add a splitter operation here.
    // needed splitter(combined_field, internal_fields_);
    // for current use case, we copy the combined field
    // into application internal fields
    std::visit(
      [this](const auto& combined_field) {
        for (auto& field : coupled_fields_) {
          std::visit(
            [&](auto& internal_field) {
              constexpr bool can_copy = std::is_same_v<
                typename std::remove_reference_t<
                  std::remove_cv_t<decltype(combined_field)>>::value_type,
                typename std::remove_reference_t<
                  std::remove_cv_t<decltype(internal_field)>>::value_type>;
              if constexpr (can_copy) {
                copy_field(combined_field, internal_field);
              } else {
                interpolate_field(combined_field, internal_field);
              }
            },
            field.get().GetInternalField());
        }
      },
      combined_field_);
    for (auto& field : coupled_fields_) {
      field.get().SyncInternalToNative();
      field.get().Send();
    }
  };

private:
  std::vector<std::reference_wrapper<ConvertibleCoupledField>> coupled_fields_;
  std::vector<std::reference_wrapper<InternalField>> internal_fields_;
  InternalField& combined_field_;
};

class CouplerServer
{
public:
  CouplerServer(std::string name, MPI_Comm comm, redev::Partition partition,
                Omega_h::Mesh& mesh)
    : name_(std::move(name)),
      mpi_comm_(comm),
      redev_({comm, std::move(partition), ProcessType::Server}),
      internal_mesh_(mesh)
  {
  }
  Application* AddApplication(
    std::string name, std::string path = "",
    redev::TransportType transport_type = redev::TransportType::BP4,
    adios2::Params params = {{"Streaming", "On"}, {"OpenTimeoutSecs", "400"}})
  {
    auto key = path + name;
    auto [it, inserted] = applications_.template try_emplace(
      key, std::move(name), redev_, mpi_comm_, redev_, internal_mesh_, std::move(params),
      transport_type, std::move(path));
    if (!inserted) {
      std::cerr << "Application with name " << name << "already exists!\n";
      std::terminate();
    }
    return &(it->second);
  }

  // here we take a string, not string_view since we need to search map
  void ScatterFields(const std::string& name)
  {
    detail::find_or_error(name, scatter_operations_).Run();
  }
  // here we take a string, not string_view since we need to search map
  void GatherFields(const std::string& name)
  {
    detail::find_or_error(name, gather_operations_).Run();
  }
  template <typename CombinedFieldT = Real>
  [[nodiscard]]
  GatherOperation* AddGatherFieldsOp(
    const std::string& name,
    std::vector<std::reference_wrapper<ConvertibleCoupledField>> gather_fields,
    const std::string& internal_field_name, CombinerFunction func,
    Omega_h::Read<Omega_h::I8> mask = {}, std::string global_id_name = "")
  {
    static constexpr int search_nx = 10;
    static constexpr int search_ny = 10;

    auto& combined = detail::find_or_create_internal_field<CombinedFieldT>(
      internal_field_name, internal_fields_, internal_mesh_, mask,
      std::move(global_id_name), search_nx, search_ny);
    auto [it, inserted] = gather_operations_.template try_emplace(
      name, std::move(gather_fields), combined, std::move(func));
    if (!inserted) {
      std::cerr << "GatherOperation with this name" << name
                << "already exists!\n";
      std::terminate();
    }
    return &(it->second);
  }
 // template <typename CombinedFieldT = Real>
 // [[nodiscard]]
 // GatherOperation* AddGatherFieldsOp(
 //   const std::string& name, const std::vector<std::string>& fields_to_gather,
 //   const std::string& internal_field_name, CombinerFunction func,
 //   Omega_h::Read<Omega_h::I8> mask = {}, std::string global_id_name = "")
 // {
 //   auto gather_fields = detail::find_many_or_error(fields_to_gather, fields_);
 //   return AddGatherFieldsOp(name, std::move(gather_fields), internal_field_name,
 //                    std::forward<CombinerFunction>(func), std::move(mask),
 //                    std::move(global_id_name));
 // }
  template <typename CombinedFieldT = Real>
  [[nodiscard]]
  ScatterOperation* AddScatterFieldsOp(
    const std::string& name, const std::string& internal_field_name,
    std::vector<std::reference_wrapper<ConvertibleCoupledField>> scatter_fields,
    Omega_h::Read<Omega_h::I8> mask = {}, std::string global_id_name = "")
  {
    static constexpr int search_nx = 10;
    static constexpr int search_ny = 10;

    auto& combined = detail::find_or_create_internal_field<CombinedFieldT>(
      internal_field_name, internal_fields_, internal_mesh_, mask,
      std::move(global_id_name), search_nx, search_ny);
    auto [it, inserted] = scatter_operations_.template try_emplace(
      name, std::move(scatter_fields), combined);

    if (!inserted) {
      std::cerr << "Scatter with this name" << name << "already exists!\n";
      std::terminate();
    }
    return &(it->second);
  }
 // template <typename CombinedFieldT = Real>
 // [[nodiscard]]
 // ScatterOperation* AddScatterFieldsOp(
 //   const std::string& name, const std::string& internal_field_name,
 //   const std::vector<std::string>& fields_to_scatter,
 //   Omega_h::Read<Omega_h::I8> mask = {}, std::string global_id_name = "")
 // {
 //   auto scatter_fields =
 //     detail::find_many_or_error(fields_to_scatter, fields_);
 //   return AddScatterFieldsOp(name, internal_field_name,
 //                      std::move(scatter_fields), std::move(mask),
 //                      std::move(global_id_name));
 // }
  [[nodiscard]] const redev::Partition& GetPartition() const noexcept
  {
    return redev_.GetPartition();
  }
  [[nodiscard]] Omega_h::Mesh& GetMesh() noexcept { return internal_mesh_; }

  [[nodiscard]] const auto& GetInternalFields() const noexcept
  {
    return internal_fields_;
  }

  // TODO: consider an "advanced api" wrapper of some sort and protect direct
  // access to the internal fields with passkey idom or some other way. User
  // could get unexpected behavior if they mess with the internal field map.
  /// This function should not be used directly it is experimental for Philip
  /// to do Benesh development. Expect that it will be removed in the future.
  [[nodiscard]] auto& GetInternalFields() noexcept { return internal_fields_; }

private:
  std::string name_;
  MPI_Comm mpi_comm_;
  redev::Redev redev_;
  // xgc_coupler owns internal fields since both gather/scatter ops use these
  // these internal fields correspond to the "Combined" fields
  std::map<std::string, InternalField> internal_fields_;
  // gather and scatter operations have reference to internal fields
  std::map<std::string, ScatterOperation> scatter_operations_;
  std::map<std::string, GatherOperation> gather_operations_;
  std::map<std::string, Application> applications_;
  Omega_h::Mesh& internal_mesh_;
};
} // namespace wdmcpl
#endif // WDM_COUPLING_SERVER_H
