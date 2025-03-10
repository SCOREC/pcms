#ifndef PCMS_COUPLING_SERVER_H
#define PCMS_COUPLING_SERVER_H
#include "pcms/common.h"
#include "pcms/field_communicator.h"
#include "pcms/adapter/omega_h/omega_h_field.h"
#include "pcms/profile.h"
#include <map>
#include <typeinfo>

namespace pcms
{
// TODO: come up with better name for this...Don't like CoupledFieldServer
// because it's necessarily tied to the Server of the xgc_coupler
class ConvertibleCoupledField
{
public:
  template <typename FieldAdapterT, typename CommT>
  ConvertibleCoupledField(const std::string& name, FieldAdapterT field_adapter,
                          FieldCommunicator<CommT> field_comm,
                          Omega_h::Mesh& internal_mesh,
                          Omega_h::Read<Omega_h::I8> internal_field_mask = {})
    : internal_field_{OmegaHField<typename FieldAdapterT::value_type,
                                  InternalCoordinateElement>(
        name + ".__internal__", internal_mesh, internal_field_mask, "", 10, 10, field_adapter.GetEntityType())}
  {
    PCMS_FUNCTION_TIMER;
    coupled_field_ = std::make_unique<CoupledFieldModel<FieldAdapterT, CommT>>(
      std::move(field_adapter), std::move(field_comm));
  }
  template <typename FieldAdapterT>
  ConvertibleCoupledField(const std::string& name, FieldAdapterT field_adapter,
                          MPI_Comm mpi_comm, redev::Redev& redev,
                          redev::Channel& channel, Omega_h::Mesh& internal_mesh,
                          Omega_h::Read<Omega_h::I8> internal_field_mask)
    : internal_field_{OmegaHField<typename FieldAdapterT::value_type,
                                  InternalCoordinateElement>(
        name + ".__internal__", internal_mesh, internal_field_mask, "", 10, 10, field_adapter.GetEntityType())}
  {
    PCMS_FUNCTION_TIMER;
    coupled_field_ =
      std::make_unique<CoupledFieldModel<FieldAdapterT, FieldAdapterT>>(
        name, std::move(field_adapter), mpi_comm, redev, channel);
  }

  void Send(Mode mode = Mode::Synchronous)
  {
    PCMS_FUNCTION_TIMER;
    coupled_field_->Send(mode);
  }
  void Receive()
  {
    PCMS_FUNCTION_TIMER;
    coupled_field_->Receive();
  }
  [[nodiscard]] InternalField& GetInternalField() noexcept
  {
    PCMS_FUNCTION_TIMER;
    return internal_field_;
  }
  [[nodiscard]] const InternalField& GetInternalField() const noexcept
  {
    PCMS_FUNCTION_TIMER;
    return internal_field_;
  }
  template <typename T>
  [[nodiscard]] T* GetFieldAdapter() const
  {
    PCMS_FUNCTION_TIMER;
    if (typeid(T) == coupled_field_->GetFieldAdapterType()) {
      auto* adapter = coupled_field_->GetFieldAdapter();
      return reinterpret_cast<T*>(adapter);
    }
    std::cerr << "Requested type does not match field adapter type\n";
    std::abort();
  }
  struct CoupledFieldConcept
  {
    virtual void Send(Mode) = 0;
    virtual void Receive() = 0;
    [[nodiscard]] virtual const std::type_info& GetFieldAdapterType()
      const noexcept = 0;
    [[nodiscard]] virtual void* GetFieldAdapter() noexcept = 0;
    virtual ~CoupledFieldConcept() = default;
  };
  template <typename FieldAdapterT, typename CommT>
  struct CoupledFieldModel final : CoupledFieldConcept
  {
    using value_type = typename FieldAdapterT::value_type;

    CoupledFieldModel(FieldAdapterT&& field_adapter,
                      FieldCommunicator<CommT>&& comm)
      : field_adapter_(std::move(field_adapter)),
        comm_(std::move(comm)),
        type_info_(typeid(FieldAdapterT))
    {
      PCMS_FUNCTION_TIMER;
    }
    CoupledFieldModel(const std::string& name, FieldAdapterT&& field_adapter,
                      MPI_Comm mpi_comm, redev::Redev& redev,
                      redev::Channel& channel)
      : field_adapter_(std::move(field_adapter)),
        comm_(FieldCommunicator<FieldAdapterT>(name, mpi_comm, redev, channel,
                                               field_adapter_)),
        type_info_(typeid(FieldAdapterT))
    {
      PCMS_FUNCTION_TIMER;
    }
    void Send(Mode mode) final
    {
      PCMS_FUNCTION_TIMER;
      comm_.Send(mode);
    };
    void Receive() final
    {
      PCMS_FUNCTION_TIMER;
      comm_.Receive();
    };
    virtual const std::type_info& GetFieldAdapterType() const noexcept
    {
      return type_info_;
    }
    virtual void* GetFieldAdapter() noexcept
    {
      return reinterpret_cast<void*>(&field_adapter_);
    };

    FieldAdapterT field_adapter_;
    FieldCommunicator<CommT> comm_;
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
  Application(std::string name, redev::Redev& rdv, MPI_Comm comm,
              redev::Redev& redev, Omega_h::Mesh& internal_mesh,
              adios2::Params params, redev::TransportType transport_type,
              std::string path)
    : mpi_comm_(comm),
      redev_(redev),
      channel_{rdv.CreateAdiosChannel(std::move(name), std::move(params),
                                      transport_type, std::move(path))},
      internal_mesh_{internal_mesh}
  {
    PCMS_FUNCTION_TIMER;
  }
  // FIXME should take a file path for the parameters, not take adios2 params.
  // These fields are supposed to be agnostic to adios2...
  template <typename FieldAdapterT>
  ConvertibleCoupledField* AddField(
    std::string name, FieldAdapterT&& field_adapter,
    Omega_h::Read<Omega_h::I8> internal_field_mask = {})
  {
    PCMS_FUNCTION_TIMER;
    auto [it, inserted] = fields_.template try_emplace(
      name, name, std::forward<FieldAdapterT>(field_adapter), mpi_comm_, redev_,
      channel_, internal_mesh_,
      internal_field_mask);
    if (!inserted) {
      std::cerr << "OHField with this name" << name << "already exists!\n";
      std::terminate();
    }
    return &(it->second);
  }
  void SendField(const std::string& name, Mode mode = Mode::Synchronous)
  {
    PCMS_FUNCTION_TIMER;
    PCMS_ALWAYS_ASSERT(InSendPhase());
    detail::find_or_error(name, fields_).Send(mode);
  };
  void ReceiveField(const std::string& name)
  {
    PCMS_FUNCTION_TIMER;
    PCMS_ALWAYS_ASSERT(InReceivePhase());
    detail::find_or_error(name, fields_).Receive();
  };
  [[nodiscard]] bool InSendPhase() const noexcept
  {
    PCMS_FUNCTION_TIMER;
    return channel_.InSendCommunicationPhase();
  }
  [[nodiscard]] bool InReceivePhase() const noexcept
  {
    PCMS_FUNCTION_TIMER;
    return channel_.InReceiveCommunicationPhase();
  }
  void BeginSendPhase()
  {
    PCMS_FUNCTION_TIMER;
    channel_.BeginSendCommunicationPhase();
  }
  void EndSendPhase()
  {
    PCMS_FUNCTION_TIMER;
    channel_.EndSendCommunicationPhase();
  }
  void BeginReceivePhase()
  {
    PCMS_FUNCTION_TIMER;
    channel_.BeginReceiveCommunicationPhase();
  }
  void EndReceivePhase()
  {
    PCMS_FUNCTION_TIMER;
    channel_.EndReceiveCommunicationPhase();
  }

  template <typename Func, typename... Args>
  auto SendPhase(const Func& func, Args&&... args)
  {
    PCMS_FUNCTION_TIMER;
    return channel_.SendPhase(func, std::forward<Args>(args)...);
  }
  template <typename Func, typename... Args>
  auto ReceivePhase(const Func& func, Args&&... args)
  {
    PCMS_FUNCTION_TIMER;
    return channel_.ReceivePhase(func, std::forward<Args>(args)...);
  }

private:
  MPI_Comm mpi_comm_;
  redev::Redev& redev_;
  redev::Channel channel_;
  // map is used rather than unordered_map because we give pointers to the
  // internal data and rehash of unordered_map can cause pointer invalidation.
  // map is less cache friendly, but pointers are not invalidated.
  std::map<std::string, ConvertibleCoupledField> fields_;
  Omega_h::Mesh& internal_mesh_;
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
    PCMS_FUNCTION_TIMER;
  }
  Application* AddApplication(
    std::string name, std::string path = "",
    redev::TransportType transport_type = redev::TransportType::BP4,
    adios2::Params params = {{"Streaming", "On"}, {"OpenTimeoutSecs", "400"}})
  {
    PCMS_FUNCTION_TIMER;
    auto key = path + name;
    auto [it, inserted] = applications_.template try_emplace(
      key, std::move(name), redev_, mpi_comm_, redev_, internal_mesh_,
      std::move(params), transport_type, std::move(path));
    if (!inserted) {
      std::cerr << "Application with name " << name << "already exists!\n";
      std::terminate();
    }
    return &(it->second);
  }

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
  std::map<std::string, Application> applications_;
  Omega_h::Mesh& internal_mesh_;
};
} // namespace pcms
#endif // PCMS_COUPLING_SERVER_H
