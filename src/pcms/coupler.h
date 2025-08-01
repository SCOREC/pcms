#ifndef PCMS_COUPLER_H
#define PCMS_COUPLER_H
#include "pcms/common.h"
#include "pcms/field_communicator.h"
#include "pcms/adapter/omega_h/omega_h_field.h"
#include "pcms/profile.h"

namespace pcms
{

class CoupledField
{
public:
  template <typename FieldAdapterT>
  CoupledField(const std::string& name, FieldAdapterT field_adapter,
               MPI_Comm mpi_comm, redev::Redev& redev, redev::Channel& channel,
               bool participates = true)
  {
    PCMS_FUNCTION_TIMER;
    MPI_Comm mpi_comm_subset = MPI_COMM_NULL;
    PCMS_ALWAYS_ASSERT((mpi_comm == MPI_COMM_NULL) ? (participates == false)
                                                   : true);
    if (mpi_comm != MPI_COMM_NULL) {
      int rank = -1;
      MPI_Comm_rank(mpi_comm, &rank);
      MPI_Comm_split(mpi_comm, participates ? 0 : MPI_UNDEFINED, rank,
                     &mpi_comm_subset);
    }
    coupled_field_ =
      std::make_unique<CoupledFieldModel<FieldAdapterT, FieldAdapterT>>(
        name, std::move(field_adapter), mpi_comm_subset, redev, channel,
        participates);
  }

  void Send(Mode mode = Mode::Synchronous)
  {
    PCMS_FUNCTION_TIMER;
    coupled_field_->Send(mode);
  }
  void Receive(Mode mode = Mode::Synchronous)
  {
    PCMS_FUNCTION_TIMER;
    coupled_field_->Receive(mode);
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
    virtual void Receive(Mode) = 0;
    [[nodiscard]] virtual const std::type_info& GetFieldAdapterType()
      const noexcept = 0;
    [[nodiscard]] virtual void* GetFieldAdapter() noexcept = 0;
    virtual ~CoupledFieldConcept() = default;
  };
  template <typename FieldAdapterT, typename CommT>
  struct CoupledFieldModel final : CoupledFieldConcept
  {
    using value_type = typename FieldAdapterT::value_type;

    CoupledFieldModel(const std::string& name, FieldAdapterT&& field_adapter,
                      MPI_Comm mpi_comm_subset, redev::Redev& redev,
                      redev::Channel& channel, bool participates)
      : mpi_comm_subset_(mpi_comm_subset),
        field_adapter_(std::move(field_adapter)),
        comm_(FieldCommunicator<CommT>(name, mpi_comm_subset_, redev, channel,
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
    void Receive(Mode mode) final
    {
      PCMS_FUNCTION_TIMER;
      comm_.Receive(mode);
    };
    virtual const std::type_info& GetFieldAdapterType() const noexcept
    {
      return type_info_;
    }
    virtual void* GetFieldAdapter() noexcept
    {
      return reinterpret_cast<void*>(&field_adapter_);
    };
    ~CoupledFieldModel()
    {
      PCMS_FUNCTION_TIMER;
      if (mpi_comm_subset_ != MPI_COMM_NULL)
        MPI_Comm_free(&mpi_comm_subset_);
    }

    MPI_Comm mpi_comm_subset_;
    FieldAdapterT field_adapter_;
    FieldCommunicator<CommT> comm_;
    const std::type_info& type_info_;
  };

private:
  std::unique_ptr<CoupledFieldConcept> coupled_field_;
};

class Application
{
public:
  Application(std::string name, MPI_Comm comm,
              redev::Redev& redev, adios2::Params params,
              redev::TransportType transport_type,
              std::string path)
    : mpi_comm_(comm),
      redev_(redev),
      channel_{redev_.CreateAdiosChannel(std::move(name), std::move(params),
                                      transport_type, std::move(path))}
  {
    PCMS_FUNCTION_TIMER;
  }
  // FIXME should take a file path for the parameters, not take adios2 params.
  // These fields are supposed to be agnostic to adios2...
  template <typename FieldAdapterT>
  CoupledField* AddField(
    std::string name, FieldAdapterT&& field_adapter, bool participates = true)
  {
    PCMS_FUNCTION_TIMER;
    auto [it, inserted] = fields_.template try_emplace(
      name, name, std::forward<FieldAdapterT>(field_adapter), mpi_comm_, redev_,
      channel_, participates);
    if (!inserted) {
      std::cerr << "Field with this name" << name << "already exists!\n";
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
  void ReceiveField(const std::string& name, Mode mode = Mode::Synchronous)
  {
    PCMS_FUNCTION_TIMER;
    PCMS_ALWAYS_ASSERT(InReceivePhase());
    detail::find_or_error(name, fields_).Receive(mode);
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
  std::map<std::string, CoupledField> fields_;
};

class Coupler
{
private:
  redev::Redev SetUpRedev(bool isServer, redev::Partition partition) {
    if (isServer)
      return redev::Redev(mpi_comm_, std::move(partition), ProcessType::Server);
    else
      return redev::Redev(mpi_comm_);
  }
public:
  Coupler(std::string name, MPI_Comm comm, bool isServer, redev::Partition partition)
    : name_(std::move(name)),
      mpi_comm_(comm),
      redev_(SetUpRedev(isServer, std::move(partition)))
  {
    PCMS_FUNCTION_TIMER;
  }
  Application* AddApplication(
    std::string name, std::string path = "",
    redev::TransportType transport_type = redev::TransportType::BP4,
    adios2::Params params = {{"Streaming", "On"}, {"OpenTimeoutSecs", "60"}})
  {
    PCMS_FUNCTION_TIMER;
    auto key = path + name;
    auto [it, inserted] = applications_.template try_emplace(
      key, std::move(name), mpi_comm_, redev_,
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

private:
  std::string name_;
  MPI_Comm mpi_comm_;
  redev::Redev redev_;
  // gather and scatter operations have reference to internal fields
  std::map<std::string, Application> applications_;
};

} // namespace pcms

#endif // PCMS_COUPLER_H
