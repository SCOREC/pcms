#ifndef COUPLER2_H_
#define COUPLER2_H_

#include "field.h"
#include "field_communicator2.h"
#include "field_layout.h"
#include "field_layout_communicator.h"
#include "pcms/field_layout.h"
#include "pcms/field_layout_communicator.h"
#include "pcms/field_communicator2.h"
#include "pcms/utility/assert.h"
#include "pcms/utility/common.h"
#include "pcms/utility/profile.h"
#include <memory>

namespace pcms
{

// to avoid having any redev:: types in the user interface
using ProcessType = redev::ProcessType;

class Application2
{
public:
  Application2(std::string name, MPI_Comm comm, redev::Redev& redev,
               adios2::Params params, redev::TransportType transport_type,
               std::string path)
    : mpi_comm_(comm),
      redev_(redev),
      channel_{redev_.CreateAdiosChannel(std::move(name), std::move(params),
                                         transport_type, std::move(path))}
  {
    PCMS_FUNCTION_TIMER;
  }

  const FieldLayout& AddLayout(std::string name,
                               std::unique_ptr<FieldLayout> layout);

  // FIXME should take a file path for the parameters, not take adios2 params.
  // These fields are supposed to be agnostic to adios2...
  void AddField(std::string name, OwnedFieldPtr field,
                bool participates = true);

  void SendField(const std::string& name,
                 redev::Mode mode = redev::Mode::Synchronous)
  {
    PCMS_FUNCTION_TIMER;
    PCMS_ALWAYS_ASSERT(InSendPhase());
    FieldCommunicator2Ptr& communicator =
      detail::find_or_error(name, field_communicators_);
    std::visit(
      [mode](auto& field_communicator) { field_communicator->Send(mode); },
      communicator);
  };
  void ReceiveField(const std::string& name,
                    redev::Mode mode = redev::Mode::Synchronous)
  {
    PCMS_FUNCTION_TIMER;
    PCMS_ALWAYS_ASSERT(InReceivePhase());
    std::visit(
      [mode](auto& field_communicator) { field_communicator->Receive(); },
      detail::find_or_error(name, field_communicators_));
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
  FieldLayoutCommunicator& GetLayoutCommunicator(const FieldLayout& layout);

  MPI_Comm mpi_comm_;
  redev::Redev& redev_;
  redev::Channel channel_;
  std::vector<std::unique_ptr<FieldLayout>> layouts_;
  std::vector<OwnedFieldPtr> fields_;
  // map is used rather than unordered_map because we give pointers to the
  // internal data and rehash of unordered_map can cause pointer invalidation.
  // map is less cache friendly, but pointers are not invalidated.
  std::map<std::string, FieldCommunicator2Ptr> field_communicators_;
  std::map<const FieldLayout*, std::unique_ptr<FieldLayoutCommunicator>>
    field_layout_communicators_;
};

class Coupler2
{
private:
  redev::Redev SetUpRedev(bool isServer, redev::Partition partition)
  {
    if (isServer)
      return redev::Redev(mpi_comm_, std::move(partition), ProcessType::Server);
    else
      return redev::Redev(mpi_comm_);
  }

public:
  Coupler2(std::string name, MPI_Comm comm, bool isServer,
           redev::Partition partition)
    : name_(std::move(name)),
      mpi_comm_(comm),
      redev_(SetUpRedev(isServer, std::move(partition)))
  {
    PCMS_FUNCTION_TIMER;
  }
  Application2* AddApplication(
    std::string name, std::string path = "",
    redev::TransportType transport_type = redev::TransportType::BP4,
    adios2::Params params = {{"Streaming", "On"}, {"OpenTimeoutSecs", "60"}})
  {
    PCMS_FUNCTION_TIMER;
    auto key = path + name;
    auto [it, inserted] = applications_.try_emplace(
      key, std::move(name), mpi_comm_, redev_, std::move(params),
      transport_type, std::move(path));
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
  std::map<std::string, Application2> applications_;
};

} // namespace pcms

#endif // COUPLER2_H_