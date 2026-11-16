#ifndef PCMS_GLOBAL_COMMUNICATOR_H
#define PCMS_GLOBAL_COMMUNICATOR_H
#include <redev.h>
#include <pcms/utility/profile.h>
#include <pcms/utility/assert.h>
namespace pcms
{
using redev::Mode;
template <typename T>
struct GlobalCommunicator
{
  using value_type = T;

public:
  GlobalCommunicator(std::string name, MPI_Comm mpi_comm,
                     redev::Channel& channel)
    : mpi_comm(mpi_comm), channel_(channel), name_(std::move(name))
  {
    PCMS_FUNCTION_TIMER;
    comm_ = channel_.CreateComm<T>(name_, mpi_comm, redev::CommType::Global);
  }
  GlobalCommunicator(const GlobalCommunicator&) = delete;
  GlobalCommunicator& operator=(const GlobalCommunicator&) = delete;
  GlobalCommunicator(GlobalCommunicator&&) = default;
  GlobalCommunicator& operator=(GlobalCommunicator&&) = default;

  void Send(Rank1View<T, HostMemorySpace> msg, std::string VarName,
            Mode mode = Mode::Synchronous)
  {
    PCMS_FUNCTION_TIMER;
    PCMS_ALWAYS_ASSERT(channel_.InSendCommunicationPhase());
    auto msg_size = static_cast<std::size_t>(msg.extent(0));
    comm_.SetCommParams(VarName, msg_size);
    comm_.Send(msg.data_handle(), mode);
  }
  void Receive(Rank1View<T, HostMemorySpace> destination, std::string VarName,
               Mode mode = Mode::Synchronous)
  {
    PCMS_FUNCTION_TIMER;
    PCMS_ALWAYS_ASSERT(channel_.InReceiveCommunicationPhase());
    auto msg_size = static_cast<std::size_t>(destination.extent(0));
    comm_.SetCommParams(VarName, msg_size);
    comm_.Recv(destination.data_handle(), msg_size, mode);
  }

private:
  MPI_Comm mpi_comm;
  redev::Channel& channel_;
  std::string name_;
  redev::BidirectionalComm<T> comm_;
};
} // namespace pcms
#endif // PCMS_GLOBAL_COMMUNICATOR_H