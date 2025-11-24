#ifndef PCMS_GLOBAL_COMMUNICATOR_H
#define PCMS_GLOBAL_COMMUNICATOR_H
#endif // PCMS_GLOBAL_COMMUNICATOR_H

#include <redev.h>

namespace pcms
{
  using redev::Mode;
  template <typename T>
  struct GlobalCommunicator
  {
    using value_type = T;
    public:
    GlobalCommunicator(std::string name, MPI_Comm mpi_comm, redev::Channel& channel)
      : mpi_comm(mpi_comm),
        channel_(channel),
        name_(std::move(name))
    {
      PCMS_FUNCTION_TIMER;
      comm_ = channel_.CreateComm<T>(name_, mpi_comm, redev::CommType::Global );
    }
    GlobalCommunicator(const GlobalCommunicator&) = delete;
    GlobalCommunicator& operator=(const GlobalCommunicator&) = delete;
    GlobalCommunicator(GlobalCommunicator&&)= default;
    GlobalCommunicator& operator=(GlobalCommunicator&&) = default;

    void Send(T* msg, size_t msg_size, std::string VarName, Mode mode = Mode::Synchronous)
    {
      PCMS_FUNCTION_TIMER;
      PCMS_ALWAYS_ASSERT(channel_.InSendCommunicationPhase());
      comm_.SetCommParams( VarName, msg_size);
      comm_.Send(msg, VarName, mode);
    }
    std::vector<T> Receive(std::string VarName, size_t msg_size, Mode mode = Mode::Synchronous)
    {
      PCMS_FUNCTION_TIMER;
      PCMS_ALWAYS_ASSERT(channel_.InReceiveCommunicationPhase());
      comm_.SetCommParams(VarName, msg_size);
      auto data = comm_.Recv(mode);
      return data;
    }
  private:
    MPI_Comm mpi_comm;
    redev::Channel& channel_;
    std::string name_;
    redev::BidirectionalComm<T> comm_;
  };
}