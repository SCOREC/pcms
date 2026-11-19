#ifndef PCMS_COUPLER_DETAIL_GLOBAL_DATA_INTERFACE_HPP
#define PCMS_COUPLER_DETAIL_GLOBAL_DATA_INTERFACE_HPP

#include "pcms/coupler/coupler_types.hpp"
#include "pcms/coupler/detail/global_communicator.hpp"
#include "pcms/utility/arrays.hpp"
#include "pcms/utility/detail/profile.hpp"
#include "pcms/utility/memory_spaces.hpp"
#include <mpi.h>
#include <redev.h>
#include <string>
#include <variant>

namespace pcms::detail
{

// Storage of the global data of an ApplicationComm (AddData): a host array
// sent and received whole.
template <typename T>
class GlobalDataInterface
{
public:
  GlobalDataInterface(std::string name,
                      Rank1View<T, pcms::HostMemorySpace> data,
                      MPI_Comm mpi_comm, redev::Channel& channel)
    : data_(data),
      variable_name_(name),
      mpi_comm_(mpi_comm),
      comm_(name, mpi_comm_, channel)
  {
    PCMS_FUNCTION_TIMER;
  }

  void Send(Mode mode = Mode::Synchronous)
  {
    PCMS_FUNCTION_TIMER;
    comm_.Send(data_, variable_name_, mode);
  }

  void Receive(Mode mode = Mode::Synchronous)
  {
    PCMS_FUNCTION_TIMER;
    comm_.Receive(data_, variable_name_, mode);
  }

private:
  Rank1View<T, pcms::HostMemorySpace> data_;
  std::string variable_name_;
  MPI_Comm mpi_comm_;
  GlobalCommunicator<T> comm_;
};

using GlobalDataVariant =
  std::variant<GlobalDataInterface<int8_t>, GlobalDataInterface<int32_t>,
               GlobalDataInterface<int64_t>, GlobalDataInterface<float>,
               GlobalDataInterface<double>>;

} // namespace pcms::detail

#endif // PCMS_COUPLER_DETAIL_GLOBAL_DATA_INTERFACE_HPP
