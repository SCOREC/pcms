#ifndef COUPLER2_H_
#define COUPLER2_H_

#include "pcms/coupler/coupler_types.h"
#include "pcms/field/field.h"
#include "pcms/field/field_layout.h"
#include "pcms/coupler/field_layout_communicator.h"
#include "pcms/coupler/field_communicator.hpp"
#include "pcms/coupler/field_exchange_planner.h"
#include "pcms/coupler/overlap_mask.h"
#include "pcms/utility/assert.h"
#include "pcms/utility/common.h"
#include "pcms/utility/profile.h"
#include "pcms/coupler/global_communicator.h"
#include <memory>

namespace pcms
{
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

  void Send(redev::Mode mode = redev::Mode::Synchronous)
  {
    PCMS_FUNCTION_TIMER;
    comm_.Send(data_, variable_name_, mode);
  }

  void Receive(redev::Mode mode = redev::Mode::Synchronous)
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

class ApplicationComm;

template <typename T>
class FieldHandle
{
public:
  FieldHandle(ApplicationComm* app, std::string name)
    : app_(app), name_(std::move(name))
  {
  }

  void Send(redev::Mode mode = redev::Mode::Synchronous) const;
  void Receive(redev::Mode mode = redev::Mode::Synchronous) const;
  [[nodiscard]] Field<T>& GetField() const;

private:
  ApplicationComm* app_;
  std::string name_;
};
template <typename T>
class DataHandle
{
public:
  DataHandle(ApplicationComm* app, std::string name)
    : app_(app), name_(std::move(name))
  {
  }
  void Send(redev::Mode mode = redev::Mode::Synchronous) const;
  void Receive(redev::Mode mode = redev::Mode::Synchronous) const;

  [[nodiscard]] const std::string& GetName() const noexcept { return name_; }

private:
  ApplicationComm* app_;
  std::string name_;
};

// FunctionHandle is a FieldHandle that additionally carries the function space,
// so it can be evaluated / used to build transfer operators. Only AddFunction
// produces one; CreateTransfer accepts FunctionHandle (not FieldHandle), which
// makes passing a comm-only field to a transfer a compile error.
template <typename T>
class FunctionHandle : public FieldHandle<T>
{
public:
  FunctionHandle(ApplicationComm* app, std::string name,
                 std::shared_ptr<const FunctionSpace> space)
    : FieldHandle<T>(app, std::move(name)), space_(std::move(space))
  {
  }

  [[nodiscard]] const FunctionSpace& GetSpace() const noexcept
  {
    return *space_;
  }
  [[nodiscard]] std::shared_ptr<const FunctionSpace> GetSpacePtr()
    const noexcept
  {
    return space_;
  }

private:
  std::shared_ptr<const FunctionSpace> space_;
};

class ApplicationComm
{
public:
  ApplicationComm(std::string name, MPI_Comm comm, redev::Redev& redev,
                  adios2::Params params, redev::TransportType transport_type,
                  std::string path)
    : mpi_comm_(comm),
      redev_(redev),
      channel_{redev_.CreateAdiosChannel(std::move(name), std::move(params),
                                         transport_type, std::move(path))}
  {
    PCMS_FUNCTION_TIMER;
  }

  // Set the overlap mask for a field's layout by name. Must be set before the
  // AddField/AddFunction that first creates that layout's communicator.
  void SetLayoutOverlapMask(const std::string& layout_name,
                            std::unique_ptr<OverlapMask> overlap_mask);

  // Register a comm-only field, identified by its own name. The layout's
  // communicator is created on first use and shared by any later
  // field/function whose layout has the same name.
  template <typename T>
  FieldHandle<T> AddField(Field<T>&& field, bool participates = true);

  template <typename T>
  FieldHandle<T> AddField(Field<T>&& field,
                          std::unique_ptr<FieldSerializer<T>> serializer,
                          bool participates = true);

  // Register a transferable field: like AddField but retains the function space
  // (via the stored Function) and returns a FunctionHandle usable in transfers.
  template <typename T>
  FieldHandle<T> AddField(std::string name, Field<T>&& field,
                          std::unique_ptr<FieldSerializer<T>> serializer,
                          bool participates = true);
  // Registers a reference to user-owned data. The Application does not take
  // ownership; it uses the registered reference for subsequent send and receive
  // operations.
  template <typename T>
  DataHandle<T> AddData(std::string name, Rank1View<T, HostMemorySpace> data,
                        MPI_Comm mpi_comm);
  template <typename T>
  FunctionHandle<T> AddFunction(Function<T>&& function,
                                bool participates = true);

  template <typename T>
  FunctionHandle<T> AddFunction(Function<T>&& function,
                                std::unique_ptr<FieldSerializer<T>> serializer,
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
  void SendData(const std::string& name,
                redev::Mode mode = redev::Mode::Synchronous)
  {
    PCMS_FUNCTION_TIMER;
    PCMS_ALWAYS_ASSERT(InSendPhase());

    std::visit([mode](auto& data_interface) { data_interface.Send(mode); },
               detail::find_or_error(name, global_data_interfaces_));
  }

  void ReceiveData(const std::string& name,
                   redev::Mode mode = redev::Mode::Synchronous)
  {
    PCMS_FUNCTION_TIMER;
    PCMS_ALWAYS_ASSERT(InReceivePhase());

    std::visit([mode](auto& data_interface) { data_interface.Receive(mode); },
               detail::find_or_error(name, global_data_interfaces_));
  }
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

  [[nodiscard]] std::size_t GetLayoutCommunicatorCount() const noexcept
  {
    return field_layout_communicators_.size();
  }

  template <typename T>
  [[nodiscard]] Field<T>& GetField(const std::string& name);

private:
  // Returns the communicator for `layout`, creating it (MPI split + overlap
  // mask lookup + planner) on first use and reusing it for any later
  // field/function whose layout has the same name.
  FieldLayoutCommunicator& GetOrCreateLayoutCommunicator(
    const FieldLayout& layout, bool participates);

  template <typename T>
  void RegisterFieldCommunicator(Field<T>& field_obj,
                                 std::unique_ptr<FieldSerializer<T>> serializer,
                                 bool participates);

  MPI_Comm mpi_comm_;
  redev::Redev& redev_;
  redev::Channel channel_;
  std::map<std::string, FieldVariant> fields_;
  std::map<std::string, FunctionVariant> functions_;
  std::vector<std::unique_ptr<FieldLayoutCommunicator>>
    owned_field_layout_communicators_;
  // map is used rather than unordered_map because we give pointers to the
  // internal data and rehash of unordered_map can cause pointer invalidation.
  // map is less cache friendly, but pointers are not invalidated.
  std::map<std::string, FieldCommunicator2Ptr> field_communicators_;
  std::map<std::string, std::unique_ptr<FieldLayoutCommunicator>>
    field_layout_communicators_;
  std::map<std::string, std::unique_ptr<OverlapMask>> layout_overlap_masks_;
  std::map<std::string, GlobalDataVariant> global_data_interfaces_;
};

template <typename T>
DataHandle<T> ApplicationComm::AddData(std::string name,
                                       Rank1View<T, HostMemorySpace> data,
                                       MPI_Comm mpi_comm)
{
  PCMS_FUNCTION_TIMER;
  auto [it, inserted] = global_data_interfaces_.try_emplace(
    name, std::in_place_type<GlobalDataInterface<T>>, name, data, mpi_comm,
    channel_);
  if (!inserted) {
    throw pcms_error("Global data interface with this name already exists");
  }
  return DataHandle<T>{this, std::move(name)};
}

template <typename T>
void DataHandle<T>::Send(redev::Mode mode) const
{
  PCMS_ALWAYS_ASSERT(app_ != nullptr);

  app_->SendData(name_, mode);
}
template <typename T>
void DataHandle<T>::Receive(redev::Mode mode) const
{
  PCMS_ALWAYS_ASSERT(app_ != nullptr);

  app_->ReceiveData(name_, mode);
}

class CouplerComm
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
  CouplerComm(std::string name, MPI_Comm comm, bool isServer,
              redev::Partition partition)
    : name_(std::move(name)),
      mpi_comm_(comm),
      redev_(SetUpRedev(isServer, std::move(partition)))
  {
    PCMS_FUNCTION_TIMER;
  }
  ApplicationComm* AddApplication(
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
  std::map<std::string, ApplicationComm> applications_;
};

} // namespace pcms

template <typename T>
void pcms::FieldHandle<T>::Send(redev::Mode mode) const
{
  PCMS_ALWAYS_ASSERT(app_ != nullptr);
  app_->SendField(name_, mode);
}

template <typename T>
void pcms::FieldHandle<T>::Receive(redev::Mode mode) const
{
  PCMS_ALWAYS_ASSERT(app_ != nullptr);
  app_->ReceiveField(name_, mode);
}

template <typename T>
pcms::Field<T>& pcms::FieldHandle<T>::GetField() const
{
  PCMS_ALWAYS_ASSERT(app_ != nullptr);
  return app_->GetField<T>(name_);
}

template <typename T>
pcms::Field<T>& pcms::ApplicationComm::GetField(const std::string& name)
{
  auto field_it = fields_.find(name);
  if (field_it != fields_.end()) {
    auto* field = std::get_if<Field<T>>(&field_it->second);
    if (field == nullptr) {
      throw pcms_error("Field stored with different type than requested");
    }
    return *field;
  }
  auto fn_it = functions_.find(name);
  if (fn_it != functions_.end()) {
    auto* function = std::get_if<Function<T>>(&fn_it->second);
    if (function == nullptr) {
      throw pcms_error("Field stored with different type than requested");
    }
    return *function;
  }
  throw pcms_error("Field '" + name + "' not found");
}

template <typename T>
void pcms::ApplicationComm::RegisterFieldCommunicator(
  Field<T>& field_obj, std::unique_ptr<FieldSerializer<T>> serializer,
  bool participates)
{
  const std::string& name = field_obj.GetName();
  FieldLayoutCommunicator& layout_communicator =
    GetOrCreateLayoutCommunicator(field_obj.GetLayout(), participates);
  FieldCommunicator2Ptr field_communicator =
    std::make_unique<FieldCommunicator<T>>(name, layout_communicator, field_obj,
                                           std::move(serializer));
  auto [it, inserted] =
    field_communicators_.emplace(name, std::move(field_communicator));
  if (!inserted) {
    throw pcms_error("Field with this name already exists");
  }
}

template <typename T>
pcms::FieldHandle<T> pcms::ApplicationComm::AddField(Field<T>&& field,
                                                     bool participates)
{
  return AddField(std::move(field), std::make_unique<FieldSerializer<T>>(),
                  participates);
}

template <typename T>
pcms::FieldHandle<T> pcms::ApplicationComm::AddField(
  Field<T>&& field, std::unique_ptr<FieldSerializer<T>> serializer,
  bool participates)
{
  PCMS_FUNCTION_TIMER;
  std::string name = field.GetName();
  if (name.empty()) {
    throw pcms_error("Application::AddField: field has no name");
  }
  auto [field_it, field_inserted] = fields_.emplace(name, std::move(field));
  if (!field_inserted) {
    throw pcms_error("Field with this name already exists");
  }
  try {
    RegisterFieldCommunicator<T>(std::get<Field<T>>(field_it->second),
                                 std::move(serializer), participates);
  } catch (...) {
    fields_.erase(field_it);
    throw;
  }
  return FieldHandle<T>{this, std::move(name)};
}

template <typename T>
pcms::FunctionHandle<T> pcms::ApplicationComm::AddFunction(
  Function<T>&& function, bool participates)
{
  return AddFunction(std::move(function),
                     std::make_unique<FieldSerializer<T>>(), participates);
}

template <typename T>
pcms::FunctionHandle<T> pcms::ApplicationComm::AddFunction(
  Function<T>&& function, std::unique_ptr<FieldSerializer<T>> serializer,
  bool participates)
{
  PCMS_FUNCTION_TIMER;
  std::string name = function.GetName();
  if (name.empty()) {
    throw pcms_error("Application::AddFunction: function has no name");
  }
  auto space = function.GetSpacePtr();
  auto [fn_it, fn_inserted] = functions_.emplace(name, std::move(function));
  if (!fn_inserted) {
    throw pcms_error("Field with this name already exists");
  }
  try {
    RegisterFieldCommunicator<T>(std::get<Function<T>>(fn_it->second),
                                 std::move(serializer), participates);
  } catch (...) {
    functions_.erase(fn_it);
    throw;
  }
  return FunctionHandle<T>{this, std::move(name), std::move(space)};
}

#endif // COUPLER2_H_
