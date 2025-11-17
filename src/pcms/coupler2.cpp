#include "coupler2.h"

namespace pcms {

FieldLayoutCommunicator& Application2::GetLayoutCommunicator(
  std::string name, const FieldLayout& layout)
{
  PCMS_FUNCTION_TIMER;
  auto ptr = std::make_unique<FieldLayoutCommunicator>(name, mpi_comm_, redev_,
                                                       channel_, layout);
  auto [it, _] = field_layout_communicators_.try_emplace(&layout, std::move(ptr));
  return *it->second;
}

void Application2::AddField(std::string name, FieldPtr field, bool participates)
{
  PCMS_FUNCTION_TIMER;
  FieldCommunicator2Ptr field_communicator = std::visit(
    [this, name](auto* field_ptr) -> FieldCommunicator2Ptr {
      using T = std::remove_pointer_t<decltype(field_ptr)>::value_type;
      FieldLayoutCommunicator& layout_communicator =
        GetLayoutCommunicator(name, field_ptr->GetLayout());
      return std::make_unique<FieldCommunicator2<T>>(layout_communicator,
                                                     *field_ptr);
    },
    field);

  auto [it, inserted] = field_communicators_.insert_or_assign(name, std::move(field_communicator));
  if (!inserted) {
    std::cerr << "Field with this name" << name << "already exists!\n";
    std::terminate();
  }
}

} // namespace pcms
