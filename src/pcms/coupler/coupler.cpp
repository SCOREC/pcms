#include "pcms/coupler/coupler.hpp"

namespace pcms
{

FieldLayoutCommunicator& ApplicationComm::GetOrCreateLayoutCommunicator(
  const FieldLayout& layout, bool participates)
{
  PCMS_FUNCTION_TIMER;
  const std::string& name = layout.GetName();
  if (name.empty()) {
    throw pcms_error(
      "Application: the field's layout has no name; name it at construction "
      "(From*(..., layout_name)) before registering it with a coupler");
  }

  // Fields whose layouts share a name share one GID-exchange communicator.
  auto it = field_layout_communicators_.find(name);
  if (it != field_layout_communicators_.end()) {
    return *it->second;
  }

  MPI_Comm mpi_comm_subset = MPI_COMM_NULL;
  bool own_mpi_comm = false;
  PCMS_ALWAYS_ASSERT((mpi_comm_ == MPI_COMM_NULL) ? (!participates) : true);
  if (mpi_comm_ != MPI_COMM_NULL) {
    int rank = -1;
    MPI_Comm_rank(mpi_comm_, &rank);
    MPI_Comm_split(mpi_comm_, participates ? 0 : MPI_UNDEFINED, rank,
                   &mpi_comm_subset);
    own_mpi_comm = true;
  }

  const OverlapMask* overlap_mask = nullptr;
  auto mask_it = layout_overlap_masks_.find(name);
  if (mask_it != layout_overlap_masks_.end()) {
    overlap_mask = mask_it->second.get();
  }

  auto [it2, inserted] = field_layout_communicators_.emplace(
    name, std::make_unique<FieldLayoutCommunicator>(
            name, mpi_comm_subset, redev_, channel_, layout,
            std::make_unique<GenericFieldExchangePlanner>(), own_mpi_comm,
            overlap_mask));
  return *it2->second;
}

void ApplicationComm::SetLayoutOverlapMask(
  const std::string& layout_name, std::unique_ptr<OverlapMask> overlap_mask)
{
  layout_overlap_masks_[layout_name] = std::move(overlap_mask);
}

} // namespace pcms
