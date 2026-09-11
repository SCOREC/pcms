#include "pcms/coupler/field_layout_communicator.h"

namespace pcms
{

FieldLayoutCommunicator::FieldLayoutCommunicator(
  const std::string& name, MPI_Comm mpi_comm, redev::Redev& redev,
  redev::Channel& channel, const FieldLayout& layout, bool own_mpi_comm)
  : FieldLayoutCommunicator(name, mpi_comm, redev, channel, layout,
                            std::make_unique<GenericFieldExchangePlanner>(),
                            own_mpi_comm)
{
}

FieldLayoutCommunicator::FieldLayoutCommunicator(
  const std::string& name, MPI_Comm mpi_comm, redev::Redev& redev,
  redev::Channel& channel, const FieldLayout& layout,
  std::unique_ptr<FieldExchangePlanner> planner, bool own_mpi_comm,
  const OverlapMask* overlap_mask)
  : mpi_comm_(mpi_comm),
    channel_(channel),
    layout_(layout),
    name_(name),
    redev_(redev),
    planner_(std::move(planner)),
    overlap_mask_(overlap_mask ? std::make_unique<OverlapMask>(*overlap_mask)
                               : std::make_unique<OverlapMask>(
                                   layout.GetNumOwnedDofHolder())),
    own_mpi_comm_(own_mpi_comm)
{
  gid_comm_ = channel.CreateComm<GO>(name_ + "_gids", mpi_comm_);
  if (mpi_comm != MPI_COMM_NULL) {
    UpdateLayout();
  } else {
    UpdateLayoutNull();
  }
}

Rank1View<const pcms::LO, pcms::HostMemorySpace>
FieldLayoutCommunicator::GetPermutationArray() const
{
  return make_const_array_view(plan_.permutation);
}

const std::string& FieldLayoutCommunicator::GetName() const
{
  return name_;
}

const FieldLayout& FieldLayoutCommunicator::GetLayout() const
{
  return layout_;
}

size_t FieldLayoutCommunicator::GetMsgSize() const
{
  return plan_.msg_size;
}

redev::Channel& FieldLayoutCommunicator::GetChannel() const
{
  return channel_;
}

MPI_Comm& FieldLayoutCommunicator::GetMPIComm()
{
  return mpi_comm_;
}

void FieldLayoutCommunicator::UpdateLayout()
{
  PCMS_FUNCTION_TIMER;
  if (redev_.GetProcessType() == redev::ProcessType::Client) {
    // overlap_mask_ is always non-null (created in constructor if not provided)
    plan_ = planner_->BuildExchangePlan(
      layout_, redev::Partition{redev_.GetPartition()}, overlap_mask_.get());
    redev::LOs gid_offsets(plan_.offsets.size());
    for (size_t i = 0; i < plan_.offsets.size(); ++i) {
      gid_offsets[i] = plan_.offsets[i] + i * ent_offsets_len;
    }
    gid_comm_.SetOutMessageLayout(plan_.dest_ranks, gid_offsets);
    std::vector<GO> gid_message(plan_.msg_size +
                                plan_.dest_ranks.size() * ent_offsets_len);
    planner_->FillGidMessage(
      layout_, plan_,
      Rank1View<GO, HostMemorySpace>(gid_message.data(), gid_message.size()));

    channel_.BeginSendCommunicationPhase();
    gid_comm_.Send(gid_message.data());
    channel_.EndSendCommunicationPhase();
  } else {
    channel_.BeginReceiveCommunicationPhase();
    auto recv_gids = gid_comm_.Recv();
    channel_.EndReceiveCommunicationPhase();

    int rank, nproc;
    MPI_Comm_rank(mpi_comm_, &rank);
    MPI_Comm_size(mpi_comm_, &nproc);

    const auto in_message_layout = gid_comm_.GetInMessageLayout();
    GlobalIDView<HostMemorySpace> recv_gids_view(recv_gids.data(),
                                                 recv_gids.size());
    plan_ = planner_->BuildReceivePlan(layout_, recv_gids_view, rank, nproc,
                                       in_message_layout);
  }
}

void FieldLayoutCommunicator::UpdateLayoutNull()
{
  PCMS_FUNCTION_TIMER;
  if (redev_.GetProcessType() == redev::ProcessType::Client) {
    channel_.BeginSendCommunicationPhase();
    channel_.EndSendCommunicationPhase();
  } else {
    channel_.BeginReceiveCommunicationPhase();
    channel_.EndReceiveCommunicationPhase();
  }
}

void FieldLayoutCommunicator::SetOverlapMask(std::unique_ptr<OverlapMask> mask)
{
  overlap_mask_ = std::move(mask);
}

FieldLayoutCommunicator::~FieldLayoutCommunicator()
{
  if (own_mpi_comm_ && mpi_comm_ != MPI_COMM_NULL) {
    MPI_Comm_free(&mpi_comm_);
  }
}

} // namespace pcms
