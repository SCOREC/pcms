#include "omega_h_field_communicator.h"
#include "pcms/assert.h"

namespace pcms {

    namespace {

    }

OmegaHFieldCommunicator::OmegaHFieldCommunicator(std::string name,
                                                 MPI_Comm mpi_comm,
                                                 redev::Redev& redev,
                                                 redev::Channel& channel,
                                                 OmegaHField2& field)
  : mpi_comm_(mpi_comm),
    channel_(channel),
    comm_buffer_{},
    message_permutation_{},
    buffer_size_needs_update_{true},
    field_(field),
    name_{std::move(name)},
    redev_(redev)
{
  comm_ = channel.CreateComm<Real>(name_, mpi_comm_);
  gid_comm_ = channel.CreateComm<GO>(name_ + "_gids", mpi_comm_);
  if (mpi_comm != MPI_COMM_NULL) {
    UpdateLayout();
  } else {
    UpdateLayoutNull();
  }
}

void OmegaHFieldCommunicator::Send(redev::Mode mode)
{
  PCMS_ALWAYS_ASSERT(channel_.InSendCommunicationPhase());
  auto n = field_.Serialize({}, {});
  REDEV_ALWAYS_ASSERT(comm_buffer_.size() == static_cast<size_t>(n));
  auto buffer = make_array_view(comm_buffer_);
  field_.Serialize(buffer, make_const_array_view(message_permutation_));
  comm_.Send(buffer.data_handle(), mode);
}

void OmegaHFieldCommunicator::Receive()
{
  PCMS_ALWAYS_ASSERT(channel_.InReceiveCommunicationPhase());
  // Current implementation requires that Receive is always called in Sync
  // mode because we make an immediate call to deserialize after a call to
  // receive.
  auto data = comm_.Recv(redev::Mode::Synchronous);
  field_.Deserialize(make_const_array_view(data),
                     make_const_array_view(message_permutation_));
}

void OmegaHFieldCommunicator::UpdateLayout()
{
  auto gids = field_.GetGids();
  if (redev_.GetProcessType() == redev::ProcessType::Client) {
    const ReversePartitionMap reverse_partition =
      field_.GetLayout().GetReversePartitionMap(Partition{redev_.GetPartition()});
    auto out_message = ConstructOutMessage(reverse_partition);
    comm_.SetOutMessageLayout(out_message.dest, out_message.offset);
    gid_comm_.SetOutMessageLayout(out_message.dest, out_message.offset);
    message_permutation_ = ConstructPermutation(reverse_partition);
    // use permutation array to send the gids
    std::vector<pcms::GO> gid_msgs(gids.size());
    REDEV_ALWAYS_ASSERT(gids.size() == message_permutation_.size());
    for (size_t i = 0; i < gids.size(); ++i) {
      gid_msgs[message_permutation_[i]] = gids[i];
    }
    channel_.BeginSendCommunicationPhase();
    gid_comm_.Send(gid_msgs.data());
    channel_.EndSendCommunicationPhase();
  } else {
    channel_.BeginReceiveCommunicationPhase();
    auto recv_gids = gid_comm_.Recv();
    channel_.EndReceiveCommunicationPhase();
    int rank, nproc;
    MPI_Comm_rank(mpi_comm_, &rank);
    MPI_Comm_size(mpi_comm_, &nproc);
    // we require that the layout for the gids and the message are the same
    const auto in_message_layout = gid_comm_.GetInMessageLayout();
    auto out_message = ConstructOutMessage(rank, nproc, in_message_layout);
    comm_.SetOutMessageLayout(out_message.dest, out_message.offset);
    // construct server permutation array
    // Verify that there are no duplicate entries in the received
    // data. Duplicate data indicates that sender is not sending data from
    // only the owned rank
    REDEV_ALWAYS_ASSERT(!HasDuplicates(recv_gids));
    message_permutation_ = ConstructPermutation(gids, recv_gids);
  }
  comm_buffer_.resize(message_permutation_.size());
}

void OmegaHFieldCommunicator::UpdateLayoutNull() {
    if (redev_.GetProcessType() == redev::ProcessType::Client) {
      channel_.BeginSendCommunicationPhase();
      channel_.EndSendCommunicationPhase();
    } else {
      channel_.BeginReceiveCommunicationPhase();
      channel_.EndReceiveCommunicationPhase();
    }
}
}
