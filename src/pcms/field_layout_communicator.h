#ifndef FIELD_LAYOUT_COMMUNICATOR_H_
#define FIELD_LAYOUT_COMMUNICATOR_H_

#include "field_layout.h"
#include "pcms/field_layout.h"
#include "pcms/field.h"
#include "pcms/utility/profile.h"
#include "pcms/utility/assert.h"
#include "pcms/utility/inclusive_scan.h"
#include "pcms/utility/arrays.h"
#include <redev.h>
#include <memory>

namespace pcms
{

namespace field_layout_communicator
{
struct OutMsg
{
  redev::LOs dest;
  redev::LOs offset;
};

// reverse partition is a map that has the partition rank as a key
// and the values are an vector where each entry is the index into
// the array of data to send
OutMsg ConstructOutMessage(const ReversePartitionMap2& reverse_partition);

size_t count_entries(const ReversePartitionMap2& reverse_partition);

// note this function can be parallelized by making use of the offsets
redev::LOs ConstructPermutation(const ReversePartitionMap2& reverse_partition,
                                size_t num_entries, int* length);

/**
 *
 * @param local_gids local gids are the mesh GIDs in local mesh iteration order
 * @param received_gids received GIDs are the GIDS in the order of the incomming
 * message1
 * @return permutation array such that GIDS(Permutation[i]) = msgs
 */
redev::LOs ConstructPermutation(GlobalIDView<HostMemorySpace> local_gids,
                                GlobalIDView<HostMemorySpace> received_msg,
                                EntOffsetsArray ent_offsets);

OutMsg ConstructOutMessage(int rank, int nproc,
                           const redev::InMessageLayout& in);

template <typename ItBegin, typename ItEnd>
bool HasDuplicates(ItBegin begin, ItEnd end)
{
  PCMS_FUNCTION_TIMER;
  std::sort(begin, end);
  auto it = std::adjacent_find(begin, end);
  return it != end;
}

template <typename T>
bool IsValid(std::vector<T> recv_msg)
{
  auto gids = recv_msg.begin() + 4;
  for (int i = 0; i < 4; ++i) {
    auto begin = gids + recv_msg[i];
    auto end = i + 1 < 4 ? gids + recv_msg[i + 1] : recv_msg.end();
    if (HasDuplicates(begin, end)) {
      return false;
    }
  }

  return true;
}
} // namespace field_layout_communicator

class FieldLayoutCommunicator
{
public:
  FieldLayoutCommunicator(std::string name, MPI_Comm mpi_comm,
                          redev::Redev& redev, redev::Channel& channel,
                          const FieldLayout& layout)
    : mpi_comm_(mpi_comm),
      channel_(channel),
      message_permutation_{},
      buffer_size_needs_update_{true},
      layout_(layout),
      name_{std::move(name)},
      redev_(redev)
  {
    gid_comm_ = channel.CreateComm<GO>(name_ + "_gids", mpi_comm_);
    if (mpi_comm != MPI_COMM_NULL) {
      UpdateLayout();
    } else {
      UpdateLayoutNull();
    }
  }

  Rank1View<const pcms::LO, pcms::HostMemorySpace> GetPermutationArray() const
  {
    return make_const_array_view(message_permutation_);
  }

  const std::string& GetName() const { return name_; }

  const FieldLayout& GetLayout() const { return layout_; }

  size_t GetMsgSize() const { return msg_size_; }

  redev::Channel& GetChannel() { return channel_; }

  MPI_Comm& GetMPIComm() { return mpi_comm_; }

  template <typename T>
  void SetOutMessageLayout(redev::BidirectionalComm<T>& comm)
  {
    comm.SetOutMessageLayout(out_msg_.dest, out_msg_.offset);
  }

  void UpdateLayout()
  {
    namespace flc = field_layout_communicator;

    PCMS_FUNCTION_TIMER;
    auto gids = layout_.GetGids();
    auto owned = layout_.GetOwned();
    auto ent_offsets = layout_.GetEntOffsets();
    if (redev_.GetProcessType() == redev::ProcessType::Client) {
      const ReversePartitionMap2 reverse_partition =
        layout_.GetReversePartitionMap(redev::Partition{redev_.GetPartition()});
      out_msg_ = flc::ConstructOutMessage(reverse_partition);
      gid_comm_.SetOutMessageLayout(out_msg_.dest, out_msg_.offset);
      int length;
      message_permutation_ =
        flc::ConstructPermutation(reverse_partition, gids.size(), &length);
      // use permutation array to send the gids
      msg_size_ = length;
      std::vector<pcms::GO> msg(msg_size_);
      for (size_t i = 0; i < gids.size(); ++i) {
        if (owned[i])
          msg[message_permutation_[i]] = gids[i];
      }
      for (auto& rank : reverse_partition) {
        size_t i_offsets =
          message_permutation_[rank.second.indices[0]] - ent_offsets_len;
        for (int i = 0; i < rank.second.ent_offsets.size(); ++i) {
          msg[i_offsets + i] = rank.second.ent_offsets[i];
        }
      }

      channel_.BeginSendCommunicationPhase();
      gid_comm_.Send(msg.data());
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
      out_msg_ = flc::ConstructOutMessage(rank, nproc, in_message_layout);
      // construct server permutation array
      // Verify that there are no duplicate entries in the received
      // data. Duplicate data indicates that sender is not sending data from
      // only the owned rank
      // REDEV_ALWAYS_ASSERT(IsValid(recv_gids));
      GlobalIDView<HostMemorySpace> recv_gids_view(recv_gids.data(),
                                                   recv_gids.size());

      message_permutation_ =
        flc::ConstructPermutation(gids, recv_gids_view, ent_offsets);
      msg_size_ = recv_gids.size();
    }
  }

  void UpdateLayoutNull()
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

private:
  MPI_Comm mpi_comm_;
  redev::Channel& channel_;
  std::vector<pcms::LO> message_permutation_;
  redev::BidirectionalComm<GO> gid_comm_;
  bool buffer_size_needs_update_;
  field_layout_communicator::OutMsg out_msg_;
  const FieldLayout& layout_;
  redev::Redev& redev_;
  std::string name_;
  size_t msg_size_;
};
} // namespace pcms

#endif // FIELD_LAYOUT_COMMUNICATOR_H_
