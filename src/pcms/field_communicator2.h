#ifndef FIELD_COMMUNICATOR2_H_
#define FIELD_COMMUNICATOR2_H_

#include "pcms/field_layout.h"
#include "pcms/field.h"
#include "pcms/profile.h"
#include "pcms/assert.h"
#include "pcms/inclusive_scan.h"
#include "pcms/arrays.h"
#include <redev.h>

namespace pcms {

namespace
{
struct OutMsg
{
  redev::LOs dest;
  redev::LOs offset;
};

// reverse partition is a map that has the partition rank as a key
// and the values are an vector where each entry is the index into
// the array of data to send
OutMsg ConstructOutMessage(const ReversePartitionMap& reverse_partition)
{
  PCMS_FUNCTION_TIMER;
  OutMsg out;
  redev::LOs counts;
  counts.reserve(reverse_partition.size());
  out.dest.clear();
  out.dest.reserve(reverse_partition.size());
  // number of entries for each rank
  for (auto& rank : reverse_partition) {
    out.dest.push_back(rank.first);
    counts.push_back(rank.second.size());
  }
  out.offset.resize(counts.size() + 1);
  out.offset[0] = 0;
  pcms::inclusive_scan(counts.begin(), counts.end(),
                         std::next(out.offset.begin(), 1));
  return out;
}

size_t count_entries(const ReversePartitionMap& reverse_partition)
{
  PCMS_FUNCTION_TIMER;
  size_t num_entries = 0;
  for (const auto& v : reverse_partition) {
    num_entries += v.second.size();
  }
  return num_entries;
}

// note this function can be parallelized by making use of the offsets
redev::LOs ConstructPermutation(const ReversePartitionMap& reverse_partition)
{
  PCMS_FUNCTION_TIMER;
  auto num_entries = count_entries(reverse_partition);
  redev::LOs permutation(num_entries);
  LO entry = 0;
  for (auto& rank : reverse_partition) {
    for (auto& idx : rank.second) {
      PCMS_ALWAYS_ASSERT(idx < num_entries);
      permutation[idx] = entry++;
    }
  }
  return permutation;
}

/**
 *
 * @param local_gids local gids are the mesh GIDs in local mesh iteration order
 * @param received_gids received GIDs are the GIDS in the order of the incomming
 * message1
 * @return permutation array such that GIDS(Permutation[i]) = msgs
 */
redev::LOs ConstructPermutation(GlobalIDView<HostMemorySpace> local_gids,
                                GlobalIDView<HostMemorySpace> received_gids)
{
  PCMS_FUNCTION_TIMER;
  REDEV_ALWAYS_ASSERT(local_gids.size() == received_gids.size());
  REDEV_ALWAYS_ASSERT(std::is_permutation(
    local_gids.data_handle(), local_gids.data_handle() + local_gids.size(),
    received_gids.data_handle()));
  std::map<pcms::GO, pcms::LO> global_to_local_ids;
  for (size_t i = 0; i < local_gids.size(); ++i) {
    global_to_local_ids[local_gids[i]] = i;
  }
  redev::LOs permutation;
  permutation.reserve(local_gids.size());
  for (size_t i = 0; i < received_gids.size(); ++i) {
    permutation.push_back(global_to_local_ids[received_gids[i]]);
  }
  return permutation;
}

OutMsg ConstructOutMessage(int rank, int nproc,
                           const redev::InMessageLayout& in)
{
  PCMS_FUNCTION_TIMER;
  REDEV_ALWAYS_ASSERT(!in.srcRanks.empty());
  // auto nAppProcs =
  // Omega_h::divide_no_remainder(in.srcRanks.size(),static_cast<size_t>(nproc));
  auto nAppProcs = in.srcRanks.size() / static_cast<size_t>(nproc);
  // build dest and offsets arrays from incoming message metadata
  redev::LOs senderDeg(nAppProcs);
  for (size_t i = 0; i < nAppProcs - 1; i++) {
    senderDeg[i] =
      in.srcRanks[(i + 1) * nproc + rank] - in.srcRanks[i * nproc + rank];
  }
  const auto totInMsgs = in.offset[rank + 1] - in.offset[rank];
  senderDeg[nAppProcs - 1] =
    totInMsgs - in.srcRanks[(nAppProcs - 1) * nproc + rank];
  OutMsg out;
  for (size_t i = 0; i < nAppProcs; i++) {
    if (senderDeg[i] > 0) {
      out.dest.push_back(i);
    }
  }
  redev::GO sum = 0;
  for (auto deg : senderDeg) { // exscan over values > 0
    if (deg > 0) {
      out.offset.push_back(sum);
      sum += deg;
    }
  }
  out.offset.push_back(sum);
  return out;
}

template <typename T>
bool HasDuplicates(std::vector<T> v)
{
  PCMS_FUNCTION_TIMER;
  std::sort(v.begin(), v.end());
  auto it = std::adjacent_find(v.begin(), v.end());
  return it != v.end();
}
} // namespace

template <typename T>
class FieldCommunicator2
{
public:
  FieldCommunicator2(std::string name, MPI_Comm mpi_comm, redev::Redev& redev,
                     redev::Channel& channel, FieldLayout& layout)
    : mpi_comm_(mpi_comm),
      channel_(channel),
      comm_buffer_{},
      message_permutation_{},
      buffer_size_needs_update_{true},
      layout_(layout),
      name_{std::move(name)},
      redev_(redev)
  {
    comm_ = channel.CreateComm<T>(name_, mpi_comm_);
    gid_comm_ = channel.CreateComm<GO>(name_ + "_gids", mpi_comm_);
    if (mpi_comm != MPI_COMM_NULL) {
      UpdateLayout();
    } else {
      UpdateLayoutNull();
    }
  }

  void Send(const FieldT<T>& field, redev::Mode mode = redev::Mode::Synchronous)
  {
    PCMS_FUNCTION_TIMER;
    PCMS_ALWAYS_ASSERT(channel_.InSendCommunicationPhase());
    PCMS_ALWAYS_ASSERT(&field.GetLayout() == &layout_);
    auto n = field.Serialize({}, {});
    REDEV_ALWAYS_ASSERT(comm_buffer_.size() == static_cast<size_t>(n));
    auto buffer = make_array_view(comm_buffer_);
    field.Serialize(buffer, make_const_array_view(message_permutation_));
    comm_.Send(buffer.data_handle(), mode);
  }

  void Receive(FieldT<T>& field)
  {
    PCMS_FUNCTION_TIMER;
    PCMS_ALWAYS_ASSERT(channel_.InReceiveCommunicationPhase());
    PCMS_ALWAYS_ASSERT(&field.GetLayout() == &layout_);
    // Current implementation requires that Receive is always called in Sync
    // mode because we make an immediate call to deserialize after a call to
    // receive.
    auto data = comm_.Recv(redev::Mode::Synchronous);
    field.Deserialize(make_const_array_view(data),
                       make_const_array_view(message_permutation_));
  }

  void UpdateLayout()
  {
    PCMS_FUNCTION_TIMER;
    auto gids = layout_.GetGids();
    if (redev_.GetProcessType() == redev::ProcessType::Client) {
      const ReversePartitionMap reverse_partition =
        layout_.GetReversePartitionMap(
          redev::Partition{redev_.GetPartition()});
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
      GlobalIDView<HostMemorySpace> recv_gids_view(recv_gids.data(), recv_gids.size());
      message_permutation_ = ConstructPermutation(gids, recv_gids_view);
    }
    comm_buffer_.resize(message_permutation_.size());
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
  std::vector<T> comm_buffer_;
  std::vector<pcms::LO> message_permutation_;
  redev::BidirectionalComm<T> comm_;
  redev::BidirectionalComm<GO> gid_comm_;
  bool buffer_size_needs_update_;
  FieldLayout& layout_;
  redev::Redev& redev_;
  std::string name_;
};
} // namespace pcms

#endif // FIELD_COMMUNICATOR2_H_
