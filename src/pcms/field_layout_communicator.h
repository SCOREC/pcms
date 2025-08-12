#ifndef FIELD_LAYOUT_COMMUNICATOR_H_
#define FIELD_LAYOUT_COMMUNICATOR_H_

#include "field_layout.h"
#include "pcms/field_layout.h"
#include "pcms/field.h"
#include "pcms/profile.h"
#include "pcms/assert.h"
#include "pcms/inclusive_scan.h"
#include "pcms/arrays.h"
#include <redev.h>

namespace pcms {

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
OutMsg ConstructOutMessage(const ReversePartitionMap2& reverse_partition)
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
    counts.push_back(rank.second.indices.size() + rank.second.ent_offsets.size());
  }
  out.offset.resize(counts.size() + 1);
  out.offset[0] = 0;
  pcms::inclusive_scan(counts.begin(), counts.end(),
                         std::next(out.offset.begin(), 1));
  return out;
}

size_t count_entries(const ReversePartitionMap2& reverse_partition)
{
  PCMS_FUNCTION_TIMER;
  size_t num_entries = 0;
  for (const auto& v : reverse_partition) {
    num_entries += v.second.indices.size();
  }
  return num_entries;
}

// note this function can be parallelized by making use of the offsets
redev::LOs ConstructPermutation(const ReversePartitionMap2& reverse_partition)
{
  PCMS_FUNCTION_TIMER;
  auto num_entries = count_entries(reverse_partition);
  redev::LOs permutation(num_entries);
  LO entry = 0;
  for (auto& rank : reverse_partition) {
    entry += 5;

    for (int e = 0; e < rank.second.ent_offsets.size() - 1; ++e) {
      int start = rank.second.ent_offsets[e];
      int end = rank.second.ent_offsets[e + 1];

      for (int i = start; i < end; ++i) {
        LO index = rank.second.indices[i];
        PCMS_ALWAYS_ASSERT(index < permutation.size());
        permutation[index] = entry++;
      }
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
                                GlobalIDView<HostMemorySpace> received_msg,
                                std::array<size_t, 5> ent_offsets)
{
  PCMS_FUNCTION_TIMER;
  std::array<std::map<pcms::GO, pcms::LO>, 4> gid_to_buffer_index;
  size_t offset = 0;
  while (true) {
    GlobalIDView<HostMemorySpace> received_offsets(received_msg.data_handle() + offset, 5);
    int length = received_offsets[received_offsets.size() - 1];
    GlobalIDView<HostMemorySpace> received_gids(
      received_msg.data_handle() + offset + 5, length);

    PCMS_ALWAYS_ASSERT(offset + 5 + length - 1 < received_msg.size());

    for (int e = 0; e < received_offsets.size() - 1; ++e) {
      size_t start = received_offsets[e];
      size_t end = received_offsets[e + 1];

      for (int i = start; i < end; ++i) {
        gid_to_buffer_index[e][received_gids[i]] = offset + 5 + i;
      }
    }

    offset += length + 5;
    if (offset >= received_msg.size()) break;
  }

  LO local_index = 0;
  redev::LOs permutation;
  permutation.reserve(local_gids.size());
  for (int e = 0; e < ent_offsets.size() - 1; ++e) {
    size_t start = ent_offsets[e];
    size_t end = ent_offsets[e + 1];

    for (int i = start; i < end; ++i) {
      permutation.push_back(gid_to_buffer_index[e][local_gids[i]]);
    }
  }

  REDEV_ALWAYS_ASSERT(permutation.size() == local_gids.size());
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
} // namespace

template <typename T>
class FieldLayoutCommunicator
{
public:
  FieldLayoutCommunicator(std::string name, MPI_Comm mpi_comm,
                          redev::Redev& redev, redev::Channel& channel,
                          FieldLayout& layout)
    : mpi_comm_(mpi_comm),
      channel_(channel),
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

  void Send(T* msg, redev::Mode mode = redev::Mode::Synchronous)
  {
    PCMS_FUNCTION_TIMER;
    comm_.Send(msg, mode);
  }

  std::vector<T> Recv(redev::Mode mode = redev::Mode::Synchronous)
  {
    PCMS_FUNCTION_TIMER;
    return comm_.Recv(mode);
  }

  Rank1View<const pcms::LO, pcms::HostMemorySpace> GetPermutationArray() const
  {
    return make_const_array_view(message_permutation_);
  }

  const redev::Channel& GetChannel() const { return channel_; }

  const FieldLayout& GetLayout() const { return layout_; }

  size_t GetMsgSize() const { return msg_size_; }

  void UpdateLayout()
  {
    using namespace field_layout_communicator;

    PCMS_FUNCTION_TIMER;
    auto gids = layout_.GetGids();
    auto ent_offsets = layout_.GetEntOffsets();
    if (redev_.GetProcessType() == redev::ProcessType::Client) {
      const ReversePartitionMap2 reverse_partition =
        layout_.GetReversePartitionMap(
          redev::Partition{redev_.GetPartition()});
      auto out_message = ConstructOutMessage(reverse_partition);
      comm_.SetOutMessageLayout(out_message.dest, out_message.offset);
      gid_comm_.SetOutMessageLayout(out_message.dest, out_message.offset);
      message_permutation_ = ConstructPermutation(reverse_partition);
      // use permutation array to send the gids
      msg_size_ = gids.size() + reverse_partition.size() * 5;
      std::vector<pcms::GO> msg(msg_size_);
      for (size_t i = 0; i < gids.size(); ++i) {
        msg[message_permutation_[i]] = gids[i];
      }
      for (auto& rank : reverse_partition) {
        size_t i_offsets = message_permutation_[rank.second.indices[0]] - 5;
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
      auto out_message = ConstructOutMessage(rank, nproc, in_message_layout);
      comm_.SetOutMessageLayout(out_message.dest, out_message.offset);
      // construct server permutation array
      // Verify that there are no duplicate entries in the received
      // data. Duplicate data indicates that sender is not sending data from
      // only the owned rank
      // REDEV_ALWAYS_ASSERT(IsValid(recv_gids));
      GlobalIDView<HostMemorySpace> recv_gids_view(recv_gids.data(),
                                                   recv_gids.size());

      message_permutation_ = ConstructPermutation(gids, recv_gids_view, ent_offsets);
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
  redev::BidirectionalComm<T> comm_;
  redev::BidirectionalComm<GO> gid_comm_;
  bool buffer_size_needs_update_;
  FieldLayout& layout_;
  redev::Redev& redev_;
  std::string name_;
  size_t msg_size_;
};
} // namespace pcms

#endif // FIELD_LAYOUT_COMMUNICATOR_H_
