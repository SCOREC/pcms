#ifndef WDM_COUPLING_FIELD_COMMUNICATOR_H
#define WDM_COUPLING_FIELD_COMMUNICATOR_H
#include <redev.h>
#include "wdmcpl/field.h"
namespace wdmcpl
{

namespace detail
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
  std::inclusive_scan(counts.begin(), counts.end(),
                      std::next(out.offset.begin(), 1));
  return out;
}
// note this function can be parallelized by making use of the offsets
redev::LOs ConstructPermutation(const ReversePartitionMap& reverse_partition)
{

  auto num_entries = std::transform_reduce(
    reverse_partition.begin(), reverse_partition.end(), 0, std::plus<>(),
    [](const std::pair<const LO, std::vector<LO>>& v) {
      return v.second.size();
    });
  redev::LOs permutation(num_entries);
  LO entry = 0;
  for (auto& rank : reverse_partition) {
    for (auto& idx : rank.second) {
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
redev::LOs ConstructPermutation(const std::vector<wdmcpl::GO>& local_gids,
                                const std::vector<wdmcpl::GO>& received_gids)
{
  REDEV_ALWAYS_ASSERT(local_gids.size() == received_gids.size());
  std::map<wdmcpl::GO, wdmcpl::LO> global_to_local_ids;
  for (size_t i = 0; i < local_gids.size(); ++i) {
    global_to_local_ids[local_gids[i]] = i;
  }
  redev::LOs permutation;
  permutation.reserve(local_gids.size());
  for (auto gid : received_gids) {
    permutation.push_back(global_to_local_ids[gid]);
  }
  return permutation;
}
OutMsg ConstructOutMessage(int rank, int nproc,
                           const redev::InMessageLayout& in)
{

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
  std::sort(v.begin(), v.end());
  auto it = std::adjacent_find(v.begin(), v.end());
  return it != v.end();
}
} // namespace detail
template <typename FieldAdapterT>
struct FieldCommunicator
{
  using T = typename FieldAdapterT::value_type;

public:
  FieldCommunicator(
    std::string name, redev::Redev& redev, MPI_Comm mpi_comm,
    FieldAdapterT& field_adapter,
    redev::TransportType transport_type = redev::TransportType::BP4,
    adios2::Params params = adios2::Params{{"Streaming", "On"},
                                           {"OpenTimeoutSecs", "30"}})
    : mpi_comm_(mpi_comm),
      comm_buffer_{},
      message_permutation_{},
      buffer_size_needs_update_{true},
      field_adapter_(field_adapter),
      redev_{redev},
      name_{std::move(name)}
  {
    std::string transport_name = name_;
    comm_ = redev_.CreateAdiosClient<T>(transport_name, params, transport_type);
    // set up GID comm to do setup phase and get the
    // FIXME: use  one directional comm instead of the adios bidirectional
    // comm
    transport_name = transport_name.append("_gids");
    gid_comm_ = redev_.CreateAdiosClient<wdmcpl::GO>(transport_name, params,
                                                     transport_type);
    UpdateLayout();
  }

  FieldCommunicator(const FieldCommunicator&) = delete;
  FieldCommunicator(FieldCommunicator&&) = default;
  FieldCommunicator& operator=(const FieldCommunicator&) = delete;
  FieldCommunicator& operator=(FieldCommunicator&&) = default;

  void Send()
  {
    auto n = field_adapter_.Serialize({}, {});
    REDEV_ALWAYS_ASSERT(comm_buffer_.size() == static_cast<size_t>(n));

    auto buffer = ScalarArrayView<typename decltype(comm_buffer_)::value_type,
                                  HostMemorySpace>(comm_buffer_.data(),
                                                   comm_buffer_.size());
    const auto permutation =
      ScalarArrayView<const typename decltype(message_permutation_)::value_type,
                      HostMemorySpace>(message_permutation_.data(),
                                       message_permutation_.size());

    field_adapter_.Serialize(buffer, permutation);
    comm_.Send(buffer.data_handle());
  }
  void Receive()
  {
    auto data = comm_.Recv();
    auto buffer = ScalarArrayView<T, HostMemorySpace>(data.data(), data.size());
    static_assert(std::is_same_v<T, typename decltype(data)::value_type>);
    auto permutation =
      ScalarArrayView<const typename decltype(message_permutation_)::value_type,
                      HostMemorySpace>(message_permutation_.data(),
                                       message_permutation_.size());
    // load data into the field based on user specified function/functor
    field_adapter_.Deserialize(buffer, permutation);
  }

  /** update the permutation array and buffer sizes upon mesh change
   * @WARNING this function mut be called on *both* the client and server
   * after any modifications on the client
   */
  void UpdateLayout()
  {
    auto gids = field_adapter_.GetGids();
    if (redev_.GetProcessType() == redev::ProcessType::Client) {
      const ReversePartitionMap reverse_partition =
        field_adapter_.GetReversePartitionMap(redev_.GetPartition());
      auto out_message = detail::ConstructOutMessage(reverse_partition);
      comm_.SetOutMessageLayout(out_message.dest, out_message.offset);
      gid_comm_.SetOutMessageLayout(out_message.dest, out_message.offset);
      message_permutation_ = detail::ConstructPermutation(reverse_partition);
      // use permutation array to send the gids
      std::vector<wdmcpl::GO> gid_msgs(gids.size());
      REDEV_ALWAYS_ASSERT(gids.size() == message_permutation_.size());
      for (size_t i = 0; i < gids.size(); ++i) {
        gid_msgs[message_permutation_[i]] = gids[i];
      }
      gid_comm_.Send(gid_msgs.data());
    } else {
      auto recv_gids = gid_comm_.Recv();
      int rank, nproc;
      MPI_Comm_rank(mpi_comm_, &rank);
      MPI_Comm_size(mpi_comm_, &nproc);
      // we require that the layout for the gids and the message are the same
      const auto in_message_layout = gid_comm_.GetInMessageLayout();
      auto out_message = detail::ConstructOutMessage(rank, nproc, in_message_layout);
      comm_.SetOutMessageLayout(out_message.dest, out_message.offset);
      // construct server permutation array
      // Verify that there are no duplicate entries in the received
      // data. Duplicate data indicates that sender is not sending data from
      // only the owned rank
      REDEV_ALWAYS_ASSERT(!detail::HasDuplicates(recv_gids));
      message_permutation_ = detail::ConstructPermutation(gids, recv_gids);
    }
    comm_buffer_.resize(message_permutation_.size());
  }

private:
  MPI_Comm mpi_comm_;
  std::vector<T> comm_buffer_;
  std::vector<wdmcpl::LO> message_permutation_;
  redev::BidirectionalComm<T> comm_;
  redev::BidirectionalComm<GO> gid_comm_;
  bool buffer_size_needs_update_;
  // Stored functions used for updated field
  // info/serialization/deserialization
  FieldAdapterT& field_adapter_;
  redev::Redev& redev_;
  std::string name_;
};
template <>
struct FieldCommunicator<void>
{
  void Send() {}
  void Receive() {}
};
} // namespace wdmcpl

#endif // WDM_COUPLING_FIELD_COMMUNICATOR_H
