#ifndef WDMCPL_H_
#define WDMCPL_H_
#include <mpi.h>
#include <redev.h>
#include "wdmcpl/coordinate_systems.h"
#include <unordered_map>
#include "wdmcpl/assert.h"
#include "wdmcpl/external/span.h"
#include "wdmcpl/types.h"
#include <variant>
#include <numeric>
#include <functional>

namespace wdmcpl
{
using ProcessType = redev::ProcessType;
/**
 * Key: result of partition object i.e. rank that the data is sent to on host
 * Value: Vector of local index (ordered)
 */
using ReversePartitionMap = std::map<wdmcpl::LO, std::vector<wdmcpl::LO>>;

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
    reverse_partition.begin(), reverse_partition.end(), 0, std::plus<LO>(),
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
  for (wdmcpl::LO i = 0; i < local_gids.size(); ++i) {
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

  REDEV_ALWAYS_ASSERT(in.srcRanks.size() > 0);
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

/**
 * serializer the serializer is an invocable that takes a Field* and a
 * span of one of the wdmcpl data types. It is used in a two pass algorithm,
 * so the implimenter must check the size of the buffer that is passed in. If
 * the buffer has a size of 0, routine must not write any data into the
 * buffer. The serializer invocable must return the number of entries to be
 * serialized. Note, this algorithm is not guaranteed to call the first
 * counting pass.
 * */
template <typename T>
using SerializerFunction = std::function<int(std::string_view, nonstd::span<T>,
                                             nonstd::span<const wdmcpl::LO>)>;
template <typename T>
using DeserializerFunction = std::function<void(
  std::string_view, nonstd::span<const T>, nonstd::span<const wdmcpl::LO>)>;
using ReversePartitionMapFunction =
  std::function<ReversePartitionMap(std::string_view, const redev::Partition&)>;
using GlobalIDFunction =
  std::function<std::vector<wdmcpl::GO>(std::string_view)>;

class FieldCommunicator
{
public:
  virtual void Send() = 0;
  virtual void Receive() = 0;
  virtual ~FieldCommunicator() = default;
};

// Plan is to have AdiosFieldCommunicator and BaneshFieldCommunicator
// Note: This class specifically does not take reference to an instance of the
// coupler to avoid tight coupling that would otherwise ensue.
//
// TODO question for @cwsmith
// Would the communicator rather take some type of Field object either
// inheritance or type erased that user could inherit from? Nice thing about the
// function objects is that someone can pass in a lambda if they want quick and
// dirty.
template <typename T>
class FieldCommunicatorT : public FieldCommunicator
{
public:
  FieldCommunicatorT(
    std::string name, redev::Redev& redev, MPI_Comm mpi_comm,
    GlobalIDFunction global_id_function,
    ReversePartitionMapFunction reverse_partition_map_function,
    SerializerFunction<T> serializer, DeserializerFunction<T> deserializer,
    std::string_view name_prefix = "",
    redev::TransportType transport_type = redev::TransportType::BP4,
    adios2::Params params = adios2::Params{{"Streaming", "On"},
                                           {"OpenTimeoutSecs", "12"}})
    : mpi_comm_(mpi_comm),
      comm_buffer_{},
      message_permutation_{},
      buffer_size_needs_update_{true},
      reverse_partition_map_function_{
        std::move(reverse_partition_map_function)},
      global_id_function_{std::move(global_id_function)},
      serializer_{std::move(serializer)},
      deserializer_{std::move(deserializer)},
      redev_{redev},
      name_{std::move(name)}

  {

    std::string transport_name = std::string(name_prefix);
    transport_name.append("_").append(name_);
    comm_ = redev_.CreateAdiosClient<T>(transport_name, params, transport_type);
    // set up GID comm to do setup phase and get the
    // FIXME: use  one directional comm instead of the adios bidirectional comm
    transport_name = transport_name.append("_gids");
    gid_comm_ = redev_.CreateAdiosClient<wdmcpl::GO>(transport_name,
                                                    params, transport_type);
    UpdateLayout();
  }

  FieldCommunicatorT(const FieldCommunicatorT&) = delete;
  FieldCommunicatorT(FieldCommunicatorT&&) = default;
  FieldCommunicatorT& operator=(const FieldCommunicatorT&) = delete;
  FieldCommunicatorT& operator=(FieldCommunicatorT&&) = default;

  void Send() override { Send(serializer_); }
  void Receive() override { Receive(deserializer_); }
  /*
   *  Provide templated Send/Recv so that lambda/functor with local data could
   *  be used if desired
   */
  template <typename Serializer>
  void Send(Serializer&& serializer)
  {
    auto n = serializer(name_, {}, {});
    if(comm_buffer_.size() != n) {
      UpdateLayout();
      REDEV_ALWAYS_ASSERT(comm_buffer_.size() == n);
    }

    auto buffer = nonstd::span(comm_buffer_);
    const auto permutation = nonstd::span<const typename decltype(message_permutation_)::value_type>(message_permutation_);
    serializer(name_, buffer, permutation);
    comm_.Send(buffer.data());
  }
  template <typename Deserializer>
  void Receive(Deserializer&& deserializer)
  {
    auto data = comm_.Recv();
    const auto buffer =
      nonstd::span<const T>(data);
    static_assert(std::is_same_v<T, typename decltype(data)::value_type>);
    const auto permutation = nonstd::span<const typename decltype(message_permutation_)::value_type>(message_permutation_);
    // load data into the field based on user specified function/functor
    deserializer(name_, buffer, permutation);
  }

  /** update the permutation array and buffer sizes upon mesh change
   * @WARNING this function mut be called on *both* the client and server after
   * any modifications on the client
   */
  void UpdateLayout() {
    auto gids = global_id_function_(name_);
    if (redev_.GetProcessType() == redev::ProcessType::Client) {
      const ReversePartitionMap reverse_partition =
        reverse_partition_map_function_(name_, redev_.GetPartition());

      auto out_message = ConstructOutMessage(reverse_partition);
      comm_.SetOutMessageLayout(out_message.dest, out_message.offset);
      gid_comm_.SetOutMessageLayout(out_message.dest, out_message.offset);
      message_permutation_ = ConstructPermutation(reverse_partition);
      // use permutation array to send the gids
      std::vector<wdmcpl::GO> gid_msgs(gids.size());
      REDEV_ALWAYS_ASSERT(gids.size() == message_permutation_.size());
      for (int i = 0; i < gids.size(); ++i) {
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

private:
  MPI_Comm mpi_comm_;
  std::vector<T> comm_buffer_;
  std::vector<wdmcpl::LO> message_permutation_;
  redev::BidirectionalComm<T> comm_;
  redev::BidirectionalComm<T> gid_comm_;
  bool buffer_size_needs_update_;
  // Stored functions used for updated field info/serialization/deserialization
  ReversePartitionMapFunction reverse_partition_map_function_;
  GlobalIDFunction global_id_function_;
  SerializerFunction<T> serializer_;
  DeserializerFunction<T> deserializer_;
  redev::Redev& redev_;
  std::string name_;
};

class Application {
public:
  Application(std::string name, redev::Redev& redev, MPI_Comm mpi_comm) :
 name_(std::move(name)), redev_(redev), mpi_comm_(mpi_comm) {}
  // Copy is deleted to prevent user from not grabbing a reference. This class
  // is intended to be used with data stored in the coupler
  Application(const Application&) = delete;
  Application& operator=(const Application& ) = delete;
  // moving can/should be enabled. Errors due to reference to redev_. Can probably
  // wrap that in reference_wrapper to solve the compiler errors.
  Application(Application&&) = delete;
  Application& operator=(Application&&) = delete;
  template <typename T>
  FieldCommunicatorT<T>& AddField(std::string field_name, GlobalIDFunction global_id_function,
                                  ReversePartitionMapFunction reverse_partition_map_function,
                                  SerializerFunction<T> serializer,
                                  DeserializerFunction<T> deserializer)
  {
    auto [it, inserted] = fields_.template try_emplace(
      field_name,
      std::make_unique<FieldCommunicatorT<T>>(field_name,
                                              redev_,
                                              mpi_comm_,
                                              global_id_function,
                                              reverse_partition_map_function,
                                              serializer, deserializer, name_
      ));
    if (!inserted) {
      std::cerr << "Field with this name" << field_name << "already exists!\n";
      std::exit(EXIT_FAILURE);
    }
    REDEV_ALWAYS_ASSERT(it->second != nullptr);
    return static_cast<FieldCommunicatorT<T>&>(*(it->second));
  }
private:
  redev::Redev& redev_;
  MPI_Comm mpi_comm_;
  std::string name_;

  std::unordered_map<std::string, std::unique_ptr<FieldCommunicator>> fields_;
};
// Coupler needs to have both a standalone mesh definitions to setup rdv comms
// and a list of fields
// in the server it also needs sets of fields that will be combined
template <ProcessType PT>
class Coupler
{
public:
  Coupler(std::string name, MPI_Comm comm,
          redev::Partition partition)
    : name_(std::move(name)),
      mpi_comm_(comm),
      redev_({comm, std::move(partition), PT})
  {
  }
  Application& AddApplication(std::string_view name) {
    std::string prefixed_name = name_;
    prefixed_name.append("_").append(name);
    auto [it, inserted] = applications_.template try_emplace(std::string(name),
                                                             prefixed_name, redev_, mpi_comm_);
    if (!inserted) {
      std::cerr << "Application with name" << name << "already exists!\n";
      std::exit(EXIT_FAILURE);
    }
    return it->second;
  }

  // register_field sets up a rdv::BidirectionalComm on the given mesh rdv
  // object for the field
  // GlobalIDFunc gives a list of global IDs in mesh iteration order
  // rank_count_func gives the ReversePartitionMap on the clients
  // TODO: decide...function objects can become regular types if field sizes are
  // static
  // FIXME: rank_count_func is only needed on client

private:
  std::string name_;
  static constexpr ProcessType process_type_{PT};
  MPI_Comm mpi_comm_;
  redev::Redev redev_;
  std::unordered_map<std::string, Application> applications_;
};
} // namespace wdmcpl

#endif
