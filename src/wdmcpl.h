#ifndef WDMCPL_H_
#define WDMCPL_H_
#include <mpi.h>
#include <redev.h>
#include "wdmcpl/coordinate_systems.h"
#include <unordered_map>
#include "wdmcpl/assert.h"
#include "wdmcpl/types.h"
#include <variant>
#include <numeric>
#include <functional>
#include "wdmcpl/field.h"
#include "wdmcpl/memory_spaces.h"
#include "wdmcpl/transfer_field.h"
#include "wdmcpl/omega_h_field.h"

namespace wdmcpl
{
using ProcessType = redev::ProcessType;
class GatherOperation;
class ScatterOperation;
namespace detail
{
using InternalCoordinateElement = Real;

// internal field can only be one of the types supported by Omega_h
// The coordinate element for all internal fields is the same since
// all internal fields are on the same mesh
using InternalField =
  std::variant<OmegaHField<Omega_h::I8, InternalCoordinateElement>,
               OmegaHField<Omega_h::I32, InternalCoordinateElement>,
               OmegaHField<Omega_h::I64, InternalCoordinateElement>,
               OmegaHField<Omega_h::Real, InternalCoordinateElement>>;

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

template <typename FieldShimT>
struct FieldCommunicator
{
  using T = typename FieldShimT::value_type;

public:
  FieldCommunicator(
    std::string name, redev::Redev& redev, MPI_Comm mpi_comm,
    FieldShimT& field_shim,
    redev::TransportType transport_type = redev::TransportType::BP4,
    adios2::Params params = adios2::Params{{"Streaming", "On"},
                                           {"OpenTimeoutSecs", "30"}})
    : mpi_comm_(mpi_comm),
      comm_buffer_{},
      message_permutation_{},
      buffer_size_needs_update_{true},
      redev_{redev},
      name_{std::move(name)},
      field_shim_(field_shim)
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
    auto n = field_shim_.Serialize({}, {});
    REDEV_ALWAYS_ASSERT(comm_buffer_.size() == n);

    auto buffer = ScalarArrayView<typename decltype(comm_buffer_)::value_type,
                                  HostMemorySpace>(comm_buffer_.data(),
                                                   comm_buffer_.size());
    const auto permutation =
      ScalarArrayView<const typename decltype(message_permutation_)::value_type,
                      HostMemorySpace>(message_permutation_.data(),
                                       message_permutation_.size());

    field_shim_.Serialize(buffer, permutation);
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
    field_shim_.Deserialize(buffer, permutation);
  }

  /** update the permutation array and buffer sizes upon mesh change
   * @WARNING this function mut be called on *both* the client and server
   * after any modifications on the client
   */
  void UpdateLayout()
  {
    auto gids = field_shim_.GetGids();
    if (redev_.GetProcessType() == redev::ProcessType::Client) {
      const ReversePartitionMap reverse_partition =
        field_shim_.GetReversePartitionMap(redev_.GetPartition());
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
  // Stored functions used for updated field
  // info/serialization/deserialization
  FieldShimT& field_shim_;
  redev::Redev& redev_;
  std::string name_;
};
template <>
struct FieldCommunicator<void>
{
  void Send() {}
  void Receive() {}
};
// helper function for dealing with field maps
template <typename T, typename U>
auto& find_or_error(const std::string& name,
                    const std::unordered_map<T, U>& map)
{
  auto it = map.find(name);
  WDMCPL_ALWAYS_ASSERT(it != map.end());
  return it->second;
}

template <typename T, typename U>
auto& find_or_error(const std::string& name, std::unordered_map<T, U>& map)
{
  auto it = map.find(name);
  WDMCPL_ALWAYS_ASSERT(it != map.end());
  return it->second;
}
template <typename T, typename U>
auto find_many_or_error(const std::vector<T>& keys,
                        const std::unordered_map<T, U>& map)
{

  std::vector<std::reference_wrapper<U>> results;
  results.reserve(keys.size());
  std::transform(keys.begin(), keys.end(), std::back_inserter(results),
                 [&map](const std::string& key) {
                   return std::ref(detail::find_or_error(key, map));
                 });
  return results;
}
template <typename T, typename U>
auto find_many_or_error(const std::vector<T>& keys,
                        std::unordered_map<T, U>& map)
{

  std::vector<std::reference_wrapper<U>> results;
  results.reserve(keys.size());
  std::transform(keys.begin(), keys.end(), std::back_inserter(results),
                 [&map](const std::string& key) {
                   return std::ref(detail::find_or_error(key, map));
                 });
  return results;
}
template <typename T, typename... Args>
auto& find_or_create_internal_field(
  const std::string& key,
  std::unordered_map<std::string, InternalField>& internal_fields,
  Args&&... args)
{
  auto [it, inserted] = internal_fields.try_emplace(
    key, std::in_place_type<OmegaHField<T, InternalCoordinateElement>>, key,
    std::forward<Args>(args)...);
  WDMCPL_ALWAYS_ASSERT(
    (std::holds_alternative<OmegaHField<T, InternalCoordinateElement>>(
      it->second)));
  return it->second;
}
} // namespace detail

struct TransferOptions
{
  FieldTransferMethod transfer_method;
  FieldEvaluationMethod evaluation_method;
};

// TODO: come up with better name for this...Don't like CoupledFieldServer
// because it's necessarily tied to the Server of the coupler
class ConvertibleCoupledField
{
public:
  template <typename FieldShimT, typename CommT>
  ConvertibleCoupledField(const std::string& name, FieldShimT field_shim,
                          detail::FieldCommunicator<CommT> field_comm,
                          Omega_h::Mesh& internal_mesh,
                          TransferOptions native_to_internal,
                          TransferOptions internal_to_native,
                          Omega_h::Read<Omega_h::I8> internal_field_mask = {})
    : internal_field_{OmegaHField<typename FieldShimT::value_type,
                                  detail::InternalCoordinateElement>(
        name + ".__internal__", internal_mesh, internal_field_mask)}
  {
    coupled_field_ = std::make_unique<CoupledFieldModel<FieldShimT, CommT>>(
      std::move(field_shim), std::move(field_comm),
      std::move(native_to_internal), std::move(internal_to_native));
  }
  template <typename FieldShimT>
  ConvertibleCoupledField(const std::string& name, FieldShimT field_shim,
                          redev::Redev& redev, MPI_Comm mpi_comm,
                          Omega_h::Mesh& internal_mesh,
                          TransferOptions native_to_internal,
                          TransferOptions internal_to_native,
                          Omega_h::Read<Omega_h::I8> internal_field_mask = {})
    : internal_field_{OmegaHField<typename FieldShimT::value_type,
                                  detail::InternalCoordinateElement>(
        name + ".__internal__", internal_mesh, internal_field_mask)}
  {
    coupled_field_ =
      std::make_unique<CoupledFieldModel<FieldShimT, FieldShimT>>(
        name, std::move(field_shim), redev, mpi_comm,
        std::move(native_to_internal), std::move(internal_to_native));
  }

  void Send() { coupled_field_->Send(); }
  void Receive() { coupled_field_->Receive(); }
  void SyncNativeToInternal()
  {
    coupled_field_->SyncNativeToInternal(internal_field_);
  }
  void SyncInternalToNative()
  {
    coupled_field_->SyncInternalToNative(internal_field_);
  }
  [[nodiscard]] detail::InternalField& GetInternalField() noexcept
  {
    return internal_field_;
  }
  [[nodiscard]] const detail::InternalField& GetInternalField() const noexcept
  {
    return internal_field_;
  }
  struct CoupledFieldConcept
  {
    virtual void Send() = 0;
    virtual void Receive() = 0;
    virtual void SyncNativeToInternal(detail::InternalField&) = 0;
    virtual void SyncInternalToNative(const detail::InternalField&) = 0;
    virtual ~CoupledFieldConcept() = default;
  };
  template <typename FieldShimT, typename CommT>
  struct CoupledFieldModel final : CoupledFieldConcept
  {
    using value_type = typename FieldShimT::value_type;

    CoupledFieldModel(FieldShimT&& field_shim,
                      detail::FieldCommunicator<CommT>&& comm,
                      TransferOptions&& native_to_internal,
                      TransferOptions&& internal_to_native)
      : field_shim_(std::move(field_shim)),
        comm_(std::move(comm)),
        native_to_internal_(std::move(native_to_internal)),
        internal_to_native_(std::move(internal_to_native))
    {
    }
    CoupledFieldModel(const std::string& name, FieldShimT&& field_shim,
                      redev::Redev& redev, MPI_Comm mpi_comm,
                      TransferOptions&& native_to_internal,
                      TransferOptions&& internal_to_native)
      : field_shim_(std::move(field_shim)),
        comm_(detail::FieldCommunicator<FieldShimT>(name, redev, mpi_comm,
                                                    field_shim_)),
        native_to_internal_(std::move(native_to_internal)),
        internal_to_native_(std::move(internal_to_native))
    {
    }
    void Send() { comm_.Send(); };
    void Receive() { comm_.Receive(); };
    void SyncNativeToInternal(detail::InternalField& internal_field)
    {
      field_shim_.ToOmegaH(
        std::get<OmegaHField<typename FieldShimT::value_type,
                             detail::InternalCoordinateElement>>(
          internal_field),
        native_to_internal_.transfer_method,
        native_to_internal_.evaluation_method);
    };
    void SyncInternalToNative(const detail::InternalField& internal_field)
    {
      field_shim_.FromOmegaH(
        std::get<OmegaHField<typename FieldShimT::value_type,
                             detail::InternalCoordinateElement>>(
          internal_field),
        internal_to_native_.transfer_method,
        internal_to_native_.evaluation_method);
    };

    FieldShimT field_shim_;
    detail::FieldCommunicator<CommT> comm_;
    TransferOptions native_to_internal_;
    TransferOptions internal_to_native_;
  };

private:
  std::unique_ptr<CoupledFieldConcept> coupled_field_;
  // even though we know the type of the internal field,
  // we store it as the InternalField variant since this avoids any copies
  // This comes at the cost of a slightly larger type with need to use the get<>
  // function
  detail::InternalField internal_field_;
};
class CoupledField
{
public:
  template <typename FieldShimT, typename CommT>
  CoupledField(const std::string& name, FieldShimT field_shim,
               detail::FieldCommunicator<CommT> field_comm,
               Omega_h::Mesh& internal_mesh, TransferOptions native_to_internal,
               TransferOptions internal_to_native)
  {
    coupled_field_ = std::make_unique<CoupledFieldModel<FieldShimT, CommT>>(
      std::move(field_shim), std::move(field_comm));
  }
  template <typename FieldShimT>
  CoupledField(const std::string& name, FieldShimT field_shim,
               redev::Redev& redev, MPI_Comm mpi_comm)
  {

    coupled_field_ =
      std::make_unique<CoupledFieldModel<FieldShimT, FieldShimT>>(
        name, std::move(field_shim), redev, mpi_comm);
  }

  void Send() { coupled_field_->Send(); }
  void Receive() { coupled_field_->Receive(); }
  struct CoupledFieldConcept
  {
    virtual void Send() = 0;
    virtual void Receive() = 0;
    virtual ~CoupledFieldConcept() = default;
  };
  template <typename FieldShimT, typename CommT>
  struct CoupledFieldModel final : CoupledFieldConcept
  {
    using value_type = typename FieldShimT::value_type;

    CoupledFieldModel(FieldShimT&& field_shim,
                      detail::FieldCommunicator<CommT>&& comm)
      : field_shim_(std::move(field_shim)), comm_(std::move(comm))
    {
    }
    CoupledFieldModel(const std::string& name, FieldShimT&& field_shim,
                      redev::Redev& redev, MPI_Comm mpi_comm)
      : field_shim_(std::move(field_shim)),
        comm_(
          detail::FieldCommunicator<CommT>(name, redev, mpi_comm, field_shim_))
    {
    }
    void Send() { comm_.Send(); };
    void Receive() { comm_.Receive(); };

    FieldShimT field_shim_;
    detail::FieldCommunicator<CommT> comm_;
  };

private:
  std::unique_ptr<CoupledFieldConcept> coupled_field_;
};

using CombinerFunction = std::function<void(
  nonstd::span<const std::reference_wrapper<detail::InternalField>>,
  detail::InternalField&)>;
class GatherOperation
{
public:
  GatherOperation(std::vector<std::reference_wrapper<ConvertibleCoupledField>>
                    fields_to_gather,
                  detail::InternalField& combined_field,
                  CombinerFunction combiner)
    : coupled_fields_(std::move(fields_to_gather)),
      combined_field_(combined_field),
      combiner_(std::move(combiner))
  {
    internal_fields_.reserve(coupled_fields_.size());
    std::transform(coupled_fields_.begin(), coupled_fields_.end(),
                   std::back_inserter(internal_fields_),
                   [](ConvertibleCoupledField& fld) {
                     return std::ref(fld.GetInternalField());
                   });
  }
  void Run() const
  {
    for (auto& field : coupled_fields_) {
      field.get().Receive();
      field.get().SyncNativeToInternal();
    }
    combiner_(internal_fields_, combined_field_);
  };

private:
  std::vector<std::reference_wrapper<ConvertibleCoupledField>> coupled_fields_;
  std::vector<std::reference_wrapper<detail::InternalField>> internal_fields_;
  detail::InternalField& combined_field_;
  CombinerFunction combiner_;
};
class ScatterOperation
{
public:
  ScatterOperation(std::vector<std::reference_wrapper<ConvertibleCoupledField>>
                     fields_to_scatter,
                   detail::InternalField& combined_field)
    : coupled_fields_(std::move(fields_to_scatter)),
      combined_field_{combined_field}
  {

    internal_fields_.reserve(coupled_fields_.size());
    std::transform(begin(coupled_fields_), end(coupled_fields_),
                   std::back_inserter(internal_fields_),
                   [](ConvertibleCoupledField& fld) {
                     return std::ref(fld.GetInternalField());
                   });
  }
  void Run() const
  {
    // possible we may need to add a splitter operation here. will evaluate if
    // needed splitter(combined_field, internal_fields_);
    for (auto& field : coupled_fields_) {
      field.get().SyncInternalToNative();
      field.get().Send();
    }
  };

private:
  std::vector<std::reference_wrapper<ConvertibleCoupledField>> coupled_fields_;
  std::vector<std::reference_wrapper<detail::InternalField>> internal_fields_;
  detail::InternalField& combined_field_;
};

class CouplerServer
{
public:
  CouplerServer(std::string name, MPI_Comm comm, redev::Partition partition,
                Omega_h::Mesh& mesh)
    : name_(std::move(name)),
      mpi_comm_(comm),
      redev_({comm, std::move(partition), ProcessType::Server}),
      internal_mesh_(mesh)
  {
  }
  template <typename FieldShimT>
  ConvertibleCoupledField* AddField(
    std::string unique_name, FieldShimT field_shim,
    FieldTransferMethod to_field_transfer_method,
    FieldEvaluationMethod to_field_eval_method,
    FieldTransferMethod from_field_transfer_method,
    FieldEvaluationMethod from_field_eval_method,
    Omega_h::Read<Omega_h::I8> internal_field_mask = {})
  {
    auto [it, inserted] = fields_.template try_emplace(
      unique_name, unique_name, std::move(field_shim), redev_, mpi_comm_,
      internal_mesh_,
      TransferOptions{to_field_transfer_method, to_field_eval_method},
      TransferOptions{from_field_transfer_method, from_field_eval_method},
      internal_field_mask);
    if (!inserted) {
      std::cerr << "OHField with this name" << unique_name
                << "already exists!\n";
      std::terminate();
    }
    return &(it->second);
  }

  // here we take a string, not string_view since we need to search map
  void ScatterFields(const std::string& name)
  {
    detail::find_or_error(name, scatter_operations_).Run();
  }
  // here we take a string, not string_view since we need to search map
  void GatherFields(const std::string& name)
  {
    detail::find_or_error(name, gather_operations_).Run();
  }
  void SendField(const std::string& name)
  {
    detail::find_or_error(name, fields_).Send();
  };
  void ReceiveField(const std::string& name)
  {
    detail::find_or_error(name, fields_).Receive();
  };

  // TODO: refactor searchnx/seachny into search input parameter struct
  template <typename CombinedFieldT = Real>
  GatherOperation* AddGatherFieldsOp(
    const std::string& name, const std::vector<std::string>& fields_to_gather,
    const std::string& internal_field_name, CombinerFunction func,
    Omega_h::Read<Omega_h::I8> mask = {})
  {
    auto gather_fields = detail::find_many_or_error(fields_to_gather, fields_);
    static constexpr int seach_nx = 10;
    static constexpr int seach_ny = 10;

    auto& combined = detail::find_or_create_internal_field<CombinedFieldT>(
      internal_field_name, internal_fields_, internal_mesh_,mask, seach_nx, seach_ny);
    auto [it, inserted] = gather_operations_.template try_emplace(
      name, std::move(gather_fields), combined, std::move(func));
    if (!inserted) {
      std::cerr << "GatherOperation with this name" << name
                << "already exists!\n";
      std::terminate();
    }
    return &(it->second);
  }
  template <typename CombinedFieldT = Real>
  ScatterOperation* AddScatterFieldsOp(
    const std::string& name, const std::string& internal_field_name,
    const std::vector<std::string>& fields_to_scatter,
    Omega_h::Read<Omega_h::I8> mask = {})
  {
    auto scatter_fields =
      detail::find_many_or_error(fields_to_scatter, fields_);
    static constexpr int seach_nx = 10;
    static constexpr int seach_ny = 10;

    auto& combined = detail::find_or_create_internal_field<CombinedFieldT>(
      internal_field_name, internal_fields_, internal_mesh_, mask, seach_nx,
      seach_ny);
    auto [it, inserted] = scatter_operations_.template try_emplace(
      name, std::move(scatter_fields), combined);

    if (!inserted) {
      std::cerr << "Scatter with this name" << name << "already exists!\n";
      std::terminate();
    }
    return &(it->second);
  }
  [[nodiscard]] const redev::Partition& GetPartition() const
  {
    return redev_.GetPartition();
  }
  [[nodiscard]] Omega_h::Mesh& GetMesh() { return internal_mesh_; }

private:
  std::string name_;
  MPI_Comm mpi_comm_;
  redev::Redev redev_;
  std::unordered_map<std::string, ConvertibleCoupledField> fields_;
  // coupler owns internal fields since both gather/scatter ops use these
  std::unordered_map<std::string, detail::InternalField> internal_fields_;
  // gather and scatter operations have reference to internal fields
  std::unordered_map<std::string, ScatterOperation> scatter_operations_;
  std::unordered_map<std::string, GatherOperation> gather_operations_;
  Omega_h::Mesh& internal_mesh_;
};
class CouplerClient
{
public:
  CouplerClient(std::string name, MPI_Comm comm, redev::Partition partition)
    : name_(std::move(name)),
      mpi_comm_(comm),
      redev_({comm, std::move(partition), ProcessType::Client})
  {
  }
  template <typename FieldShimT>
  CoupledField* AddField(std::string unique_name, FieldShimT field_shim)
  {
    auto [it, inserted] = fields_.template try_emplace(
      unique_name, unique_name, std::move(field_shim), redev_, mpi_comm_);
    if (!inserted) {
      std::cerr << "OHField with this name" << unique_name
                << "already exists!\n";
      std::terminate();
    }
    return &(it->second);
  }

  void SendField(const std::string& name)
  {
    detail::find_or_error(name, fields_).Send();
  };
  void ReceiveField(const std::string& name)
  {
    detail::find_or_error(name, fields_).Receive();
  };
  [[nodiscard]] const redev::Partition& GetPartition() const
  {
    return redev_.GetPartition();
  }

private:
  std::string name_;
  MPI_Comm mpi_comm_;
  redev::Redev redev_;
  std::unordered_map<std::string, CoupledField> fields_;
};
} // namespace wdmcpl

#endif
