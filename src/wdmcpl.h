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
namespace detail
{
using InternalCoordinateElement = Real;

// internal field can only be one of the types supported by Omega_h
// The coordinate element for all internal fields is the same since
// all internal fields are on the same mesh
using InternalField =
  std::variant<std::monostate,
               OmegaHField<Omega_h::I8, InternalCoordinateElement>,
               OmegaHField<Omega_h::I32, InternalCoordinateElement>,
               OmegaHField<Omega_h::I64, InternalCoordinateElement>,
               OmegaHField<Omega_h::Real, InternalCoordinateElement>>;
} // namespace detail
using ProcessType = redev::ProcessType;

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

class FieldCommunicator
{
public:
  template <typename FieldShimT>
  FieldCommunicator(
    std::string name, redev::Redev& redev, MPI_Comm mpi_comm,
    FieldShimT& field_shim,
    redev::TransportType transport_type = redev::TransportType::BP4,
    adios2::Params params = adios2::Params{{"Streaming", "On"},
                                           {"OpenTimeoutSecs", "12"}})
  {
    field_communicator_ = std::make_unique<FieldCommunicatorModel<FieldShimT>>(
      std::move(name), redev, mpi_comm, field_shim, transport_type,
      std::move(params));
  }
  // default constructor creates field communicator that doesn't send or receive
  FieldCommunicator()
  {
    field_communicator_ = std::make_unique<FieldCommunicatorModel<void>>();
  }
  void Send() const { field_communicator_->Send(); }
  void Receive() const { field_communicator_->Receive(); }

private:
  struct FieldCommunicatorConcept
  {
  public:
    virtual void Send() = 0;
    virtual void Receive() = 0;
    virtual ~FieldCommunicatorConcept() = default;
  };
  template <typename FieldShimT>
  struct FieldCommunicatorModel final : public FieldCommunicatorConcept
  {
    using T = typename FieldShimT::value_type;

  public:
    FieldCommunicatorModel(
      std::string name, redev::Redev& redev, MPI_Comm mpi_comm,
      FieldShimT& field_shim,
      redev::TransportType transport_type = redev::TransportType::BP4,
      adios2::Params params = adios2::Params{{"Streaming", "On"},
                                             {"OpenTimeoutSecs", "12"}})
      : mpi_comm_(mpi_comm),
        comm_buffer_{},
        message_permutation_{},
        buffer_size_needs_update_{true},
        redev_{redev},
        name_{std::move(name)},
        field_shim_(std::move(field_shim))
    {
      std::string transport_name = name_;
      comm_ =
        redev_.CreateAdiosClient<T>(transport_name, params, transport_type);
      // set up GID comm to do setup phase and get the
      // FIXME: use  one directional comm instead of the adios bidirectional
      // comm
      transport_name = transport_name.append("_gids");
      gid_comm_ = redev_.CreateAdiosClient<wdmcpl::GO>(transport_name, params,
                                                       transport_type);
      UpdateLayout();
    }

    FieldCommunicatorModel(const FieldCommunicatorModel&) = delete;
    FieldCommunicatorModel(FieldCommunicatorModel&&) = default;
    FieldCommunicatorModel& operator=(const FieldCommunicatorModel&) = delete;
    FieldCommunicatorModel& operator=(FieldCommunicatorModel&&) = default;

    void Send() override
    {
      auto n = field_shim_.Serialize({}, {});
      if (comm_buffer_.size() != n) {
        UpdateLayout();
        REDEV_ALWAYS_ASSERT(comm_buffer_.size() == n);
      }

      auto buffer = ScalarArrayView<typename decltype(comm_buffer_)::value_type,
                                    HostMemorySpace>(comm_buffer_.data(),
                                                     comm_buffer_.size());
      const auto permutation = ScalarArrayView<
        const typename decltype(message_permutation_)::value_type,
        HostMemorySpace>(message_permutation_.data(),
                         message_permutation_.size());

      field_shim_.Serialize(buffer, permutation);
      comm_.Send(buffer.data_handle());
    }
    void Receive() override
    {
      auto data = comm_.Recv();
      auto buffer =
        ScalarArrayView<T, HostMemorySpace>(data.data(), data.size());
      static_assert(std::is_same_v<T, typename decltype(data)::value_type>);
      auto permutation = ScalarArrayView<
        const typename decltype(message_permutation_)::value_type,
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
  struct FieldCommunicatorModel<void> final : public FieldCommunicatorConcept
  {
    void Send() {}
    void Receive() {}
  };

  std::unique_ptr<FieldCommunicatorConcept> field_communicator_;
};

struct FieldTransferOptions
{
  FieldTransferMethod transfer_method;
  FieldEvaluationMethod evaluation_method;
};

namespace detail
{

// helper function for dealing with field maps
template <typename T, typename U>
auto& find_or_error(const std::string& name,
                    const std::unordered_map<T, U>& map)
{
  auto it = map.find(name);
  WDMCPL_ALWAYS_ASSERT(it != map.end());
  return it->second;
}

template <typename T>
struct SourceFieldIdentity
{
  using type = T;
};
template <typename T>
struct TargetFieldIdentity
{
  using type = T;
};
struct TransferOperation
{
  template <typename SourceField, typename TargetField>
  TransferOperation(SourceFieldIdentity<SourceField>,
                    TargetFieldIdentity<TargetField>, FieldTransferMethod ftm,
                    FieldEvaluationMethod fem)
  {
    transfer_operation_ =
      std::make_unique<TransferOperationModel<SourceField, TargetField>>(ftm,
                                                                         fem);
  }
  template <typename SourceField, typename TargetField>
  TransferOperation(const SourceField& source, TargetField& target,
                    FieldTransferMethod ftm, FieldEvaluationMethod fem)
  {
    auto model_ptr =
      std::make_unique<TransferOperationModel<SourceField, TargetField>>(ftm,
                                                                         fem);
    // directly set ptrs to avoid creation of type erased field type
    model_ptr->source_ = &source;
    model_ptr->target_ = &target;
    transfer_operation_ = std::move(model_ptr);
  }
  void Run() { transfer_operation_->Run(); }
  void SetTargetField(Field& f) { transfer_operation_->SetTargetField(f); }
  void SetSourceField(const Field& f)
  {
    transfer_operation_->SetSourceField(f);
  }

private:
  struct TransferOperationConcept
  {
    virtual void Run() = 0;
    virtual void SetTargetField(Field&) = 0;
    virtual void SetSourceField(const Field&) = 0;
    virtual ~TransferOperationConcept() = default;
  };

  template <typename SourceField, typename TargetField>
  struct TransferOperationModel final : TransferOperationConcept
  {
    TransferOperationModel(FieldTransferMethod ftm, FieldEvaluationMethod fem)
      : source_(nullptr),
        target_(nullptr),
        transfer_method_(ftm),
        evaluation_method_(fem)
    {
    }

    void SetTargetField(Field& field)
    {
      target_ = field_cast<TargetField>(&field);
    }
    void SetSourceField(const Field& field)
    {
      source_ = field_cast<SourceField>(&field);
    }

    void Run()
    {
      // source and target fields must be set with SetSourceField and
      // SetTargetField
      WDMCPL_ALWAYS_ASSERT(source_ != nullptr && target_ != nullptr);
      transfer_field(*source_, *target_, transfer_method_, evaluation_method_);
    }
    // nonowning pointers to source/target field
    const SourceField* source_;
    TargetField* target_;

    FieldTransferMethod transfer_method_;
    FieldEvaluationMethod evaluation_method_;
  };
  std::unique_ptr<TransferOperationConcept> transfer_operation_;
};

} // namespace detail

class CoupledField
{
public:
  class GatherScatterPasskey
  {
  private:
    friend class GatherOperation;
    friend class ScatterOperation;
    GatherScatterPasskey() = default;
    GatherScatterPasskey(const GatherScatterPasskey&) = default;
    GatherScatterPasskey(GatherScatterPasskey&&) = delete;
    GatherScatterPasskey& operator=(const GatherScatterPasskey&) = delete;
    GatherScatterPasskey& operator=(GatherScatterPasskey&&) = delete;
  };
  // only coupler creates coupled field??
  // don't really want "internal" field or TransferOperation in interface
  CoupledField(FieldCommunicator comm, detail::InternalField internal_field,
               detail::TransferOperation native_to_internal,
               detail::TransferOperation internal_to_native)
    : comm_(std::move(comm)),
      native_to_internal_(std::move(native_to_internal)),
      internal_to_native_(std::move(internal_to_native)),
      internal_field_(std::move(internal_field))
  {
  }
  template <typename FieldShimT>
  CoupledField(std::string name, redev::Redev&, MPI_Comm mpi_comm,
               FieldShimT field_shim, FieldTransferMethod,
               FieldEvaluationMethod, FieldTransferMethod,
               FieldEvaluationMethod)
  {
    // 1. Create field communicator
  }
  void Send() const { comm_.Send(); }
  void Receive() const { comm_.Receive(); }
  void SyncNativeToInternal(GatherScatterPasskey) { native_to_internal_.Run(); }
  void SyncInternalToNative(GatherScatterPasskey) { internal_to_native_.Run(); }
  [[nodiscard]] detail::InternalField& GetInternalField(
    GatherScatterPasskey) noexcept
  {
    return internal_field_;
  }
  [[nodiscard]] const detail::InternalField& GetInternalField(
    GatherScatterPasskey) const noexcept
  {
    return internal_field_;
  }

private:
  // FIXME need to own Field shim here!
  FieldCommunicator comm_;
  detail::TransferOperation native_to_internal_;
  detail::TransferOperation internal_to_native_;
  detail::InternalField internal_field_;
};

class GatherOperation
{
public:
  using Combiner = std::function<void(
    const std::vector<std::reference_wrapper<detail::InternalField>>&,
    detail::InternalField&)>;
  GatherOperation(
    std::vector<std::reference_wrapper<CoupledField>> fields_to_gather,
    detail::InternalField& combined_field, Combiner combiner)
    : coupled_fields_(std::move(fields_to_gather)),
      combined_field_{combined_field},
      combiner_(std::move(combiner))
  {

    internal_fields_.reserve(coupled_fields_.size());
    std::transform(begin(coupled_fields_), end(coupled_fields_),
                   std::back_inserter(internal_fields_), [](CoupledField& fld) {
                     return std::ref(fld.GetInternalField(
                       CoupledField::GatherScatterPasskey()));
                   });
  }
  void Run() const
  {
    for (auto& field : coupled_fields_) {
      field.get().Receive();
      field.get().SyncNativeToInternal(CoupledField::GatherScatterPasskey());
    }
    combiner_(internal_fields_, combined_field_);
  };

private:
  std::vector<std::reference_wrapper<CoupledField>> coupled_fields_;
  std::vector<std::reference_wrapper<detail::InternalField>> internal_fields_;
  detail::InternalField& combined_field_;
  Combiner combiner_;
};
class ScatterOperation
{
public:
  ScatterOperation(
    std::vector<std::reference_wrapper<CoupledField>> fields_to_scatter,
    detail::InternalField& combined_field)
    : coupled_fields_(std::move(fields_to_scatter)),
      combined_field_{combined_field}
  {

    internal_fields_.reserve(coupled_fields_.size());
    std::transform(begin(coupled_fields_), end(coupled_fields_),
                   std::back_inserter(internal_fields_), [](CoupledField& fld) {
                     return std::ref(fld.GetInternalField(
                       CoupledField::GatherScatterPasskey()));
                   });
  }
  void Run() const
  {
    // possible we may need to add a splitter operation here. will evaluate if
    // needed splitter(combined_field, internal_fields_);
    for (auto& field : coupled_fields_) {
      field.get().SyncInternalToNative(CoupledField::GatherScatterPasskey());
      field.get().Send();
    }
  };

private:
  std::vector<std::reference_wrapper<CoupledField>> coupled_fields_;
  std::vector<std::reference_wrapper<detail::InternalField>> internal_fields_;
  detail::InternalField& combined_field_;
};

// TODO: refactor into ClientCoupler, ServerCoupler, BaseCoupler or Coupler<PT>
// with base with shared
template <ProcessType PT>
class Coupler
{
public:
  Coupler(std::string name, MPI_Comm comm, redev::Partition partition)
    : name_(std::move(name)),
      mpi_comm_(comm),
      redev_({comm, std::move(partition), PT})
  {
  }
  template <typename FieldShimT>
  CoupledField* AddField(std::string unique_name, FieldShimT field_shim)
  {
    auto [it, inserted] = fields_.template try_emplace(
      unique_name, unique_name, redev_, mpi_comm_, std::move(field_shim));
    if (!inserted) {
      std::cerr << "OHField with this name" << unique_name
                << "already exists!\n";
      std::exit(EXIT_FAILURE);
    }
    REDEV_ALWAYS_ASSERT(it->second != nullptr);
    return it->second;
  }
  template <typename FieldShimT>
  CoupledField* AddField(std::string unique_name, FieldShimT field_shim,
                         FieldTransferMethod to_field_transfer_method,
                         FieldEvaluationMethod to_field_eval_method,
                         FieldTransferMethod from_field_transfer_method,
                         FieldEvaluationMethod from_field_eval_method)
  {
    auto [it, inserted] = fields_.template try_emplace(
      unique_name, unique_name, redev_, mpi_comm_, std::move(field_shim),
      to_field_transfer_method, to_field_eval_method,
      from_field_transfer_method, from_field_eval_method);
    if (!inserted) {
      std::cerr << "OHField with this name" << unique_name
                << "already exists!\n";
      std::exit(EXIT_FAILURE);
    }
    REDEV_ALWAYS_ASSERT(it->second != nullptr);
    return it->second;
  }

  // here we take a string, not string_view since we need to search map
  void ScatterFields(const std::string& name)
  {
    /*
    auto& scatter_operation = detail::find_or_error(name, scatter_operations_);
    auto& internal_field = scatter_operation.GetInternalField();
    auto& fields = scatter_operation.GetFields();

    for (auto& field : fields) {
      auto& field_name = field.GetName();
      // using invoke here to clarify that we are calling the returned function
      std::invoke(find_or_error(field_name, scatter_transfer_));
      auto& field_comm = detail::find_or_error(field_name, fields_);
      field_comm->Send();
    }
     */
    detail::find_or_error(name, scatter_operations_).Run();
  }
  // here we take a string, not string_view since we need to search map
  void GatherFields(const std::string& name)
  {
    detail::find_or_error(name, gather_operations_).Run();
    // gather op needs to store a copy of the appropriate field transfer op.
    // Upon construction of gather op the field_transfer_op is set for each
    // field since we only then have
    // the pointer to the internal field (it's created in CreateGatherOp).
    // note that we end up taking a copy of the TransferOp since a single field
    // may get transferred onto multiple internal fields!
    // gather_op.Receive();
    ////for field in fields:
    // gather_op.NativeToInternal();
    // gather_op.Combine();
    //  or just gather_op.Run();
  }
  void SendField(const std::string& name)
  {
    detail::find_or_error(name, fields_).Send();
  };
  void ReceiveField(const std::string& name)
  {
    detail::find_or_error(name, fields_).Receive();
  };

private:
  std::string name_;
  static constexpr ProcessType process_type_{PT};
  MPI_Comm mpi_comm_;
  redev::Redev redev_;
  std::unordered_map<std::string, CoupledField> fields_;
  // coupler owns internal fields since both gather/scatter ops use these
  std::unordered_map<std::string, detail::InternalField> internal_fields_;
  // gather and scatter operations have reference to internal fields
  std::unordered_map<std::string, GatherOperation> scatter_operations_;
  std::unordered_map<std::string, ScatterOperation> gather_operations_;
};
} // namespace wdmcpl

#endif
