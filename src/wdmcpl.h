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

#include "wdmcpl/server.h"
#include "wdmcpl/client.h"
#include "wdmcpl/field_communicator.h"


namespace wdmcpl
{
using ProcessType = redev::ProcessType;
class GatherOperation;
class ScatterOperation;

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
  template <typename FieldAdapterT, typename CommT>
  ConvertibleCoupledField(const std::string& name, FieldAdapterT field_adapter,
                          FieldCommunicator<CommT> field_comm,
                          Omega_h::Mesh& internal_mesh,
                          TransferOptions native_to_internal,
                          TransferOptions internal_to_native,
                          Omega_h::Read<Omega_h::I8> internal_field_mask = {})
    : internal_field_{OmegaHField<typename FieldAdapterT::value_type,
                                  InternalCoordinateElement>(
        name + ".__internal__", internal_mesh, internal_field_mask)}
  {
    coupled_field_ = std::make_unique<CoupledFieldModel<FieldAdapterT, CommT>>(
      std::move(field_adapter), std::move(field_comm),
      std::move(native_to_internal), std::move(internal_to_native));
  }
  template <typename FieldAdapterT>
  ConvertibleCoupledField(const std::string& name, FieldAdapterT field_adapter,
                          redev::Redev& redev, MPI_Comm mpi_comm,
                          Omega_h::Mesh& internal_mesh,
                          TransferOptions native_to_internal,
                          TransferOptions internal_to_native,
                          Omega_h::Read<Omega_h::I8> internal_field_mask = {})
    : internal_field_{OmegaHField<typename FieldAdapterT::value_type,
                                  InternalCoordinateElement>(
        name + ".__internal__", internal_mesh, internal_field_mask)}
  {
    coupled_field_ =
      std::make_unique<CoupledFieldModel<FieldAdapterT, FieldAdapterT>>(
        name, std::move(field_adapter), redev, mpi_comm,
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
  [[nodiscard]] InternalField& GetInternalField() noexcept
  {
    return internal_field_;
  }
  [[nodiscard]] const InternalField& GetInternalField() const noexcept
  {
    return internal_field_;
  }
  struct CoupledFieldConcept
  {
    virtual void Send() = 0;
    virtual void Receive() = 0;
    virtual void SyncNativeToInternal(InternalField&) = 0;
    virtual void SyncInternalToNative(const InternalField&) = 0;
    virtual ~CoupledFieldConcept() = default;
  };
  template <typename FieldAdapterT, typename CommT>
  struct CoupledFieldModel final : CoupledFieldConcept
  {
    using value_type = typename FieldAdapterT::value_type;

    CoupledFieldModel(FieldAdapterT&& field_adapter,
                      FieldCommunicator<CommT>&& comm,
                      TransferOptions&& native_to_internal,
                      TransferOptions&& internal_to_native)
      : field_adapter_(std::move(field_adapter)),
        comm_(std::move(comm)),
        native_to_internal_(std::move(native_to_internal)),
        internal_to_native_(std::move(internal_to_native))
    {
    }
    CoupledFieldModel(const std::string& name, FieldAdapterT&& field_adapter,
                      redev::Redev& redev, MPI_Comm mpi_comm,
                      TransferOptions&& native_to_internal,
                      TransferOptions&& internal_to_native)
      : field_adapter_(std::move(field_adapter)),
        comm_(FieldCommunicator<FieldAdapterT>(name, redev, mpi_comm,
                                                    field_adapter_)),
        native_to_internal_(std::move(native_to_internal)),
        internal_to_native_(std::move(internal_to_native))
    {
    }
    void Send() final { comm_.Send(); };
    void Receive() final { comm_.Receive(); };
    void SyncNativeToInternal(InternalField& internal_field) final
    {
      ConvertFieldAdapterToOmegaH(field_adapter_, internal_field,
                                  native_to_internal_.transfer_method,
                                  native_to_internal_.evaluation_method);
    };
    void SyncInternalToNative(const InternalField& internal_field) final
    {
      ConvertOmegaHToFieldAdapter(internal_field, field_adapter_,
                                  internal_to_native_.transfer_method,
                                  internal_to_native_.evaluation_method);
    };

    FieldAdapterT field_adapter_;
    FieldCommunicator<CommT> comm_;
    TransferOptions native_to_internal_;
    TransferOptions internal_to_native_;
  };

private:
  std::unique_ptr<CoupledFieldConcept> coupled_field_;
  // even though we know the type of the internal field,
  // we store it as the InternalField variant since this avoids any copies
  // This comes at the cost of a slightly larger type with need to use the get<>
  // function
  InternalField internal_field_;
};
class CoupledField
{
public:
  template <typename FieldAdapterT, typename CommT>
  CoupledField(const std::string& name, FieldAdapterT field_adapter,
               FieldCommunicator<CommT> field_comm,
               Omega_h::Mesh& internal_mesh, TransferOptions native_to_internal,
               TransferOptions internal_to_native)
  {
    coupled_field_ = std::make_unique<CoupledFieldModel<FieldAdapterT, CommT>>(
      std::move(field_adapter), std::move(field_comm));
  }
  template <typename FieldAdapterT>
  CoupledField(const std::string& name, FieldAdapterT field_adapter,
               redev::Redev& redev, MPI_Comm mpi_comm)
  {

    coupled_field_ =
      std::make_unique<CoupledFieldModel<FieldAdapterT, FieldAdapterT>>(
        name, std::move(field_adapter), redev, mpi_comm);
  }

  void Send() { coupled_field_->Send(); }
  void Receive() { coupled_field_->Receive(); }
  struct CoupledFieldConcept
  {
    virtual void Send() = 0;
    virtual void Receive() = 0;
    virtual ~CoupledFieldConcept() = default;
  };
  template <typename FieldAdapterT, typename CommT>
  struct CoupledFieldModel final : CoupledFieldConcept
  {
    using value_type = typename FieldAdapterT::value_type;

    CoupledFieldModel(FieldAdapterT&& field_adapter,
                      FieldCommunicator<CommT>&& comm)
      : field_adapter_(std::move(field_adapter)), comm_(std::move(comm))
    {
    }
    CoupledFieldModel(const std::string& name, FieldAdapterT&& field_adapter,
                      redev::Redev& redev, MPI_Comm mpi_comm)
      : field_adapter_(std::move(field_adapter)),
        comm_(
          FieldCommunicator<CommT>(name, redev, mpi_comm,
                                               field_adapter_))
    {
    }
    void Send() final { comm_.Send(); };
    void Receive() final { comm_.Receive(); };

    FieldAdapterT field_adapter_;
    FieldCommunicator<CommT> comm_;
  };

private:
  std::unique_ptr<CoupledFieldConcept> coupled_field_;
};

using CombinerFunction = std::function<void(
  nonstd::span<const std::reference_wrapper<InternalField>>,
  InternalField&)>;
class GatherOperation
{
public:
  GatherOperation(std::vector<std::reference_wrapper<ConvertibleCoupledField>>
                    fields_to_gather,
                  InternalField& combined_field,
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
  std::vector<std::reference_wrapper<InternalField>> internal_fields_;
  InternalField& combined_field_;
  CombinerFunction combiner_;
};
class ScatterOperation
{
public:
  ScatterOperation(std::vector<std::reference_wrapper<ConvertibleCoupledField>>
                     fields_to_scatter,
                   InternalField& combined_field)
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
  std::vector<std::reference_wrapper<InternalField>> internal_fields_;
  InternalField& combined_field_;
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
  template <typename FieldAdapterT>
  ConvertibleCoupledField* AddField(
    std::string unique_name, FieldAdapterT field_adapter,
    FieldTransferMethod to_field_transfer_method,
    FieldEvaluationMethod to_field_eval_method,
    FieldTransferMethod from_field_transfer_method,
    FieldEvaluationMethod from_field_eval_method,
    Omega_h::Read<Omega_h::I8> internal_field_mask = {})
  {
    auto [it, inserted] = fields_.template try_emplace(
      unique_name, unique_name, std::move(field_adapter), redev_, mpi_comm_,
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
      internal_field_name, internal_fields_, internal_mesh_, mask, seach_nx,
      seach_ny);
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
  std::unordered_map<std::string, InternalField> internal_fields_;
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
  template <typename FieldAdapterT>
  CoupledField* AddField(std::string unique_name, FieldAdapterT field_adapter)
  {
    auto [it, inserted] = fields_.template try_emplace(
      unique_name, unique_name, std::move(field_adapter), redev_, mpi_comm_);
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
