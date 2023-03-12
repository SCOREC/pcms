#ifndef WDM_COUPLING_CLIENT_H
#define WDM_COUPLING_CLIENT_H
#include "wdmcpl/common.h"
#include "wdmcpl/field_communicator.h"
namespace wdmcpl
{

class CoupledField
{
public:
  template <typename FieldAdapterT>
  CoupledField(const std::string& name, FieldAdapterT field_adapter,
               MPI_Comm mpi_comm, redev::Redev& redev, redev::BidirectionalChannel& channel)
  {

    coupled_field_ =
      std::make_unique<CoupledFieldModel<FieldAdapterT, FieldAdapterT>>(
        name, std::move(field_adapter), mpi_comm, redev, channel);
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

    CoupledFieldModel(const std::string& name, FieldAdapterT&& field_adapter,
                      MPI_Comm mpi_comm,
                      redev::Redev& redev,
                      redev::BidirectionalChannel& channel)
      : field_adapter_(std::move(field_adapter)),
        comm_(FieldCommunicator<CommT>(name, mpi_comm, redev, channel, field_adapter_))
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
class CouplerClient
{
public:
  CouplerClient(std::string name, MPI_Comm comm,
                redev::TransportType transport_type = redev::TransportType::BP4,
                adios2::Params params = {{"Streaming", "On"}, {"OpenTimeoutSecs", "400"}},
                std::string path="")
    : name_(std::move(name)), mpi_comm_(comm), redev_(comm),
      channel_{redev_.CreateAdiosChannel(name_, std::move(params),
                                  transport_type, std::move(path))}
  {
  }

  [[nodiscard]] const redev::Partition& GetPartition() const
  {
    return redev_.GetPartition();
  }

  template <typename FieldAdapterT>
  CoupledField* AddField(
    std::string name, FieldAdapterT field_adapter,
    redev::TransportType transport_type = redev::TransportType::BP4,
    adios2::Params params = {{"Streaming", "On"}, {"OpenTimeoutSecs", "400"}})
  {
    auto [it, inserted] = fields_.template try_emplace(name, name, std::move(field_adapter),
                                                       mpi_comm_, redev_, channel_);
    if (!inserted) {
      std::cerr << "OHField with this name" << name << "already exists!\n";
      std::terminate();
    }
    return &(it->second);
  }

  // take a string& since map cannot be searched with string_view
  // (heterogeneous lookup)
  void SendField(const std::string& name)
  {
    detail::find_or_error(name, fields_).Send();
  };
  // take a string& since map cannot be searched with string_view
  // (heterogeneous lookup)
  void ReceiveField(const std::string& name)
  {
    detail::find_or_error(name, fields_).Receive();
  };

private:
  std::string name_;
  MPI_Comm mpi_comm_;
  redev::Redev redev_;
  // map rather than unordered_map is necessary to avoid iterator invalidation.
  // This is important because we pass pointers to the fields out of this class
  std::map<std::string, CoupledField> fields_;
  redev::BidirectionalChannel channel_;

};
} // namespace wdmcpl

#endif // WDM_COUPLING_CLIENT_H
