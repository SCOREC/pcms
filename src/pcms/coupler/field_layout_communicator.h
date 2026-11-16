#ifndef FIELD_LAYOUT_COMMUNICATOR_H_
#define FIELD_LAYOUT_COMMUNICATOR_H_

#include "pcms/field/field_layout.h"
#include "field_exchange_planner.h"
#include "overlap_mask.h"
#include "pcms/utility/profile.h"
#include "pcms/utility/arrays.h"
#include <redev.h>
#include <memory>

namespace pcms
{

class FieldLayoutCommunicator
{
public:
  FieldLayoutCommunicator(const std::string& name, MPI_Comm mpi_comm,
                          redev::Redev& redev, redev::Channel& channel,
                          const FieldLayout& layout, bool own_mpi_comm = false);

  FieldLayoutCommunicator(const std::string& name, MPI_Comm mpi_comm,
                          redev::Redev& redev, redev::Channel& channel,
                          const FieldLayout& layout,
                          std::unique_ptr<FieldExchangePlanner> planner,
                          bool own_mpi_comm = false,
                          const OverlapMask* overlap_mask = nullptr);

  Rank1View<const pcms::LO, pcms::HostMemorySpace> GetPermutationArray() const;

  const std::string& GetName() const;

  const FieldLayout& GetLayout() const;

  size_t GetMsgSize() const;

  redev::Channel& GetChannel() const;

  MPI_Comm& GetMPIComm();

  template <typename T>
  void SetOutMessageLayout(redev::BidirectionalComm<T>& comm)
  {
    // plan_.offsets is a per-holder CSR layout. Each holder carries
    // num_components field values, so the field message's offsets are scaled
    // by num_components while the destination ranks are unchanged.
    const int num_comp = layout_.GetNumComponents();
    if (num_comp == 1) {
      comm.SetOutMessageLayout(plan_.dest_ranks, plan_.offsets);
      return;
    }
    redev::LOs scaled_offsets(plan_.offsets.size());
    for (size_t i = 0; i < plan_.offsets.size(); ++i) {
      scaled_offsets[i] = plan_.offsets[i] * num_comp;
    }
    comm.SetOutMessageLayout(plan_.dest_ranks, scaled_offsets);
  }

  void UpdateLayout();

  void UpdateLayoutNull();

  void SetOverlapMask(std::unique_ptr<OverlapMask> mask);

  ~FieldLayoutCommunicator();

private:
  MPI_Comm mpi_comm_;
  redev::Channel& channel_;
  redev::BidirectionalComm<GO> gid_comm_;
  ExchangePlan plan_;
  const FieldLayout& layout_;
  std::string name_;
  redev::Redev& redev_;
  std::unique_ptr<FieldExchangePlanner> planner_;
  std::unique_ptr<OverlapMask> overlap_mask_;
  bool own_mpi_comm_ = false;
};
} // namespace pcms

#endif // FIELD_LAYOUT_COMMUNICATOR_H_
