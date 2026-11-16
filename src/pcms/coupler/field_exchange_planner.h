#ifndef PCMS_FIELD_EXCHANGE_PLANNER_H
#define PCMS_FIELD_EXCHANGE_PLANNER_H

#include "pcms/field/field_layout.h"
#include "pcms/utility/arrays.h"
#include "overlap_mask.h"
#include <redev.h>
#include <vector>

namespace pcms
{

struct ExchangePlan
{
  redev::LOs dest_ranks;
  // Field-payload offsets and permutation. GID-message headers are inserted
  // only while exchanging the layout and are not represented in this plan.
  redev::LOs offsets;
  std::vector<LO> permutation;
  size_t msg_size = 0;
};

class FieldExchangePlanner
{
public:
  virtual ExchangePlan BuildExchangePlan(
    const FieldLayout& layout, const redev::Partition& partition,
    const OverlapMask* overlap_mask = nullptr) const = 0;

  virtual ExchangePlan BuildReceivePlan(
    const FieldLayout& layout, GlobalIDView<HostMemorySpace> received_gids,
    int rank, int nproc,
    const redev::InMessageLayout& in_message_layout) const = 0;

  virtual void FillGidMessage(
    const FieldLayout& layout, const ExchangePlan& plan,
    Rank1View<GO, HostMemorySpace> gid_message) const = 0;

  virtual ~FieldExchangePlanner() noexcept = default;
};

class GenericFieldExchangePlanner : public FieldExchangePlanner
{
public:
  ExchangePlan BuildExchangePlan(
    const FieldLayout& layout, const redev::Partition& partition,
    const OverlapMask* overlap_mask = nullptr) const override;

  ExchangePlan BuildReceivePlan(
    const FieldLayout& layout, GlobalIDView<HostMemorySpace> received_gids,
    int rank, int nproc,
    const redev::InMessageLayout& in_message_layout) const override;

  void FillGidMessage(
    const FieldLayout& layout, const ExchangePlan& plan,
    Rank1View<GO, HostMemorySpace> gid_message) const override;
};

} // namespace pcms
#endif // PCMS_FIELD_EXCHANGE_PLANNER_H
