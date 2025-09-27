#ifndef PCMS_PARTITION_H
#define PCMS_PARTITION_H
#include "pcms/common.h"
#include "pcms/profile.h"

namespace pcms
{
struct GetRank
{
  GetRank(const LO id, const LO dim, const std::array<Real, 3>& coord)
    : id_(id), dim_(dim), coord_(coord)
  {
    PCMS_FUNCTION_TIMER;
  }
  auto operator()(const redev::ClassPtn& ptn) const
  {
    PCMS_FUNCTION_TIMER;
    const auto ent = redev::ClassPtn::ModelEnt({dim_, id_});
    return ptn.GetRank(ent);
  }
  auto operator()(const redev::RCBPtn& ptn)
  {
    PCMS_FUNCTION_TIMER;
    return ptn.GetRank(coord_);
  }
  LO id_;
  LO dim_;
  std::array<Real, 3> coord_;
};

struct Partition
{
  Partition(const redev::Partition& partition) : partition_(partition)
  {
    PCMS_FUNCTION_TIMER;
  }

  auto GetDr(const LO id, const LO dim,
             const std::array<Real, 3> coord = {}) const
  {
    return std::visit(GetRank{id, dim, coord}, partition_);
  }

  redev::Partition partition_;
};
} // namespace pcms

#endif //PCMS_PARTITION_H
