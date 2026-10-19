#ifndef PCMS_PARTITION_H
#define PCMS_PARTITION_H
#include "pcms/utility/common.h"
#include "pcms/utility/profile.h"
#include "pcms/utility/types.h"
#include <redev.h>
#include <redev_partition.h>

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
  // PointSearchPtn maps a point to a *list* of candidate ranks (a point may
  // fall in the overlap region of several ranks' grid cells), unlike
  // ClassPtn/RCBPtn which map to a single owning rank. GetDr() is a
  // single-rank query, so this returns the first candidate as a
  // representative owner; callers that need every candidate should use
  // ptn.GetCandidateRanks() directly instead of GetDr().
  auto operator()(const redev::PointSearchPtn& ptn)
  {
    PCMS_FUNCTION_TIMER;
    const auto& candidates = ptn.GetCandidateRanks(coord_);
    return candidates.empty() ? static_cast<LO>(-1) : candidates.front();
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

#endif // PCMS_PARTITION_H
