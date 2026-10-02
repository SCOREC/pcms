#ifndef PCMS_DISTRIBUTED_EVALUATION_PARTITION_H
#define PCMS_DISTRIBUTED_EVALUATION_PARTITION_H
#include <array>
#include <vector>
#include "pcms/field/distributed_partition.h"
#include "pcms/utility/types.h"

namespace pcms
{

/**
 * @brief Routing-hint abstraction for the distributed-field evaluation server:
 * "which rank(s) might own this point?".
 *
 * It is the extension point for supplying the server's routing partition from a
 * different source without the rest of the protocol caring which source that
 * is; the only implementation today is BoundingBoxEvaluationPartition.
 */
class EvaluationPartitionQuery
{
public:
  /**
   * @brief Appends the candidate ranks for point to candidate_ranks.
   *
   * Does not clear candidate_ranks first, so callers may accumulate across
   * several points.
   */
  virtual void GetCandidateRanks(const std::array<Real, 3>& point,
                                 std::vector<int>& candidate_ranks) const = 0;

  virtual ~EvaluationPartitionQuery() noexcept = default;
};

/**
 * @brief Routes by wrapping a DistributedBoundingBoxPartition.
 *
 * May return more than one candidate; see DistributedBoundingBoxPartition's own
 * documentation for why that is expected and handled by the caller.
 */
class BoundingBoxEvaluationPartition : public EvaluationPartitionQuery
{
public:
  explicit BoundingBoxEvaluationPartition(
    DistributedBoundingBoxPartition partition)
    : partition_(std::move(partition))
  {
  }

  void GetCandidateRanks(const std::array<Real, 3>& point,
                         std::vector<int>& candidate_ranks) const override
  {
    partition_.GetCandidateRanks(point, candidate_ranks);
  }

private:
  DistributedBoundingBoxPartition partition_;
};

} // namespace pcms

#endif // PCMS_DISTRIBUTED_EVALUATION_PARTITION_H
