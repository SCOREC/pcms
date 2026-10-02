#ifndef PCMS_DISTRIBUTED_PARTITION_H
#define PCMS_DISTRIBUTED_PARTITION_H
#include <array>
#include <vector>
#include "pcms/utility/types.h"

namespace pcms
{

class FieldLayout;

/**
 * @brief Coarse spatial routing index answering "which ranks might own this
 * point?" for a distributed field.
 *
 * Boxes are padded axis-aligned bounds of each rank's actual owned region,
 * derived from its owned DOF-holder coordinates (see LocalOwnedBounds), so they
 * never under-cover a rank's real data (no false negatives) but can overlap or
 * over-cover near irregular partition boundaries (false positives are
 * possible). It is a hint, not a ground-truth owner map: callers must verify
 * every candidate with a local, authoritative point-location test and treat a
 * point matched by no candidate as genuinely outside the global domain.
 *
 * The partition is built once, from bounds a field job published, and reused
 * across many queries. It carries no MPI communicator of its own.
 */
class DistributedBoundingBoxPartition
{
public:
  /**
   * @brief Builds a partition from already-gathered per-rank bounds.
   *
   * Has no MPI communication of its own, so it can be used by a process (e.g.
   * a standalone evaluation server) that received every rank's box over some
   * other channel (such as redev) instead of owning a local FieldLayout to
   * gather from itself. mins/maxs/owns_nothing are flat, rank-major arrays as
   * produced by LocalOwnedBounds (owns_nothing entries are 0/1).
   */
  static DistributedBoundingBoxPartition FromGatheredBounds(
    int dim, Real padding, std::vector<Real> mins, std::vector<Real> maxs,
    std::vector<char> owns_nothing);

  /**
   * @brief Computes this rank's own bounding box of its owned DOF-holder
   * coordinates, with no MPI communication of its own.
   *
   * This is the per-rank information a field job ships to a standalone
   * evaluation server so that the server can assemble a routing partition
   * without a mesh or a partition file of its own (see
   * pcms::DistributedEvaluationChannel::PublishOwnedBounds). min_out/max_out
   * use the first `dim` entries and are only meaningful when owns_nothing_out
   * is false; on a rank that owns nothing they are zeroed.
   */
  static void LocalOwnedBounds(const FieldLayout& layout, int dim,
                               std::array<Real, 3>& min_out,
                               std::array<Real, 3>& max_out,
                               bool& owns_nothing_out);

  /**
   * @brief Appends the ranks whose padded bounding box contains point.
   *
   * Only the first dim coordinates of point are used; the rest are ignored.
   * Ranks that own nothing are never candidates. Does not clear
   * candidate_ranks first, so callers may accumulate across several points.
   */
  void GetCandidateRanks(const std::array<Real, 3>& point,
                         std::vector<int>& candidate_ranks) const;

private:
  DistributedBoundingBoxPartition(int dim, Real padding,
                                  std::vector<Real> mins,
                                  std::vector<Real> maxs,
                                  std::vector<char> owns_nothing);

  int num_ranks_ = 0;
  int dim_ = 0;
  Real padding_;
  // Flattened [rank * dim_ + d] bounds of each rank's owned region.
  std::vector<Real> mins_;
  std::vector<Real> maxs_;
  // owns_nothing_[rank] is true for a rank with no owned DOF holders, whose
  // (otherwise undefined) box must never match any point.
  std::vector<char> owns_nothing_;
};

} // namespace pcms

#endif // PCMS_DISTRIBUTED_PARTITION_H
