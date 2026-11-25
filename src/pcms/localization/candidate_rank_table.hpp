#ifndef PCMS_LOCALIZATION_CANDIDATE_RANK_TABLE_HPP
#define PCMS_LOCALIZATION_CANDIDATE_RANK_TABLE_HPP

#include "pcms/localization/overlap.hpp"
#include "pcms/utility/bounding_box.hpp"
#include "pcms/utility/mpi_utilities.hpp"
#include "pcms/utility/types.hpp"
#include "pcms/utility/uniform_grid.hpp"

#include <mpi.h>

#include <optional>
#include <vector>

namespace pcms
{

/**
 * \brief Which ranks may own the element containing a query point.
 *
 * Make sure if a point is contained by an element that `rank` owns, then
 * `rank` is among the candidates for that point. The converse does not hold: a
 * candidate list is a superset of the ranks that can actually localize the
 * point, so the caller must still search and merge the answers.
 */
template <int Dim>
class CandidateRankTable
{
public:
  static constexpr int dim = Dim;
  using GridT = UniformGrid<Dim>;

  /**
   * Publishes this rank's box and gathers the others'. Collective.
   *
   * \param comm The communicator of the decomposition.
   * \param grid The global search grid; must be identical on every rank.
   * \param local_bounds A box containing this rank's part of the domain, or
   *        none when the rank holds nothing (it is then never a candidate).
   */
  CandidateRankTable(MPI_Comm comm, const GridT& grid,
                     const std::optional<AABBox<Dim>>& local_bounds)
    : grid_(grid)
  {
    int nranks = 0;
    MPI_Comm_size(comm, &nranks);

    // Every rank needs every other rank's box, so this is an all-ranks gather
    // of a fixed 2*Dim values per rank; a negative half width marks a rank
    // with no box. The box is grown by one cell, so a point the local search
    // accepts just outside this rank's part still routes here.
    // \see pcms/utility/mpi_utilities.hpp
    Kokkos::Array<Real, 2 * static_cast<size_t>(Dim)> local{};
    for (int d = 0; d < Dim; ++d) {
      const Real cell_width = grid_.edge_length[d] / grid_.divisions[d];
      local[static_cast<size_t>(d)] =
        local_bounds ? local_bounds->center[d] : 0;
      local[static_cast<size_t>(Dim) + d] =
        local_bounds ? local_bounds->half_width[d] + cell_width : -1;
    }
    const auto gathered = AllGatherAllRanks(comm, local.data(), 2 * Dim);

    boxes_.assign(static_cast<size_t>(nranks), std::nullopt);
    for (int r = 0; r < nranks; ++r) {
      const auto* values = &gathered[static_cast<size_t>(r) * 2 * Dim];
      if (values[Dim] < 0) {
        continue;
      }
      AABBox<Dim> box{};
      for (int d = 0; d < Dim; ++d) {
        box.center[d] = values[d];
        box.half_width[d] = values[Dim + d];
      }
      boxes_[static_cast<size_t>(r)] = box;
    }
  }

  /** The grid a point's cell is taken in. Identical on every rank. */
  [[nodiscard]] const GridT& GetGrid() const noexcept { return grid_; }

  /** Number of ranks in the decomposition, i.e. the size of the table. */
  [[nodiscard]] int NumRanks() const noexcept
  {
    return static_cast<int>(boxes_.size());
  }

  /** The grown box `rank` published, or none when it holds nothing. */
  [[nodiscard]] const std::optional<AABBox<Dim>>& BoundsForRank(int rank) const
  {
    return boxes_.at(static_cast<size_t>(rank));
  }

  /** Ranks that may own the element containing `point` (`Dim` coordinates). */
  template <typename P>
  [[nodiscard]] std::vector<int> CandidatesForPoint(const P& point) const
  {
    const auto cell = grid_.GetCellBBOX(grid_.ClosestCellID(point));
    std::vector<int> candidates;
    for (size_t r = 0; r < boxes_.size(); ++r) {
      if (boxes_[r] && Overlaps(cell, *boxes_[r])) {
        candidates.push_back(static_cast<int>(r));
      }
    }
    return candidates;
  }

  /** True when at least one rank may own the element containing `point`. */
  template <typename P>
  [[nodiscard]] bool HasCandidates(const P& point) const
  {
    return !CandidatesForPoint(point).empty();
  }

private:
  GridT grid_{};
  std::vector<std::optional<AABBox<Dim>>> boxes_;
};

} // namespace pcms

#endif // PCMS_LOCALIZATION_CANDIDATE_RANK_TABLE_HPP
