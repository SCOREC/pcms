#ifndef PCMS_LOCALIZATION_CANDIDATE_RANK_TABLE_HPP
#define PCMS_LOCALIZATION_CANDIDATE_RANK_TABLE_HPP

#include "pcms/localization/overlap.h"
#include "pcms/utility/bounding_box.h"
#include "pcms/utility/mpi_utilities.hpp"
#include "pcms/utility/types.h"
#include "pcms/utility/uniform_grid.h"

#include <Omega_h_array.hpp>
#include <Omega_h_bbox.hpp>
#include <Omega_h_mesh.hpp>

#include <stdexcept>
#include <vector>

namespace pcms
{

/**
 * \brief Which ranks may own the element containing a query point.
 *
 * The search structure is implicit. `UniformGrid` is a handful of scalars, so
 * every rank holds the entire global structure without copying anything
 * per cell, and the only thing exchanged is one cell range per rank.
 *
 * Make sure if a point is contained by an element that `rank` owns, then
 * `rank` is among the candidates for that point. The converse does not hold: a
 * candidate list is a superset of the ranks that can actually localize the
 * point, so the caller must still search and merge the answers.
 *
 */
template <int Dim>
class CandidateRankTable
{
public:
  static constexpr int dim = Dim;
  using GridT = UniformGrid<Dim>;
  using CellRange = pcms::CellRange<Dim>;

  /**
   * Build on the shared routing grid.
   *
   * The grid is derived from the global element count and the globally reduced
   * bounding box, so every rank constructs the same one.
   *
   * \param mesh The mesh whose elements are being located.
   */
  explicit CandidateRankTable(Omega_h::Mesh& mesh)
    : CandidateRankTable(mesh, SharedRoutingGrid(mesh))
  {
  }

  /**
   * Build on a caller-supplied grid.
   *
   * \param mesh The mesh whose elements are being located.
   * \param grid The global search grid; must be identical on every rank.
   */
  CandidateRankTable(Omega_h::Mesh& mesh, const GridT& grid) : grid_(grid)
  {
    if (mesh.dim() != Dim) {
      throw std::invalid_argument(
        "CandidateRankTable: mesh dimension does not match the table");
    }

    auto comm = mesh.comm();
    const int nranks = comm->size();
    local_range_ = LocalCellRange(mesh);

    // Every rank needs every other rank's range, so this is an all-ranks gather
    // of a fixed 2*Dim values per rank. `Omega_h::Comm::allgather` cannot be
    // used: it only exchanges with the ranks of the mesh's communication graph.
    // \see pcms/utility/mpi_utilities.hpp
    Kokkos::Array<LO, 2 * static_cast<size_t>(Dim)> local{};
    for (int d = 0; d < Dim; ++d) {
      local[static_cast<size_t>(d)] = local_range_.lo[d];
      local[static_cast<size_t>(Dim) + d] = local_range_.hi[d];
    }
    const auto gathered =
      AllGatherAllRanks(comm->get_impl(), local.data(), 2 * Dim);

    ranges_.assign(static_cast<size_t>(nranks), CellRange{});
    for (int r = 0; r < nranks; ++r) {
      for (int d = 0; d < Dim; ++d) {
        ranges_[static_cast<size_t>(r)].lo[d] =
          gathered[static_cast<size_t>(r) * 2 * Dim + d];
        ranges_[static_cast<size_t>(r)].hi[d] =
          gathered[static_cast<size_t>(r) * 2 * Dim + Dim + d];
      }
    }
  }

  /** The shared routing grid for `mesh`, identical on every rank. */
  [[nodiscard]] static GridT SharedRoutingGrid(Omega_h::Mesh& mesh,
                                               LO min_divisions = 10)
  {
    return CreateUniformGridFromMesh<Dim>(
      mesh, DivisionsPerAxisForGlobalMesh<Dim>(mesh, min_divisions));
  }

  /**
   * The routing grid for `mesh` with `divisions` cells per axis, or the
   * default resolution (`SharedRoutingGrid`) when `divisions` is 0.
   *
   * `divisions` must be the same on every rank, or cell ids would name
   * different regions of space on different ranks.
   */
  [[nodiscard]] static GridT RoutingGrid(Omega_h::Mesh& mesh, LO divisions)
  {
    if (divisions > 0) {
      return CreateUniformGridFromMesh<Dim>(mesh, divisions);
    }
    return SharedRoutingGrid(mesh);
  }

  /** The grid a cell ID is interpreted in. Identical on every rank. */
  [[nodiscard]] const GridT& GetGrid() const noexcept { return grid_; }

  /** Number of ranks in the decomposition, i.e. the size of the table. */
  [[nodiscard]] int NumRanks() const noexcept
  {
    return static_cast<int>(ranges_.size());
  }

  /** The range this rank published. Empty when it holds no geometry. */
  [[nodiscard]] const CellRange& LocalRange() const noexcept
  {
    return local_range_;
  }

  /** The range `rank` published. */
  [[nodiscard]] const CellRange& RangeForRank(int rank) const
  {
    return ranges_.at(static_cast<size_t>(rank));
  }

  /** All published ranges, indexed by rank. */
  [[nodiscard]] const std::vector<CellRange>& Ranges() const noexcept
  {
    return ranges_;
  }

  /**
   * Ranks that may own the element containing a point in `cell`.
   *
   * Empty when no rank can own the cell.
   */
  [[nodiscard]] std::vector<int> CandidatesForCell(LO cell) const
  {
    const auto axis_indices = grid_.AxisIndicesFromCellIndex(cell);
    std::vector<int> candidates;
    candidates.reserve(ranges_.size());
    for (size_t r = 0; r < ranges_.size(); ++r) {
      if (ranges_[r].Contains(axis_indices)) {
        candidates.push_back(static_cast<int>(r));
      }
    }
    return candidates;
  }

  /** Ranks that may own the element containing `point`. */
  [[nodiscard]] std::vector<int> CandidatesForPoint(
    const Omega_h::Vector<Dim>& point) const
  {
    return CandidatesForCell(grid_.ClosestCellID(point));
  }

  /** True when at least one rank may own the element containing `point`. */
  [[nodiscard]] bool HasCandidates(const Omega_h::Vector<Dim>& point) const
  {
    return !CandidatesForPoint(point).empty();
  }

private:
  // Cell range covering this rank's local geometry, plus one cell on every
  // side. Empty when the rank holds no coordinates.
  [[nodiscard]] CellRange LocalCellRange(Omega_h::Mesh& mesh) const
  {
    // A rank with no local coordinates owns no elements either, so there is
    // nothing for a query to find and the range stays empty.
    if (mesh.nverts() == 0) {
      return CellRange::Empty();
    }

    const auto bbox = Omega_h::find_bounding_box<Dim>(mesh.coords());

    AABBox<Dim> box{};
    for (int d = 0; d < Dim; ++d) {
      box.center[d] = (bbox.min[d] + bbox.max[d]) / 2.0;
      box.half_width[d] = (bbox.max[d] - bbox.min[d]) / 2.0;
    }
    // Cells that touch the box count
    return candidate_cells<Dim, TriangleVsBbox>(box, grid_);
  }

  GridT grid_{};
  CellRange local_range_{};
  std::vector<CellRange> ranges_;
};

} // namespace pcms

#endif // PCMS_LOCALIZATION_CANDIDATE_RANK_TABLE_HPP