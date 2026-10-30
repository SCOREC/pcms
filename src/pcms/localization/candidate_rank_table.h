#ifndef PCMS_LOCALIZATION_CANDIDATE_RANK_TABLE_H
#define PCMS_LOCALIZATION_CANDIDATE_RANK_TABLE_H

#include "pcms/localization/overlap.h"
#include "pcms/utility/bounding_box.h"
#include "pcms/utility/mpi_utilities.h"
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
 * A distributed point search cannot tell from a point's coordinates which rank
 * owns it -- that is the whole problem. This table answers the routing question
 * without a search: every rank publishes the range of global search cells that
 * its own geometry falls in, and a point is then routed to the ranks whose
 * published range contains the point's cell.
 *
 * The search structure is implicit. `UniformGrid` is a handful of scalars, so
 * every rank holds the entire global structure without materializing anything
 * per cell, and the only thing exchanged is one cell range per rank. Memory is
 * therefore O(ranks) per rank and does not grow with the grid resolution, which
 * is what lets the caller choose a dense grid for selectivity without a
 * corresponding memory cost.
 *
 * Make sure if a point is contained by an element that `rank` owns, then
 * `rank` is among the candidates for that point. The converse does not hold: a
 * candidate list is a superset of the ranks that can actually localize the
 * point, so the caller must still search and merge the answers.
 *
 * Two properties of the grid make the guarantee work:
 *
 * - **The grid must be identical on every rank.** A linear cell ID only names
 *   the same region of space everywhere if every rank derived the same grid.
 *   `CreateUniformGridFromMesh` reduces the bounding box globally, so that half
 *   is automatic, but the resolution must come from a global quantity too --
 *   `DivisionsPerAxisForGlobalMesh`, not `DivisionsPerAxisForMesh(nelems())`,
 *   since `nelems()` counts ghost copies and so varies with the partition. The
 *   convenience constructor uses the global form; the explicit-grid constructor
 *   exists so a caller that already owns a shared grid uses that one.
 *
 * - **The published range must be widened by the containment tolerance.** A
 *   point search accepts a point up to some distance outside an element (the
 *   vertex/edge tolerance bands in 2D, the facet plane band in 3D), so a point
 *   legitimately contained by an element can lie slightly outside that element
 *   and outside the unpadded bounding box that publishes it. Without the
 *   widening such a point would be routed to no rank at all and silently
 *   treated as out of bounds. `tolerance_reach` is that distance.
 *
 * The bounding box is taken over the rank's *local* coordinates, which include
 * ghost vertices. That is a superset of the owned region -- an owned element
 * may reference a vertex it does not own -- so it can only ever add candidate
 * ranks, never miss one. Keeping it simple this way also makes a rank with no
 * local geometry fall out naturally: it publishes an empty range and is never a
 * candidate.
 *
 * Queries are host-side and intended to be called once per point, on the subset
 * of points that failed to localize locally, after the search results have been
 * copied back from the device. A lookup scans the published ranges, so it is
 * O(ranks) per point; the work is bounded by, and much smaller than, the search
 * it replaces.
 *
 * Setup uses a gather rather than a scatter: every rank needs every other
 * rank's range, so this is an `allgather` of O(ranks) values, not an all-to-all
 * with per-rank payloads.
 */
template <int Dim>
class CandidateRankTable
{
public:
  static constexpr int dim = Dim;
  using GridT = UniformGrid<Dim>;
  // Qualified: the target is the namespace-scope `pcms::CellRange`, which an
  // unqualified spelling would shadow with this very alias.
  using CellRange = pcms::CellRange<Dim>;

  /**
   * Build on the shared routing grid.
   *
   * The grid is derived from the global element count and the globally reduced
   * bounding box, so every rank constructs the same one.
   *
   * \param mesh The mesh whose elements are being located.
   * \param tolerance_reach How far outside an element the point search may
   * still accept a point; the published ranges are widened by this much. \see
   * PointLocalizationSearch::PointSearchTolerances
   */
  explicit CandidateRankTable(Omega_h::Mesh& mesh, Real tolerance_reach = 0.0)
    : CandidateRankTable(mesh, SharedRoutingGrid(mesh), tolerance_reach)
  {
  }

  /**
   * Build on a caller-supplied grid.
   *
   * Use this when the caller already holds the shared grid, so that the grid
   * used to route and the grid used for anything else are provably the same.
   *
   * \param mesh The mesh whose elements are being located.
   * \param grid The global search grid; must be identical on every rank.
   * \param tolerance_reach \see the convenience constructor.
   */
  CandidateRankTable(Omega_h::Mesh& mesh, const GridT& grid,
                     Real tolerance_reach = 0.0)
    : grid_(grid)
  {
    if (mesh.dim() != Dim) {
      throw std::invalid_argument(
        "CandidateRankTable: mesh dimension does not match the table");
    }

    auto comm = mesh.comm();
    const int nranks = comm->size();
    local_range_ = LocalCellRange(mesh, tolerance_reach);

    // Every rank needs every other rank's range, so this is an all-ranks gather
    // of a fixed 2*Dim values per rank. `Omega_h::Comm::allgather` cannot be
    // used: it only exchanges with the ranks of the mesh's communication graph.
    // \see pcms/utility/mpi_utilities.h
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

  /// The shared routing grid for `mesh`, identical on every rank.
  [[nodiscard]] static GridT SharedRoutingGrid(Omega_h::Mesh& mesh,
                                               LO min_divisions = 10)
  {
    return CreateUniformGridFromMesh<Dim>(
      mesh, DivisionsPerAxisForGlobalMesh<Dim>(mesh, min_divisions));
  }

  /// The routing grid for `mesh` with `divisions` cells per axis, or the
  /// default resolution (`SharedRoutingGrid`) when `divisions` is 0.
  ///
  /// `divisions` must be the same on every rank, or cell ids would name
  /// different regions of space on different ranks.
  [[nodiscard]] static GridT RoutingGrid(Omega_h::Mesh& mesh, LO divisions)
  {
    if (divisions > 0) {
      return CreateUniformGridFromMesh<Dim>(mesh, divisions);
    }
    return SharedRoutingGrid(mesh);
  }

  /// The grid a cell ID is interpreted in. Identical on every rank.
  [[nodiscard]] const GridT& GetGrid() const noexcept { return grid_; }

  /// Number of ranks in the decomposition, i.e. the size of the table.
  [[nodiscard]] int NumRanks() const noexcept
  {
    return static_cast<int>(ranges_.size());
  }

  /// The range this rank published. Empty when it holds no geometry.
  [[nodiscard]] const CellRange& LocalRange() const noexcept
  {
    return local_range_;
  }

  /// The range `rank` published.
  [[nodiscard]] const CellRange& RangeForRank(int rank) const
  {
    return ranges_.at(static_cast<size_t>(rank));
  }

  /// All published ranges, indexed by rank.
  [[nodiscard]] const std::vector<CellRange>& Ranges() const noexcept
  {
    return ranges_;
  }

  /// Ranks that may own the element containing a point in `cell`, ascending.
  ///
  /// Empty when no rank can own the cell: a cell in a hole of a non-convex
  /// domain, or outside the decomposition entirely. That is a legitimate answer
  /// and the caller treats it as out of bounds.
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

  /// Ranks that may own the element containing `point`, ascending.
  ///
  /// The point is mapped to its cell first, clamped to the grid, so a point
  /// outside the domain is routed by its nearest boundary cell rather than
  /// being rejected here. That is deliberate: a point just outside the domain
  /// can still be accepted by a boundary element within the containment
  /// tolerance, and it is that tolerance the published ranges were widened by.
  [[nodiscard]] std::vector<int> CandidatesForPoint(
    const Omega_h::Vector<Dim>& point) const
  {
    return CandidatesForCell(grid_.ClosestCellID(point));
  }

  /// True when at least one rank may own the element containing `point`.
  [[nodiscard]] bool HasCandidates(const Omega_h::Vector<Dim>& point) const
  {
    return !CandidatesForPoint(point).empty();
  }

private:
  /// Cell range covering this rank's local geometry, widened by the containment
  /// tolerance. Empty when the rank holds no coordinates.
  [[nodiscard]] CellRange LocalCellRange(Omega_h::Mesh& mesh,
                                         Real tolerance_reach) const
  {
    // A rank with no local coordinates owns no elements either, so there is
    // nothing for a query to find and the range stays empty. This also keeps
    // `find_bounding_box` away from a rank whose reduced bounding box would
    // come back inverted: it starts from plus/minus the largest representable
    // value, so an empty point set yields max < min.
    if (mesh.nverts() == 0) {
      return CellRange::Empty();
    }

    // Local coordinates, not owned ones: an owned element can reference a
    // vertex it does not own, and including those vertices can only widen the
    // range, which is safe. Never allreduces, so this stays per-rank.
    const auto bbox = Omega_h::find_bounding_box<Dim>(mesh.coords());

    AABBox<Dim> box{};
    for (int d = 0; d < Dim; ++d) {
      box.center[d] = (bbox.min[d] + bbox.max[d]) / 2.0;
      // Widen by the tolerance so that a point the search accepts just outside
      // this rank's geometry still routes here. \see the class comment.
      box.half_width[d] = (bbox.max[d] - bbox.min[d]) / 2.0 + tolerance_reach;
    }
    // Cells that touch the box count, which is what `candidate_cells` encodes:
    // a vertex sitting on a cell boundary makes the touched cell a candidate.
    return candidate_cells<Dim, TriangleVsBbox>(box, grid_);
  }

  GridT grid_{};
  CellRange local_range_{};
  std::vector<CellRange> ranges_;
};

} // namespace pcms

#endif // PCMS_LOCALIZATION_CANDIDATE_RANK_TABLE_H