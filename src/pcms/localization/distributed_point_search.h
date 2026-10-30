#ifndef PCMS_LOCALIZATION_DISTRIBUTED_POINT_SEARCH_H
#define PCMS_LOCALIZATION_DISTRIBUTED_POINT_SEARCH_H

#include "pcms/localization/candidate_rank_table.h"
#include "pcms/localization/point_search.h"
#include "pcms/utility/arrays.h"
#include "pcms/utility/memory_spaces.h"
#include "pcms/utility/types.h"

#include <Kokkos_Core.hpp>
#include <Omega_h_mesh.hpp>

#include <cstddef>
#include <vector>

namespace pcms
{

/**
 * \brief Locate points on a mesh that is distributed over several ranks.
 *
 * The algorithm, for each query point:
 *  1. Search locally. A point inside an element this rank owns is done.
 *  2. Otherwise, look the point up in the global search grid
 *     (`CandidateRankTable`) to get the ranks that may own it.
 *  3. Send the point to those ranks; each searches it locally.
 *  4. Receive the answers. If several ranks own an element containing the
 *     point (it sits on a partition interface), the lowest rank wins.
 *
 * A rank answers only for elements it owns, never for ghosts, because ghost
 * data belongs to another rank.
 *
 * This class only locates points. Evaluating a field at them is the job of
 * `DistributedPointEvaluator`, which uses the result to know which rank to ask.
 */
template <int Dim>
class DistributedPointSearch
{
public:
  using LocalSearch = PointLocalizationSearch<Dim>;

  /// What this rank could determine about a point on its own.
  enum class LocalStatus
  {
    /// Contained by an element this rank owns, so this rank can answer. The
    /// search is finished for this point.
    LOCATED_HERE,
    /// Contained by an element this rank holds only as a ghost. The owner has
    /// to answer; routing must not send the point back to this rank.
    LOCATED_ELSEWHERE,
    /// Not contained by anything this rank holds. The owner, if there is one,
    /// is somewhere else.
    NOT_FOUND_HERE
  };

  struct LocalClassification
  {
    LocalStatus status = LocalStatus::NOT_FOUND_HERE;
    /// Local element id: the owning element for `LOCATED_HERE`, the ghost that
    /// contains the point for `LOCATED_ELSEWHERE`, `-1` otherwise.
    LO element_id = -1;
    /// Barycentric coordinates in `element_id`, zero when not found.
    Omega_h::Vector<Dim + 1> barycentric{};
  };

  /// The outcome of a distributed locate for one point.
  enum class ResultStatus
  {
    /// Some rank owns an element containing the point. `owner_rank`,
    /// `element_id` and `barycentric` describe it.
    LOCATED,
    /// No rank could place the point: it is outside the mesh, or outside every
    /// cell any rank published. The caller decides whether that is an error or
    /// a fill.
    NOT_FOUND
  };

  struct Result
  {
    ResultStatus status = ResultStatus::NOT_FOUND;
    /// Rank that owns the containing element, `-1` when not located.
    int owner_rank = -1;
    /// Global id of the containing element on `owner_rank`, `-1` when the mesh
    /// carries no global ids or the point was not located.
    GO global_element_id = -1;
    /// Local element id **on `owner_rank`**. Only meaningful together with
    /// `owner_rank`; it is not an index into this rank's mesh.
    LO element_id = -1;
    /// Barycentric coordinates inside `element_id`.
    Omega_h::Vector<Dim + 1> barycentric{};
  };

  /**
   * \param mesh The mesh to locate points on.
   * \param local_search A search over the same mesh. Must outlive this object.
   * \param tolerance_reach How far outside an element the local search may
   * still accept a point; the routing ranges are widened by this much.
   * \param grid_divisions Cells per axis of the global search grid. 0 (the
   * default) derives it from the global element count, about one element per
   * cell. Must be the same on every rank; the constructor checks this and
   * throws on every rank if not.
   *
   * Constructing this publishes this rank's cell range and gathers the others',
   * so the constructor is collective.
   */
  DistributedPointSearch(Omega_h::Mesh& mesh, const LocalSearch& local_search,
                         Real tolerance_reach = 0.0, LO grid_divisions = 0);

  /**
   * \brief Locate a batch of points, searching remote ranks where necessary.
   *
   * Collective: every rank must call it, because the exchange is a collective.
   * Points that this rank owns are answered locally at no communication cost;
   * only the rest are routed, and only to the ranks whose published cell range
   * covers the point. That keeps the common case -- a rank asking about points
   * in its own partition -- free of messages.
   *
   * A point outside the routing grid costs nothing at all. The grid spans the
   * global mesh bounding box and every element lies inside that box, so no rank
   * can contain such a point: it is reported `NOT_FOUND` without being routed.
   * The box is widened by the search tolerance, so a point the local search
   * would still accept just outside a boundary face is not rejected here.
   *
   * When several candidate ranks report a containing element -- a point on a
   * partition interface is contained by an element on each side -- the answer
   * from the **lowest rank id** wins. Any of them is geometrically right; the
   * rule just makes the outcome deterministic and identical on every rank.
   *
   * To evaluate at the same points repeatedly, keep the result rather than
   * locating again: `DistributedRoutingPlan` is built from it once and reused.
   *
   * \param points Query coordinates, one per row.
   * \return One result per row, in the same order.
   */
  [[nodiscard]] std::vector<Result> Locate(
    Kokkos::View<const Real* [Dim]> points);

  /// How many points the most recent `Locate` sent to another rank.
  ///
  /// Diagnostics. Zero means every point was answered here, either because this
  /// rank owns it or because it lies outside the mesh.
  [[nodiscard]] std::size_t LastRoutedPointCount() const noexcept
  {
    return last_routed_;
  }

  /**
   * \brief Classify a batch by what this rank can answer by itself.
   *
   * A point is `LOCATED_HERE` only when the search reports an element that this
   * rank owns. A point the search resolves to a ghost is `LOCATED_ELSEWHERE`
   * rather than `LOCATED_HERE` even though the ghost geometrically contains it:
   * the ghost's data belongs to another rank, so the owner has to answer. That
   * happens for a sizeable fraction of points near a partition interface --
   * roughly a quarter of element-centroid queries on a two-rank mesh -- and for
   * points on shared vertices or edges, where the reported element is resolved
   * through local adjacency and can land on either side.
   *
   * \param points Query coordinates, one per row.
   * \return One classification per row, in the same order.
   */
  [[nodiscard]] std::vector<LocalClassification> ClassifyLocally(
    Kokkos::View<const Real* [Dim]> points) const;

  /// The mesh this search locates points on.
  [[nodiscard]] Omega_h::Mesh& GetMesh() const noexcept { return mesh_; }

  /// The local search the distributed search is built on.
  [[nodiscard]] const LocalSearch& GetLocalSearch() const noexcept
  {
    return local_search_;
  }

  /// True when this rank holds ghost elements, i.e. when ownership can differ
  /// from "every local element".
  [[nodiscard]] bool HasGhostLayer() const noexcept { return has_ghost_layer_; }

  /// Whether this rank owns the local element `element_id`.
  [[nodiscard]] bool OwnsLocalElement(LO element_id) const
  {
    return element_id >= 0 &&
           static_cast<size_t>(element_id) < owned_elements_.size() &&
           owned_elements_[static_cast<size_t>(element_id)] != 0;
  }

  /// The routing table used to pick candidate ranks.
  [[nodiscard]] const CandidateRankTable<Dim>& GetRoutingTable() const noexcept
  {
    return routing_;
  }

  /// This rank's id in the mesh's communicator.
  [[nodiscard]] int Rank() const noexcept { return rank_; }

private:
  /// Global id of a local element, or -1 when the mesh carries no global ids.
  [[nodiscard]] GO GlobalElementId(LO element_id) const
  {
    if (global_ids_.empty() || element_id < 0 ||
        static_cast<size_t>(element_id) >= global_ids_.size()) {
      return -1;
    }
    return global_ids_[static_cast<size_t>(element_id)];
  }

  Omega_h::Mesh& mesh_;
  const LocalSearch& local_search_;
  CandidateRankTable<Dim> routing_;
  int rank_ = 0;
  /// Per-local-element ownership, indexed by local element id. All true when
  /// the mesh has no ghost layer, since then every local element is owned.
  std::vector<char> owned_elements_;
  /// Per-local-element global ids, empty when the mesh has no "global" tag.
  std::vector<GO> global_ids_;
  /// How far outside its geometry the local search still accepts a point.
  Real tolerance_reach_ = 0;
  bool has_ghost_layer_ = false;
  std::size_t last_routed_ = 0;
};

extern template class DistributedPointSearch<2>;
extern template class DistributedPointSearch<3>;

} // namespace pcms

#endif // PCMS_LOCALIZATION_DISTRIBUTED_POINT_SEARCH_H