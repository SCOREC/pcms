#ifndef PCMS_FIELD_EVALUATOR_DISTRIBUTED_POINT_EVALUATOR_H
#define PCMS_FIELD_EVALUATOR_DISTRIBUTED_POINT_EVALUATOR_H

#include "pcms/field/evaluation_request.h"
#include "pcms/field/field.h"
#include "pcms/field/function_space.h"
#include "pcms/field/out_of_bounds_policy.h"
#include "pcms/field/point_evaluator.h"
#include "pcms/localization/distributed_point_search.h"
#include "pcms/localization/point_search.h"
#include "pcms/utility/arrays.h"
#include "pcms/utility/assert.h"
#include "pcms/utility/memory_spaces.h"
#include "pcms/utility/mpi_utilities.h"
#include "pcms/utility/uniform_grid.h"

#include <Kokkos_Core.hpp>
#include <Omega_h_mesh.hpp>

#include <algorithm>
#include <limits>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace pcms
{

namespace detail
{

/// The concrete local search for a spatial dimension.
template <int Dim>
struct GridSearchOf;
template <>
struct GridSearchOf<2>
{
  using type = GridPointSearch2D;
};
template <>
struct GridSearchOf<3>
{
  using type = GridPointSearch3D;
};

/// A search over the whole decomposition, so every rank builds the same grid.
///
/// The local evaluator factories build their search from the local element
/// count, which differs per rank. Routing needs the opposite: a grid every rank
/// agrees on, because the rank a point is routed to is read off that grid.
template <int Dim>
typename GridSearchOf<Dim>::type MakeGlobalGridSearch(Omega_h::Mesh& mesh)
{
  const auto n = DivisionsPerAxisForGlobalMesh<Dim>(mesh);
  if constexpr (Dim == 3) {
    return typename GridSearchOf<Dim>::type(mesh, n, n, n);
  } else {
    return typename GridSearchOf<Dim>::type(mesh, n, n);
  }
}

/// The reading an owner writes for a point its own evaluator could not place.
///
/// Owner-side evaluation has to say "no answer" inside the payload, because the
/// reply is sent before the requester can be told anything else, and because
/// the local evaluator API has no way to report an unplaced point other than
/// through the value it wrote. The out-of-bounds policy supplies the marker: an
/// owner evaluates with FILL and a fill value nothing else can produce.
///
/// An exception would be the obvious alternative and is not available: an owner
/// evaluates while its peers are blocked waiting for its reply, so throwing
/// there would leave the other ranks in the exchange, not raise an error.
///
/// The cost is that a field value equal to the marker is read as "no answer".
/// NaN cannot be a coefficient of a meaningful field, and the integral marker
/// is the type's lowest value, so both are safe in practice; it is the reason
/// the confidence is in the test rather than in the type system.
template <typename T>
[[nodiscard]] inline T UnansweredMarker()
{
  if constexpr (std::is_floating_point_v<T>) {
    return std::numeric_limits<T>::quiet_NaN();
  } else {
    return std::numeric_limits<T>::lowest();
  }
}

/// Whether `value` is the marker. Device-callable; `marker` is passed in
/// because `std::numeric_limits` is host-only.
template <typename T>
[[nodiscard]] KOKKOS_INLINE_FUNCTION bool IsUnanswered(T value, T marker)
{
  if constexpr (std::is_floating_point_v<T>) {
    (void)marker;
    return value != value; // NaN, without <cmath>'s overload set
  } else {
    return value == marker;
  }
}

using DeviceExec = DeviceMemorySpace::execution_space;

/// Set every entry of `values` to `fill`.
template <typename T, typename Values>
void FillAll(Values values, T fill)
{
  const int nc = static_cast<int>(values.extent(1));
  Kokkos::parallel_for(
    "DistributedPointEvaluator::Fill",
    Kokkos::RangePolicy<DeviceExec>(0, values.extent(0)), KOKKOS_LAMBDA(LO i) {
      for (int c = 0; c < nc; ++c) {
        values(i, c) = fill;
      }
    });
}

/// Write row `k` of `rows` to row `targets(k)` of `values`, skipping a row that
/// holds the no-answer marker. Returns how many rows were skipped.
template <typename T, typename Values>
LO ScatterAnswered(Values values, Kokkos::View<T**, DeviceMemorySpace> rows,
                   Kokkos::View<const LO*, DeviceMemorySpace> targets)
{
  const int nc = static_cast<int>(values.extent(1));
  const T marker = UnansweredMarker<T>();
  LO skipped = 0;
  Kokkos::parallel_reduce(
    "DistributedPointEvaluator::Scatter",
    Kokkos::RangePolicy<DeviceExec>(0, targets.extent(0)),
    KOKKOS_LAMBDA(LO k, LO & n_skipped) {
      for (int c = 0; c < nc; ++c) {
        if (IsUnanswered(rows(k, c), marker)) {
          ++n_skipped;
          return;
        }
      }
      const LO i = targets(k);
      for (int c = 0; c < nc; ++c) {
        values(i, c) = rows(k, c);
      }
    },
    skipped);
  return skipped;
}

/// `counts` times `factor`, for sending `factor` scalars per point.
inline std::vector<int> Scaled(const std::vector<int>& counts, int factor)
{
  std::vector<int> out(counts.size());
  for (size_t r = 0; r < counts.size(); ++r) {
    out[r] = counts[r] * factor;
  }
  return out;
}

} // namespace detail

/**
 * \brief Evaluate a field at points that are distributed over ranks.
 *
 * Any rank can evaluate the field at any point, including points outside its
 * own partition, and every rank gets the same answer.
 *
 * Construction (collective, runs once because the query points are fixed):
 *  1. `DistributedPointSearch::Locate` finds the owning rank of every point:
 *     local search first, then the global search grid gives candidate ranks,
 *     which search the point locally; the lowest owning rank wins.
 *  2. Each point owned by another rank is sent to its owner, once.
 *  3. Each rank builds two ordinary local evaluators through
 *     `FunctionSpace::CreatePointEvaluator`: one for its own points, and one
 *     for the points other ranks sent it.
 *
 * `Evaluate` (collective): each rank evaluates both batches locally, and one
 * `MPI_Alltoallv` sends the values for the incoming batch back to the ranks
 * that asked. Locations are never searched again; values are recomputed on
 * every call, from the coefficients as they are now.
 *
 * Nothing here is specific to a field type: a backend gets distributed
 * evaluation as soon as it can evaluate locally.
 *
 * A rank only ever evaluates points it owns, so a ghost element's (possibly
 * stale) coefficients never enter an answer. If an owner's local evaluator
 * still cannot place a point -- its own search disagreeing with the routing
 * search -- it writes `detail::UnansweredMarker<T>()` and the point is treated
 * as not found.
 *
 * \tparam T Value type. Each backend accepts `CreatePointEvaluator<T>` for
 *         whichever types it supports.
 * \tparam Dim Spatial dimension of the mesh, 2 or 3.
 */
template <typename T, int Dim,
          typename LayoutPolicy =
            detail::default_layout_for_memory_space_t<DeviceMemorySpace>>
class DistributedPointEvaluator : public PointEvaluator<T, LayoutPolicy>
{
public:
  /**
   * \param mesh The mesh points are located on, and whose ownership decides who
   *        answers for a point. Must outlive this object.
   * \param space The field's function space. Every local evaluation goes
   *        through it.
   * \param points Query coordinates, `[num_points][Dim]`, fixed for the life of
   *        the evaluator.
   * \param policy What to do with a point no rank can locate.
   * \param tolerance_reach How far outside an element the search may still
   *        accept a point.
   * \param grid_divisions Cells per axis of the global search grid; 0 derives
   *        it from the mesh. Must be the same on every rank.
   */
  DistributedPointEvaluator(
    Omega_h::Mesh& mesh, std::shared_ptr<const FunctionSpace> space,
    Kokkos::View<const Real* [Dim], DeviceMemorySpace> points,
    OutOfBoundsPolicy policy = {}, Real tolerance_reach = 1e-12,
    LO grid_divisions = 0)
    : mesh_(mesh), space_(std::move(space)), policy_(policy)
  {
    if (!space_) {
      throw std::invalid_argument(
        "DistributedPointEvaluator: function space is null");
    }
    if (mesh_.dim() != Dim) {
      throw std::invalid_argument(
        "DistributedPointEvaluator: mesh dimension does not match");
    }
    n_components_ = space_->GetLayout()->GetNumComponents();
    n_points_ = static_cast<LO>(points.extent(0));
    Setup(points, tolerance_reach, grid_divisions);
  }

  /// Evaluate the field at every query point. Collective.
  ///
  /// \param field Supplies the coefficients. Must be compatible with this
  ///        evaluator's function space.
  /// \param values Output, `[num_query_points][num_components]`.
  void Evaluate(
    const Field<T>& field,
    Rank2View<T, DeviceMemorySpace, LayoutPolicy> values) const override;

  /// How many query points get their value from another rank.
  [[nodiscard]] std::size_t LastRoutedPointCount() const noexcept
  {
    return remote_points_.extent(0);
  }

  /// The space whose local evaluator produces the values.
  [[nodiscard]] const FunctionSpace& GetFunctionSpace() const noexcept
  {
    return *space_;
  }

private:
  void Setup(Kokkos::View<const Real* [Dim], DeviceMemorySpace> points,
             Real tolerance_reach, LO grid_divisions);

  /// A local evaluator over `coords`, or null when there are none. `coords`
  /// must stay alive for as long as the evaluator.
  [[nodiscard]] std::unique_ptr<PointEvaluator<T>> MakeLocalEvaluator(
    Kokkos::View<Real**, DeviceMemorySpace> coords) const;

  /// Run `evaluator` over its `n` points, on the device.
  [[nodiscard]] Kokkos::View<T**, DeviceMemorySpace> EvaluateLocally(
    const PointEvaluator<T>* evaluator, LO n, const Field<T>& field) const;

  Omega_h::Mesh& mesh_;
  std::shared_ptr<const FunctionSpace> space_;
  OutOfBoundsPolicy policy_;
  int n_components_ = 0;
  LO n_points_ = 0;

  // Query points this rank owns (indices into the query set), their
  // coordinates, and the evaluator over them.
  Kokkos::View<LO*, DeviceMemorySpace> own_points_;
  Kokkos::View<Real**, DeviceMemorySpace> own_coords_;
  std::unique_ptr<PointEvaluator<T>> own_evaluator_;

  // Query points another rank owns, grouped by owner rank in ascending order,
  // and how many go to each rank.
  Kokkos::View<LO*, DeviceMemorySpace> remote_points_;
  std::vector<int> remote_counts_;

  // Points other ranks sent this rank to evaluate, grouped by requesting rank,
  // and how many came from each rank.
  std::vector<int> incoming_counts_;
  Kokkos::View<Real**, DeviceMemorySpace> incoming_coords_;
  std::unique_ptr<PointEvaluator<T>> incoming_evaluator_;
};

template <typename T, int Dim, typename LayoutPolicy>
void DistributedPointEvaluator<T, Dim, LayoutPolicy>::Setup(
  Kokkos::View<const Real* [Dim], DeviceMemorySpace> points,
  Real tolerance_reach, LO grid_divisions)
{
  // ---- 1. find the owning rank of every point ----
  // The search is only needed here: once the owners are known, `Evaluate`
  // never searches again.
  const auto local_search = detail::MakeGlobalGridSearch<Dim>(mesh_);
  DistributedPointSearch<Dim> search(mesh_, local_search, tolerance_reach,
                                     grid_divisions);
  const auto located = search.Locate(points);

  // Split the query indices by owner. The search result is host data, so this
  // is host work on indices only; coordinates stay on the device.
  const int rank = search.Rank();
  const int nranks = mesh_.comm()->size();
  std::vector<std::vector<LO>> by_owner(static_cast<size_t>(nranks));
  std::vector<LO> own;
  for (LO i = 0; i < n_points_; ++i) {
    const auto& r = located[static_cast<size_t>(i)];
    if (r.status != DistributedPointSearch<Dim>::ResultStatus::LOCATED) {
      continue; // in neither list; `Evaluate` fills it
    }
    if (r.owner_rank == rank) {
      own.push_back(i);
    } else {
      by_owner[static_cast<size_t>(r.owner_rank)].push_back(i);
    }
  }
  std::vector<LO> remote;
  remote_counts_.assign(static_cast<size_t>(nranks), 0);
  for (int r = 0; r < nranks; ++r) {
    const auto& points_of_r = by_owner[static_cast<size_t>(r)];
    remote.insert(remote.end(), points_of_r.begin(), points_of_r.end());
    remote_counts_[static_cast<size_t>(r)] =
      static_cast<int>(points_of_r.size());
  }
  own_points_ = CopyToDevice(own);
  remote_points_ = CopyToDevice(remote);

  // ---- 2. send each remote point to its owner, once ----
  // Only the remote points' coordinates leave the device, because MPI needs
  // them on the host.
  const auto send =
    CopyToHostRowMajor(GatherRows(points, remote_points_));
  std::vector<int> incoming_scalar_counts;
  const auto incoming =
    ExchangeRecords(mesh_.comm()->get_impl(), detail::Scaled(remote_counts_, Dim),
                    send, incoming_scalar_counts);
  incoming_counts_.resize(incoming_scalar_counts.size());
  for (size_t r = 0; r < incoming_counts_.size(); ++r) {
    incoming_counts_[r] = incoming_scalar_counts[r] / Dim;
  }

  // ---- 3. one local evaluator per batch, built once ----
  own_coords_ = GatherRows(points, own_points_); // device only
  incoming_coords_ = CopyRowMajorToDevice(incoming, Dim);
  own_evaluator_ = MakeLocalEvaluator(own_coords_);
  incoming_evaluator_ = MakeLocalEvaluator(incoming_coords_);
}

template <typename T, int Dim, typename LayoutPolicy>
std::unique_ptr<PointEvaluator<T>>
DistributedPointEvaluator<T, Dim, LayoutPolicy>::MakeLocalEvaluator(
  Kokkos::View<Real**, DeviceMemorySpace> coords) const
{
  if (coords.extent(0) == 0) {
    return nullptr;
  }
  // Only owned points reach this evaluator, so one it cannot place is a
  // disagreement between its search and the routing search. FILL with a
  // marker lets `Evaluate` treat such a point as not found.
  OutOfBoundsPolicy fallback;
  fallback.mode = OutOfBoundsMode::FILL;
  fallback.fill_value = static_cast<Real>(detail::UnansweredMarker<T>());

  const CoordinateView<DeviceMemorySpace> view(space_->GetCoordinateSystem(),
                                               MakeConstRank2View(coords));
  return space_->CreatePointEvaluator<T>(
    EvaluationRequest::FromCoordinates(view, fallback));
}

template <typename T, int Dim, typename LayoutPolicy>
Kokkos::View<T**, DeviceMemorySpace>
DistributedPointEvaluator<T, Dim, LayoutPolicy>::EvaluateLocally(
  const PointEvaluator<T>* evaluator, LO n, const Field<T>& field) const
{
  Kokkos::View<T**, DeviceMemorySpace> values("local values", n,
                                              n_components_);
  if (evaluator != nullptr) {
    evaluator->Evaluate(field, MakeRank2View(values));
  }
  return values;
}

template <typename T, int Dim, typename LayoutPolicy>
void DistributedPointEvaluator<T, Dim, LayoutPolicy>::Evaluate(
  const Field<T>& field,
  Rank2View<T, DeviceMemorySpace, LayoutPolicy> values) const
{
  PCMS_ALWAYS_ASSERT(values.extent(0) == static_cast<size_t>(n_points_));
  PCMS_ALWAYS_ASSERT(values.extent(1) == static_cast<size_t>(n_components_));
  const int nc = n_components_;
  const auto n_own = static_cast<LO>(own_points_.extent(0));
  const auto n_remote = static_cast<LO>(remote_points_.extent(0));

  // Owner side: evaluate this rank's own points and the points sent to it.
  const auto own_values = EvaluateLocally(own_evaluator_.get(), n_own, field);
  const auto incoming_values = EvaluateLocally(
    incoming_evaluator_.get(), static_cast<LO>(incoming_coords_.extent(0)),
    field);

  // Route the incoming values back to the ranks that asked: the only values
  // that cross to the host and back, because MPI takes host buffers. One scalar per component per point.
  const auto remote_host = ExchangeRecordsWithKnownCounts(
    mesh_.comm()->get_impl(), detail::Scaled(incoming_counts_, nc),
    CopyToHostRowMajor(incoming_values), detail::Scaled(remote_counts_, nc));
  const auto remote_values = CopyRowMajorToDevice(remote_host, nc);

  // Assemble on the device. Every point starts at the policy's fill value, and
  // a point is overwritten only when its owner answered for it.
  detail::FillAll(values, static_cast<T>(policy_.fill_value));
  LO missing = n_points_ - n_own - n_remote; // located by no rank
  missing += detail::ScatterAnswered(values, own_values, own_points_);
  missing += detail::ScatterAnswered(values, remote_values, remote_points_);

  // Thrown after the exchange, so every rank has completed it.
  if (policy_.mode == OutOfBoundsMode::ERROR && missing > 0) {
    throw pcms_error(
      "DistributedPointEvaluator: a query point lies outside the mesh and the "
      "out-of-bounds policy is ERROR");
  }
}

/**
 * \brief Make a distributed evaluator for `space`, dispatching on the mesh
 * dimension. Collective.
 *
 * \param points Query coordinates, `[num_points][mesh_dim]`; copied into the
 *        layout the evaluator uses.
 * \param grid_divisions Cells per axis of the global search grid; 0 derives it
 *        from the mesh. Must be the same on every rank.
 *
 * The returned evaluator is an ordinary `PointEvaluator<T>`, so it can be used
 * wherever a serial one can.
 */
template <typename T>
[[nodiscard]] std::unique_ptr<PointEvaluator<T>> MakeDistributedPointEvaluator(
  Omega_h::Mesh& mesh, std::shared_ptr<const FunctionSpace> space,
  Kokkos::View<const Real**, DeviceMemorySpace> points,
  OutOfBoundsPolicy policy = {}, Real tolerance_reach = 1e-12,
  LO grid_divisions = 0)
{
  const int mesh_dim = mesh.dim();
  const auto n_points = static_cast<LO>(points.extent(0));
  if (mesh_dim != 2 && mesh_dim != 3) {
    throw std::invalid_argument(
      "MakeDistributedPointEvaluator: only 2D and 3D meshes are supported");
  }
  if (points.extent(1) != static_cast<std::size_t>(mesh_dim)) {
    throw std::invalid_argument(
      "MakeDistributedPointEvaluator: query coordinates are not [num_points x "
      "mesh_dim]");
  }

  if (mesh_dim == 2) {
    Kokkos::View<Real* [2], DeviceMemorySpace> coords("distributed coords",
                                                      n_points);
    Kokkos::parallel_for(
      "copy_distributed_coords",
      Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, n_points),
      KOKKOS_LAMBDA(LO i) {
        coords(i, 0) = points(i, 0);
        coords(i, 1) = points(i, 1);
      });
    return std::make_unique<DistributedPointEvaluator<T, 2>>(
      mesh, std::move(space), coords, policy, tolerance_reach, grid_divisions);
  }
  Kokkos::View<Real* [3], DeviceMemorySpace> coords("distributed coords",
                                                    n_points);
  Kokkos::parallel_for(
    "copy_distributed_coords",
    Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, n_points),
    KOKKOS_LAMBDA(LO i) {
      coords(i, 0) = points(i, 0);
      coords(i, 1) = points(i, 1);
      coords(i, 2) = points(i, 2);
    });
  return std::make_unique<DistributedPointEvaluator<T, 3>>(
    mesh, std::move(space), coords, policy, tolerance_reach, grid_divisions);
}

/// How many query points of a distributed evaluator get their value from
/// another rank, or 0 when `evaluator` is not one.
///
/// Diagnostics for callers that hold only the `PointEvaluator` interface, which
/// is the interface a distributed evaluator is meant to be used through.
template <typename T>
[[nodiscard]] std::size_t LastRoutedPointCount(
  const PointEvaluator<T>& evaluator)
{
  if (auto* d2 =
        dynamic_cast<const DistributedPointEvaluator<T, 2>*>(&evaluator)) {
    return d2->LastRoutedPointCount();
  }
  if (auto* d3 =
        dynamic_cast<const DistributedPointEvaluator<T, 3>*>(&evaluator)) {
    return d3->LastRoutedPointCount();
  }
  return 0;
}

} // namespace pcms

#endif // PCMS_FIELD_EVALUATOR_DISTRIBUTED_POINT_EVALUATOR_H