#ifndef PCMS_FIELD_EVALUATOR_DISTRIBUTED_POINT_EVALUATOR_H
#define PCMS_FIELD_EVALUATOR_DISTRIBUTED_POINT_EVALUATOR_H

#include "pcms/field/evaluation_request.h"
#include "pcms/field/field.h"
#include "pcms/field/function_space.h"
#include "pcms/field/out_of_bounds_policy.h"
#include "pcms/field/point_evaluator.h"
#include "pcms/field/value_view.hpp"
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

/// Write row `k` of `rows` to row `targets(k)` of `values`.
template <typename T, typename Values>
void ScatterRows(Values values, Kokkos::View<T**, DeviceMemorySpace> rows,
                 Kokkos::View<const LO*, DeviceMemorySpace> targets)
{
  const int nc = static_cast<int>(values.extent(1));
  Kokkos::parallel_for(
    "DistributedPointEvaluator::Scatter",
    Kokkos::RangePolicy<DeviceExec>(0, targets.extent(0)), KOKKOS_LAMBDA(LO k) {
      const LO i = targets(k);
      for (int c = 0; c < nc; ++c) {
        values(i, c) = rows(k, c);
      }
    });
}

/// Whether each of the first `n` rows appears in `filled`, as host flags.
inline std::vector<char> FilledFlags(
  Kokkos::View<const LO*, DeviceMemorySpace> filled, LO n)
{
  std::vector<char> flags(static_cast<size_t>(n), 0);
  const auto filled_h =
    Kokkos::create_mirror_view_and_copy(HostMemorySpace{}, filled);
  for (size_t k = 0; k < filled_h.extent(0); ++k) {
    flags[static_cast<size_t>(filled_h(k))] = 1;
  }
  return flags;
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
 * Precondition: the field's ghost DOF values are synchronized
 * (`Field::SynchronizeGhosts`). A rank evaluates only points in elements it
 * owns, but an owned element can have vertices another rank owns, and their
 * values on this rank are ghost copies. `FieldCommunicator::Receive`
 * synchronizes after receiving; a field written some other way, with only its
 * owned values set, must be synchronized before `Evaluate`.
 *
 * A point no rank can locate, or one its owner's local evaluator still cannot
 * place (that evaluator's search disagreeing with the routing search), is an
 * unanswered point. Which points those are is settled at construction: under
 * FILL they get the fill value and are reported by `FilledPoints()`; under
 * ERROR, `Evaluate` throws after its exchange, so no rank is left blocked.
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
    ValueView<T, DeviceMemorySpace, LayoutPolicy> values) const override;

  /// The query points that get the fill value, ascending. Empty unless the
  /// policy is FILL. Fixed at construction.
  [[nodiscard]] Kokkos::View<const LO*, DeviceMemorySpace> FilledPoints()
    const override
  {
    if (policy_.mode != OutOfBoundsMode::FILL) {
      return {};
    }
    return unanswered_points_;
  }

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

  // Query points nothing answers for, ascending: never located, or located
  // but not placed by the owner's local evaluator.
  Kokkos::View<LO*, DeviceMemorySpace> unanswered_points_;
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

  // ---- 4. which points nothing answers for ----
  // Each owner tells the requester which of its points its local evaluator
  // filled, once, so `Evaluate` never has to decode that from values.
  const auto n_incoming = static_cast<LO>(incoming_coords_.extent(0));
  const auto incoming_filled =
    incoming_evaluator_ ? detail::FilledFlags(incoming_evaluator_->FilledPoints(),
                                              n_incoming)
                        : std::vector<char>(static_cast<size_t>(n_incoming), 0);
  const auto remote_filled = ExchangeRecordsWithKnownCounts(
    mesh_.comm()->get_impl(), incoming_counts_, incoming_filled,
    remote_counts_);
  const auto own_filled =
    own_evaluator_
      ? detail::FilledFlags(own_evaluator_->FilledPoints(),
                            static_cast<LO>(own.size()))
      : std::vector<char>{};

  std::vector<char> answered(static_cast<size_t>(n_points_), 0);
  for (size_t k = 0; k < own.size(); ++k) {
    answered[static_cast<size_t>(own[k])] = own_filled[k] == 0;
  }
  for (size_t k = 0; k < remote.size(); ++k) {
    answered[static_cast<size_t>(remote[k])] = remote_filled[k] == 0;
  }
  std::vector<LO> unanswered;
  for (LO i = 0; i < n_points_; ++i) {
    if (answered[static_cast<size_t>(i)] == 0) {
      unanswered.push_back(i);
    }
  }
  unanswered_points_ = CopyToDevice(unanswered);
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
  // disagreement between its search and the routing search. It must not throw
  // or guess: FILL marks it, and its `FilledPoints()` reports it at setup.
  OutOfBoundsPolicy fallback;
  fallback.mode = OutOfBoundsMode::FILL;
  fallback.fill_value = policy_.fill_value;

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
    evaluator->Evaluate(field, ValueView<T, DeviceMemorySpace>(
                                 field.GetData().GetValueBasis(),
                                 MakeRank2View(values)));
  }
  return values;
}

template <typename T, int Dim, typename LayoutPolicy>
void DistributedPointEvaluator<T, Dim, LayoutPolicy>::Evaluate(
  const Field<T>& field,
  ValueView<T, DeviceMemorySpace, LayoutPolicy> out) const
{
  this->CheckEvaluateWriteTag("DistributedPointEvaluator::Evaluate", field,
                              out);
  const auto values = out.GetValues();
  PCMS_ALWAYS_ASSERT(values.extent(0) == static_cast<size_t>(n_points_));
  PCMS_ALWAYS_ASSERT(values.extent(1) == static_cast<size_t>(n_components_));
  const int nc = n_components_;
  const auto n_own = static_cast<LO>(own_points_.extent(0));

  // Owner side: evaluate this rank's own points and the points sent to it.
  const auto own_values = EvaluateLocally(own_evaluator_.get(), n_own, field);
  const auto incoming_values = EvaluateLocally(
    incoming_evaluator_.get(), static_cast<LO>(incoming_coords_.extent(0)),
    field);

  // Route the incoming values back to the ranks that asked: the only values
  // that cross to the host and back, because MPI takes host buffers.
  const auto remote_host = ExchangeRecordsWithKnownCounts(
    mesh_.comm()->get_impl(), detail::Scaled(incoming_counts_, nc),
    CopyToHostRowMajor(incoming_values), detail::Scaled(remote_counts_, nc));
  const auto remote_values = CopyRowMajorToDevice(remote_host, nc);

  // Assemble on the device. A point no rank located keeps the fill value; a
  // point an owner could not place already holds it, from the owner's FILL.
  detail::FillAll(values, static_cast<T>(policy_.fill_value));
  detail::ScatterRows(values, own_values, own_points_);
  detail::ScatterRows(values, remote_values, remote_points_);

  // Thrown after the exchange, so every rank has completed it.
  if (policy_.mode == OutOfBoundsMode::ERROR &&
      unanswered_points_.extent(0) > 0) {
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