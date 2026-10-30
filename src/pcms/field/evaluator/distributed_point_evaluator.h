#ifndef PCMS_FIELD_EVALUATOR_DISTRIBUTED_POINT_EVALUATOR_H
#define PCMS_FIELD_EVALUATOR_DISTRIBUTED_POINT_EVALUATOR_H

#include "pcms/field/coordinate_system.hpp"
#include "pcms/field/coordinate_view.hpp"
#include "pcms/field/evaluation_request.h"
#include "pcms/field/evaluator/omega_h_lagrange.h"
#include "pcms/field/field.h"
#include "pcms/field/function_space.h"
#include "pcms/field/out_of_bounds_policy.h"
#include "pcms/field/point_evaluator.h"
#include "pcms/field/point_evaluator_factory.hpp"
#include "pcms/field/value_view.hpp"
#include "pcms/localization/distributed_point_search.h"
#include "pcms/localization/distributed_routing_plan.h"
#include "pcms/localization/point_search.h"
#include "pcms/utility/arrays.h"
#include "pcms/utility/assert.h"
#include "pcms/utility/memory_spaces.h"
#include "pcms/utility/mpi_utilities.h"

#include <Kokkos_Core.hpp>
#include <Omega_h_mesh.hpp>

#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

namespace pcms
{

namespace detail
{

using DeviceExec = DeviceMemorySpace::execution_space;
using DeviceCoords =
  Rank2View<const Real, DeviceMemorySpace,
            default_layout_for_memory_space_t<DeviceMemorySpace>>;

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

/// Rows `rows(k)` of a coordinate view, as a new device matrix.
inline Kokkos::View<Real**, DeviceMemorySpace> GatherCoordinateRows(
  DeviceCoords src, Kokkos::View<const LO*, DeviceMemorySpace> rows)
{
  const int dim = static_cast<int>(src.extent(1));
  Kokkos::View<Real**, DeviceMemorySpace> out("gathered coords", rows.extent(0),
                                              dim);
  Kokkos::parallel_for(
    "DistributedPointEvaluator::GatherCoords",
    Kokkos::RangePolicy<DeviceExec>(0, rows.extent(0)), KOKKOS_LAMBDA(LO k) {
      for (int d = 0; d < dim; ++d) {
        out(k, d) = src(rows(k), d);
      }
    });
  return out;
}

/// A coordinate view copied into the fixed-dimension layout the search takes.
template <int Dim>
Kokkos::View<Real* [Dim], DeviceMemorySpace> ToSearchCoords(DeviceCoords src) {
  const auto n = static_cast<LO>(src.extent(0));
  Kokkos::View<Real* [Dim], DeviceMemorySpace> out("search coords", n);
  Kokkos::parallel_for(
    "DistributedPointEvaluator::ToSearchCoords",
    Kokkos::RangePolicy<DeviceExec>(0, n), KOKKOS_LAMBDA(LO i) {
      for (int d = 0; d < Dim; ++d) {
        out(i, d) = src(i, d);
      }
    });
  return out;
}

/// Whether each of the first `n` rows appears in `rows`, as host flags.
inline std::vector<char> FlagsOf(
  Kokkos::View<const LO*, DeviceMemorySpace> rows, LO n)
{
  std::vector<char> flags(static_cast<size_t>(n), 0);
  const auto rows_h =
    Kokkos::create_mirror_view_and_copy(HostMemorySpace{}, rows);
  for (size_t k = 0; k < rows_h.extent(0); ++k) {
    flags[static_cast<size_t>(rows_h(k))] = 1;
  }
  return flags;
}

/// `view` as a host vector.
inline std::vector<LO> ToHostVector(
  Kokkos::View<const LO*, DeviceMemorySpace> view)
{
  const auto h = Kokkos::create_mirror_view_and_copy(HostMemorySpace{}, view);
  return std::vector<LO>(h.data(), h.data() + h.extent(0));
}

/// Throws on every rank if `error` is non-empty on any rank, so a failure on
/// one rank never leaves the others blocked in the next collective.
inline void ThrowIfAnyRankFailed(const Omega_h::CommPtr& comm,
                                 const std::string& context,
                                 const std::string& error)
{
  const int failed =
    comm->allreduce(error.empty() ? 0 : 1, Omega_h_Op::OMEGA_H_MAX);
  if (failed != 0) {
    throw pcms_error(
      context + ": " +
      (error.empty() ? std::string("failed on another rank") : error));
  }
}

inline void CheckSameCoordinateSystem(
  const char* context, const CoordinateView<DeviceMemorySpace>& coords,
  const std::shared_ptr<const CoordinateSystem>& expected)
{
  if (!SameCoordinateSystem(coords.GetCoordinateSystem(), expected)) {
    throw pcms_error(std::string(context) + ": query coordinate system '" +
                     std::string(coords.GetCoordinateSystem()->Kind()) +
                     "' is not the factory's coordinate system '" +
                     std::string(expected->Kind()) + "'");
  }
}

} // namespace detail

/**
 * \brief Evaluate a field at points that are distributed over ranks.
 *
 * Built by `DistributedPointEvaluatorFactory`, from a routing plan and an
 * `EvaluationRequest`. Any rank can evaluate the field at any point, including
 * points outside its own partition.
 *
 * Construction (collective): each rank builds two ordinary local evaluators
 * through `FunctionSpace::CreatePointEvaluator`, one for the points it owns and
 * one for the points other ranks sent it. Their requests are derived from the
 * caller's:
 *  - The query layout is passed on only when this rank owns every one of its
 *    query points, because only then do the local points match the layout.
 *    Otherwise the local request carries coordinates only.
 *  - The local evaluators always use FILL. The caller's policy is applied here,
 *    after the exchange, so that ERROR throws on every rank together.
 *
 * `Evaluate` (collective): each rank evaluates both batches locally, and one
 * `MPI_Alltoallv` sends the values for the incoming batch back to the ranks
 * that asked. Points are never searched again.
 *
 * Precondition: the field's ghost DOF values are synchronized
 * (`Field::SynchronizeGhosts`). A rank evaluates only points in elements it
 * owns, but an owned element can have vertices another rank owns, and their
 * values on this rank are ghost copies.
 *
 * A point no rank can locate, or one its owner's local evaluator still cannot
 * place, is an unanswered point. Which points those are is settled at
 * construction: under FILL they get the fill value and are reported by
 * `FilledPoints()`; under ERROR, `Evaluate` throws after its exchange.
 */
template <typename T,
          typename LayoutPolicy =
            detail::default_layout_for_memory_space_t<DeviceMemorySpace>>
class DistributedPointEvaluator : public PointEvaluator<T, LayoutPolicy>
{
public:
  /// Collective. Prefer `DistributedPointEvaluatorFactory`, which builds the
  /// plan and checks that it matches the request.
  ///
  /// \param space The field's function space. Every local evaluation goes
  ///        through it.
  /// \param plan Routing for exactly the request's query points.
  /// \param request Query points, policy and optional query layout.
  DistributedPointEvaluator(std::shared_ptr<const FunctionSpace> space,
                            std::shared_ptr<const DistributedRoutingPlan> plan,
                            const EvaluationRequest& request)
    : space_(std::move(space)), plan_(std::move(plan)), policy_(request.policy)
  {
    if (!space_ || !plan_) {
      throw std::invalid_argument(
        "DistributedPointEvaluator: function space and plan must not be null");
    }
    n_components_ = space_->GetLayout()->GetNumComponents();
    Setup(request);
  }

  /// Evaluate the field at every query point. Collective.
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
    return plan_->RoutedPointCount();
  }

  [[nodiscard]] const DistributedRoutingPlan& GetRoutingPlan() const noexcept
  {
    return *plan_;
  }

  [[nodiscard]] const FunctionSpace& GetFunctionSpace() const noexcept
  {
    return *space_;
  }

private:
  void Setup(const EvaluationRequest& request);

  /// Run `evaluator` over its `n` points, on the device.
  [[nodiscard]] Kokkos::View<T**, DeviceMemorySpace> EvaluateLocally(
    const PointEvaluator<T>* evaluator, LO n, const Field<T>& field) const;

  std::shared_ptr<const FunctionSpace> space_;
  std::shared_ptr<const DistributedRoutingPlan> plan_;
  OutOfBoundsPolicy policy_;
  int n_components_ = 0;

  // Own points' coordinates, when they are a subset of the request's.
  Kokkos::View<Real**, DeviceMemorySpace> own_coords_;
  std::unique_ptr<PointEvaluator<T>> own_evaluator_;
  std::unique_ptr<PointEvaluator<T>> incoming_evaluator_;

  // Query points nothing answers for, ascending: never located, or located
  // but not placed by the owner's local evaluator.
  Kokkos::View<LO*, DeviceMemorySpace> unanswered_points_;
};

template <typename T, typename LayoutPolicy>
void DistributedPointEvaluator<T, LayoutPolicy>::Setup(
  const EvaluationRequest& request)
{
  const auto& plan = *plan_;
  const auto n_own = static_cast<LO>(plan.own_points.extent(0));
  const auto n_incoming = static_cast<LO>(plan.IncomingPointCount());

  // Points this evaluator cannot place were already routed to their owner, so
  // a local evaluator that still cannot place one must not throw or guess:
  // FILL marks it, and its `FilledPoints()` reports it below.
  OutOfBoundsPolicy fallback;
  fallback.mode = OutOfBoundsMode::FILL;
  fallback.fill_value = policy_.fill_value;
  const auto& system = space_->GetCoordinateSystem();

  // Local evaluator creation is local, but a failure must reach every rank
  // before the exchange below.
  std::string error;
  try {
    if (n_own > 0) {
      if (plan.AllOwnedHere()) {
        // The local points are the request's points, so the query layout still
        // describes them and the backend may use it.
        EvaluationRequest own_request = request;
        own_request.policy = fallback;
        own_evaluator_ = space_->CreatePointEvaluator<T>(own_request);
      } else {
        own_coords_ = detail::GatherCoordinateRows(request.coords.GetValues(),
                                                   plan.own_points);
        own_evaluator_ =
          space_->CreatePointEvaluator<T>(EvaluationRequest::FromCoordinates(
            CoordinateView<DeviceMemorySpace>(system,
                                              MakeConstRank2View(own_coords_)),
            fallback));
      }
    }
    if (n_incoming > 0) {
      incoming_evaluator_ =
        space_->CreatePointEvaluator<T>(EvaluationRequest::FromCoordinates(
          CoordinateView<DeviceMemorySpace>(
            system, MakeConstRank2View(plan.incoming_coords)),
          fallback));
    }
  } catch (const std::exception& e) {
    error = e.what();
  }
  detail::ThrowIfAnyRankFailed(plan.comm, "DistributedPointEvaluator", error);

  // Each owner tells the requester which of its points its local evaluator
  // filled, once, so `Evaluate` never has to decode that from values.
  const auto incoming_filled =
    incoming_evaluator_
      ? detail::FlagsOf(incoming_evaluator_->FilledPoints(), n_incoming)
      : std::vector<char>(static_cast<size_t>(n_incoming), 0);
  const auto remote_filled =
    ExchangeRecordsWithKnownCounts(plan.comm->get_impl(), plan.incoming_counts,
                                   incoming_filled, plan.remote_counts);
  const auto own_filled =
    own_evaluator_ ? detail::FlagsOf(own_evaluator_->FilledPoints(), n_own)
                   : std::vector<char>{};

  // Unlocated points stay unanswered.
  std::vector<char> answered(static_cast<size_t>(plan.n_points), 0);
  const auto own = detail::ToHostVector(plan.own_points);
  const auto remote = detail::ToHostVector(plan.remote_points);
  for (size_t k = 0; k < own.size(); ++k) {
    answered[static_cast<size_t>(own[k])] = own_filled[k] == 0;
  }
  for (size_t k = 0; k < remote.size(); ++k) {
    answered[static_cast<size_t>(remote[k])] = remote_filled[k] == 0;
  }
  std::vector<LO> unanswered;
  for (LO i = 0; i < plan.n_points; ++i) {
    if (answered[static_cast<size_t>(i)] == 0) {
      unanswered.push_back(i);
    }
  }
  unanswered_points_ = CopyToDevice(unanswered);
}

template <typename T, typename LayoutPolicy>
Kokkos::View<T**, DeviceMemorySpace>
DistributedPointEvaluator<T, LayoutPolicy>::EvaluateLocally(
  const PointEvaluator<T>* evaluator, LO n, const Field<T>& field) const
{
  Kokkos::View<T**, DeviceMemorySpace> values("local values", n, n_components_);
  if (evaluator != nullptr) {
    evaluator->Evaluate(
      field, ValueView<T, DeviceMemorySpace>(field.GetData().GetValueBasis(),
                                             MakeRank2View(values)));
  }
  return values;
}

template <typename T, typename LayoutPolicy>
void DistributedPointEvaluator<T, LayoutPolicy>::Evaluate(
  const Field<T>& field,
  ValueView<T, DeviceMemorySpace, LayoutPolicy> out) const
{
  this->CheckEvaluateWriteTag("DistributedPointEvaluator::Evaluate", field,
                              out);
  const auto& plan = *plan_;
  const auto values = out.GetValues();
  PCMS_ALWAYS_ASSERT(values.extent(0) == static_cast<size_t>(plan.n_points));
  PCMS_ALWAYS_ASSERT(values.extent(1) == static_cast<size_t>(n_components_));
  const int nc = n_components_;

  // Owner side: evaluate this rank's own points and the points sent to it.
  const auto own_values = EvaluateLocally(
    own_evaluator_.get(), static_cast<LO>(plan.own_points.extent(0)), field);
  const auto incoming_values =
    EvaluateLocally(incoming_evaluator_.get(),
                    static_cast<LO>(plan.IncomingPointCount()), field);

  // Route the incoming values back to the ranks that asked: the only values
  // that cross to the host and back, because MPI takes host buffers.
  const auto remote_host = ExchangeRecordsWithKnownCounts(
    plan.comm->get_impl(), detail::Scaled(plan.incoming_counts, nc),
    CopyToHostRowMajor(incoming_values),
    detail::Scaled(plan.remote_counts, nc));
  const auto remote_values = CopyRowMajorToDevice(remote_host, nc);

  // Assemble on the device. A point no rank located keeps the fill value; a
  // point an owner could not place already holds it, from the owner's FILL.
  detail::FillAll(values, static_cast<T>(policy_.fill_value));
  detail::ScatterRows(values, own_values, plan.own_points);
  detail::ScatterRows(values, remote_values, plan.remote_points);

  // Thrown after the exchange, so every rank has completed it.
  if (policy_.mode == OutOfBoundsMode::ERROR &&
      unanswered_points_.extent(0) > 0) {
    throw pcms_error(
      "DistributedPointEvaluator: a query point lies outside the mesh and the "
      "out-of-bounds policy is ERROR");
  }
}

/**
 * \brief Creates distributed point evaluators for a function space on a
 * partitioned mesh.
 *
 * The distributed counterpart of a function space's own evaluator factory, with
 * the same `CreatePointEvaluator(EvaluationRequest)` call. Where a serial
 * factory turns a request into a localization hint, this one turns it into a
 * `DistributedRoutingPlan`, then builds local evaluators for the points each
 * rank owns.
 *
 * Unlike a serial factory, **construction and every `CreatePointEvaluator` /
 * `BuildRoutingPlan` call are collective**: every rank of the mesh's
 * communicator must make the same call, with its own points (possibly none).
 *
 * The search is built once, here, and reused for every request. A plan can be
 * built once and shared by several evaluators: it depends only on the mesh
 * partition and the points, so evaluators for other function spaces on the
 * same mesh can reuse it through `CreatePointEvaluator(plan, request)`.
 */
class DistributedPointEvaluatorFactory : public PointEvaluatorFactory
{
public:
  /**
   * \param mesh The partitioned mesh. Its ownership decides who answers for a
   *        point. Must outlive this factory.
   * \param space The function space evaluators are created for. Its layout
   *        must live on `mesh`.
   * \param tolerance_reach How far outside an element the search may still
   *        accept a point.
   * \param grid_divisions Cells per axis of the global routing grid; 0 derives
   *        it from the mesh. Must be the same on every rank.
   */
  DistributedPointEvaluatorFactory(Omega_h::Mesh& mesh,
                                   std::shared_ptr<const FunctionSpace> space,
                                   Real tolerance_reach = 1e-12,
                                   LO grid_divisions = 0)
    : mesh_(mesh), space_(std::move(space))
  {
    if (!space_) {
      throw std::invalid_argument(
        "DistributedPointEvaluatorFactory: function space is null");
    }
    // The same local search the Omega_h Lagrange factory builds, so for that
    // backend the routing and the owner's local evaluator place every point
    // in the same element.
    local_search_ =
      std::make_unique<GridPointSearchVariant>(detail::MakeSearch(mesh_));
    if (mesh_.dim() == 2) {
      search_ = std::make_unique<DistributedPointSearch<2>>(
        mesh_, std::get<GridPointSearch2D>(*local_search_), tolerance_reach,
        grid_divisions);
    } else {
      search_ = std::make_unique<DistributedPointSearch<3>>(
        mesh_, std::get<GridPointSearch3D>(*local_search_), tolerance_reach,
        grid_divisions);
    }
  }

  [[nodiscard]] const std::shared_ptr<const CoordinateSystem>&
  GetCoordinateSystem() const override
  {
    return space_->GetCoordinateSystem();
  }

  /// Locate `coords` over the decomposition and route them. Collective.
  [[nodiscard]] std::shared_ptr<const DistributedRoutingPlan> BuildRoutingPlan(
    const CoordinateView<DeviceMemorySpace>& coords) const
  {
    detail::CheckSameCoordinateSystem(
      "DistributedPointEvaluatorFactory::BuildRoutingPlan", coords,
      GetCoordinateSystem());
    const auto values = coords.GetValues();
    if (values.extent(1) != static_cast<std::size_t>(mesh_.dim())) {
      throw std::invalid_argument(
        "DistributedPointEvaluatorFactory: query coordinates are not "
        "[num_points x mesh_dim]");
    }
    return std::visit(
      [&](auto& search) {
        constexpr int Dim = std::decay_t<decltype(*search)>::LocalSearch::DIM;
        return std::make_shared<const DistributedRoutingPlan>(
          pcms::BuildRoutingPlan<Dim>(*search,
                                      detail::ToSearchCoords<Dim>(values)));
      },
      search_);
  }

  using PointEvaluatorFactory::CreatePointEvaluator;

  /// An evaluator for `request`'s points, routed by an existing `plan`.
  /// Collective. `plan` must have been built from the same points.
  template <typename T>
  [[nodiscard]] std::unique_ptr<PointEvaluator<T>> CreatePointEvaluator(
    std::shared_ptr<const DistributedRoutingPlan> plan,
    const EvaluationRequest& request) const
  {
    detail::CheckSameCoordinateSystem(
      "DistributedPointEvaluatorFactory::CreatePointEvaluator", request.coords,
      GetCoordinateSystem());
    CheckRequest(request);
    if (!plan) {
      throw std::invalid_argument(
        "DistributedPointEvaluatorFactory: routing plan is null");
    }
    const auto values = request.coords.GetValues();
    if (plan->comm.get() != mesh_.comm().get() ||
        plan->n_points != static_cast<LO>(values.extent(0)) ||
        plan->dim != static_cast<int>(values.extent(1))) {
      throw std::invalid_argument(
        "DistributedPointEvaluatorFactory: the routing plan was not built for "
        "these query points on this mesh");
    }
    return std::make_unique<DistributedPointEvaluator<T>>(
      space_, std::move(plan), request);
  }

  [[nodiscard]] const FunctionSpace& GetFunctionSpace() const noexcept
  {
    return *space_;
  }

  [[nodiscard]] Omega_h::Mesh& GetMesh() const noexcept { return mesh_; }

protected:
  PointEvaluatorVariant CreatePointEvaluatorImpl(
    Type value_type, const EvaluationRequest& request) const override
  {
    CheckRequest(request);
    auto plan = BuildRoutingPlan(request.coords);
    switch (value_type) {
      case Type::Real: return CreatePointEvaluator<Real>(plan, request);
      case Type::LO: return CreatePointEvaluator<LO>(plan, request);
      case Type::GO: return CreatePointEvaluator<GO>(plan, request);
      case Type::Int8: return CreatePointEvaluator<int8_t>(plan, request);
      case Type::Float: return CreatePointEvaluator<float>(plan, request);
    }
    throw pcms_error(
      "DistributedPointEvaluatorFactory: unsupported value type");
  }

private:
  static void CheckRequest(const EvaluationRequest& request)
  {
    if (request.policy.mode == OutOfBoundsMode::NEAREST_BOUNDARY) {
      throw pcms_error("DistributedPointEvaluatorFactory: NEAREST_BOUNDARY is "
                       "not supported");
    }
  }

  Omega_h::Mesh& mesh_;
  std::shared_ptr<const FunctionSpace> space_;
  // Declared before `search_`, which holds a reference to it.
  std::unique_ptr<GridPointSearchVariant> local_search_;
  std::variant<std::unique_ptr<DistributedPointSearch<2>>,
               std::unique_ptr<DistributedPointSearch<3>>>
    search_;
};

/// How many query points of a distributed evaluator get their value from
/// another rank, or 0 when `evaluator` is not one.
template <typename T>
[[nodiscard]] std::size_t LastRoutedPointCount(
  const PointEvaluator<T>& evaluator)
{
  if (auto* d = dynamic_cast<const DistributedPointEvaluator<T>*>(&evaluator)) {
    return d->LastRoutedPointCount();
  }
  return 0;
}

} // namespace pcms

#endif // PCMS_FIELD_EVALUATOR_DISTRIBUTED_POINT_EVALUATOR_H
