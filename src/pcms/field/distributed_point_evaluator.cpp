#include "pcms/field/distributed_evaluator.hpp"
#include "pcms/field/field.hpp"
#include "pcms/field/point_evaluator.hpp"
#include "pcms/field/value_view.hpp"
#include "pcms/geometry/distributed_binding.hpp"
#include "pcms/geometry/distributed_routing_plan.hpp"
#include "pcms/utility/detail/arrays.hpp"
#include "pcms/utility/detail/types.hpp"
#include "pcms/utility/memory_spaces.hpp"
#include "pcms/utility/mpi_utilities.hpp"
#include "pcms/utility/types.hpp"

#include <Kokkos_Core.hpp>

#include <memory>

namespace pcms
{

namespace detail
{

// Write row `k` of `rows` to row `targets(k)` of `values`.
template <typename T, typename Values>
void ScatterRows(Values values, Kokkos::View<T**, DeviceMemorySpace> rows,
                 Kokkos::View<const LO*, DeviceMemorySpace> targets)
{
  const int nc = static_cast<int>(values.extent(1));
  Kokkos::parallel_for(
    "DistributedPointEvaluator::Scatter",
    Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0,
                                                            targets.extent(0)),
    KOKKOS_LAMBDA(LO k) {
      const LO i = targets(k);
      for (int c = 0; c < nc; ++c) {
        values(i, c) = rows(k, c);
      }
    });
}

} // namespace detail

namespace
{

/**
 * \brief Evaluates a field at points distributed over ranks.
 *
 * Each rank evaluates the field at its own points and at the points routed to
 * it, with evaluators of the local evaluator source, and sends the values of
 * the routed points back to the ranks that asked. Evaluate is collective. The
 * field's ghost values must be synchronized.
 */
template <typename T>
class DistributedPointEvaluator final : public PointEvaluator<T>
{
public:
  /**
   * \param bound bound points of a DistributedGeometry
   * \param local_source evaluator source on the local geometry
   * \param op the operator the evaluator applies
   */
  DistributedPointEvaluator(std::shared_ptr<const DistributedBinding> bound,
                            const PointEvaluatorFactory& local_source,
                            DifferentialOperator op)
    : PointEvaluator<T>(bound, op),
      own_(local_source.CreatePointEvaluator<T>(bound->OwnBound(), op)),
      incoming_(
        local_source.CreatePointEvaluator<T>(bound->IncomingBound(), op))
  {
  }

  [[nodiscard]] ValueBasis OutputBasis(const ValueBasis& stored) const override
  {
    return own_->OutputBasis(stored);
  }

protected:
  void EvaluateImpl(const Field<T>& field,
                    ValueView<T, DeviceMemorySpace> out) const override
  {
    const auto& plan =
      static_cast<const DistributedBinding&>(this->Bound()).GetRoutingPlan();
    const auto values = out.GetValues();
    const int nc = static_cast<int>(values.extent(1));

    const auto own_values = EvaluateLocally(*own_, field, out.GetBasis(), nc);
    const auto incoming_values =
      EvaluateLocally(*incoming_, field, out.GetBasis(), nc);

    // Send the incoming points' values back to the ranks that asked.
    const auto remote_host = ExchangeRecordsWithKnownCounts(
      plan.comm, detail::Scaled(plan.incoming_counts, nc),
      pcms::detail::CopyToHostRowMajor(incoming_values),
      detail::Scaled(plan.remote_counts, nc));
    const auto remote_values =
      pcms::detail::CopyRowMajorToDevice(remote_host, nc);

    detail::ScatterRows(values, own_values, plan.own_points);
    detail::ScatterRows(values, remote_values, plan.remote_points);
  }

private:
  static Kokkos::View<T**, DeviceMemorySpace> EvaluateLocally(
    const PointEvaluator<T>& evaluator, const Field<T>& field,
    const ValueBasis& basis, int nc)
  {
    Kokkos::View<T**, DeviceMemorySpace> values("local values",
                                                evaluator.NumPoints(), nc);
    evaluator.Evaluate(
      field, ValueView<T, DeviceMemorySpace>(basis, MakeRank2View(values)));
    return values;
  }

  std::unique_ptr<PointEvaluator<T>> own_;
  std::unique_ptr<PointEvaluator<T>> incoming_;
};

/**
 * \brief Evaluator source over a DistributedGeometry: wraps an evaluator source
 * on its local geometry (e.g. a space), which evaluates the points each rank
 * owns.
 */
class DistributedPointEvaluatorFactory final : public PointEvaluatorFactory
{
public:
  /**
   * \param geometry the distributed geometry
   * \param local_source evaluator source on geometry's local geometry
   */
  DistributedPointEvaluatorFactory(
    std::shared_ptr<const DistributedGeometry> geometry,
    std::shared_ptr<const PointEvaluatorFactory> local_source)
    : geometry_(std::move(geometry)), local_source_(std::move(local_source))
  {
    if (geometry_ == nullptr || local_source_ == nullptr) {
      throw pcms_error("DistributedPointEvaluatorFactory: the geometry and "
                       "the local source must not be null");
    }
    if (local_source_->GetGeometry() != geometry_->GetLocalGeometry()) {
      throw pcms_error("DistributedPointEvaluatorFactory: the local source "
                       "is not on the geometry's local geometry");
    }
  }

  [[nodiscard]] std::shared_ptr<const Geometry> GetGeometry()
    const noexcept override
  {
    return geometry_;
  }

  [[nodiscard]] bool Supports(DifferentialOperator op) const noexcept override
  {
    return local_source_->Supports(op);
  }

protected:
  PointEvaluatorVariant CreatePointEvaluatorImpl(
    Type value_type, std::shared_ptr<const BoundPoints> bound,
    DifferentialOperator op) const override
  {
    auto binding =
      std::dynamic_pointer_cast<const DistributedBinding>(std::move(bound));
    if (binding == nullptr) {
      throw pcms_error("DistributedPointEvaluatorFactory: bound points of an "
                       "unexpected kind for this geometry");
    }
    return pcms::detail::apply_to_type(
      value_type, [&](auto tag) -> PointEvaluatorVariant {
        using T = typename decltype(tag)::type;
        return std::unique_ptr<PointEvaluator<T>>(
          std::make_unique<DistributedPointEvaluator<T>>(binding,
                                                         *local_source_, op));
      });
  }

private:
  std::shared_ptr<const DistributedGeometry> geometry_;
  std::shared_ptr<const PointEvaluatorFactory> local_source_;
};

} // namespace

std::shared_ptr<const PointEvaluatorFactory> CreateDistributedEvaluator(
  std::shared_ptr<const DistributedGeometry> geometry,
  std::shared_ptr<const PointEvaluatorFactory> local_source)
{
  return std::make_shared<DistributedPointEvaluatorFactory>(
    std::move(geometry), std::move(local_source));
}

} // namespace pcms
