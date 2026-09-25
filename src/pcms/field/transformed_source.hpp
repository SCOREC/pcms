#ifndef PCMS_FIELD_TRANSFORMED_SOURCE_HPP
#define PCMS_FIELD_TRANSFORMED_SOURCE_HPP

#include "pcms/field/coordinate_map.hpp"
#include "pcms/field/coordinate_system.hpp"
#include "pcms/field/evaluation_request.h"
#include "pcms/field/field.h"
#include "pcms/field/field_data.h"
#include "pcms/field/out_of_bounds_policy.h"
#include "pcms/field/point_evaluator.h"
#include "pcms/field/point_evaluator_factory.hpp"
#include "pcms/field/point_status.hpp"
#include "pcms/field/value_view.hpp"
#include "pcms/utility/arrays.h"
#include "pcms/utility/assert.h"
#include "pcms/utility/memory_spaces.h"
#include "pcms/utility/profile.h"
#include "pcms/utility/types.h"
#include <Kokkos_Core.hpp>
#include <memory>
#include <type_traits>
#include <utility>

namespace pcms
{

namespace detail
{

/// Throws pcms_error naming `context` if any status is not Valid.
void ThrowIfAnyPointInvalid(const PointStatusView& status, const char* context);

/// Sets every component of `values` row i to `fill_value` where all
/// components of `pre_transform` row i equal `fill_value`.
void RestoreFillRows(
  const Kokkos::View<Real**, DeviceMemorySpace>& pre_transform,
  Rank2View<Real, DeviceMemorySpace> values, Real fill_value);

} // namespace detail

/// Point evaluator over mapped query points that re-expresses evaluated
/// components in the bound map's output basis.
template <typename T>
class TransformedPointEvaluator final : public PointEvaluator<T>
{
public:
  /// @param bound map bound to the query points; owns the mapped coordinates
  ///        `inner` was built over, so it is declared first and outlives it
  /// @param inner evaluator bound to `bound`'s mapped points
  /// @param policy out-of-bounds policy `inner` was built with
  TransformedPointEvaluator(std::unique_ptr<BoundCoordinateMap> bound,
                            std::unique_ptr<PointEvaluator<T>> inner,
                            OutOfBoundsPolicy policy)
    : bound_(std::move(bound)), inner_(std::move(inner)), policy_(policy)
  {
  }

  void Evaluate(const Field<T>& field,
                ValueView<T, DeviceMemorySpace> values) const override
  {
    PCMS_FUNCTION_TIMER;
    const ValueBasis inner_out =
      inner_->OutputBasis(field.GetData().GetValueBasis());
    if (inner_out.Rank() == 0) {
      inner_->Evaluate(field, values);
      return;
    }
    if constexpr (std::is_same_v<T, Real>) {
      if (scratch_.extent(0) != values.extent(0) ||
          scratch_.extent(1) != values.extent(1)) {
        scratch_ = Kokkos::View<Real**, DeviceMemorySpace>(
          "transformed_point_evaluator_scratch", values.extent(0),
          values.extent(1));
      }
      inner_->Evaluate(field, ValueView<Real, DeviceMemorySpace>(
                                inner_out, MakeRank2View(scratch_)));
      bound_->TransformValues(ValueView<const Real, DeviceMemorySpace>(
                                inner_out, MakeConstRank2View(scratch_)),
                              values);
      if (policy_.mode == OutOfBoundsMode::FILL) {
        detail::RestoreFillRows(scratch_, values.GetValues(),
                                policy_.fill_value);
      }
    } else {
      throw pcms_error("TransformedPointEvaluator::Evaluate: basis "
                       "transformations require T == Real");
    }
  }

  [[nodiscard]] ValueBasis OutputBasis(const ValueBasis& stored) const override
  {
    const ValueBasis inner_out = inner_->OutputBasis(stored);
    if (inner_out.Rank() == 0) {
      return inner_out;
    }
    return bound_->OutputBasis(inner_out);
  }

private:
  std::unique_ptr<BoundCoordinateMap> bound_;
  std::unique_ptr<PointEvaluator<T>> inner_;
  OutOfBoundsPolicy policy_;
  mutable Kokkos::View<Real**, DeviceMemorySpace> scratch_;
};

/// Presents a point-evaluator factory as one in another coordinate system:
/// query points are mapped into the source's system and evaluated components
/// are re-expressed at the query points.
class TransformedSource final : public PointEvaluatorFactory
{
public:
  /// @param source factory whose evaluators do the evaluation; co-owned, so
  ///        nesting requires the outer view to hold a shared_ptr to the inner
  /// @param to_source map from this object's coordinate system to `source`'s
  TransformedSource(std::shared_ptr<const PointEvaluatorFactory> source,
                    std::shared_ptr<const CoordinateMap> to_source);

  [[nodiscard]] const std::shared_ptr<const CoordinateSystem>&
  GetCoordinateSystem() const override
  {
    return system_;
  }

protected:
  PointEvaluatorVariant CreatePointEvaluatorImpl(
    Type value_type, const EvaluationRequest& request) const override;

private:
  template <typename T>
  std::unique_ptr<PointEvaluator<T>> Create(
    const EvaluationRequest& request) const;

  std::shared_ptr<const PointEvaluatorFactory> source_;
  std::shared_ptr<const CoordinateMap> to_source_;
  std::shared_ptr<const CoordinateSystem> system_;
};

} // namespace pcms

#endif // PCMS_FIELD_TRANSFORMED_SOURCE_HPP
