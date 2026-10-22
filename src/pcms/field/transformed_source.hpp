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

/// Sets every component of the listed rows of `values` to `fill_value`.
void FillRows(Rank2View<Real, DeviceMemorySpace> values,
              const Kokkos::View<const LO*, DeviceMemorySpace>& rows,
              Real fill_value);

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
      detail::CheckWrittenValueBasis("TransformedPointEvaluator::Evaluate",
                                     values.GetBasis(),
                                     bound_->OutputBasis(inner_out));
      const ValueView<Real, DeviceMemorySpace> source_values(
        inner_out, values.GetValues());
      inner_->Evaluate(field, source_values);
      bound_->TransformValues(source_values);
      // Fill rows hold no field value to re-express, so the rotation must
      // not change them.
      detail::FillRows(values.GetValues(), inner_->FilledPoints(),
                       policy_.fill_value);
    } else {
      throw pcms_error("TransformedPointEvaluator::Evaluate: basis "
                       "transformations require T == Real");
    }
  }

  [[nodiscard]] Kokkos::View<const LO*, DeviceMemorySpace> FilledPoints()
    const override
  {
    return inner_->FilledPoints();
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
