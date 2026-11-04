#ifndef PCMS_FIELD_DETAIL_TRANSFORMED_POINT_EVALUATOR_HPP
#define PCMS_FIELD_DETAIL_TRANSFORMED_POINT_EVALUATOR_HPP

#include "pcms/field/coordinate_map.hpp"
#include "pcms/field/field.hpp"
#include "pcms/field/field_data.hpp"
#include "pcms/field/point_evaluator.hpp"
#include "pcms/field/value_view.hpp"
#include "pcms/geometry/bound_source.hpp"
#include "pcms/utility/assert.hpp"
#include "pcms/utility/memory_spaces.hpp"
#include "pcms/utility/detail/profile.hpp"
#include <memory>
#include <type_traits>
#include <utility>

namespace pcms::field
{

/// Point evaluator at the points of a bound source with a coordinate map: the
/// inner evaluator evaluates at the mapped points, and this one re-expresses
/// the evaluated components in the bound map's output basis.
/// PointEvaluatorFactory::CreatePointEvaluator creates it for every bound
/// source whose GetBoundMap() is set.
///
/// Implementer API: not covered by the public-API stability promise.
template <typename T>
class TransformedPointEvaluator final : public PointEvaluator<T>
{
public:
  /// @param bound bound source of the query points; GetBoundMap() not null
  /// @param inner evaluator created from `bound` that evaluates at the mapped
  ///        points without re-expressing the values; not null
  TransformedPointEvaluator(std::shared_ptr<const BoundSource> bound,
                            std::unique_ptr<PointEvaluator<T>> inner)
    : PointEvaluator<T>(bound, inner->Operator(), inner->AcceptedLayout()),
      map_(bound->GetBoundMap()),
      inner_(std::move(inner))
  {
    PCMS_ALWAYS_ASSERT(map_ != nullptr);
  }

  // Public because CUDA forbids extended lambdas in non-public member
  // functions. User should call Evaluate
  // TODO: add passkey or other mechanism to avoid user calling this
  void EvaluateImpl(const Field<T>& field,
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
      const ValueView<Real, DeviceMemorySpace> source_values(
        inner_out, values.GetValues());
      inner_->Evaluate(field, source_values);
      map_->TransformValues(source_values);
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
    return map_->OutputBasis(inner_out);
  }

private:
  std::shared_ptr<const BoundCoordinateMap> map_;
  std::unique_ptr<PointEvaluator<T>> inner_;
};

} // namespace pcms::field

#endif // PCMS_FIELD_DETAIL_TRANSFORMED_POINT_EVALUATOR_HPP
