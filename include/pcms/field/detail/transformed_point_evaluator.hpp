#ifndef PCMS_FIELD_DETAIL_TRANSFORMED_POINT_EVALUATOR_HPP
#define PCMS_FIELD_DETAIL_TRANSFORMED_POINT_EVALUATOR_HPP

#include "pcms/field/field.hpp"
#include "pcms/field/field_data.hpp"
#include "pcms/field/detail/transform_chain.hpp"
#include "pcms/field/point_evaluator.hpp"
#include "pcms/field/value_transform.hpp"
#include "pcms/field/value_view.hpp"
#include "pcms/geometry/bound_points.hpp"
#include "pcms/utility/assert.hpp"
#include "pcms/utility/memory_spaces.hpp"
#include "pcms/utility/detail/profile.hpp"
#include <functional>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>

namespace pcms::detail
{

/// The value transforms, in order, that re-express values in a basis.
using ValueTransforms =
  std::function<std::vector<const ValueTransform*>(const ValueBasis&)>;

/// Point evaluator that re-expresses the values of an inner evaluator through
/// the value transforms of their basis (those of a bind with a coordinate
/// map), through a TransformChain.
/// PointEvaluatorFactory::CreatePointEvaluator creates it.
template <typename T>
class TransformedPointEvaluator final : public PointEvaluator<T>
{
public:
  /// @param bound BoundPoints of the query points
  /// @param inner evaluator created from `bound` that writes the values
  ///        before the transforms; not null
  /// @param transforms the transforms of the values of `inner`
  TransformedPointEvaluator(std::shared_ptr<const BoundPoints> bound,
                            std::unique_ptr<PointEvaluator<T>> inner,
                            ValueTransforms transforms)
    : PointEvaluator<T>(std::move(bound), inner->Operator(),
                        inner->AcceptedLayout()),
      inner_(std::move(inner)),
      transforms_(std::move(transforms))
  {
  }

  // Public because CUDA forbids extended lambdas in non-public member
  // functions. User should call Evaluate
  // TODO: add passkey or other mechanism to avoid user calling this
  void EvaluateImpl(const Field<T>& field,
                    ValueView<T, DeviceMemorySpace> values) const override
  {
    PCMS_FUNCTION_TIMER;
    const auto& chain = Chain(field.GetData().GetValueBasis());
    if (chain.Empty()) {
      inner_->Evaluate(field, values);
      return;
    }
    if constexpr (std::is_same_v<T, Real>) {
      const auto in = chain.InputStorage(values.GetValues());
      inner_->Evaluate(
        field, ValueView<Real, DeviceMemorySpace>(chain.InputBasis(), in));
      chain.Apply(in, values.GetValues());
    } else {
      throw pcms_error("TransformedPointEvaluator::Evaluate: basis "
                       "transformations require T == Real");
    }
  }

  /// @return whether the inner evaluator supports the transpose
  [[nodiscard]] bool SupportsTranspose() const noexcept override
  {
    return inner_->SupportsTranspose();
  }

  // Public for the same reason as EvaluateImpl.
  void EvaluateTransposeImpl(Rank2View<const T, DeviceMemorySpace> values,
                             Field<T>& functional) const override
  {
    PCMS_FUNCTION_TIMER;
    const auto& chain = Chain(functional.GetData().GetValueBasis());
    if constexpr (std::is_same_v<T, Real>) {
      inner_->EvaluateTranspose(
        ValueView<const Real, DeviceMemorySpace>(chain.InputBasis(),
                                                 chain.ApplyTranspose(values)),
        functional);
    } else {
      if (!chain.Empty()) {
        throw pcms_error("TransformedPointEvaluator::EvaluateTranspose: basis "
                         "transformations require T == Real");
      }
      inner_->EvaluateTranspose(
        ValueView<const T, DeviceMemorySpace>(chain.InputBasis(), values),
        functional);
    }
  }

  [[nodiscard]] ValueBasis OutputBasis(const ValueBasis& stored) const override
  {
    const auto input = inner_->OutputBasis(stored);
    return TransformChain(input, transforms_(input)).OutputBasis();
  }

private:
  // the chain of the inner evaluator's basis for `stored`
  const TransformChain& Chain(const ValueBasis& stored) const
  {
    const auto input = inner_->OutputBasis(stored);
    chain_.Reset(input, transforms_(input));
    return chain_;
  }

  std::unique_ptr<PointEvaluator<T>> inner_;
  ValueTransforms transforms_;
  mutable TransformChain chain_;
};

} // namespace pcms::detail

#endif // PCMS_FIELD_DETAIL_TRANSFORMED_POINT_EVALUATOR_HPP
