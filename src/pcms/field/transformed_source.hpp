#ifndef PCMS_FIELD_TRANSFORMED_SOURCE_HPP
#define PCMS_FIELD_TRANSFORMED_SOURCE_HPP

#include "pcms/field/basis_transformation.hpp"
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
/// components in a target value basis.
template <typename T>
class TransformedPointEvaluator final : public PointEvaluator<T>
{
public:
  /// @param inner evaluator bound to `mapped`'s points
  /// @param mapped query points mapped into the source's coordinate system
  /// @param law transformation bound at the query points; null for rank 0
  /// @param source_basis basis the source field must be stored in
  /// @param target_basis basis written when `law` is non-null
  /// @param policy out-of-bounds policy `inner` was built with
  TransformedPointEvaluator(std::unique_ptr<PointEvaluator<T>> inner,
                            MappedPoints mapped,
                            std::unique_ptr<BoundBasisTransformation> law,
                            ValueBasis source_basis, ValueBasis target_basis,
                            OutOfBoundsPolicy policy)
    : inner_(std::move(inner)),
      mapped_(std::move(mapped)),
      law_(std::move(law)),
      source_basis_(std::move(source_basis)),
      target_basis_(std::move(target_basis)),
      policy_(policy)
  {
  }

  void Evaluate(const Field<T>& field,
                ValueView<T, DeviceMemorySpace> values) const override
  {
    PCMS_FUNCTION_TIMER;
    if (law_ == nullptr) {
      inner_->Evaluate(field, values);
      return;
    }
    if constexpr (std::is_same_v<T, Real>) {
      if (!SameValueBasis(inner_->OutputBasis(field.GetData().GetValueBasis()),
                          source_basis_)) {
        throw pcms_error("TransformedPointEvaluator::Evaluate: the field's "
                         "stored basis is not the source basis this "
                         "evaluator was built for");
      }
      detail::CheckWrittenValueBasis("TransformedPointEvaluator::Evaluate",
                                     values.GetBasis(), target_basis_);
      if (scratch_.extent(0) != values.extent(0) ||
          scratch_.extent(1) != values.extent(1)) {
        scratch_ = Kokkos::View<Real**, DeviceMemorySpace>(
          "transformed_point_evaluator_scratch", values.extent(0),
          values.extent(1));
      }
      inner_->Evaluate(field, ValueView<Real, DeviceMemorySpace>(
                                source_basis_, MakeRank2View(scratch_)));
      law_->Apply(ValueView<const Real, DeviceMemorySpace>(
                    source_basis_, MakeConstRank2View(scratch_)),
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
    if (law_ == nullptr) {
      return stored;
    }
    if (!SameValueBasis(inner_->OutputBasis(stored), source_basis_)) {
      throw pcms_error("TransformedPointEvaluator::OutputBasis: the stored "
                       "basis is not the source basis this evaluator was "
                       "built for");
    }
    return target_basis_;
  }

private:
  std::unique_ptr<PointEvaluator<T>> inner_;
  MappedPoints mapped_;
  std::unique_ptr<BoundBasisTransformation> law_;
  ValueBasis source_basis_;
  ValueBasis target_basis_;
  OutOfBoundsPolicy policy_;
  mutable Kokkos::View<Real**, DeviceMemorySpace> scratch_;
};

/// Presents a point-evaluator factory as one in another coordinate system:
/// query points are mapped into the source's system and evaluated components
/// are re-expressed in the target basis at the query points.
class TransformedSource final : public PointEvaluatorFactory
{
public:
  /// @param source factory whose evaluators do the evaluation; must outlive
  ///        this object
  /// @param to_source map from this object's coordinate system to `source`'s
  /// @param source_basis basis fields evaluated through this object are
  ///        stored in; its system must be the map's target system for rank > 0
  /// @param target_basis basis evaluators created here write; its system must
  ///        be the map's source system for rank > 0, with `source_basis`'s
  ///        rank and variance
  TransformedSource(const PointEvaluatorFactory& source,
                    std::shared_ptr<const CoordinateMap> to_source,
                    ValueBasis source_basis, ValueBasis target_basis);

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

  const PointEvaluatorFactory* source_;
  std::shared_ptr<const CoordinateMap> to_source_;
  std::shared_ptr<const CoordinateSystem> system_;
  ValueBasis source_basis_;
  ValueBasis target_basis_;
};

} // namespace pcms

#endif // PCMS_FIELD_TRANSFORMED_SOURCE_HPP
