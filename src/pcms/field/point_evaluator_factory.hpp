#ifndef PCMS_FIELD_POINT_EVALUATOR_FACTORY_HPP
#define PCMS_FIELD_POINT_EVALUATOR_FACTORY_HPP

#include "pcms/field/differential_operator.hpp"
#include "pcms/field/field_factory.h"
#include "pcms/field/point_evaluator.h"
#include "pcms/field/transformed_point_evaluator.hpp"
#include "pcms/geometry/bound_source.hpp"
#include "pcms/geometry/geometry.hpp"
#include "pcms/utility/assert.h"
#include "pcms/utility/types.h"
#include <initializer_list>
#include <memory>
#include <set>
#include <string>
#include <variant>

namespace pcms
{

/// @brief Source of point evaluators: an evaluator source on a geometry.
///
/// The geometry's Bind computes the per-point data for a set of query points
/// once; CreatePointEvaluator creates evaluators of any value type and
/// supported operator from the result.
class PointEvaluatorFactory
{
public:
  /// @return the geometry whose bound sources this source accepts
  [[nodiscard]] virtual std::shared_ptr<const Geometry> GetGeometry()
    const noexcept = 0;

  /// @return whether CreatePointEvaluator accepts `op`
  [[nodiscard]] virtual bool Supports(
    DifferentialOperator op) const noexcept = 0;

  /// @brief Creates an evaluator at the points of a bound source.
  ///
  /// If the bind had a coordinate map, the evaluator evaluates at the mapped
  /// points and re-expresses the values through the bound map
  /// (TransformedPointEvaluator).
  /// @tparam T field value type
  /// @param bound a bound source of this source's geometry (a result of
  ///        GetGeometry()->Bind); bound sources of other geometries are
  ///        refused, even on the same mesh
  /// @param op the operator the evaluator applies
  /// @return evaluator that writes OutputBasis of the field's stored basis
  /// @throws pcms_error if the source has no geometry (it does not evaluate
  ///         at points), `bound` is null or belongs to another geometry, or
  ///         the source does not support `op`
  template <typename T>
  [[nodiscard]] std::unique_ptr<PointEvaluator<T>> CreatePointEvaluator(
    std::shared_ptr<const BoundSource> bound,
    DifferentialOperator op = DifferentialOperator::Value) const;

  virtual ~PointEvaluatorFactory() noexcept = default;

protected:
  /// @brief CreatePointEvaluator of `source` without the bound map's value
  ///        transformation: the evaluator evaluates at the located points and
  ///        writes the source's basis.
  ///
  /// A source that forwards to another source (a space to its backend source,
  /// an operator source to the source it wraps) creates the inner evaluator
  /// with this, so that its own CreatePointEvaluator transforms the values
  /// once.
  /// @throws pcms_error as CreatePointEvaluator
  template <typename T>
  [[nodiscard]] static std::unique_ptr<PointEvaluator<T>>
  CreateUntransformedPointEvaluator(const PointEvaluatorFactory& source,
                                    std::shared_ptr<const BoundSource> bound,
                                    DifferentialOperator op);

  /// @brief Creates an evaluator from a bound source.
  /// @param value_type field value type
  /// @param bound bound source of GetGeometry()
  /// @param op the operator the evaluator applies; Supports(op) is true
  /// @return evaluator for `value_type`
  /// @throws pcms_error for an unsupported value type
  virtual PointEvaluatorVariant CreatePointEvaluatorImpl(
    Type value_type, std::shared_ptr<const BoundSource> bound,
    DifferentialOperator op) const = 0;
};

/// @brief Source that evaluates fields of one value type from bound sources
///        of one kind.
///
/// Implements CreatePointEvaluatorImpl: it checks the value type and passes
/// the bound source on as a `Bound`. Supports reports the operators given to
/// the constructor.
/// @tparam T value type of the evaluated fields
/// @tparam Bound bound-source kind the geometry's Bind creates
template <typename T, typename Bound>
class TypedPointEvaluatorFactory : public PointEvaluatorFactory
{
public:
  [[nodiscard]] std::shared_ptr<const Geometry> GetGeometry()
    const noexcept final
  {
    return geometry_;
  }

  [[nodiscard]] bool Supports(DifferentialOperator op) const noexcept final
  {
    return supported_.count(op) != 0;
  }

protected:
  /// @param geometry geometry whose bound sources are `Bound`s
  /// @param supported operators CreateEvaluator implements; Value is always
  ///        included
  explicit TypedPointEvaluatorFactory(
    std::shared_ptr<const Geometry> geometry,
    std::initializer_list<DifferentialOperator> supported = {})
    : geometry_(std::move(geometry)), supported_{DifferentialOperator::Value}
  {
    supported_.insert(supported.begin(), supported.end());
  }

  /// @brief Creates an evaluator from a bound source of the geometry.
  /// @param bound bound source of GetGeometry()
  /// @param op one of the operators given to the constructor
  /// @return the evaluator
  virtual std::unique_ptr<PointEvaluator<T>> CreateEvaluator(
    std::shared_ptr<const Bound> bound, DifferentialOperator op) const = 0;

private:
  PointEvaluatorVariant CreatePointEvaluatorImpl(
    Type value_type, std::shared_ptr<const BoundSource> bound,
    DifferentialOperator op) const final
  {
    if (value_type != TypeEnumFromType<T>()) {
      throw pcms_error("CreatePointEvaluator: this source only supports "
                       "fields of its own value type");
    }
    auto typed = std::dynamic_pointer_cast<const Bound>(std::move(bound));
    if (typed == nullptr) {
      // a programming error: the geometry check passed
      throw pcms_error("CreatePointEvaluator: bound source of an unexpected "
                       "kind for this geometry");
    }
    return CreateEvaluator(std::move(typed), op);
  }

  std::shared_ptr<const Geometry> geometry_;
  std::set<DifferentialOperator> supported_;
};

template <typename T>
std::unique_ptr<PointEvaluator<T>>
PointEvaluatorFactory::CreateUntransformedPointEvaluator(
  const PointEvaluatorFactory& source, std::shared_ptr<const BoundSource> bound,
  DifferentialOperator op)
{
  static_assert(is_supported_field_type_v<T>,
                "T is not a supported field type");
  const auto geometry = source.GetGeometry();
  if (geometry == nullptr) {
    throw pcms_error(
      "CreatePointEvaluator: this source does not evaluate fields at points");
  }
  if (bound == nullptr) {
    throw pcms_error("CreatePointEvaluator: bound source must not be null");
  }
  if (bound->GetGeometry().get() != geometry.get()) {
    throw pcms_error("CreatePointEvaluator: the bound source belongs to "
                     "another geometry; bind this source's geometry to the "
                     "points");
  }
  if (!source.Supports(op)) {
    throw pcms_error("CreatePointEvaluator: " + std::string(ToString(op)) +
                     " is not supported by this source");
  }
  return std::get<std::unique_ptr<PointEvaluator<T>>>(
    source.CreatePointEvaluatorImpl(TypeEnumFromType<T>(), std::move(bound),
                                    op));
}

template <typename T>
std::unique_ptr<PointEvaluator<T>> PointEvaluatorFactory::CreatePointEvaluator(
  std::shared_ptr<const BoundSource> bound, DifferentialOperator op) const
{
  auto evaluator = CreateUntransformedPointEvaluator<T>(*this, bound, op);
  if (bound->GetBoundMap() == nullptr) {
    return evaluator;
  }
  return std::make_unique<TransformedPointEvaluator<T>>(std::move(bound),
                                                        std::move(evaluator));
}

} // namespace pcms

#endif // PCMS_FIELD_POINT_EVALUATOR_FACTORY_HPP
