#ifndef PCMS_FIELD_DETAIL_POINT_EVALUATOR_FACTORY_HPP
#define PCMS_FIELD_DETAIL_POINT_EVALUATOR_FACTORY_HPP

#include "pcms/field/point_evaluator_factory.hpp"
#include <memory>
#include <set>

namespace pcms::detail
{

/// @brief Source that evaluates fields of one value type from BoundPoints
///        of one kind.
///
/// Implements CreatePointEvaluatorImpl: it checks the value type and passes
/// the BoundPoints on as a `Bound`. Supports reports the operators given to
/// the constructor.
/// @tparam T value type of the evaluated fields
/// @tparam Bound BoundPoints kind the geometry's Bind creates
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
  /// @param geometry geometry whose BoundPoints are `Bound`s
  /// @param supported operators CreateEvaluator implements; Value is always
  ///        included
  explicit TypedPointEvaluatorFactory(
    std::shared_ptr<const Geometry> geometry,
    std::initializer_list<DifferentialOperator> supported = {})
    : geometry_(std::move(geometry)), supported_{DifferentialOperator::Value}
  {
    supported_.insert(supported.begin(), supported.end());
  }

  /// @brief Creates an evaluator from a BoundPoints of the geometry.
  /// @param bound BoundPoints of GetGeometry()
  /// @param op one of the operators given to the constructor
  /// @return the evaluator
  virtual std::unique_ptr<PointEvaluator<T>> CreateEvaluator(
    std::shared_ptr<const Bound> bound, DifferentialOperator op) const = 0;

private:
  PointEvaluatorVariant CreatePointEvaluatorImpl(
    Type value_type, std::shared_ptr<const BoundPoints> bound,
    DifferentialOperator op) const final
  {
    if (value_type != detail::TypeEnumFromType<T>()) {
      throw pcms_error("CreatePointEvaluator: this source only supports "
                       "fields of its own value type");
    }
    auto typed = std::dynamic_pointer_cast<const Bound>(std::move(bound));
    if (typed == nullptr) {
      // a programming error: the geometry check passed
      throw pcms_error("CreatePointEvaluator: BoundPoints of an unexpected "
                       "kind for this geometry");
    }
    return CreateEvaluator(std::move(typed), op);
  }

  std::shared_ptr<const Geometry> geometry_;
  std::set<DifferentialOperator> supported_;
};

} // namespace pcms::detail

#endif // PCMS_FIELD_DETAIL_POINT_EVALUATOR_FACTORY_HPP
