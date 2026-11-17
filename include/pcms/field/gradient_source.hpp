#ifndef PCMS_FIELD_GRADIENT_SOURCE_HPP
#define PCMS_FIELD_GRADIENT_SOURCE_HPP

#include "pcms/field/differential_operator.hpp"
#include "pcms/field/point_evaluator_factory.hpp"
#include <memory>

namespace pcms
{

/// @brief Source whose evaluators apply `op` to the fields of `source`.
///
/// It has the geometry of `source`, so it accepts the bound sources `source`
/// accepts. Only Value is supported on the returned source.
/// @param op operator to apply; not Value
/// @param source source whose evaluators compute `op`; not null
/// @return source whose evaluators write the basis `op` gives for the field's
///         stored basis
/// @throws pcms_error if `source` is null, `op` is Value, or `source` does
///         not support `op`
///
/// Implementer API: not covered by the public-API stability promise.
[[nodiscard]] std::shared_ptr<const PointEvaluatorFactory> Apply(
  DifferentialOperator op, std::shared_ptr<const PointEvaluatorFactory> source);

/// @brief Apply(DifferentialOperator::Gradient, source).
///
/// Implementer API: not covered by the public-API stability promise.
[[nodiscard]] std::shared_ptr<const PointEvaluatorFactory> Gradient(
  std::shared_ptr<const PointEvaluatorFactory> source);

/// @brief Apply(DifferentialOperator::Divergence, source).
///
/// Implementer API: not covered by the public-API stability promise.
[[nodiscard]] std::shared_ptr<const PointEvaluatorFactory> Divergence(
  std::shared_ptr<const PointEvaluatorFactory> source);

/// @brief Apply(DifferentialOperator::Curl, source).
///
/// Implementer API: not covered by the public-API stability promise.
[[nodiscard]] std::shared_ptr<const PointEvaluatorFactory> Curl(
  std::shared_ptr<const PointEvaluatorFactory> source);

} // namespace pcms

#endif // PCMS_FIELD_GRADIENT_SOURCE_HPP
