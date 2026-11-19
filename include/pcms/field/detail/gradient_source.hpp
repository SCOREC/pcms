#ifndef PCMS_FIELD_DETAIL_GRADIENT_SOURCE_HPP
#define PCMS_FIELD_DETAIL_GRADIENT_SOURCE_HPP

#include "pcms/field/differential_operator.hpp"
#include "pcms/field/point_evaluator_factory.hpp"
#include <memory>

namespace pcms::detail
{

/// @brief Source whose evaluators apply `op` to the fields of `source`.
///
/// It has the geometry of `source`, so it accepts the BoundPoints `source`
/// accepts. Only Value is supported on the returned source.
/// @param op operator to apply; not Value
/// @param source source whose evaluators compute `op`; not null
/// @return source whose evaluators write the basis `op` gives for the field's
///         stored basis
/// @throws pcms_error if `source` is null, `op` is Value, or `source` does
///         not support `op`
[[nodiscard]] std::shared_ptr<const PointEvaluatorFactory> Apply(
  DifferentialOperator op, std::shared_ptr<const PointEvaluatorFactory> source);

} // namespace pcms::detail

#endif // PCMS_FIELD_DETAIL_GRADIENT_SOURCE_HPP
