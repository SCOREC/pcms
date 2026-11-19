#ifndef PCMS_FIELD_GRADIENT_SOURCE_HPP
#define PCMS_FIELD_GRADIENT_SOURCE_HPP

#include "pcms/field/point_evaluator_factory.hpp"
#include <memory>

namespace pcms
{

// Derivative sources: evaluator sources whose evaluators write a derivative of
// the fields of `source`, e.g. as the source of an interpolator or projection
// (CreateInterpolator(Gradient(space), target_space)). They have the geometry
// of `source`, so they accept the BoundPoints `source` accepts.

/// @brief Source whose evaluators write the gradient of the fields of
///        `source`.
/// @throws pcms_error if `source` is null or does not support the gradient
[[nodiscard]] std::shared_ptr<const PointEvaluatorFactory> Gradient(
  std::shared_ptr<const PointEvaluatorFactory> source);

/// @brief Source whose evaluators write the divergence of the fields of
///        `source`.
/// @throws pcms_error if `source` is null or does not support the divergence
[[nodiscard]] std::shared_ptr<const PointEvaluatorFactory> Divergence(
  std::shared_ptr<const PointEvaluatorFactory> source);

/// @brief Source whose evaluators write the curl of the fields of `source`.
/// @throws pcms_error if `source` is null or does not support the curl
[[nodiscard]] std::shared_ptr<const PointEvaluatorFactory> Curl(
  std::shared_ptr<const PointEvaluatorFactory> source);

} // namespace pcms

#endif // PCMS_FIELD_GRADIENT_SOURCE_HPP
