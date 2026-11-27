// Part of pcms/field/point_evaluator.hpp, which includes it after declaring
// PointEvaluator and before defining it.
#include "pcms/field/point_evaluator.hpp"
#ifndef PCMS_FIELD_DETAIL_EVALUATE_TRANSPOSE_HPP
#define PCMS_FIELD_DETAIL_EVALUATE_TRANSPOSE_HPP

namespace pcms::detail
{

/// @brief functional = E^T values + beta functional, with the checks of
///        PointEvaluator::EvaluateTranspose; beta == 0 does not read the
///        functional (collective)
template <typename T, typename LayoutPolicy>
void EvaluateTranspose(
  const PointEvaluator<T, LayoutPolicy>& evaluator,
  ValueView<const T, DeviceMemorySpace, LayoutPolicy> values,
  Field<T>& functional, T beta);

} // namespace pcms::detail

#endif // PCMS_FIELD_DETAIL_EVALUATE_TRANSPOSE_HPP
