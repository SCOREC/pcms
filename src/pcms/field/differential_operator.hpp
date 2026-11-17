#ifndef PCMS_FIELD_DIFFERENTIAL_OPERATOR_HPP
#define PCMS_FIELD_DIFFERENTIAL_OPERATOR_HPP

#include <string_view>

namespace pcms
{

/// @brief The quantity a point evaluator computes from a field.
enum class DifferentialOperator
{
  Value,      ///< the field's values
  Gradient,   ///< first partial derivatives: one covariant index added
  Divergence, ///< contraction of the gradient: one index removed
  Curl,       ///< rank preserved; odd parity
  Hessian     ///< second partial derivatives: two covariant indices added
};

[[nodiscard]] constexpr std::string_view ToString(
  DifferentialOperator op) noexcept
{
  switch (op) {
    case DifferentialOperator::Value: return "Value";
    case DifferentialOperator::Gradient: return "Gradient";
    case DifferentialOperator::Divergence: return "Divergence";
    case DifferentialOperator::Curl: return "Curl";
    case DifferentialOperator::Hessian: return "Hessian";
  }
  return "?";
}

} // namespace pcms

#endif // PCMS_FIELD_DIFFERENTIAL_OPERATOR_HPP
