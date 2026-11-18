#include "pcms/field/gradient_source.hpp"
#include "pcms/utility/assert.hpp"
#include <string>

namespace pcms::detail
{
namespace
{

/// @brief Source whose evaluators apply one operator to the fields of the
///        wrapped source.
class OperatorSource final : public PointEvaluatorFactory
{
public:
  /// @param source source whose evaluators compute `op`; not null
  /// @param op operator to apply
  OperatorSource(std::shared_ptr<const PointEvaluatorFactory> source,
                 DifferentialOperator op)
    : source_(std::move(source)), op_(op)
  {
  }

  [[nodiscard]] std::shared_ptr<const Geometry> GetGeometry()
    const noexcept override
  {
    return source_->GetGeometry();
  }

  [[nodiscard]] bool Supports(DifferentialOperator op) const noexcept override
  {
    return op == DifferentialOperator::Value && source_->Supports(op_);
  }

protected:
  PointEvaluatorVariant CreatePointEvaluatorImpl(
    Type value_type, std::shared_ptr<const BoundSource> bound,
    DifferentialOperator op) const override
  {
    PCMS_ALWAYS_ASSERT(op == DifferentialOperator::Value); // Supports gated it
    return apply_to_type(value_type, [&](auto type) -> PointEvaluatorVariant {
      return CreateUntransformedPointEvaluator<typename decltype(type)::type>(
        *source_, std::move(bound), op_);
    });
  }

private:
  std::shared_ptr<const PointEvaluatorFactory> source_;
  DifferentialOperator op_;
};

} // namespace

std::shared_ptr<const PointEvaluatorFactory> Apply(
  DifferentialOperator op, std::shared_ptr<const PointEvaluatorFactory> source)
{
  if (source == nullptr) {
    throw pcms_error("Apply: source must not be null");
  }
  if (op == DifferentialOperator::Value) {
    throw pcms_error("Apply: Value is not an operator source");
  }
  if (!source->Supports(op)) {
    throw pcms_error("Apply: " + std::string(ToString(op)) +
                     " is not supported by the wrapped source");
  }
  return std::make_shared<OperatorSource>(std::move(source), op);
}

std::shared_ptr<const PointEvaluatorFactory> Gradient(
  std::shared_ptr<const PointEvaluatorFactory> source)
{
  return Apply(DifferentialOperator::Gradient, std::move(source));
}

std::shared_ptr<const PointEvaluatorFactory> Divergence(
  std::shared_ptr<const PointEvaluatorFactory> source)
{
  return Apply(DifferentialOperator::Divergence, std::move(source));
}

std::shared_ptr<const PointEvaluatorFactory> Curl(
  std::shared_ptr<const PointEvaluatorFactory> source)
{
  return Apply(DifferentialOperator::Curl, std::move(source));
}

} // namespace pcms::detail
