#ifndef PCMS_FIELD_FUNCTION_SPACE_QUADRATURE_SPACE_HPP
#define PCMS_FIELD_FUNCTION_SPACE_QUADRATURE_SPACE_HPP

#include "pcms/field/detail/simple_field_data.hpp"
#include "pcms/field/function_space.hpp"
#include <functional>
#include <memory>

namespace pcms::detail
{

/// A space with a natural quadrature a function makes.
class QuadratureSpace final : public FunctionSpace
{
public:
  QuadratureSpace(std::shared_ptr<const FieldLayout> layout,
                  std::shared_ptr<const PointEvaluatorFactory> evaluator,
                  std::shared_ptr<const Geometry> geometry,
                  std::function<Quadrature()> make_quadrature)
    : FunctionSpace(layout, std::move(evaluator), std::move(geometry)),
      layout_(std::move(layout)),
      make_quadrature_(std::move(make_quadrature))
  {
  }

  [[nodiscard]] Quadrature CreateQuadrature() const override
  {
    return make_quadrature_();
  }

protected:
  FieldDataVariant CreateFieldData(Type type, ValueBasis basis) const override
  {
    if (type != Type::Real) {
      throw pcms_error("this space has Real fields only");
    }
    return std::make_unique<SimpleFieldData<Real>>(layout_, std::move(basis));
  }

private:
  std::shared_ptr<const FieldLayout> layout_;
  std::function<Quadrature()> make_quadrature_;
};

} // namespace pcms::detail

#endif // PCMS_FIELD_FUNCTION_SPACE_QUADRATURE_SPACE_HPP
