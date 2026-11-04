#ifndef PCMS_FIELD_FUNCTION_SPACE_SPLINE_HPP
#define PCMS_FIELD_FUNCTION_SPACE_SPLINE_HPP

#include "pcms/field/field.hpp"
#include "pcms/field/field_data.hpp"
#include "pcms/field/field_evaluator_factory.hpp"
#include "pcms/field/field_layout.hpp"
#include "pcms/field/value_view.hpp"
#include "pcms/field/function_space.hpp"
#include "pcms/field/out_of_bounds_policy.hpp"
#include "pcms/field/point_evaluator.hpp"
#include "pcms/utility/arrays.hpp"
#include "pcms/utility/eqdsk.hpp"
#include "pcms/utility/memory_spaces.hpp"
#include "pcms/utility/uniform_grid.hpp"

#include <memory>

namespace pcms::field
{

class SplineFunctionSpace : public FunctionSpace
{
  // Passkey: only members (the From* factories) can create a Key, so only they
  // can construct a space, while std::make_shared can still call the ctor.
  struct Key
  {
    explicit Key() = default;
  };

public:
  SplineFunctionSpace(
    Key, std::shared_ptr<const FieldLayout> layout,
    std::shared_ptr<FieldEvaluatorFactory<Real>> evaluator_factory) noexcept;

  [[nodiscard]] static std::shared_ptr<SplineFunctionSpace> FromUniformGrid(
    const UniformGrid<2>& grid,
    std::shared_ptr<const CoordinateSystem> mesh_coordinate_system);

  [[nodiscard]] std::shared_ptr<const FieldLayout> GetLayout()
    const noexcept override;

protected:
  [[nodiscard]] FieldVariant CreateFieldImpl(Type storage_type,
                                             ValueBasis basis) const override;

  [[nodiscard]] FieldVariant CreateFieldImpl(
    FieldDataVariant data) const override;

  [[nodiscard]] PointEvaluatorVariant CreatePointEvaluatorImpl(
    Type value_type, const EvaluationRequest& request) const override;

private:
  std::shared_ptr<const FieldLayout> layout_;
  std::shared_ptr<FieldEvaluatorFactory<Real>> evaluator_factory_;
};

} // namespace pcms::field

#endif // PCMS_FIELD_FUNCTION_SPACE_SPLINE_HPP
