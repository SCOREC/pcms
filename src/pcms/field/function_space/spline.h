#ifndef PCMS_SPLINE_FUNCTION_SPACE_H
#define PCMS_SPLINE_FUNCTION_SPACE_H

#include "pcms/field/field.h"
#include "pcms/field/field_data.h"
#include "pcms/field/field_evaluator_factory.h"
#include "pcms/field/field_layout.h"
#include "pcms/field/field_metadata.h"
#include "pcms/field/function_space.h"
#include "pcms/field/out_of_bounds_policy.h"
#include "pcms/field/point_evaluator.h"
#include "pcms/utility/arrays.h"
#include "pcms/utility/eqdsk.h"
#include "pcms/utility/memory_spaces.h"
#include "pcms/utility/uniform_grid.h"

#include <memory>

namespace pcms
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
  [[nodiscard]] FieldVariant CreateFieldImpl(Type value_type, FieldMetadata metadata) const override;

  [[nodiscard]] FieldVariant CreateFieldImpl(
    FieldDataVariant data) const override;

  [[nodiscard]] PointEvaluatorVariant CreatePointEvaluatorImpl(
    Type value_type, const EvaluationRequest& request) const override;

private:
  std::shared_ptr<const FieldLayout> layout_;
  std::shared_ptr<FieldEvaluatorFactory<Real>> evaluator_factory_;
};

} // namespace pcms

#endif // PCMS_SPLINE_FUNCTION_SPACE_H
