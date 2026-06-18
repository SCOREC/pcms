#include "pcms/field/function_space/spline.h"

#include "pcms/field/data/simple.h"
#include "pcms/field/evaluator/uniform_grid_spline.h"
#include "pcms/field/layout/uniform_grid.h"
#include "pcms/utility/common.h"

#include <stdexcept>

namespace pcms
{

SplineFunctionSpace::SplineFunctionSpace(
  std::shared_ptr<const FieldLayout> layout,
  std::shared_ptr<FieldEvaluatorFactory<Real>> evaluator_factory) noexcept
  : layout_(std::move(layout)), evaluator_factory_(std::move(evaluator_factory))
{
}

SplineFunctionSpace SplineFunctionSpace::FromUniformGrid(
  const UniformGrid<2>& grid, CoordinateSystem coordinate_system)
{
  auto layout =
    std::make_shared<UniformGridFieldLayout<2>>(grid, 1, coordinate_system, 1);
  auto evaluator_factory =
    std::make_shared<UniformGridSplineEvaluatorFactory2D>(layout);
  return SplineFunctionSpace(layout, std::move(evaluator_factory));
}

std::shared_ptr<const FieldLayout> SplineFunctionSpace::GetLayout()
  const noexcept
{
  return layout_;
}

CoordinateSystem SplineFunctionSpace::GetCoordinateSystem() const noexcept
{
  return evaluator_factory_->GetCoordinateSystem();
}

FieldVariant SplineFunctionSpace::CreateFieldImpl(Type value_type,
                                                  FieldMetadata metadata) const
{
  return apply_to_type(value_type, [&](auto tag) -> FieldVariant {
    using T = typename decltype(tag)::type;
    if constexpr (!std::is_same_v<T, double>) {
      throw pcms_error("SplineFunctionSpace: only double (Real) is supported");
    } else {
      return WrapField<double>(
        layout_, std::make_unique<SimpleFieldData<double>>(layout_, metadata));
    }
  });
}

FieldVariant SplineFunctionSpace::CreateFieldImpl(FieldDataVariant data) const
{
  if (!std::holds_alternative<std::unique_ptr<FieldData<double>>>(data)) {
    throw pcms_error("SplineFunctionSpace: only double (Real) is supported");
  }
  auto fd = std::move(std::get<std::unique_ptr<FieldData<double>>>(data));
  PCMS_ALWAYS_ASSERT(fd != nullptr);
  if (dynamic_cast<const SimpleFieldData<double>*>(fd.get()) == nullptr) {
    throw pcms_error(
      "SplineFunctionSpace::CreateField: requires SimpleFieldData<double>");
  }
  if (fd->GetDOFHolderDataHost().size() !=
      detail::ExpectedFlatFieldDataSize(*layout_)) {
    throw pcms_error(
      "SplineFunctionSpace::CreateField: field data size does not match "
      "layout");
  }
  return WrapField<double>(layout_, std::move(fd));
}

PointEvaluatorVariant SplineFunctionSpace::CreatePointEvaluatorImpl(
  Type value_type, const EvaluationRequest& request) const
{
  if (value_type != Type::Real) {
    throw pcms_error(
      "SplineFunctionSpace: point evaluation only supports double (Real)");
  }
  if (!evaluator_factory_) {
    throw pcms_error(
      "SplineFunctionSpace::CreatePointEvaluatorImpl: evaluator construction "
      "is not available for this backend");
  }
  return evaluator_factory_->CreatePointEvaluator(request);
}

} // namespace pcms
