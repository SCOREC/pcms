#include "pcms/field/function_spaces.hpp"

#include "pcms/field/simple_field_data.hpp"
#include "pcms/backends/uniform_grid/geometry.hpp"
#include "pcms/backends/uniform_grid/spline_evaluator.hpp"
#include "pcms/backends/uniform_grid/layout.hpp"
#include "pcms/utility/detail/common.hpp"

namespace pcms
{
using backend::uniform_grid::UniformGridFieldLayout;
using backend::uniform_grid::UniformGridGeometry;
using backend::uniform_grid::UniformGridSplineEvaluatorFactory2D;
} // namespace pcms

namespace pcms
{

std::shared_ptr<const FunctionSpace> CreateSplineSpace(
  std::shared_ptr<const Geometry> geometry)
{
  if (geometry == nullptr) {
    throw pcms_error("CreateSplineSpace: the geometry must not be null");
  }
  auto grid_geometry =
    std::dynamic_pointer_cast<const UniformGridGeometry<2>>(geometry);
  if (grid_geometry == nullptr) {
    throw pcms_error(
      "CreateSplineSpace: only 2D uniform-grid geometries are supported");
  }
  auto layout =
    std::make_shared<UniformGridFieldLayout<2>>(grid_geometry, 1, 1);
  auto evaluator_factory =
    std::make_shared<UniformGridSplineEvaluatorFactory2D>(
      std::move(grid_geometry), layout);
  return std::make_shared<FunctionSpace>(
    layout,
    [layout](Type t, ValueBasis basis) -> FieldDataVariant {
      return apply_to_type(t, [&](auto tag) -> FieldDataVariant {
        using T = typename decltype(tag)::type;
        if constexpr (!std::is_same_v<T, double>) {
          throw pcms_error("CreateSplineSpace: only double (Real) is "
                           "supported");
        } else {
          return std::make_unique<SimpleFieldData<double>>(layout,
                                                           std::move(basis));
        }
      });
    },
    std::move(evaluator_factory));
}

} // namespace pcms
