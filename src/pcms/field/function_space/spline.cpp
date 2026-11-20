#include "pcms/field/function_spaces.hpp"

#include "pcms/field/detail/simple_field_data.hpp"
#include "pcms/backends/uniform_grid/uniform_grid_geometry.hpp"
#include "pcms/backends/uniform_grid/spline_evaluator.hpp"
#include "pcms/backends/uniform_grid/layout.hpp"
#include "pcms/utility/detail/common.hpp"

namespace pcms
{
using uniform_grid::UniformGridFieldLayout;
using uniform_grid::UniformGridGeometry;
using uniform_grid::UniformGridSplineEvaluatorFactory;
using uniform_grid::VisitUniformGrid;
} // namespace pcms

namespace pcms
{

namespace
{

template <unsigned Dim>
std::shared_ptr<const FunctionSpace> CreateUniformGridSplineSpace(
  std::shared_ptr<const UniformGridGeometry<Dim>> grid_geometry,
  int num_components)
{
  auto layout = std::make_shared<UniformGridFieldLayout<Dim>>(
    grid_geometry, num_components, 1);
  auto evaluator_factory =
    std::make_shared<UniformGridSplineEvaluatorFactory<Dim>>(
      std::move(grid_geometry), layout);
  return std::make_shared<FunctionSpace>(layout, evaluator_factory,
                                         evaluator_factory->GetGeometry());
}

} // namespace

std::shared_ptr<const FunctionSpace> CreateSplineSpace(
  std::shared_ptr<const Geometry> geometry, int num_components)
{
  if (geometry == nullptr) {
    throw pcms_error("CreateSplineSpace: the geometry must not be null");
  }
  if (auto space =
        VisitUniformGrid(geometry,
                         [num_components]<unsigned Dim>(
                           std::shared_ptr<const UniformGridGeometry<Dim>> grid)
                           -> std::shared_ptr<const FunctionSpace> {
                           if constexpr (Dim == 3) {
                             return nullptr;
                           } else {
                             return CreateUniformGridSplineSpace(
                               std::move(grid), num_components);
                           }
                         })) {
    return space;
  }
  throw pcms_error(
    "CreateSplineSpace: only 1D and 2D uniform-grid geometries are supported");
}

} // namespace pcms
