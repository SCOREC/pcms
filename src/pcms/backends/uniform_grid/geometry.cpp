#include "pcms/backends/uniform_grid/uniform_grid_geometry.hpp"
#include "pcms/utility/detail/profile.hpp"

namespace pcms::uniform_grid
{

template <unsigned Dim>
std::shared_ptr<BoundPoints> UniformGridGeometry<Dim>::BindImpl(
  const BindRequest& request) const
{
  PCMS_FUNCTION_TIMER;
  return std::make_shared<pcms::detail::UniformGridBinding<Dim>>(
    shared_from_this(), request, grid_);
}

template class UniformGridGeometry<1>;
template class UniformGridGeometry<2>;
template class UniformGridGeometry<3>;

template <unsigned Dim>
std::shared_ptr<const Geometry> CreateGeometry(
  UniformGrid<Dim> grid, std::shared_ptr<const CoordinateSystem> system)
{
  return std::make_shared<UniformGridGeometry<Dim>>(grid, std::move(system));
}

template std::shared_ptr<const Geometry> CreateGeometry(
  UniformGrid<1>, std::shared_ptr<const CoordinateSystem>);
template std::shared_ptr<const Geometry> CreateGeometry(
  UniformGrid<2>, std::shared_ptr<const CoordinateSystem>);
template std::shared_ptr<const Geometry> CreateGeometry(
  UniformGrid<3>, std::shared_ptr<const CoordinateSystem>);

} // namespace pcms::uniform_grid
