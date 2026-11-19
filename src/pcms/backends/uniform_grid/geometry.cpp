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

template class UniformGridGeometry<2>;
template class UniformGridGeometry<3>;

std::shared_ptr<const Geometry> CreateGeometry(
  UniformGrid<2> grid, std::shared_ptr<const CoordinateSystem> system)
{
  return std::make_shared<UniformGridGeometry<2>>(grid, std::move(system));
}

std::shared_ptr<const Geometry> CreateGeometry(
  UniformGrid<3> grid, std::shared_ptr<const CoordinateSystem> system)
{
  return std::make_shared<UniformGridGeometry<3>>(grid, std::move(system));
}

} // namespace pcms::uniform_grid
