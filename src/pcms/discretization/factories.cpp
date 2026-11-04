#include "pcms/discretization/discretization.hpp"
#include "pcms/discretization/empty.hpp"
#include "pcms/backends/uniform_grid/discretization.hpp"

namespace pcms
{

std::shared_ptr<const Discretization> CreateUniformGridDiscretization(
  const UniformGrid<2>& grid)
{
  return std::make_shared<backend::uniform_grid::UniformGridDiscretization<2>>(
    grid);
}

std::shared_ptr<const Discretization> CreateUniformGridDiscretization(
  const UniformGrid<3>& grid)
{
  return std::make_shared<backend::uniform_grid::UniformGridDiscretization<3>>(
    grid);
}

std::shared_ptr<const Discretization> CreateEmptyDiscretization()
{
  return std::make_shared<EmptyDiscretization>();
}

} // namespace pcms
