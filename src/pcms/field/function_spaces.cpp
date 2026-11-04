#include "pcms/field/function_spaces.hpp"
#include "pcms/field/layout/empty.hpp"

namespace pcms
{

std::shared_ptr<const FieldLayout> CreateEmptyLayout(
  std::shared_ptr<const CoordinateSystem> coordinate_system, int num_components)
{
  return std::make_shared<EmptyFieldLayout>(std::move(coordinate_system),
                                            num_components);
}

} // namespace pcms
