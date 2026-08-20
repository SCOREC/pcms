#include "pcms/field/coordinate_system.hpp"
#include "pcms/field/coordinate_systems/cartesian.hpp"
#include "pcms/utility/assert.h"
#include <string>

namespace pcms
{
std::shared_ptr<const CoordinateSystem> ResolveCoordinateSystem(
  std::shared_ptr<const CoordinateSystem> system, int dim)
{
  if (system == nullptr) {
    throw pcms_error(
      "ResolveCoordinateSystem: coordinate system must not be null");
  }
  // To avoid requiring the user specifying the dimension twice in
  // FunctionSpace constructors, we have a default
  // Cartesian constuctor that constructs with dimension 0. 
  if (dynamic_cast<const csys::Cartesian*>(system.get()) != nullptr &&
      system->Dimension() == 0) {
    return csys::Cartesian::Create(dim);
  }
  if (system->Dimension() != dim) {
    throw pcms_error("ResolveCoordinateSystem: coordinate system '" +
                     std::string(system->Kind()) + "' has " +
                     std::to_string(system->Dimension()) +
                     " coordinate columns but the data has " +
                     std::to_string(dim));
  }
  return system;
}

} // namespace pcms
