#include "create_field.h"
#include "adapter/omega_h/omega_h_field2.h"
#include "adapter/omega_h/omega_h_field_layout.h"

#include <utility>

namespace pcms {

std::unique_ptr<FieldLayout> CreateLagrangeLayout(
  Omega_h::Mesh& mesh, int order, int num_components,
  CoordinateSystem coordinate_system)
{

  std::array<int, 4> nodes_per_dim;

  switch (order) {
    case 1: nodes_per_dim = {1, 0, 0, 0}; break;
    case 2: nodes_per_dim = {1, 1, 0, 0}; break;
    default: throw std::runtime_error("Unimplemented order");
  }

  return std::make_unique<OmegaHFieldLayout>(mesh, nodes_per_dim,
                                             num_components, coordinate_system);
}

} // namespace pcms
