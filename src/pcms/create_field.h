#ifndef CREATE_FIELD_H_
#define CREATE_FIELD_H_

#include "adapter/omega_h/omega_h_field_layout.h"
#include "field_layout.h"
#include "field.h"
#include "coordinate_system.h"

#include <Omega_h_mesh.hpp>
#include <memory>

namespace pcms
{
std::unique_ptr<FieldLayout> CreateLagrangeLayout(
  Omega_h::Mesh& mesh, int order, int num_components,
  CoordinateSystem coordinate_system);
} // namespace pcms

#endif // CREATE_FIELD_H_
