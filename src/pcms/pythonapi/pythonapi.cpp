
#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace pcms
{

void bind_omega_h_field2(py::module& m);

void bind_transfer_field2_module(py::module& m);

void bind_omega_h_mesh_module(py::module& m);

void bind_coordinate_system_module(py::module& m);

void bind_coordinate_module(py::module& m);

void bind_create_field_module(py::module& m);

void bind_field_layout_module(py::module& m);

void bind_field_module(py::module& m);

void bind_omega_h_field_layout_module(py::module& m);

void bind_uniform_grid_field_layout_module(py::module& m);

void bind_uniform_grid_field_module(py::module& m);

void bind_mls_interpolation_module(py::module& m);

void bind_mesh_utilities_module(py::module& m);

} // namespace pcms
PYBIND11_MODULE(pcms, m)
{
  // Bind fundamental types first (coordinate systems, etc.)
  pcms::bind_coordinate_system_module(m);
  pcms::bind_coordinate_module(m);

  // Bind mesh and field infrastructure
  pcms::bind_omega_h_mesh_module(m);
  pcms::bind_field_layout_module(m);
  pcms::bind_field_module(m);
  pcms::bind_omega_h_field_layout_module(m);
  pcms::bind_uniform_grid_field_layout_module(m);
  pcms::bind_create_field_module(m);

  // Bind field types
  pcms::bind_omega_h_field2(m);
  pcms::bind_uniform_grid_field_module(m);

  // Bind field operations
  pcms::bind_transfer_field2_module(m);

  // Bind interpolator operations
  pcms::bind_mls_interpolation_module(m);

  // Bind mesh utility functions
  pcms::bind_mesh_utilities_module(m);
}