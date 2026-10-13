#include <pybind11/pybind11.h>
#include "pcms/configuration.h"
#include "numpy_view.h"
#include <stdexcept>

namespace py = pybind11;

namespace pcms
{

// Defined in bind_runtime.cpp. initialize() brings up MPI/Kokkos/PETSc only if
// they are not already running (e.g. an embedding application or another
// pybind11 module such as PyOmega_h already owns them); finalize() tears
// down only the pieces this module itself started.
void initialize();
void finalize();

void bind_coordinate_system_module(py::module& m);

void bind_value_view_module(py::module& m);

void bind_omega_h_field(py::module& m);

void bind_uniform_grid_field_layout_module(py::module& m);

void bind_create_field_module(py::module& m);

void bind_transfer_field_module(py::module& m);

void bind_mls_interpolation_module(py::module& m);

void bind_mesh_utilities_module(py::module& m);

} // namespace pcms

namespace
{
void py_atexit_finalize()
{
  pcms::finalize();
}
} // namespace

PYBIND11_MODULE(pcms, m)
{
  // Bring up MPI/Kokkos/PETSc if nothing else (e.g. PyOmega_h) has already done
  // so, and finalize whatever we own when the interpreter shuts down.
  pcms::initialize();
  if (Py_AtExit(&py_atexit_finalize) != 0) {
    throw std::runtime_error("pcms: failed to register Py_AtExit finalizer");
  }

  // Whether arrays shared with numpy can alias pcms memory.
  m.attr("device_memory_is_host_accessible") =
    pcms::is_host_accessible_v<pcms::DeviceMemorySpace>;

  // Types used as default arguments are registered before their users:
  // coordinate systems and value bases, then OutOfBoundsPolicy, then the
  // function spaces, evaluators, and transfer operators that take them.
  pcms::bind_coordinate_system_module(m);
  pcms::bind_value_view_module(m);
  pcms::bind_omega_h_field(m);
  pcms::bind_uniform_grid_field_layout_module(m);
  pcms::bind_create_field_module(m);
  pcms::bind_transfer_field_module(m);
  pcms::bind_mls_interpolation_module(m);
  pcms::bind_mesh_utilities_module(m);
}
