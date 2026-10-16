#include <pybind11/pybind11.h>
#include "pcms/configuration.h"
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

void bind_omega_h_field(py::module& m);

void bind_transfer_field_module(py::module& m);

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

  // Bind fundamental types first (coordinate systems, etc.)
  pcms::bind_coordinate_system_module(m);
  pcms::bind_coordinate_module(m);

  // bind_field_module is a no-op stub —
  // FieldT<T>/LocalizationHint/FieldDataView have been removed from the C++
  // API.
  pcms::bind_field_module(m);
  pcms::bind_uniform_grid_field_layout_module(m);
  // Bind OutOfBoundsPolicy before FunctionSpace so the default argument in
  // create_point_evaluator is a registered Python type at binding time.
  pcms::bind_omega_h_field(m);
  // bind_create_field_module registers FunctionSpace,
  // LagrangeFunctionSpace, and Field<Real>.
  pcms::bind_create_field_module(m);
  // bind_uniform_grid_field_module is a no-op stub — UniformGridField<N> has
  // been removed; use LagrangeFunctionSpace::from_uniform_grid instead.
  pcms::bind_uniform_grid_field_module(m);

  // Bind field operations such as Interpolator and Copy.
  pcms::bind_transfer_field_module(m);

  // Bind interpolator operations
  pcms::bind_mls_interpolation_module(m);

  // Bind mesh utility functions
  pcms::bind_mesh_utilities_module(m);
}
