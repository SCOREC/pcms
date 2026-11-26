#include <pybind11/pybind11.h>
#include "pcms/configuration.hpp"
#include "pcms/version.hpp"
#include "pcms/pythonapi/numpy_view.hpp"
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

void bind_coordinate_map_module(py::module& m);

void bind_point_cloud_module(py::module& m);

#ifdef PCMS_ENABLE_OMEGA_H
void bind_omega_h_module(py::module& m);
#endif

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

  m.attr("__version__") = pcms::pcms_version();
  m.attr("__git_revision__") = pcms::pcms_git_revision();

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
  pcms::bind_coordinate_map_module(m);
  pcms::bind_point_cloud_module(m);
#ifdef PCMS_ENABLE_OMEGA_H
  pcms::bind_omega_h_module(m);
#endif

  // The capabilities this build has (read-only): the keys are the
  // PCMS_ENABLE_* flags without the prefix, in lowercase; a backend's bindings
  // are the submodule of the same name.
  py::dict features;
  for (const char* name :
       {"omega_h", "omega_h_meshfields", "petsc", "polynomial_reconstruction",
        "xgc", "mfem", "coupler", "c", "fortran"}) {
    features[name] = false;
  }
#ifdef PCMS_ENABLE_OMEGA_H
  features["omega_h"] = true;
#endif
#ifdef PCMS_ENABLE_OMEGA_H_MESHFIELDS
  features["omega_h_meshfields"] = true;
#endif
#ifdef PCMS_ENABLE_PETSC
  features["petsc"] = true;
#endif
#ifdef PCMS_ENABLE_POLYNOMIAL_RECONSTRUCTION
  features["polynomial_reconstruction"] = true;
#endif
#ifdef PCMS_ENABLE_XGC
  features["xgc"] = true;
#endif
#ifdef PCMS_ENABLE_MFEM
  features["mfem"] = true;
#endif
#ifdef PCMS_ENABLE_COUPLER
  features["coupler"] = true;
#endif
#ifdef PCMS_ENABLE_C
  features["c"] = true;
#endif
#ifdef PCMS_ENABLE_FORTRAN
  features["fortran"] = true;
#endif
  m.attr("features") =
    py::module_::import("types").attr("MappingProxyType")(features);
}
