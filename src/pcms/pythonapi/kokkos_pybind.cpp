#include <pcms/capi/kokkos.h>
#include <Kokkos_Core.hpp>
#include <pybind11/pybind11.h>


namespace py = pybind11;


PYBIND11_MODULE(pcms_kokkos_module, m) {
    m.def("initialize_kokkos", &pcms_kokkos_initialize_without_args, "Initialize Kokkos");
    m.def("finalize_kokkos", &pcms_kokkos_finalize, "Finalize Kokkos");
}



