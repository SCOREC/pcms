#include <pcms/capi/kokkos.h>
#include <Kokkos_Core.hpp>
void pcms_kokkos_initialize_without_args() {
  Kokkos::initialize();
}
void pcms_kokkos_finalize() {
  Kokkos::finalize();
}
