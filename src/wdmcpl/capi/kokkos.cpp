#include <wdmcpl/capi/kokkos.h>
#include <Kokkos_Core.hpp>
void wdmcpl_kokkos_initialize_without_args() {
  Kokkos::initialize();
}
void wdmcpl_kokkos_finalize() {
  Kokkos::finalize();
}
