#define CATCH_CONFIG_RUNNER
#include <catch2/catch.hpp>

int main( int argc, char* argv[] ) {
  // initialize kokkos/redev/mpi/omega_h/whatever
  int result = Catch::Session().run( argc, argv );
  //deinitialize kokkos/redev/mpi
  return result;
}