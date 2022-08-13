#define CATCH_CONFIG_RUNNER
#include <catch2/catch.hpp>
#include <Kokkos_Core.hpp>
#include <mpi.h>


int main( int argc, char* argv[] ) {
  MPI_Init(&argc, &argv);
  Kokkos::ScopeGuard kokkos_scope();
  int result = Catch::Session().run( argc, argv );
  MPI_Finalize();
  return result;
}