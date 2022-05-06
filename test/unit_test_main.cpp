#define CATCH_CONFIG_RUNNER
#include <catch2/catch.hpp>

#include <wdmcpl/point_search.h>
#include <Omega_h_build.hpp>
#include <Kokkos_Core.hpp>
#include <mpi.h>

//Omega_h::Library omega_h_library;
int main( int argc, char* argv[] )
{
  MPI_Init(&argc, &argv);
  //omega_h_library = Omega_h::Library(&argc, &argv);
  //[[maybe_unused]] auto world = omega_h_library.world();
  Kokkos::ScopeGuard kokkos{};
  // initialize kokkos/redev/mpi/omega_h/whatever
  int result = Catch::Session().run(argc, argv);
  MPI_Finalize();
  return result;
}