#include <catch2/catch_session.hpp>
#include <Kokkos_Core.hpp>
#include <mpi.h>


int main( int argc, char* argv[] )
{
  MPI_Init(&argc, &argv);
  Kokkos::ScopeGuard kokkos{};
  int result = Catch::Session().run(argc, argv);
  MPI_Finalize();
  return result;
}
