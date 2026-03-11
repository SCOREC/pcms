#include <catch2/catch_session.hpp>
#include <Kokkos_Core.hpp>
#include <mpi.h>
#include <petsc_defaults.hpp>

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  int result = 0;
  {
    // petsc uses kokkos, so it must be initialized before petsc
    Kokkos::ScopeGuard kokkos{};
    PetscBool petsc_initialized = PETSC_FALSE;
    PetscBool petsc_initialized_by_main = PETSC_FALSE;
    PetscInitialized(&petsc_initialized);
    if (!petsc_initialized) {
      PetscInitialize(&argc, &argv, nullptr, nullptr);
      pcms::apply_petsc_defaults();
      petsc_initialized_by_main = PETSC_TRUE;
    }
    // petsc must be finalized in the kokkos ScopeGuard scope kokkos needs to
    // deallocate after petsc is done finalizing.
    result = Catch::Session().run(argc, argv);
    if (petsc_initialized_by_main) {
      PetscFinalize();
    }
  }
  MPI_Finalize();
  return result;
}
