#include <catch2/catch_session.hpp>
#include <mpi.h>
#include <Kokkos_Core.hpp>
#include <pcms/configuration.h>

#ifdef PCMS_ENABLE_PETSC
#include <petscsys.h>
#endif

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  int result = 0;
  {
    // PETSc uses Kokkos, so initialize Kokkos before PETSc when enabled.
    Kokkos::ScopeGuard kokkos{argc, argv};
#ifdef PCMS_ENABLE_PETSC
    PetscBool petsc_initialized = PETSC_FALSE;
    PetscBool petsc_initialized_by_main = PETSC_FALSE;
    PetscInitialized(&petsc_initialized);
    if (!petsc_initialized) {
      PetscInitialize(&argc, &argv, nullptr, nullptr);
      petsc_initialized_by_main = PETSC_TRUE;
    }
    // petsc must be finalized in the kokkos ScopeGuard scope kokkos needs to
    // deallocate after petsc is done finalizing.
    result = Catch::Session().run(argc, argv);
    if (petsc_initialized_by_main) {
      PetscFinalize();
    }
#else
    result = Catch::Session().run(argc, argv);
#endif
  }
  MPI_Finalize();
  return result;
}
