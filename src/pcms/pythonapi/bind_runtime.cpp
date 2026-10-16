#include "pcms/configuration.h"
#include <mpi.h>
#include "pcms/utility/assert.h"
#include <Kokkos_Core.hpp>
#if defined(PCMS_ENABLE_PETSC) && defined(PCMS_ENABLE_MESHFIELDS)
#include <petscsys.h>
#endif

namespace pcms
{

namespace
{
bool g_py_owns_mpi = false;
bool g_py_owns_kokkos = false;
bool g_py_owns_petsc = false;
} // namespace

// Initializes MPI/Kokkos/PETSc only if they are not already running (e.g.
// started by an embedding application or another module such as PyOmega_h), and
// records ownership so finalize() only tears down what this module started.
void initialize()
{
  int is_initialized = 0;
  MPI_Initialized(&is_initialized);
  if (!is_initialized) {
    if (MPI_Init(nullptr, nullptr) != MPI_SUCCESS) {
      throw pcms_error("pcms: MPI_Init failed");
    }
    g_py_owns_mpi = true;
  }
  if (!Kokkos::is_initialized()) {
    Kokkos::initialize();
    g_py_owns_kokkos = true;
  }
#if defined(PCMS_ENABLE_PETSC) && defined(PCMS_ENABLE_MESHFIELDS)
  PetscBool petsc_initialized = PETSC_FALSE;
  PetscInitialized(&petsc_initialized);
  if (!petsc_initialized) {
    if (PetscInitializeNoArguments() != 0) {
      throw pcms_error("pcms: PetscInitialize failed");
    }
    g_py_owns_petsc = true;
  }
#endif
}

void finalize()
{
#if defined(PCMS_ENABLE_PETSC) && defined(PCMS_ENABLE_MESHFIELDS)
  if (g_py_owns_petsc) {
    PetscBool petsc_finalized = PETSC_FALSE;
    PetscFinalized(&petsc_finalized);
    if (!petsc_finalized) {
      PetscFinalize();
    }
    g_py_owns_petsc = false;
  }
#endif
  if (g_py_owns_kokkos && !Kokkos::is_finalized()) {
    Kokkos::finalize();
    g_py_owns_kokkos = false;
  }
  int is_finalized = 0;
  MPI_Finalized(&is_finalized);
  if (g_py_owns_mpi && !is_finalized) {
    MPI_Finalize();
    g_py_owns_mpi = false;
  }
}

} // namespace pcms
