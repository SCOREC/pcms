#ifdef PCMS_TRANSFER_PETSC_DEFAULTS_HPP
#define PCMS_TRANSFER_PETSC_DEFAULTS_HPP
#include <petscsys.h>

namespace pcms
{

inline void apply_petsc_defaults()
{
  PetscBool initialized = PETSC_FALSE;
  PetscInitialized(&initialized);
  if (!initialized) {
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ORDER,
            "apply_petsc_defaults() must be called after PetscInitialize().");
  }

  // Set defaults (user can override by passing options before initialization)
  PetscOptionsSetValue(NULL, "-mat_type", "aijkokkos");
  PetscOptionsSetValue(NULL, "-vec_type", "kokkos");
  PetscOptionsSetValue(NULL, "-use_gpu_aware_mpi", "0");
}

} // namespace pcms
#endif
