#define XGC1 1
#include <petscversion.h>
#if PETSC_VERSION_LT(3,6,0)
#include <finclude/petscsnesdef.h>
#else
#include <petsc/finclude/petscsnesdef.h>
#endif


!< Wrapper routine for PetscOptionsClearValue
!< due to change of its interface from PETSc
!< 3.6 to 3.7
!! @param value (in) name of the option to be cleared (character)
!! @param ierr (inout) error handler (PetscErrorCode)
subroutine my_PetscOptionsClearValue(value,ierr)
  use petscsys
  implicit none
  character(*), intent(in) :: value
  PetscErrorCode, intent(inout) :: ierr
  logical :: petsc_370

#if PETSC_VERSION_LT(3,7,0)
  petsc_370=.false.
#else
  petsc_370=.true.
#endif

  if (petsc_370) then
    call PetscOptionsClearValue(PETSC_NULL_OBJECT,value,ierr)
  else
    call PetscOptionsClearValue(value,ierr)
  endif

  return

end subroutine my_PetscOptionsClearValue

