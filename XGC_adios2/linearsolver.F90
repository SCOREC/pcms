#if !defined(NO_PETSC)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     
!     Mark Adams 
!     petsc_xgc.F90
!     13 January 2006
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine petsc_init( ierr )
  !-----[--.----+----.----+----.-----------------------------------------]
  !     initialize solve linear system and solver.
  !     
  !     input:
  !     
  !     output:
  !     ierr: error code
  !     
  !     side effects:
  !     - creates objects in global data in 'linear_solver.h'
  !
  !-----[--.----+----.----+----.-----------------------------------------]

  implicit none 
#include <petscversion.h>
#if PETSC_VERSION_LT(3,6,0)
#include <finclude/petsc.h90>
#else
#include <petsc/finclude/petsc.h90>
#endif
  PetscErrorCode::ierr
  call PetscInitialize('./petsc.rc',ierr)
  !call PetscInitialize(0,ierr)
  CHKERRQ(ierr)

  ! log initialization
  !-----[--.----+----.----+----.-----------------------------------------]
  !call PetscLogEventRegister(init_log, 'Solver init     ',0,ierr)
  !call PetscLogEventRegister(solve_log, 'Linear Solve    ',0,ierr)
  !call PetscLogEventRegister(aux_log, 'Aux 1           ',0,ierr)
  !call PetscLogEventRegister(aux2_log, 'Aux 2           ',0,ierr)
  !-----[--.----+----.----+----.-----------------------------------------]
end subroutine petsc_init

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine petsc_end( ierr )
  !-----[--.----+----.----+----.-----------------------------------------]
  !     cleanup common blocks and modules
  !     
  !     input:
  !     output:
  !     ierr: error code
  !     
  !-----[--.----+----.----+----.-----------------------------------------]
  implicit none
#include <petscversion.h>
#if PETSC_VERSION_LT(3,6,0)
#include <finclude/petsc.h90>
#else
#include <petsc/finclude/petsc.h90>
#endif
  PetscErrorCode::ierr
  call PetscFinalize(ierr)
  CHKERRQ(ierr)
end subroutine petsc_end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

#endif


