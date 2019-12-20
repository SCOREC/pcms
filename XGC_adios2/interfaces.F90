module xgc_interfaces
#include <petscversion.h>
#if PETSC_VERSION_LT(3,6,0)
#include <finclude/petscdef.h>
#else
#include <petsc/finclude/petscdef.h>
#endif
  ! TS interfaces
#ifdef XGC1_EM
  interface ts_init
     ! TS init
     subroutine ts_init(a_grid,a_psn,a_bc,a_ts,ierr)
       use grid_class
       use psn_class
       use boundary_class
       use xgc_ts_module
       type(grid_type),intent(in)::a_grid
       type(psn_type)::a_psn
       type(boundary2_type),intent(in)::a_bc
       type(xgc_ts)::a_ts
       integer,intent(out)::ierr
     end subroutine ts_init
  end interface
  ! TS solve
  interface ts_solve
     subroutine ts_solve(ts, dt, n1, apar, phi, Je, ierr)
       use xgc_ts_module
       implicit none
       type(xgc_ts)::ts
       real (kind=8),dimension(ts%nnode)::n1,apar,phi,Je
       PetscErrorCode::ierr
       !integer,intent(out)::ierr
       real (kind=8),intent(in)::dt
     end subroutine ts_solve
  end interface
  !
#endif
end module xgc_interfaces
