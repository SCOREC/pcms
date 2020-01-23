#include "redef.h"
#include "petscversion.h"
#if (PETSC_VERSION_MAJOR>2)
#if (PETSC_VERSION_MAJOR>3) || (PETSC_VERSION_MINOR>1)
#include "slepcversion.h"
#endif
#endif
!>Routines for handling the EPS and ST objects used in SLEPc
!!
!!It provides a module variable eps_obj which should be used by all routines based on SLEPc.
module slepc_aux
  use par_in
  use par_mod
  use par_other, only: print_ini_msg
  use discretization
  use eigen_parameters
  use petsc_aux
  use petsc_precond
#if (SLEPC_VERSION_GE(3,8,0))
#include "slepc/finclude/slepceps.h"
#include "slepc/finclude/slepcmfn.h"
  use slepceps
  use slepcmfn
#endif
  use communications
  use petscvec
  use diagnostics, only: diag_nrg,momentum_flux
implicit none

#if (PETSC_VERSION_GE(3,8,0))
!  public :: tMat, tVec
#else
#if (PETSC_VERSION_GE(3,6,0))
#include "slepc/finclude/slepc.h"
#else
#include "finclude/slepc.h"
#if (PETSC_VERSION_MINOR==0)
#include "finclude/slepcst.h"
#include "finclude/slepceps.h"
#endif
#endif
#endif

  public:: initialize_slepc_eps, finalize_slepc_eps, eps_info, eps_obj, &
       & my_SlepcInitialize, my_SlepcFinalize, n_ch, slepc_restartable
  public :: slepc_initialized_from_outside

  private

  PetscInt :: n_ch=1
  logical:: write_mem_req=.true.,slepc_initialized_from_outside=.false.

  EPS eps_obj  !<SLEPc object to describe the eigenvalue problem
  ST st_obj  !<SLEPc object to describe the spectral transform (if any)
#if SLEPC_VERSION_GE(3,7,1)
  MFN mfn
#endif

contains

!>Returns whether PETSc/SLEPc can be initialized/finalized more than once
  logical function slepc_restartable()
    !those variables are taken from c_utils.c since slepcversion.h
    !cannot be read by Fortran
    integer :: get_slepc_version_major, get_slepc_version_minor
    integer :: get_slepc_version_patch

    slepc_restartable = .false.

    !SLEPc is restartable for version > 3.1-p4 (?)
    !\todo add some 3.0 patches which worked as well
    slepc_restartable = (get_slepc_version_major().gt.3).or.&
         &((get_slepc_version_major().eq.3).and.&
         &((get_slepc_version_minor().gt.1).or.&
         &((get_slepc_version_minor().eq.1).and.&
         &(get_slepc_version_patch().gt.4))))

  end function slepc_restartable


!>Initializes the eps_obj with the appropriate problem type, dimensions etc.
!!It also sets the corresponding spectral transform (st_obj) if needed
  subroutine initialize_slepc_eps(ev_run,n_ev_in)

    logical,intent(in):: ev_run !<switches between initialization for eigenvalue runs or time step computation
    integer,intent(in),optional::n_ev_in !<sets the number of eigenvalues to be computed
    integer::n_ev_loc, ev_n_test_loc,mpd_loc
#if (SLEPC_VERSION_GE(3,4,0))
    integer :: dummyi
#endif
    STType st_type
    MatStructure flag,flag1
#if SLEPC_VERSION_GE(3,7,1)
    FN fn
    Mat shellmfn
    PetscScalar :: alpha
    PetscReal :: tol
    integer :: ncv
#endif

    !set appropriate defaults if not user-specified
    if (ev_max_it.le.0) ev_max_it = PETSC_DEFAULT_INTEGER !recommended by J. Roman
    if (ev_n_test.le.0) ev_n_test = PETSC_DECIDE

    if(present(n_ev_in))then
       n_ev_loc = n_ev_in
       ev_n_test_loc = max(2*n_ev_in,n_ev_in+15)
    else
       n_ev_loc = n_ev
       ev_n_test_loc = ev_n_test
    endif

    call EPSCreate(PETSC_COMM_WORLD,eps_obj,globerr)

    !problem definitions
#if SLEPC_VERSION_GE(3,8,0)
    call EPSSetOperators(eps_obj,shellmat,PETSC_NULL_MAT,globerr)
#else
    call EPSSetOperators(eps_obj,shellmat,PETSC_NULL_OBJECT,globerr)
#endif
    call EPSSetProblemType(eps_obj,EPS_NHEP,globerr)

    if (.not.ev_run.and.which_ev=='smallest_real') then
       !smallest real timestep computatios have special (fixed) parameters
       call EPSSetDimensions(eps_obj,n_ev_sr,ev_n_test_sr,PETSC_DECIDE,globerr)
       call EPSSetTolerances(eps_obj,ev_prec,ev_max_it_sr,globerr)
    else
       call EPSSetDimensions(eps_obj,n_ev_loc,ev_n_test_loc,PETSC_DECIDE,globerr)
       call EPSSetTolerances(eps_obj,ev_prec,ev_max_it,globerr)
       !ask for EPS object values automatically chosen by SLEPC
       call EPSGetDimensions(eps_obj,n_ev_loc,ev_n_test_loc,mpd_loc,globerr)
    endif

    if(write_mem_req) then
       !output memory requirements
       if (mype.eq.0) then
          if(ev_run.and.(which_ev.eq.'all_lapack')) then
             write(*,"(a,i8,a,f8.2,a)") 'using LAPACK (not parallelized!), dimensionality of the problem is ',&
                  vlen*n_procs_sim, ', the corresponding matrix has ',&
                  vlen*vlen*n_procs_sim*n_procs_sim*16.0/1024**2,'MB'
          elseif (print_ini_msg) then
             write(*,"(4(a))") 'EV solver/preconditioner: ',trim(which_ev),&
                  &' / ',trim(pc_type)
             write(*,"(a,i4,a,i10,a,f8.2,a)") 'using ',ev_n_test_loc,' test vectors with length',&
                  vlen*n_procs_sim, ', i.e.',&
                  ev_n_test_loc*vlen*n_procs_sim*16.0/1024**2,'MB'
             write(*,"(a,G12.4)") 'precision: ', ev_prec
          endif
       endif
       write_mem_req=.false.
    endif

    !select eigenvalues and method
    if (.not.ev_run) then
       !slepc only used to determine timestep
       select case (which_ev)
       case('smallest_real')
          call EPSSetWhichEigenPairs(eps_obj,EPS_SMALLEST_REAL,globerr)
       case default
          call EPSSetWhichEigenPairs(eps_obj,EPS_LARGEST_MAGNITUDE,globerr)
       endselect
    else
#if ((SLEPC_VERSION_MAJOR==3)&&(SLEPC_VERSION_MINOR==8)&&(SLEPC_VERSION_SUBMINOR==0))
       call EPSSetInitialSpace(eps_obj,1,v0(1),globerr) !n_ch>1 broken in slepc-3.8.0 only
#else
#if ((PETSC_VERSION_MAJOR>3)||(PETSC_VERSION_MINOR>0))
       call EPSSetInitialSpace(eps_obj,n_ch,v0(1:n_ch),globerr)
#else
       call EPSSetInitialVector(eps_obj,v0(1),globerr)
#endif
#endif
       select case (which_ev)
       case('largest_real')
          call EPSSetWhichEigenPairs(eps_obj,EPS_LARGEST_REAL,globerr)
          call EPSSetType(eps_obj,EPSKRYLOVSCHUR,globerr)
       case('smallest_real')
          call EPSSetWhichEigenPairs(eps_obj,EPS_SMALLEST_REAL,globerr)
          call EPSSetType(eps_obj,EPSKRYLOVSCHUR,globerr)
#if (SLEPC_VERSION_GE(3,4,0))
       case('ks_as') !Krylov-Schur Arbitrary Selection
          call EPSSetArbitrarySelection(eps_obj,arb_func,dummyi,globerr)
          call EPSSetWhichEigenPairs(eps_obj,EPS_LARGEST_REAL,globerr)
          call EPSSetType(eps_obj,EPSKRYLOVSCHUR,globerr)
#endif
       case('largest_magnitude')
          call EPSSetWhichEigenPairs(eps_obj,EPS_LARGEST_MAGNITUDE,globerr)
          call EPSSetType(eps_obj,EPSKRYLOVSCHUR,globerr)
#if SLEPC_VERSION_GE(3,7,1)
       case('mfn')
          call MFNCreate(MY_MPI_COMM_WORLD,mfn,globerr)
          call MFNSetOperator(mfn,shellmat,globerr)
          ncv = 10
          call MFNSetDimensions(mfn,ncv,globerr)
          tol = 1e-6
          call MFNSetTolerances(mfn,tol,PETSC_DEFAULT_INTEGER,globerr)
          call MFNGetFN(mfn,fn,globerr)
          call FNSetType(fn,FNEXP,globerr)
          alpha = 1.0
          call FNSetScale(fn,taumfn,alpha,globerr)
          call MFNSetFromOptions(mfn,globerr)
          call MatCreateShell(MY_MPI_COMM_WORLD,vlen,vlen,vlen*n_procs_sim,&
               vlen*n_procs_sim,PETSC_NULL_INTEGER,shellmfn,globerr)
          call MatShellSetOperation(shellmfn,MATOP_MULT,matmult_mfn,globerr)
#if SLEPC_VERSION_GE(3,8,0)
          call EPSSetOperators(eps_obj,shellmfn,PETSC_NULL_MAT,globerr)
#else
          call EPSSetOperators(eps_obj,shellmfn,PETSC_NULL_OBJECT,globerr)
#endif
          call EPSSetWhichEigenPairs(eps_obj,EPS_LARGEST_MAGNITUDE,globerr)
          call EPSSetType(eps_obj,EPSKRYLOVSCHUR,globerr)
          call MatDestroy(shellmfn,globerr)
#endif
       case('harmonic')
#if ((PETSC_VERSION_MAJOR>3)||(PETSC_VERSION_MINOR>0))
          call EPSSetWhichEigenPairs(eps_obj,EPS_TARGET_MAGNITUDE,globerr)
#endif
          call EPSSetTarget(eps_obj,ev_shift,globerr)
          call EPSSetExtraction(eps_obj,EPS_HARMONIC,globerr)
          call EPSSetType(eps_obj,EPSKRYLOVSCHUR,globerr)
#if (SLEPC_VERSION_MAJOR==3) && (SLEPC_VERSION_MINOR==6)
          call EPSKrylovSchurSetLocking(eps_obj,PETSC_FALSE,globerr)
#endif
#if ((PETSC_VERSION_MAJOR>3)||(PETSC_VERSION_MINOR>0))
       case('jd')
          call EPSSetWhichEigenPairs(eps_obj,EPS_TARGET_MAGNITUDE,globerr)
          call EPSSetTarget(eps_obj,ev_shift,globerr)
          call EPSSetType(eps_obj,EPSJD,globerr)
          call EPSSetExtraction(eps_obj,EPS_HARMONIC,globerr) !EPS_HARMONIC_RIGHT?
          call EPSJDSetRestart(eps_obj,8,0,globerr)
          call EPSGetST(eps_obj,st_obj,globerr)
          call STSetType(st_obj,STPRECOND,globerr)
          call STSetShift(st_obj,ev_shift,globerr) !taken care of by EPSSetTarget (?)
          call STGetKSP(st_obj,ksp_obj,globerr)
          call KSPSetType(ksp_obj,KSPBCGS,globerr)
          if(ksp_max_it.eq.0) then
             if (pc_type.eq.'none') then
                ksp_max_it=80
             else
                ksp_max_it=5
             end if
          end if
          call KSPSetTolerances(ksp_obj,1e-16,1e-10,1e4,ksp_max_it,globerr)

       case('gd')
          call EPSSetWhichEigenPairs(eps_obj,EPS_TARGET_MAGNITUDE,globerr)
          call EPSSetTarget(eps_obj,ev_shift,globerr)
          call EPSSetExtraction(eps_obj,EPS_HARMONIC,globerr)
          call EPSSetType(eps_obj,EPSGD,globerr)
          call EPSGetST(eps_obj,st_obj,globerr)
          call STSetShift(st_obj,ev_shift,globerr)
#endif
#if (SLEPC_VERSION_GE(3,4,0))
       case('jd_as') !JD Arbitrary Selection
          call EPSSetArbitrarySelection(eps_obj,arb_func,dummyi,globerr)
          call EPSSetWhichEigenPairs(eps_obj,EPS_LARGEST_REAL,globerr)
          !?Keep the following???
          call EPSSetType(eps_obj,EPSJD,globerr)
          call EPSSetExtraction(eps_obj,EPS_HARMONIC,globerr)
          call EPSJDSetRestart(eps_obj,8,0,globerr)
          call EPSGetST(eps_obj,st_obj,globerr)
          call STSetType(st_obj,STSINVERT,globerr)
          call STSetShift(st_obj,ev_shift,globerr)
          call STGetKSP(st_obj,ksp_obj,globerr)
          call KSPSetType(ksp_obj,KSPBCGS,globerr)
          if(ksp_max_it.eq.0) then
             if (pc_type.eq.'none') then
                ksp_max_it=80
             else
                ksp_max_it=5
             end if
          end if
          call KSPSetTolerances(ksp_obj,1e-16,1e-10,1e4,ksp_max_it,globerr)
#endif
       case('shift_invert')
#if (PETSC_VERSION_GE(3,0,0))
          call EPSSetWhichEigenPairs(eps_obj,EPS_TARGET_MAGNITUDE,globerr)
          call EPSSetTarget(eps_obj,ev_shift,globerr)
          call EPSSetType(eps_obj,EPSKRYLOVSCHUR,globerr);
#else
          call EPSSetWhichEigenPairs(eps_obj,EPS_LARGEST_MAGNITUDE,globerr)
#endif
          call EPSGetST(eps_obj,st_obj,globerr)
#if (PETSC_VERSION_GE(3,0,0))
          !call EPSSetTrueResidual(eps_obj,PETSC_TRUE,globerr)
          call EPSSetTrueResidual(eps_obj,PETSC_FALSE,globerr)
          call STSetType(st_obj,STSINVERT,globerr)
          call STSetMatMode(st_obj,ST_MATMODE_SHELL,globerr)
#else
          call STSetType(st_obj,STSINV,globerr)
          call STSetMatMode(st_obj,STMATMODE_SHELL,globerr)
#endif
          call STSetFromOptions(st_obj,globerr)
          call STSetShift(st_obj,ev_shift,globerr)
          call STGetKSP(st_obj,ksp_obj,globerr)
          !call KSPSetFromOptions(ksp_obj,globerr)
          call KSPSetType(ksp_obj,KSPGMRES,globerr)
          !call STSetFromOptions(st_obj,globerr)
          call STSetUp(st_obj,globerr)
       case ('runtime')
          !do nothing; all options taken from command line
       case('all_lapack')
          call EPSSetType(eps_obj,EPSLAPACK,globerr)
#ifdef WITHSCAL
       case('shift_invert_s')
          !spectral transform is done explicitly
          call EPSSetWhichEigenPairs(eps_obj,EPS_LARGEST_MAGNITUDE,globerr)
          call EPSSetType(eps_obj,EPSKRYLOVSCHUR,globerr);
#endif
       case default
          if (mype.eq.0) print*, 'no valid which_ev selected'
          stop
       end select

       !the eps options can be overwritten at runtime!
       call EPSSetFromOptions(eps_obj,globerr)

       !setup the preconditioner
       call EPSGetST(eps_obj,st_obj,globerr)
       call STGetKSP(st_obj,ksp_obj,globerr)
       call KSPGetPC(ksp_obj,pc_obj,globerr)
       call PCSetType(pc_obj,PCNONE,globerr)
       !call KSPSetType(ksp_obj,KSPGMRES,globerr)
       !call KSPGMRESSetRestart(ksp_obj,12,globerr)
       !call KSPSetTolerances(ksp_obj,1e-5,1e-35,100000,globerr)
       if (pc_type .ne. PCNONE) then
          call STGetType(st_obj,st_type,globerr)
          !if (mype.eq.0) write(*,*) 'ST_TYPE =',trim(st_type)
          select case(st_type)
#if (PETSC_VERSION_GE(3,0,0))
          case(STPRECOND)
             call initialize_L_g
             call PCSetType(pc_obj,PCJACOBI,globerr)
             call STPrecondSetMatForPC(st_obj,L_g_mat,globerr)
             call EPSSetUp(eps_obj,globerr)
             call set_pc_obj
             call finalize_L_g
          case(STSINVERT)
#else
          case(STSINV)
#endif
            ! call EPSSetUp(eps_obj,globerr)
             call initialize_L_g
             flag = DIFFERENT_NONZERO_PATTERN
#if (PETSC_VERSION_GE(3,5,0))
             call KSPGetOperators(ksp_obj,a_mat,b_mat,globerr)
             call KSPSetOperators(ksp_obj,a_mat,L_g_mat,globerr)
#else
             call KSPGetOperators(ksp_obj,a_mat,b_mat,flag1,globerr)
             !call PetscObjectReference(a_mat,globerr)
             call KSPSetOperators(ksp_obj,a_mat,L_g_mat,flag,globerr)
#endif
!            call KSPGetPC(ksp_obj,pc_obj,globerr)
!            call PCSetOperators(pc_obj,a_mat,L_g_mat,globerr)
            !call MatDestroy(a_mat,globerr)
             call KSPSetFromOptions(ksp_obj,globerr)
             call KSPSetUp(ksp_obj,globerr)
             call set_pc_obj
             call finalize_L_g
             call EPSSetUp(eps_obj,globerr)
          case default
             if (mype.eq.0) write(*,'(2A)') 'No special treatment found for ST_TYPE=',&
                  &trim(st_type)
          end select
       end if
    end if

    !if (mype.eq.0) call EPSView(eps_obj,PETSC_VIEWER_STDOUT_SELF,globerr)
  end subroutine initialize_slepc_eps

#if SLEPC_VERSION_GE(3,7,1)
  !>Used for the shell matrix in the mfn solver
  subroutine matmult_mfn(A,x,y,ierr)
    Mat A
    Vec x,y
    PetscErrorCode ierr
    call MFNSolve(mfn,x,y,ierr)
  end subroutine matmult_mfn
#endif

  !>Deletes the eps_obj
  subroutine finalize_slepc_eps
    if (L_g_initialized) call finalize_L_g
#if SLEPC_VERSION_GE(3,7,1)
    if (which_ev.eq.'mfn') call MFNDestroy(mfn,globerr)
#endif
    call EPSDestroy(eps_obj,globerr)
  end subroutine finalize_slepc_eps

  !>Routine to get and output the number of computed eigenvalues and SLEPc iterations needed
  subroutine eps_info(ev_number)
    integer, intent(out):: ev_number !< number of eigenvalues contained in eps_obj

    call EPSGetIterationNumber(eps_obj,it_ev,globerr)
    if (mype.eq.0) write(*,"(a,i6)") 'number of iterations:',it_ev

    call EPSGetConverged(eps_obj,ev_number,globerr)
    if (mype.eq.0) write(*,"(a,i6,a)") 'calculated ',ev_number,' eigenvectors'

    if (ev_number.eq.0) then
       if (mype.eq.0) then
          write(*,"(a)") '***NO EIGENVALUES CALCULATED***'
       endif
    end if

  end subroutine eps_info

  subroutine my_SlepcInitialize(slepc_communicator,call_from_outside)
    integer,intent(in):: slepc_communicator
    logical, intent(in):: call_from_outside

    if (.not.slepc_initialized_from_outside) then
       PETSC_COMM_WORLD=slepc_communicator
       call SlepcInitialize(PETSC_NULL_CHARACTER,globerr)
    end if
    if(call_from_outside) slepc_initialized_from_outside=.true.

  end subroutine my_SlepcInitialize

  subroutine my_SlepcFinalize(call_from_outside)
    logical, intent(in):: call_from_outside

    if(call_from_outside.or.(.not.slepc_initialized_from_outside)) then
       call SlepcFinalize(globerr)
    end if
  end subroutine my_SlepcFinalize


!!!!!!!!!!!!!!!For Arbitrary Selection !!!!!!!!!!!!!!!!!!!!!!!
!Calculates the electron electromagnetic heat flux normalized to energy
  subroutine arb_func(ev,evdummy,x0,xdummy,r0,rdummy,dummyi,ierr)
#if (PETSC_VERSION_GE(3,6,0))
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscmat.h"
#else
#if !((PETSC_VERSION_MAJOR>2) && (PETSC_VERSION_MINOR>0))
#include "finclude/petsc.h"
#endif
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#endif
    Vec x0,xdummy
    PetscErrorCode ierr
    PetscScalar ev,evdummy,r0,rdummy
    logical :: dummy1,dummy2
    real, allocatable, dimension(:,:) :: nrg_data
    real :: qem_norm
    integer :: dummyi

    call pe2fo(x0,g_1)

    if(which_func_as==1) then
       if(momentum_flux) then
          allocate(nrg_data(10,ln1:ln2))
       else
          allocate(nrg_data(8,ln1:ln2))
       end if

       call diag_nrg(dummy1,dummy2,nrg_data,.false.)
       qem_norm=0.0
       if(n_spec.lt.2) then
          stop "Must have n_spec>1 for arbitrary selection with QeEM."
       end if
       !Get electron EM heat flux
       if(1.ge.ln1.and.1.le.ln2) then
          if(nrg_data(1,1).gt.0.0) then
             qem_norm=nrg_data(8,1)/nrg_data(1,1)
          else
             qem_norm=0.0
          end if
       end if
       call my_sum_to_all_real_0d(qem_norm,MY_MPI_COMM_WORLD)
       r0=cmplx(qem_norm,0)
       !write(*,*) "r0",r0,mype
       deallocate(nrg_data)
    elseif(which_func_as==2) then
       !write(*,*) 'ev',ev
       r0=ev
       !write(*,*) 'r0',r0
    else
       stop "Invalid selection for which_func_as!!!"
    end if


  end subroutine arb_func

end module slepc_aux
