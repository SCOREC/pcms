#include <petscversion.h>
#if PETSC_VERSION_LT(3,6,0)
#include <finclude/petscsnesdef.h>
#else
#include <petsc/finclude/petscsnesdef.h>
#endif
! ccccccccccccccccc subroutine petsc_solve cccccccccccccccccccccccccccccc
!
#undef __FUNCT__
#define __FUNCT__ "petsc_solve"
subroutine petsc_solve(nn,bb,bb2,xx,comm,solver,ierr)
!c     solve linear system.
!c     Solver (KSP) and matrix are in linear_solve.h
!c
!c     input:
!c       nn: global number of equations in system
!c       bb: rhs[nn]
!c       bb2: rhs[nn] -- not used
!c       comm: communicator for sysetem (ie, a poloidal plane)
!c
!c     in/output:
!c       xx: solution[nn]
!c     ierr: error code
!c
  use petscsnes
  use xgc_solver_module
  use sml_module
  use perf_monitor
  implicit none
#include <petscversion.h>
  integer::nn,comm
  PetscScalar::xx(nn),bb(nn),bb2(nn)
  type(xgc_solver)::solver
  PetscInt::i
  PetscErrorCode::ierr
  !
  PetscInt::nnodes64
  PetscScalar::norm,norm0,norm2
  Vec::xxvec,Xsub(2),Fsub(2),temp_vec
  PetscInt,parameter::ione=1
  PetscInt,parameter::itwo=2
  PetscViewer::viewer  
  KSP::innerksp,subksp(2)
  PC::spc
  
  nnodes64 = nn

  ! scatter rhs into PETSc vector
  call VecCreateSeqWithArray(PETSC_COMM_SELF,ione,nnodes64,bb,xxvec,ierr)
  call VecScatterBegin(solver%to_petsc,xxvec,solver%bVec,INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(  solver%to_petsc,xxvec,solver%bVec,INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecDestroy(xxvec,ierr)

  if( solver%n_rhs_mat>=1 ) then
     ! bb <- (I + A) * bb or I * bb
     call MatMult( solver%rhs_mat,solver%bVec, solver%xVec, ierr )
     call VecCopy( solver%xVec, solver%bVec, ierr )
  else
     stop ' solver%n_rhs_mat == 0 ???'
  end if

  if( solver%n_rhs_mat>=2 ) then
     ! add one more vector to solve
     call VecCreateSeqWithArray(PETSC_COMM_SELF,ione,nnodes64,bb2,xxvec,ierr)
     call VecDuplicate(solver%xVec,temp_vec,ierr)
     call VecScatterBegin(solver%to_petsc,xxvec,temp_vec,INSERT_VALUES,SCATTER_FORWARD,ierr)
     call VecScatterEnd(  solver%to_petsc,xxvec,temp_vec,INSERT_VALUES,SCATTER_FORWARD,ierr)
     call VecDestroy(xxvec,ierr)

     call MatMult(solver%rhs2_mat,temp_vec,solver%xVec,ierr)
     call VecAXPY(solver%bVec,1D0,solver%xVec,ierr)
     call VecDestroy(temp_vec,ierr)
  endif

  ! scatter inital guess
  call VecCreateSeqWithArray(PETSC_COMM_SELF,ione,nnodes64,xx,xxvec,ierr)
  if (solver%snes==0) then
     call VecScatterBegin(solver%to_petsc,xxvec,solver%xVec,INSERT_VALUES,SCATTER_FORWARD,ierr)
     call VecScatterEnd(  solver%to_petsc,xxvec,solver%xVec,INSERT_VALUES,SCATTER_FORWARD,ierr)
     ! 1 field Solve
     call t_startf("PETSC_SOLVE")
     call KSPSolve(solver%ksp,solver%bVec,solver%xVec,ierr)
     CHKERRQ(ierr)
     call t_stopf("PETSC_SOLVE")
  else
     call DMCompositeGetAccessArray(solver%da,solver%xVec2,itwo,PETSC_NULL_INTEGER,Xsub,ierr)
     call DMCompositeGetAccessArray(solver%da,solver%bVec2,itwo,PETSC_NULL_INTEGER,Fsub,ierr)
     ! initial guess
     call VecScatterBegin(solver%to_petsc,xxvec,Xsub(1),INSERT_VALUES,SCATTER_FORWARD,ierr)
     call VecScatterEnd(  solver%to_petsc,xxvec,Xsub(1),INSERT_VALUES,SCATTER_FORWARD,ierr)

     call VecZeroEntries(Fsub(2),ierr)
     call VecCopy(solver%bVec,Fsub(1),ierr) ! RHS 
     call MatMult(solver%Cmat,Xsub(1),Xsub(2),ierr) ! init <phi>

     ! 2 field Solve (1 or 2)
     call t_startf("PETSC_SNES_SOLVE")
     call SNESSolve(solver%snes,solver%BVec2,solver%XVec2,ierr)
     CHKERRQ(ierr)
     call t_stopf("PETSC_SNES_SOLVE")
     
     call VecCopy(Xsub(1),solver%xVec,ierr)
     if (.true.) then ! debug
        call VecNorm(Xsub(1),NORM_INFINITY,norm2,ierr)
        call MatMult(solver%Cmat,Xsub(1),Fsub(2),ierr)
        call VecAXPY(Fsub(2),-1D0,Xsub(2),ierr)
        call VecNorm(Fsub(2),NORM_INFINITY,norm,ierr)
        call VecNorm(Xsub(2),NORM_INFINITY,norm0,ierr)
        if(sml_mype==0 .and. norm0.ne.0d0) then
           write(*,'(A,ES11.4E2,A,ES11.4E2,A,ES11.4E2)') '     |C phi - <phi>| / |<phi>| = ',norm/norm0,', |<phi>|=',norm0,', |<phi>|/|phi|=',norm0/norm2
        else if(sml_mype==0) then
           write(*,'(A,ES11.4E2,A,ES11.4E2,A,ES11.4E2)') '     |C phi - <phi>| = ',norm,', |<phi>|=',norm0,', |phi|=',norm2
        end if
     end if

     if(sml_plane_index==0 .and. .false.) then ! debug
        call  PetscViewerASCIIOpen(comm, 'Phi1.m', viewer,ierr)
        call  PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB,ierr)
        call  VecView(Xsub(1),viewer,ierr)
        call  PetscViewerDestroy(viewer,ierr)
        call  PetscViewerASCIIOpen(comm, 'Phi2.m', viewer,ierr)
        call  PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB,ierr)
        call  VecView(Xsub(2),viewer,ierr)
        call  PetscViewerDestroy(viewer,ierr)
        call  PetscViewerASCIIOpen(comm, 'rho.m', viewer,ierr)
        call  PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB,ierr)
        call  VecView(Fsub(1),viewer,ierr)
        call  PetscViewerDestroy(viewer,ierr)        
     endif
     !call check_point('petscoutput')     
     call DMCompositeRestoreAccessArray(solver%da,solver%xVec2,itwo,PETSC_NULL_INTEGER,Xsub,ierr)
     call DMCompositeRestoreAccessArray(solver%da,solver%bVec2,itwo,PETSC_NULL_INTEGER,Fsub,ierr)
  end if

  ! scatter solution to xgc array
  call VecScatterBegin(solver%from_petsc,solver%xVec,xxvec,INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(  solver%from_petsc,solver%xVec,xxvec,INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecDestroy(xxvec,ierr)

  return
end subroutine petsc_solve

! ---------------------------------------------------------------------
!
!  FormJacobian - Evaluates Jacobian matrix.
!
!  Input Parameters:
!  dummy     - the SNES context
!  x         - input vector
!  solver    - solver data
!
!  Output Parameters:
!  jac      - Jacobian matrix
!  jac_prec - optionally different preconditioning matrix (not used here)
!  flag     - flag indicating matrix structure
!
#undef __FUNCT__
#define __FUNCT__ "FormJacobian"
#if (!PETSC_VERSION_LT(3,5,0))
subroutine FormJacobian(dummy,X,jac,jac_prec,solver,ierr)
#else
subroutine FormJacobian(dummy,X,jac,jac_prec,flag,solver,ierr)
#endif
  use petscsnes
  use xgc_solver_module
  implicit none
!  Input/output variables:
  SNES     dummy
  Vec      X
  Mat      jac,jac_prec
  MatStructure   flag
  type(xgc_solver)  solver
  PetscErrorCode ierr

  ! NOT doing anything now, linear preconditioner
  !  Declarations for use with local arrays:
  !Vec      Xsub(1)
  !Mat      Amat
  !PetscInt, parameter :: ione=1
  !call DMCompositeGetAccessArray(solver%da,X,ione,PETSC_NULL_INTEGER,Xsub,ierr)
  ! Compute entries for the locally owned part of the Jacobian preconditioner.
  ! call MatGetSubMatrix(jac_prec,solver%isPhi,solver%isPhi,MAT_INITIAL_MATRIX,Amat,ierr)
  ! form A, zero out and relinarize around Xsub (the 00 block)
  !call FormJacobianLocal(Xsub(1),Amat,solver,.true.,ierr)
  !call MatDestroy(Amat,ierr) ! discard our reference
  !call DMCompositeRestoreAccessArray(solver%da,X,ione,PETSC_NULL_INTEGER,Xsub,ierr)
  
  ! the rest of the matrix is not touched
  call MatAssemblyBegin(jac_prec,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(jac_prec,MAT_FINAL_ASSEMBLY,ierr)
  if (jac .ne. jac_prec) then
     call MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY,ierr)
     call MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY,ierr)
  end if

  !     Set flag to indicate that the Jacobian matrix retains an identical
  !     nonzero structure throughout all nonlinear iterations 

  flag = SAME_NONZERO_PATTERN

  !     Tell the matrix we will never add a new nonzero location to the
  !     matrix. If we do it will generate an error.
  call MatSetOption(jac_prec,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE,ierr)

  return
end subroutine FormJacobian

! ---------------------------------------------------------------------
!
!  FormFunction - Evaluates nonlinear function, F(x).
!
!  Input Parameters:
!     snes_dummy - the SNES context. not used
!     X          - input vector
!     solver     - user-defined context 
!
!  Output Parameter:
!     F          - function vector
!     ierr       - error code
!
#undef __FUNCT__
#define __FUNCT__ "FormFunction"
subroutine FormFunction(snes_dummy,X,F,solver,ierr)
  use petscsnes
  use xgc_solver_module
  use sml_module,only:sml_mype,sml_e_charge,sml_poisson_solver_type,sml_plane_index
  use eq_module
  use eq_module
  implicit none
  !  Input/output variables:
  SNES::snes_dummy
  Vec::X,F
  PetscErrorCode::ierr
  type(xgc_solver)::solver
  !
  PetscInt,parameter::itwo=2
  Vec::Xsub(itwo),Fsub(itwo),ave1,Q
  PetscInt::low,high
  PetscScalar::value,t1,t2
  PetscScalar,pointer::ave1_v(:),Q_v(:),phi_v(:),TeInv_v(:),del_v(:)
  integer::ii,ip,idxj,nloc,nloc2
  real (8)::pn,wp
#if defined(PETSC_USE_LOG)
  call PetscLogEventBegin(FormFunctionEvent,ierr)
#endif
  call DMCompositeGetAccessArray(solver%da,X,itwo,PETSC_NULL_INTEGER,Xsub,ierr);CHKERRQ(ierr)
  call DMCompositeGetAccessArray(solver%da,F,itwo,PETSC_NULL_INTEGER,Fsub,ierr);CHKERRQ(ierr)

  ! equation 1: A0 phi + n_0 (Q/<Q> - 1), Q = exp(del), del = (phi-<phi>)/Te
  if (sml_poisson_solver_type==1 .or. sml_poisson_solver_type==0) then
     ! form f0 <-- (A00 + betaM - beta M <x0>) x0
     call MatMult(    solver%Amat, Xsub(1),          Fsub(1), ierr) ! (A00 + betaM) x0
     call MatMultAdd( solver%Bmat, Xsub(2), Fsub(1), Fsub(1), ierr) ! beta M B <x0>
  else if (sml_poisson_solver_type==2) then
     call MatGetOwnershipRange(solver%Amat,low,high,ierr)
     nloc = high - low
     call VecGetArrayF90(solver%TeInv,TeInv_v,ierr)
     ! nonlinear part n_0(Q/<Q>-1); Q = exp((phi-<phi>)/Te)
     call VecDuplicate(Xsub(1),ave1,ierr)
     call VecDuplicate(Xsub(1),Q,ierr)
     !call MatMult(solver%Cmat, Xsub(1), Xsub(2), ierr)
     call MatMult(solver%BIdentMat,Xsub(2),ave1,ierr) ! <phi> for eq 1
     ! get e^(e(phi-<phi>)/Te)
     call VecGetArrayF90(ave1,ave1_v,ierr) ! realy phi_0
     call VecGetArrayReadF90(Xsub(1),phi_v,ierr) 
     call VecGetArrayF90(Q,del_v,ierr)
     do ii=1,nloc ! PETSc rows of B
        del_v(ii) = exp(TeInv_v(ii)*(phi_v(ii)-ave1_v(ii)))
        if (del_v(ii) /= del_v(ii))  print*,sml_mype,ii,TeInv_v(ii),phi_v(ii),ave1_v(ii)
     end do
     call VecRestoreArrayF90(ave1,ave1_v,ierr)
     call VecRestoreArrayReadF90(Xsub(1),phi_v,ierr) 
     call VecRestoreArrayF90(Q,del_v,ierr)
     ! <Q> into 'ave1', use fsub(2) as buffer
     call MatMult(solver%Cmat, Q, Fsub(2), ierr)
     ! add BC value 1d0 to Fsub(2)
     call MatGetOwnershipRange(solver%Cmat,low,high,ierr)
     nloc2 = high-low
     call VecGetArrayF90(solver%CConstBCVec,del_v,ierr)
     call VecGetArrayF90(Fsub(2),Q_v,ierr)
     do ii=1,nloc2
        Q_v(ii) = Q_v(ii) + del_v(ii)
     end do
     call VecRestoreArrayF90(solver%CConstBCVec,del_v,ierr)
     call VecRestoreArrayF90(Fsub(2),Q_v,ierr)
     ! scatter into phi
     call MatMult(solver%BIdentMat, Fsub(2), ave1, ierr)
     ! "Q" = Q/<Q> - 1
     call VecGetArrayF90(ave1,ave1_v,ierr)
     call VecGetArrayF90(Q,Q_v,ierr)
     do ii=1,nloc ! PETSc rows of B
        Q_v(ii) = Q_v(ii)/ave1_v(ii)-1.d0
     end do
     call VecRestoreArrayF90(ave1,ave1_v,ierr)
     call VecRestoreArrayF90(Q,Q_v,ierr)

     ! equation 1: A00 * x0 + n_0(Q/<Q> - 1)
     call MatMult(solver%FSAMass,Q,Fsub(1),ierr) 
     call MatMultAdd( solver%A0mat, Xsub(1), Fsub(1), Fsub(1), ierr)
     call VecRestoreArrayF90(solver%TeInv,TeInv_v,ierr)
     call VecDestroy(ave1,ierr)
     call VecDestroy(Q,ierr)
  else
     if (sml_mype==0) print*, '[',sml_mype,']FormFunction: Error unknown solver type ',sml_poisson_solver_type
     stop 'unknown solver type'
  end if

  ! equation 2 - linear: C phi - <phi> = 0
  call MatMult(   solver%Cmat, Xsub(1), Fsub(2), ierr)
  ! add -I <phi>
  call MatMultAdd(solver%Dmat, Xsub(2), Fsub(2), Fsub(2), ierr)
  call VecNorm(Fsub(2),NORM_INFINITY,t1,ierr)
  call VecNorm(Xsub(2),NORM_INFINITY,t2,ierr)
  if (sml_mype==0.and.t2/=0d0.and.t1/t2.gt.1.d-3) print*, 'FORMFUNCTION: Warning |C phi_0 - <phi>_0|/|<phi>_0|',t1/t2 

  call DMCompositeRestoreAccessArray(solver%da,X,itwo,PETSC_NULL_INTEGER,Xsub,ierr)
  call DMCompositeRestoreAccessArray(solver%da,F,itwo,PETSC_NULL_INTEGER,Fsub,ierr)
#if defined(PETSC_USE_LOG)
  call PetscLogEventEnd(FormFunctionEvent,ierr)
#endif
  return
end subroutine formfunction

#undef __FUNCT__
#define __FUNCT__ "write_vec"
subroutine write_vec(nn,phi,file_name,vec_name,ierr)
!c     write vector to file in matlab format
!c
!c     input:
!c       nn: global number of equations in system
!c       phi[nn]
!c       file_name: rhs[nn] -- not used
!c
!c     output:
!c       ierr: error code
!c
  use petscsnes
  implicit none
  integer,intent(in)::nn
  character (len=*),intent(in)::file_name
  character (len=*),intent(in)::vec_name
  real (kind=8),intent(in)::phi
  PetscErrorCode,intent(out)::ierr
  !
  PetscViewer::viewer
  PetscInt,parameter::ione=1
  Vec::vec
  PetscInt::n64

  n64 = nn
  call VecCreateSeqWithArray(PETSC_COMM_SELF,ione,n64,phi,vec,ierr);CHKERRQ(ierr)
  call PetscObjectSetName(vec,vec_name,ierr);CHKERRQ(ierr)
  call PetscViewerASCIIOpen(PETSC_COMM_SELF, file_name, viewer,ierr);CHKERRQ(ierr)
  call PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB,ierr);CHKERRQ(ierr)
  call VecView(vec,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  call VecDestroy(vec,ierr)
end subroutine write_vec
