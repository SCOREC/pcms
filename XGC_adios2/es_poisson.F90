#define XGC1 1
#include <petscversion.h>
#if PETSC_VERSION_LT(3,6,0)
#include <finclude/petscsnesdef.h>
#else
#include <petsc/finclude/petscsnesdef.h>
#endif
!
! init_poisson - high level init of physics and solver
!
#include "my_petscoptionsclearvalue.F90"

subroutine init_poisson(grid,psn,iflag_dummy)
  use psn_class
  use grid_class
  use sml_module
  use eq_module
  implicit none
  type(grid_type)::grid
  type(psn_type)::psn
  integer::iflag_dummy ! == 3 --obsolete
  integer::i
  real (kind=8)::psi,r,z
  PetscErrorCode::ierr

  if (iflag_dummy/=3) stop 'init_poisson: iflag\=3'

  do i=1,grid%nnode
     psi=grid%psi(i)
     r=grid%x(1,i)
     z=grid%x(2,i)
     psn%tempi_ev(i)=eq_ftn(psi,r,z,eq_tempi)  ! use real r & z value for safety
     psn%tempe_ev(i)=eq_ftn(psi,r,z,eq_tempe)
     psn%tite(i)=psn%tempi_ev(i)/psn%tempe_ev(i)
  enddo

  ! clear monitors on all planes except plane 0
  if (sml_plane_index.ne.0) then 
     call my_PetscOptionsClearValue('-fsa_snes_monitor',ierr)
     call my_PetscOptionsClearValue('-fsa_ksp_monitor',ierr)
     call my_PetscOptionsClearValue('-fsa_snes_view',ierr)
     call my_PetscOptionsClearValue('-fsa_ksp_view',ierr)
     call my_PetscOptionsClearValue('-fsa_snes_converged_reason',ierr)
     call my_PetscOptionsClearValue('-fsa_ksp_converged_reason',ierr)
     call my_PetscOptionsClearValue('-fsa_fieldsplit_phi_ksp_monitor',ierr)
     call my_PetscOptionsClearValue('-fsa_fieldsplit_phi_ksp_view',ierr)
     call my_PetscOptionsClearValue('-fsa_fieldsplit_phi_ksp_converged_reason',ierr)
     call my_PetscOptionsClearValue('-fsa_fieldsplit_lambda_upper_ksp_converged_reason',ierr)
     call my_PetscOptionsClearValue('-fsa_fieldsplit_lambda_ksp_converged_reason',ierr)
     call my_PetscOptionsClearValue('-fsa_fieldsplit_lambda_inner_ksp_converged_reason',ierr)
     call my_PetscOptionsClearValue('-fsa_fieldsplit_lambda_upper_ksp_monitor',ierr)
     call my_PetscOptionsClearValue('-fsa_fieldsplit_lambda_ksp_monitor',ierr)
     call my_PetscOptionsClearValue('-fsa_fieldsplit_lambda_inner_ksp_monitor',ierr)
     call my_PetscOptionsClearValue('-fsa_fieldsplit_lambda_upper_ksp_view',ierr)
     call my_PetscOptionsClearValue('-fsa_fieldsplit_lambda_ksp_view',ierr)
     call my_PetscOptionsClearValue('-fsa_fieldsplit_lambda_inner_ksp_view',ierr)
     call my_PetscOptionsClearValue('-help',ierr)
     call my_PetscOptionsClearValue('-info',ierr)
     call my_PetscOptionsClearValue('-vec_view',ierr)
     call my_PetscOptionsClearValue('-mat_view',ierr)
     call my_PetscOptionsClearValue('-s2_ksp_monitor',ierr)
     call my_PetscOptionsClearValue('-s2_ksp_converged_reason',ierr)
     call my_PetscOptionsClearValue('-ksp_monitor',ierr)
     call my_PetscOptionsClearValue('-ksp_converged_reason',ierr)
     call my_PetscOptionsClearValue('-s2_ksp_view',ierr)
     call my_PetscOptionsClearValue('-ksp_view',ierr)
     call my_PetscOptionsClearValue('-mat_partitioning_view',ierr)
     call my_PetscOptionsClearValue('-fsa_binva_ksp_converged_reason',ierr)
     call my_PetscOptionsClearValue('-fsa_binva_ksp_monitor',ierr)
     call my_PetscOptionsClearValue('-fsa_binva_ksp_view',ierr)
     call my_PetscOptionsClearValue('-fsa_binva_pc_gamg_verbose',ierr)
  endif

  ! init psn solver
  if (sml_iter_solver) then
     call psn_init_poisson_private(psn,grid,2) ! 2 RHS
  else 
     call psn_init_poisson_private(psn,grid,1) ! 1 RHS matrix
  endif

  if(sml_add_pot0>0) then
     !stop 'lets not allow sml_add_pot0 for now'
     call check_point('add potential offset')
     allocate( psn%add_pot0(grid%nnode))

     ! prepare additional 0 - potential
     if(sml_add_pot0==1) then
        call read_add_pot0(grid,psn)
     elseif(sml_add_pot0==2) then
!        call neo_pot0_simple(grid,psn)
        stop 'sml_add_pot0==2 is not working now'
     else
        psn%add_pot0=0D0
     endif
  endif

  call check_point('prepare gyro-averaging matrix')
#ifdef XGC1
  if(sml_plane_mype<grid%nrho*2) then
#else
  ! XGCa -->
  if(sml_plane_mype<grid%nrho) then
#endif
     call init_gyro_avg_mat(grid,psn)
  endif
end subroutine init_poisson

subroutine psn_init_poisson_private(psn,grid,nrhs)
  use petscsnes
  use psn_class
  use grid_class
  use sml_module
  use ptl_module,only:ptl_mass
  use eq_module
  implicit none
  type(grid_type)::grid
  type(psn_type)::psn
  integer::nrhs

 ! type(mat_type) :: matl, matr, matr2
  real (kind=8),allocatable :: alpha(:),beta(:),tiev(:), teev(:),den(:),b2(:), rhoi2(:)
  integer :: itr,nd(3),in_bd_node
  PetscErrorCode :: ierr
  real (kind=8) :: x_center(2),psi
  real (kind=8) :: factor,psi_in_poisson_bd,psi_out_poisson_bd
  real (kind=8), external :: gyro_radius2,psi_interpol,b_interpol
  real (kind=8), parameter :: offset=2D0

  allocate(alpha(grid%ntriangle),beta(grid%ntriangle),tiev(grid%ntriangle),teev(grid%ntriangle),den(grid%ntriangle),b2(grid%ntriangle),rhoi2(grid%ntriangle))

  ! setup n_0,b^2,Te at centroids of triangles
  do itr=1, grid%ntriangle
     nd=grid%nd(:,itr)
     x_center=(grid%x(:,nd(1))+grid%x(:,nd(2))+grid%x(:,nd(3)))/3D0
     psi=psi_interpol(x_center(1),x_center(2),0,0)
     tiev(itr)=eq_ftn(psi,x_center(1),x_center(2),eq_tempi)
     teev(itr)=eq_ftn(psi,x_center(1),x_center(2),eq_tempe)
     b2(itr)=b_interpol(x_center(1),x_center(2),0D0)**2
     den(itr)=eq_ftn(psi,x_center(1),x_center(2),eq_den)
  enddo
  rhoi2=ptl_mass(1)*tiev/sml_e_charge/b2

  ! 00 solver
  call check_point('Initialize 00 matrix solver')
  if(.NOT. sml_use_simple00 ) then
     ! **************************************
     ! initialize petsc solver
     ! **************************************
     if(.true.) then ! redundently solve 00 modes solver
        psn%solver00%comm = sml_plane_comm
        psn%solver00%mype=sml_plane_mype
        psn%solver00%totalpe=sml_plane_totalpe
     else
        psn%solver00%comm = petsc_comm_world
        psn%solver00%mype=sml_mype
        psn%solver00%totalpe=sml_totalpe
     endif
     psn%solver00%prefix=1
     psn%solver00%n_rhs_mat=nrhs ! set RHS mode (create_solver does not need this anymore!)

     ! create 00 solver
     call check_point('- 00 matrix solver memory create')
     call create_solver(psn%solver00,grid,psn%pbd0_2,ierr)
     CHKERRQ(ierr)
     ! LHS - alpha del^2 u
     alpha=ptl_mass(1)/sml_e_charge*den/b2 ! mn_0/eB^2
     psn%solver00%scale = 1d0/maxval(abs(alpha))
     alpha = alpha*psn%solver00%scale
     if (nrhs==2) then
        beta=psn%solver00%scale*den/teev
     else
        beta=0D0
     end if
     if(sml_mype==0) print *, 'psn_init_poisson_private: create Amat, first vertex in domain:',is_inside(1,psn%pbd0_2)
     !if (sml_sheath_mode/=0) stop 'psn_init_poisson_privates: sml_sheath_mode/=0'
     !call helm_matrix(psn%solver00%Amat,alpha,beta,grid,psn%pbd0_2,sml_sheath_mode/=0,psn%solver00%xgc_petsc,ierr)
     if (sml_sheath_mode/=0 .and. sml_mype==0) print *, 'psn_init_poisson_privates: sml_sheath_mode/=0  --> sml_sheath_mode==0 used instead for poisson solver'
     call helm_matrix(psn%solver00%Amat,alpha,beta,grid,psn%pbd0_2,.false.,psn%solver00%xgc_petsc,ierr)

     ! RHS mass matrix
     if(sml_mype==0) print *, 'psn_init_poisson_private: create rhs mat, first vertex in domain:',is_inside(1,psn%pbd0_2)
     alpha=0D0  ! no del term on RHS
     beta=psn%solver00%scale  ! diagonal (identity), gets ignore off proc from Amat
     call MatDuplicate(psn%solver00%Amat,MAT_DO_NOT_COPY_VALUES,psn%solver00%rhs_mat,ierr)
     call helm_matrix(psn%solver00%rhs_mat,alpha,beta,grid,psn%pbd0_2,.false.,psn%solver00%xgc_petsc,ierr)

     if (nrhs==2) then
        beta=psn%solver00%scale*den/teev  ! diagona term same as LHS, gets ignore off proc from Amat
        call MatDuplicate(psn%solver00%Amat,MAT_DO_NOT_COPY_VALUES,psn%solver00%rhs2_mat,ierr)
        call helm_matrix(psn%solver00%rhs2_mat,alpha,beta,grid,psn%pbd0_2,.false.,psn%solver00%xgc_petsc,ierr)        
        call init_simple00(grid,psn)  ! use for initial condition
     end if

     if (sml_poisson_solver_type/=0) then
        if (nrhs==2) stop 'nrhs==2 & sml_poisson_solver_type/=0'
#if (!PETSC_VERSION_LT(3,5,0))
        call create_2field_solver(grid,psn%pbd0_2,psn%solver00)
#else
        stop 'sml_poisson_solver_type/=0 & versio < 3.5'
#endif
     else
        call create_1field_solver(grid,psn%solver00)
     end if

  else ! simple 00
     call init_simple00(grid,psn)
  endif

#ifdef XGC1
  ! turb solver
  call check_point('Initialize turb solver')
  if(sml_turb_poisson) then
     ! ************************
     ! initialize petsc solver
     ! ************************
     psn%solverH%comm = sml_plane_comm
     psn%solverH%mype=sml_plane_mype
     psn%solverH%totalpe=sml_plane_totalpe
     psn%solverH%prefix=2
     psn%solverH%n_rhs_mat=1
     !
     call check_point('- turb solver memory create')
     call create_solver(psn%solverH,grid,psn%pbdH_2,ierr)
     CHKERRQ(ierr)
     if(.not. sml_use_pade) then
        ! LHS mn_0/eB^2
        alpha=ptl_mass(1)/sml_e_charge*den/b2 ! changed sign with new solver
     else
        ! LHS mn_0/eB^2 * ( 1 + Ti/Te)  - the second term comes from adiabatic response
        alpha=(1+tiev/teev)*ptl_mass(1)/sml_e_charge*den/b2 ! changed sign with new solver
     endif
     beta=den/teev
     call helm_matrix(psn%solverH%Amat,alpha,beta,grid,psn%pbdH_2,.false.,psn%solverH%xgc_petsc,ierr)
     call create_1field_solver(grid,psn%solverH)

     if(.not. sml_use_pade) then
        alpha=0D0  ! no del term on RHS
     else
        alpha=rhoi2
     endif
     beta=1D0   ! diagonal (identity)
     call MatDuplicate(psn%solverH%Amat,MAT_DO_NOT_COPY_VALUES,psn%solverH%rhs_mat,ierr)
     call helm_matrix(psn%solverH%rhs_mat,alpha,beta,grid,psn%pbdH_2,.false.,psn%solverH%xgc_petsc,ierr)
  endif
#endif

  deallocate(alpha,beta,tiev,teev,den,b2,rhoi2)

end subroutine psn_init_poisson_private
!
! f90 interfaces
!
module f90moduleinterfaces
  use petscpc
  use grid_class
  use psn_class
  use sml_module

  Interface PCSetApplicationContext
     Subroutine PCSetApplicationContext(pc,solver,ierr)
       use petscpc
       use grid_class
       use psn_class
       implicit none
       PC::pc
       type(xgc_solver)::solver
       PetscErrorCode::ierr
     End Subroutine PCSetApplicationContext
  End Interface
  
  Interface PCGetApplicationContext
     Subroutine PCGetApplicationContext(pc,solver,ierr)
       use petscpc
       use grid_class
       use psn_class
       implicit none
       PC::pc
       type(xgc_solver),pointer::solver
       PetscErrorCode::ierr
     End Subroutine PCGetApplicationContext
  End Interface
end module f90moduleinterfaces
!
! PCShermanMorrisonApply
!
subroutine PCShermanMorrisonApply(pc,xin,yout,ierr)
  use petscsys
  use grid_class
  use psn_class
  use sml_module
  use f90moduleinterfaces
  implicit none
  PetscErrorCode::ierr
  PC::pc
  Vec::xin,yout
  !
  type(xgc_solver),pointer::solver
  PetscInt,parameter::itwo=2
  Vec::Xsub(itwo),Ysub(itwo),vec
  KSP::subksp(2)
  PetscScalar,parameter::neg_one=-1d0
#if defined(PETSC_USE_LOG)
  call PetscLogStagePush(smstage,ierr)  
#endif
  call PCGetApplicationContext(pc,solver,ierr)
  call DMCompositeGetAccessArray(solver%da,xin,itwo,PETSC_NULL_INTEGER,Xsub,ierr)
  call DMCompositeGetAccessArray(solver%da,yout,itwo,PETSC_NULL_INTEGER,Ysub,ierr)
  call PCFieldSplitGetSubKSP(solver%pc,PETSC_NULL_INTEGER,subksp,ierr)
  call KSPSolve(subksp(1),Xsub(1),Ysub(1),ierr)
  call VecDuplicate(Xsub(2),vec,ierr)
  call MatMult(solver%Cmat,Ysub(1),vec,ierr) ! ignoring <> input
  call KSPSolve(subksp(2),vec,Ysub(2),ierr)
  call VecDestroy(vec,ierr)
  call MatMultAdd(solver%AinvB,Ysub(2),Ysub(1),Ysub(1),ierr)
  call VecScale(Ysub(2),neg_one,ierr)! canceled negative in SM, need for <phi> (diagnostics)
  call DMCompositeRestoreAccessArray(solver%da,xin,itwo,PETSC_NULL_INTEGER,Xsub,ierr)
  call DMCompositeRestoreAccessArray(solver%da,yout,itwo,PETSC_NULL_INTEGER,Ysub,ierr)
#if defined(PETSC_USE_LOG)
  call PetscLogStagePop(ierr)
#endif
  return
end subroutine PCShermanMorrisonApply
!
! helper function to compute nnz per row of my rows
!
subroutine getNNZ(grid,nloc,low,high,d_nnz,o_nnz,xgc_petsc,nglobal,ierr)
  use petscsnes
  use grid_class
  use sml_module
  implicit none
  type(grid_type),intent(in)::grid
  PetscInt::nloc,d_nnz(nloc),o_nnz(nloc),low,high,nglobal
  PetscInt,intent(in)::xgc_petsc(grid%nnode)
  PetscErrorCode::ierr
  !
  PetscInt::ie,ind(3),itmp,idi(1),ii,jtmp,idj(1),jj

  d_nnz = 4 ! 4/2 for edges on tiny problems
  o_nnz = 4
  do ie=1,grid%ntriangle
     ind(:) = grid%nd(:,ie)
     do itmp=1,3
        ii = xgc_petsc(ind(itmp)) + 1 ! one based, could == 0
        if (ii.gt.low.and.ii.le.high) then ! I (ii-low) touch this triangle
           do jtmp=1,3
              jj = xgc_petsc(ind(jtmp)) + 1 ! one based, could == 0
              if (ii/=jj) then
                 if (jj.gt.low.and.jj.le.high) then ! local
                    d_nnz(ii-low) = d_nnz(ii-low) + 1
                 else
                    o_nnz(ii-low) = o_nnz(ii-low) + 1 
                 end if
              end if
           end do
        end if
     end do
  end do
  
  o_nnz = o_nnz/2 ! two triangles hit each edge except outside edge
  d_nnz = d_nnz/2
  do ii=1,nloc
     if (d_nnz(ii).gt.nloc) d_nnz(ii) = nloc
     if (o_nnz(ii).gt.nglobal-nloc) o_nnz(ii) = nglobal-nloc
  end do
end subroutine getNNZ
!
! create_solver - create Amat and its data layout, and this%petscloc_xgc, this%xgc_petsc
!
subroutine create_solver(this,grid,bc,ierr)
  use petscsnes
  use petscsys
  use grid_class
  use psn_class
  use sml_module
  implicit none
  type(xgc_solver)::this
  type(grid_type),intent(in)::grid
  type(boundary2_type),intent(in)::bc
  PetscErrorCode::ierr
  !
  integer::ii
  PetscInt::nnodes64,proc,iloc,nreal
  Mat::AA
  real (kind=8),allocatable::alpha(:),beta(:)  
  MatPartitioning::part;
  IS::is,is2;
  PetscInt,pointer::xgc_proc(:)
  PetscInt,allocatable::petsc_xgc(:),d_nnz(:),o_nnz(:)
  PetscInt::npetscloc,low,high,proc_eq(0:sml_plane_totalpe)
  PetscInt,parameter::ione=1
  Vec::vec
  !
  nnodes64 = grid%nnode
  this%Amat = 0
  this%rhs_mat=0
  this%rhs2_mat=0
  ! LHS --------------
  if (allocated(this%xgc_petsc)) stop 'associated(this%xgc_petsc)'
  allocate(this%xgc_petsc(grid%nnode))  

  if(sml_mype==0) print *, 'create_solver: first vertex in domain:',is_inside(1,bc)
  
  ! get: nreal,npetscloc,petscloc_xgc,low,this%xgc_petsc
  if (.true.) then
     ! metis partitioning
     nreal=0 ! count locals with BC move to 0, set this%xgc_petsc
     do ii=1,grid%nnode
        if(is_inside(ii,bc)) then ! put BCs on proc 0
           this%xgc_petsc(ii) = nreal ! default zero based map for assembly   
           nreal = nreal + 1
        else
           this%xgc_petsc(ii) = -1 ! petsc will skip & and fem does not look anyway
        end if
     end do
     allocate(alpha(grid%ntriangle),beta(grid%ntriangle))
     ! create dummy matrix for partitioning. could be done natively
     if(sml_mype==0)write(*,1)'CREATE_SOLVER: make partitioning with ',nreal,'/',grid%nnode,' real vertices'
1       format(A,I7,A,I7,A)
     call MatCreate(this%comm,AA,ierr)
     call MatSetSizes(AA,PETSC_DECIDE,PETSC_DECIDE,nreal,nreal,ierr)
     call MatSetType(AA,MATAIJ,ierr)
     call MatSetup(AA,ierr)
     call MatGetOwnershipRange(AA,low,high,ierr) 
     npetscloc = high-low
     allocate(d_nnz(npetscloc),o_nnz(npetscloc))
     call getNNZ(grid,npetscloc,low,high,d_nnz,o_nnz,this%xgc_petsc,nreal,ierr)
     call MatSeqAIJSetPreallocation(AA,PETSC_NULL_INTEGER,d_nnz,ierr)  
     call MatMPIAIJSetPreallocation(AA,PETSC_NULL_INTEGER,d_nnz,PETSC_NULL_INTEGER,o_nnz,ierr)
     deallocate(d_nnz,o_nnz)
     call MatSetOption(AA,MAT_IGNORE_OFF_PROC_ENTRIES,PETSC_TRUE,ierr)
     call MatSetOption(AA,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,ierr)
     call MatSetup(AA,ierr)
     beta=1d0
     alpha=1d0
     call helm_matrix(AA,alpha,beta,grid,bc,.false.,this%xgc_petsc,ierr)  
     ! partition
     call MatPartitioningCreate(this%comm,part,ierr)
     call MatPartitioningSetAdjacency(part,AA,ierr)
     !call MatPartitioningSetType(part,MATPARTITIONING_PARMETIS,ierr) ! this does not get set???
     call MatPartitioningSetFromOptions(part,ierr)
     
     call MatPartitioningApply(part,is,ierr)
     ! clean up partitioning
     call MatPartitioningDestroy(part,ierr)
     call MatDestroy(AA,ierr)
     call ISGetLocalSize(is,npetscloc,ierr)
     if (npetscloc.ne.high-low) stop 'npetscloc.ne.high-low'
     call ISAllGather(is,is2,ierr) ! get global colors
     call ISDestroy(is,ierr)
     !if (sml_mype==0) call ISView(is2,PETSC_VIEWER_STDOUT_SELF,ierr)
     ! count proc sizes
     call ISGetIndicesF90(is2,xgc_proc,ierr)
     proc_eq = 0
     nreal = 0
     do ii=1,grid%nnode
        if(is_inside(ii,bc)) then
           nreal = nreal + 1
           proc = xgc_proc(nreal) + 1
           proc_eq(proc) = proc_eq(proc) + 1 ! increment, start at (1)
        end if
     end do
     !if (sml_mype==0) write(*,*) 'COUNTS proc_eq=',(proc_eq(iloc),iloc=1,sml_plane_totalpe)
     npetscloc = proc_eq(sml_plane_mype+1) ! my new local size
     ! scan
     proc_eq(0) = 0
     do ii=1,sml_plane_totalpe
        proc_eq(ii) = proc_eq(ii) + proc_eq(ii-1)
     end do
     !if (sml_mype==sml_totalpe-1)write(*,*)'STARTS proc_eq=',(proc_eq(iloc),iloc=0,sml_plane_totalpe-1)
     ! set global PETSc equation numbers
     allocate(this%petscloc_xgc(npetscloc))
     allocate(petsc_xgc(nreal))
     nreal = 0
     iloc = 0
     do ii=1,grid%nnode
        if(is_inside(ii,bc)) then
           nreal = nreal + 1
           proc = xgc_proc(nreal)
           this%xgc_petsc(ii) = proc_eq(proc)
           petsc_xgc(proc_eq(proc)+1) = ii-1
           if (proc==sml_plane_mype) then
              iloc = iloc + 1
              this%petscloc_xgc(iloc) = ii - 1
           end if
           proc_eq(proc) = proc_eq(proc) + 1 ! increment
        else 
           this%xgc_petsc(ii) = -1
        end if
     end do
     if (iloc.ne.npetscloc) stop 'iloc.ne.npetscloc'
     !if (sml_mype==sml_totalpe-1) write(*,*)'END ids proc_eq=',(proc_eq(iloc),iloc=0,sml_plane_totalpe-1)
     call ISRestoreIndicesF90(is2,xgc_proc,ierr)
     call ISDestroy(is2,ierr)
     deallocate(alpha,beta)
     ! create new matrix
     call MatCreate(this%comm,AA,ierr)
     call MatSetSizes(AA,npetscloc,npetscloc,nreal,nreal,ierr)
     call MatSetType(AA,MATAIJ,ierr)
     call MatSetup(AA,ierr)
     call MatGetOwnershipRange(AA,low,high,ierr)
     allocate(d_nnz(npetscloc),o_nnz(npetscloc))
     call getNNZ(grid,npetscloc,low,high,d_nnz,o_nnz,this%xgc_petsc,nreal,ierr)
     call MatSeqAIJSetPreallocation(AA,PETSC_NULL_INTEGER,d_nnz,ierr)  
     call MatMPIAIJSetPreallocation(AA,PETSC_NULL_INTEGER,d_nnz,PETSC_NULL_INTEGER,o_nnz,ierr)
     deallocate(d_nnz,o_nnz)
     call MatSetOption(AA,MAT_IGNORE_OFF_PROC_ENTRIES,PETSC_TRUE,ierr)
     call MatSetOption(AA,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,ierr)
     call MatSetup(AA,ierr)
     this%Amat = AA
     if (npetscloc.ne.high-low) stop 'npetscloc.ne.high-low'
     !if (sml_plane_index==0) write(*,'(A,I4,A,I4,A,I6,A,I8)') 'mype=',mype,', sml_mype=',sml_mype,', npetscloc=',npetscloc,', nreal=',nreal
  else
     if(sml_mype==0) print *, '********************** create_solver: warning simple partitionin gwith BCs not working, perhaps'
     ! simple chop partitioning
     allocate(petsc_xgc(grid%nnode))
     do ii=1,grid%nnode     
        this%xgc_petsc(ii) = ii - 1 ! default zero based map for assembly
        petsc_xgc(ii) = ii - 1 
     end do
     call MatCreate(this%comm,AA,ierr)
     call MatSetSizes(AA,PETSC_DECIDE,PETSC_DECIDE,nnodes64,nnodes64,ierr)
     call MatSetType(AA,MATAIJ,ierr)
     !call MatSeqAIJSetPreallocation(AA,maxn164,PETSC_NULL_INTEGER,ierr)  
     !call MatMPIAIJSetPreallocation(AA,maxn164,PETSC_NULL_INTEGER,maxn164/2,PETSC_NULL_INTEGER,ierr)
     stop 'create_solver: not supported'
     call MatSetOption(AA,MAT_IGNORE_OFF_PROC_ENTRIES,PETSC_TRUE,ierr)
     call MatSetOption(AA,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,ierr)
     call MatSetup(AA,ierr)
     this%Amat = AA
     call MatGetOwnershipRange(AA,low,high,ierr)
     npetscloc = high-low
     ! petscloc_xgc
     allocate(this%petscloc_xgc(npetscloc))
     iloc = 0
     do ii=low,high-1
        iloc = iloc + 1
        this%petscloc_xgc(iloc) = ii
     end do
     nreal = nnodes64
  end if
  ! vecs
  call VecCreateMPI(this%comm,npetscloc,nreal,this%xVec,ierr)
  call VecSetFromOptions(this%xVec,ierr)
  call VecDuplicate(this%xVec,this%bVec,ierr)
  ! set scatters
  call VecCreateSeq(PETSC_COMM_SELF,nnodes64,vec,ierr) ! big dumb PETSc vector
  ! forward scatter into PETSc
  call ISCreateGeneral(PETSC_COMM_SELF,npetscloc,this%petscloc_xgc,PETSC_COPY_VALUES,is,ierr) ! is: x_petscloc(n)
  call ISCreateStride(PETSC_COMM_SELF,npetscloc,low,ione,is2,ierr)
  call VecScatterCreate(vec,is,this%xVec,is2,this%to_petsc,ierr) ! forward scatter into PETSc
  CHKERRQ(ierr)
  call ISDestroy(is2,ierr)
  call ISDestroy(is,ierr)
  ! reverse scatter object
  call ISCreateGeneral(PETSC_COMM_SELF,nreal,petsc_xgc,PETSC_COPY_VALUES,is,ierr) ! is: petsc_xgc
  CHKERRQ(ierr)
  deallocate(petsc_xgc)
  !if (sml_mype==0) call ISView(is,PETSC_VIEWER_STDOUT_SELF,ierr)
  call VecScatterCreate(this%xVec,PETSC_NULL_OBJECT,vec,is,this%from_petsc,ierr) ! reverse scatter object
  call ISDestroy(is,ierr)
  ! cleanup
  call VecDestroy(vec,ierr)
  ! zero out
  this%rVec2 = 0
  this%xVec2 = 0
  this%bVec2 = 0
  this%da = 0
  this%KKTmat = 0
  this%Bmat = 0
  this%CMat = 0
  this%CConstBCVec = 0
  this%Dmat = 0
  this%A0Mat = 0
  this%BIdentMat = 0
  !this%BMassMat = 0
  this%FSAMass = 0
  this%TeInv = 0
  this%snes = 0
  this%scale = 1d0
  this%schur = 0
  this%AinvB = 0
  this%ksp = 0
  this%pc = 0
end subroutine create_solver
!
! create 1 field solver with A00 built (Amat), just setoperator & setup like 2 field
!
subroutine create_1field_solver(grid,solver)
  use petscsnes
  use grid_class
  use psn_class
  use xgc_solver_module
  use sml_module,only:sml_mype
  implicit none
  type(grid_type) :: grid
  type(xgc_solver) :: solver
  PetscErrorCode :: ierr

  ! create KSP
  call KSPCreate(solver%comm,solver%ksp,ierr)
  CHKERRQ(ierr)
  if (solver%prefix.eq.2) then ! turb
     call KSPSetOptionsPrefix(solver%ksp,'s2_',ierr)
  else ! 1 field 00 solver, no prefix
  end if
  call KSPSetFromOptions(solver%ksp,ierr)
#ifndef PETSC_ED
#if (!PETSC_VERSION_LT(3,5,0))
  call KSPSetOperators(solver%ksp,solver%Amat,solver%Amat,ierr)
#else
  call KSPSetOperators(solver%ksp,solver%Amat,solver%Amat,SAME_NONZERO_PATTERN,ierr)
#endif
#else
  call KSPSetOperators(solver%ksp,solver%Amat,solver%Amat,SAME_NONZERO_PATTERN,ierr)
#endif
  CHKERRQ(ierr)
  ! setup solver now that matrix is complete
  call KSPSetUp(solver%ksp,ierr)
  if(ierr .ne. 0) then
     print *, 'Error in setup of petsc create_1field_solver solver :',ierr
     call MPI_ABORT(solver%comm,1,ierr)
  endif
end subroutine create_1field_solver
!
! create 2 field solver with A00 built (Amat)
!
subroutine create_2field_solver(grid,bd,solver)
  use petscsnes
  use grid_class
  use psn_class
  use xgc_solver_module
  use sml_module,only:sml_mype,sml_poisson_solver_type,sml_intpl_mype,sml_plane_mype
  use eq_module
  use f90moduleinterfaces
  implicit none
  type(grid_type)::grid
  type(xgc_solver)::solver
  type(boundary2_type)::bd
  !
  PetscErrorCode :: ierr
  real (kind=8),parameter::epsilon=0D0
  integer::comm
  PetscInt::nflux,low,high,N1loc,N2loc,M,N
  PetscInt::k,j,i,ncols,nloc,ahigh,alow,anloc,ii,onnzamax
  PetscInt,allocatable::lgids(:),cols(:)
  PetscScalar,allocatable::vals(:)
  PetscScalar::dot
  Mat::matArray(4)
  KSP::innerksp,subksp(2)
  PC::spc
  Mat::n0Mass_Te
  Vec::x1Vec,x2Vec,x2Vecloc,x1Vecloc,v1,v2,Fsub(2)
  PetscViewer::viewer
  PetscScalar,pointer::parr(:)
  PetscInt,parameter::ione=1
  PetscInt,parameter::itwo=2
  PetscScalar,parameter::neg_one=-1.d0,pos_one=1d0
  real (kind=8)::x_center(2),teev,den,psi
  integer::itr,nd(3),xgc_idx
  real (kind=8),external::psi_interpol
  external FormJacobian,FormFunction,PCShermanMorrisonSetUp,PCShermanMorrisonApply
#if defined(PETSC_USE_LOG)
  PetscLogStage,save::stage
#endif
  real (kind=8)::alpha(grid%ntriangle),beta1(grid%ntriangle),beta2(grid%ntriangle)

  ! set helper vars
  comm = solver%comm
  call MatGetOwnershipRange(solver%Amat,low,high,ierr)
  CHKERRQ(ierr)
  N1loc = high - low

  ! create other matrices and vectors, DM ...
  if (sml_poisson_solver_type/=0 .and. solver%prefix==1) then
     nflux = grid%cnv_2d_00%m
     if (nflux.le.0) stop 'nflux.lt.0'
  else
     ! nflux = 0
     stop '!(sml_poisson_solver_type/=0 .and. solver%prefix==1)'
  endif
  if(sml_mype==0)  write(*,1) &
       ' FSA SOLVER: ',nflux,' flux surfaces, (doit=',sml_poisson_solver_type/=0 .and. solver%prefix==1,') make A (A_00) mat.'
1 format(A,I4,A,L,A)

#if (!PETSC_VERSION_LT(3,5,0))
#if defined(PETSC_USE_LOG)
  call PetscLogStageRegister('FSA init',stage,ierr)
  !if (sml_poisson_solver_type==1) then
  call PetscLogStageRegister('FSA Sherman-Morrison',smstage,ierr)
  !end if
  call PetscLogEventRegister('FSA FormFunction',0,FormFunctionEvent,ierr)
  call PetscLogStagePush(stage,ierr)
#endif

  ! FSAMass: mass matrix scaled with n_0 & scale
  do itr=1,grid%ntriangle
     ! setup basic values
     nd=grid%nd(:,itr)
     x_center=( grid%x(:,nd(1)) + grid%x(:,nd(2)) + grid%x(:,nd(3)) )/3D0
     psi=psi_interpol(x_center(1),x_center(2),0,0)
     teev=eq_ftn(psi,x_center(1),x_center(2),eq_tempe)
     den=eq_ftn(psi,x_center(1),x_center(2),eq_den)
     beta1(itr) =solver%scale*den
     beta2(itr) =solver%scale*den/teev
  enddo
  ! gets ignore off proc from Amat
  call MatDuplicate(solver%Amat,MAT_DO_NOT_COPY_VALUES,solver%FSAMass,ierr)
  alpha=0D0
  call helm_matrix(solver%FSAMass,alpha,beta1,grid,bd,.false.,solver%xgc_petsc,ierr)

  ! form n0Mass/Te - gets ignore off proc from Amat
  call MatDuplicate(solver%Amat,MAT_DO_NOT_COPY_VALUES,n0Mass_Te,ierr)
  !call MatDiagonalScale(n0Mass_Te,PETSC_NULL_OBJECT,solver%TeInv,ierr)
  call helm_matrix(n0Mass_Te,alpha,beta2,grid,bd,.false.,solver%xgc_petsc,ierr)

  ! form vector of 1/Te
  call VecCreate(comm,solver%TeInv,ierr)
  call VecSetSizes(solver%TeInv,N1loc,PETSC_DECIDE,ierr)
  call VecSetFromOptions(solver%TeInv,ierr) 
  do ii=low,high-1
     xgc_idx = solver%petscloc_xgc(ii-low+1) + 1 ! looking into global data structure
     if( .not.is_inside(xgc_idx,bd) ) stop 'create_2field_solver: not in'
     x_center = grid%x(:,xgc_idx)
     psi=psi_interpol(x_center(1),x_center(2),0,0)
     teev=1d0/eq_ftn(psi,x_center(1),x_center(2),eq_tempe)
     call VecSetValue(solver%TeInv,ii,teev,INSERT_VALUES,ierr)
  end do
  call VecAssemblyBegin(solver%TeInv,ierr)
  call VecAssemblyEnd(  solver%TeInv,ierr)
  
  ! store raw A_0 operator for formFunction, gets ignore off proc from Amat
  call MatDuplicate(solver%Amat,MAT_COPY_VALUES,solver%A0Mat,ierr)

  ! add linearization term to Amat
  !call add_11mat(grid,solver%Amat,bd,solver%scale,solver%xgc_petsc)
  call MatAXPY(solver%Amat,pos_one,n0Mass_Te,SAME_NONZERO_PATTERN,ierr)

  ! BIdentMat - nice to get sizes in, I know the nnz of this
  call MatCreate(comm,solver%BIdentMat,ierr)
  N2loc = nflux
  if (sml_plane_mype/=0) N2loc=0
  if (.true.) then
     call MatSetSizes(solver%BIdentMat,N1loc,N2loc,PETSC_DECIDE,nflux,ierr) ! all on proc 0
  else
     call MatSetSizes(solver%BIdentMat,N1loc,PETSC_DECIDE,PETSC_DECIDE,nflux,ierr) 
  end if
  CHKERRQ(ierr)
  call MatSetType(solver%BIdentMat,MATAIJ,ierr)
  CHKERRQ(ierr)
  call MatSeqAIJSetPreallocation(solver%BIdentMat,itwo,PETSC_NULL_INTEGER,ierr)
  CHKERRQ(ierr)
  call MatMPIAIJSetPreallocation(solver%BIdentMat,itwo,PETSC_NULL_INTEGER,itwo,PETSC_NULL_INTEGER,ierr)
  CHKERRQ(ierr)
  call MatSetOption(solver%BIdentMat,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE,ierr)
  call MatSetup(solver%BIdentMat,ierr)
  call MatGetOwnershipRangeColumn(solver%BIdentMat,low,high,ierr)
  N2loc = high-low

  ! Cmat -- short fat, need to preallocate later
  call MatCreate(comm,solver%Cmat,ierr)
  CHKERRQ(ierr)
  call MatSetSizes(solver%Cmat,N2loc,N1loc,nflux,PETSC_DECIDE,ierr)
  CHKERRQ(ierr)
  call MatSetType(solver%Cmat,MATAIJ,ierr)
  ! helper to apply constant non-homo BC
  call VecCreate(comm,solver%CConstBCVec,ierr)
  call VecSetSizes(solver%CConstBCVec,N2loc,PETSC_DECIDE,ierr)
  call VecSetFromOptions(solver%CConstBCVec,ierr) 
  ! Dmat -- -I
  call MatCreate(comm,solver%Dmat,ierr)
  CHKERRQ(ierr)
  call MatSetSizes(solver%Dmat,N2loc,N2loc,nflux,nflux,ierr)
  CHKERRQ(ierr)
  call MatSetType(solver%Dmat,MATAIJ,ierr)
  CHKERRQ(ierr)
  call MatSeqAIJSetPreallocation(solver%Dmat,ione,PETSC_NULL_INTEGER,ierr)
  CHKERRQ(ierr)
  call MatMPIAIJSetPreallocation(solver%Dmat,ione,PETSC_NULL_INTEGER,ione,PETSC_NULL_INTEGER,ierr)
  CHKERRQ(ierr)
  call MatSetOption(solver%Dmat,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE,ierr)
  call MatSetup(solver%Dmat,ierr)

  ! pbd0_2 is hard coded since nflux>0 only with solver00
  call make_12mat(grid,solver%Bmat,solver%BIdentMat,n0Mass_Te,solver,bd,solver%scale,solver%xgc_petsc)
  call MatDestroy(n0Mass_Te,ierr)
  call make_21mat(grid,solver%Cmat,solver%Bmat,solver%CConstBCVec,bd,solver%xgc_petsc)
  
  ! Identity
  do j=low,high-1
     call MatSetValues(solver%Dmat,ione,j,ione,j,neg_one,INSERT_VALUES,ierr)
     CHKERRQ(ierr)
  enddo
  call MatAssemblyBegin(solver%Dmat,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(solver%Dmat,MAT_FINAL_ASSEMBLY,ierr)

  if (.false.) then ! debug
     if (sml_intpl_mype==0) then 
        call  PetscViewerASCIIOpen(comm, 'Amat.m', viewer,ierr)
        call  PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB,ierr)
        call  MatView(solver%Amat,viewer,ierr)
        call  PetscViewerDestroy(viewer,ierr)
        call  PetscViewerASCIIOpen(comm, 'Bmat.m', viewer,ierr)
        call  PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB,ierr)
        call  MatView(solver%Bmat,viewer,ierr)
        call  PetscViewerDestroy(viewer,ierr)
        call  PetscViewerASCIIOpen(comm, 'BIdentmat.m', viewer,ierr)
        call  PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB,ierr)
        call  MatView(solver%BIdentmat,viewer,ierr)
        call  PetscViewerDestroy(viewer,ierr)
        call  PetscViewerASCIIOpen(comm, 'Cmat.m', viewer,ierr)
        call  PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB,ierr)
        call  MatView(solver%Cmat,viewer,ierr)
        call  PetscViewerDestroy(viewer,ierr)
        call  PetscViewerASCIIOpen(comm, 'Dmat.m', viewer,ierr)
        call  PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB,ierr)
        call  MatView(solver%Dmat,viewer,ierr)
        call  PetscViewerDestroy(viewer,ierr)
     endif
     call MPI_Barrier(PETSC_COMM_WORLD,ierr)
  end if

  ! phi DM
  call VecCreate(comm,x1Vec,ierr)
  call VecSetSizes(x1Vec,N1loc,PETSC_DECIDE,ierr)
  call VecSetFromOptions(x1Vec,ierr) 

  call DMShellCreate(comm,solver%daphi,ierr)
  call DMShellSetGlobalVector(solver%daphi,x1Vec,ierr)
  call DMShellSetMatrix(solver%daphi,solver%Amat,ierr)

  call VecCreate(PETSC_COMM_SELF,x1Vecloc,ierr)
  call VecSetSizes(x1Vecloc,N1loc,N1loc,ierr)
  call VecSetFromOptions(x1Vecloc,ierr) 
  call DMShellSetLocalVector(solver%daphi,x1Vecloc,ierr)
  call VecDestroy(x1Vec,ierr)
  call VecDestroy(x1Vecloc,ierr)

  ! lambda DM
  call VecCreate(comm,x2Vec,ierr)
  call VecSetSizes(x2Vec,N2loc,nflux,ierr)
  call VecSetFromOptions(x2Vec,ierr) 

  call DMShellCreate(comm,solver%dalam,ierr)
  call DMShellSetGlobalVector(solver%dalam,x2Vec,ierr)
  call DMShellSetMatrix(solver%dalam,solver%Dmat,ierr)

  call VecCreate(PETSC_COMM_SELF,x2Vecloc,ierr)
  call VecSetSizes(x2Vecloc,N2loc,N2loc,ierr)
  call VecSetFromOptions(x2Vecloc,ierr) 
  call DMShellSetLocalVector(solver%dalam,x2Vecloc,ierr)
  call VecDestroy(x2Vec,ierr)
  call VecDestroy(x2Vecloc,ierr)

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !  Create field split DM
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  call DMCompositeCreate(comm,solver%da,ierr)
  CHKERRQ(ierr)
  call DMSetOptionsPrefix(solver%daphi,'phi_',ierr)
  call DMCompositeAddDM(solver%da,solver%daphi,ierr)
  CHKERRQ(ierr)
  call DMSetOptionsPrefix(solver%dalam,'lambda_',ierr)
  call DMCompositeAddDM(solver%da,solver%dalam,ierr)
  CHKERRQ(ierr)

  matArray(1) = solver%Amat
  matArray(2) = solver%Bmat
  matArray(3) = solver%Cmat
  matArray(4) = solver%Dmat
  call MatCreateNest(comm,itwo,PETSC_NULL_OBJECT,itwo,PETSC_NULL_OBJECT,matArray,solver%KKTmat,ierr)
  call MatSetup(solver%KKTmat,ierr)
  
  !  Extract global and local vectors from DM; then duplicate for remaining
  !     vectors that are the same types
#if (PETSC_VERSION_LT(3,6,0))
  call MatGetVecs(solver%KKTmat,solver%xVec2,solver%bVec2,ierr)
#else
  call MatCreateVecs(solver%KKTmat,solver%xVec2,solver%bVec2,ierr)
#endif
  call VecDuplicate(solver%bVec2,solver%rVec2,ierr)
  call DMSetup(solver%da,ierr)

  ! Add Schur complement preconditioner to solver, make solver objects (mats)
  if(sml_mype==0) print *, 'create_2field_solver: setting up Schur complement'
  call MatGetSize(solver%Bmat,M,N,ierr)
  call MatGetOwnershipRange(solver%Cmat,low,high,ierr)
  nloc = high-low
  onnzamax = N-nloc
  call MatGetOwnershipRange(solver%Amat,alow,ahigh,ierr)
  anloc = ahigh - alow

  call VecDuplicate(solver%xVec,v1,ierr)
  call VecDuplicate(solver%xVec,v2,ierr)
  call KSPCreate(solver%comm,innerksp,ierr)
  call KSPSetOptionsPrefix(innerksp,'fsa_binva_',ierr)
  call KSPSetFromOptions(innerksp,ierr)
  call KSPSetOperators(innerksp,  solver%Amat,  solver%Amat, ierr )
  call KSPSetUp(innerksp,ierr)
  !
  call MatCreate(solver%comm,solver%schur,ierr);
  call MatSetSizes(solver%schur,nloc,nloc,PETSC_DECIDE,PETSC_DECIDE,ierr)
  call MatSetType(solver%schur,MATAIJ,ierr)
  call MatSetOptionsPrefix(solver%schur,'schur_',ierr)
  call MatSetFromOptions(solver%schur,ierr)
  call MatSeqAIJSetPreallocation(solver%schur,nloc,PETSC_NULL_INTEGER,ierr)
  call MatMPIAIJSetPreallocation(solver%schur,nloc,PETSC_NULL_INTEGER,onnzamax,PETSC_NULL_INTEGER,ierr) ! dense
  call MatSetup(solver%schur,ierr)
  
  call MatCreate(solver%comm,solver%AinvB,ierr)
  call MatSetSizes(solver%AinvB,anloc,nloc,PETSC_DECIDE,PETSC_DECIDE,ierr)
  call MatSetType(solver%AinvB,MATAIJ,ierr)
  call MatSetOptionsPrefix(solver%AinvB,'ainvb_',ierr)
  call MatSetFromOptions(solver%AinvB,ierr)
  call MatSeqAIJSetPreallocation(solver%AinvB,nloc,PETSC_NULL_INTEGER,ierr)
  call MatMPIAIJSetPreallocation(solver%AinvB,nloc,PETSC_NULL_INTEGER,onnzamax,PETSC_NULL_INTEGER,ierr)
  call MatSetup(solver%AinvB,ierr)
  
  allocate(lgids(anloc),cols(10000),vals(10000))
  do i=1,anloc
     lgids(i) = alow+i-1 ! zero based
  end do
  do k=0,N-1 ! zero based
     call MatGetColumnVector(solver%Bmat,v1,k,ierr)
     call KSPSolve(innerksp,v1,v2,ierr)
     !if (sml_poisson_solver_type==1) then
     call VecGetArrayF90(v2,parr,ierr)
     call MatSetValues(solver%AinvB,anloc,lgids,ione,k,parr,INSERT_VALUES,ierr) 
     CHKERRQ(ierr)
     call VecRestoreArrayF90(v2,parr,ierr)
     !end if
     do j=0,N-1 ! zero based
        call VecZeroEntries(v1,ierr)
        if (j.ge.low .and.j.lt.high) then
           call MatGetRow(    solver%Cmat,j,ncols,PETSC_NULL_INTEGER,PETSC_NULL_SCALAR,ierr)
           if (ncols.gt.size(cols)) then
              deallocate(cols,vals)
              allocate(cols(ncols),vals(ncols))
              if(sml_mype==0) print *, ' FSA: realloc buffers to ',ncols
           end if
           call MatRestoreRow(solver%Cmat,j,ncols,PETSC_NULL_INTEGER,PETSC_NULL_SCALAR,ierr)
           call MatGetRow(    solver%Cmat,j,ncols,cols,vals,ierr)
           if (ncols.gt.size(cols)) stop 'ncols.gt.MAX_CMAT_COLS (2)'
           call VecSetValues(v1,ncols,cols,vals,INSERT_VALUES,ierr)
           call MatRestoreRow(solver%Cmat,j,ncols,cols,vals,ierr)
        end if
        call VecAssemblyBegin(v1,ierr)
        call VecAssemblyEnd(v1,ierr)
        call VecTDot(v1,v2,dot,ierr)
        if (j.ge.low .and.j.lt.high) then
           dot = -dot
           ! D - C (A-1) B
           call MatSetValues(solver%schur,ione,j,ione,k,dot,ADD_VALUES,ierr) 
           CHKERRQ(ierr)
        end if
     end do
     if (k.ge.low.and.k.lt.high) then
        dot = -1.d0 ! should use D
        call MatSetValues(solver%schur,ione,k,ione,k,dot,ADD_VALUES,ierr) 
        CHKERRQ(ierr)
     end if
  end do
  deallocate(lgids,cols,vals)
  call MatAssemblyBegin(solver%schur,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(solver%schur,MAT_FINAL_ASSEMBLY,ierr)
  !if (sml_poisson_solver_type==1) then
  call MatAssemblyBegin(solver%AinvB,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(solver%AinvB,MAT_FINAL_ASSEMBLY,ierr)
  !end if
  call VecDestroy(v1,ierr)
  call VecDestroy(v2,ierr)
  call KSPDestroy(innerksp,ierr)
  
  !call SNESGetKSP(solver%snes,innerksp,ierr)
  !call KSPGetPC(innerksp,spc,ierr)
  !call PCSetUp(spc,ierr)
  !call PCFieldSplitGetSubKSP(spc,i,subksp,ierr)     
  !call KSPGetOperators(subksp(2),Amat,Pmat,ierr)
  !call MatComputeExplicitOperator(Amat,solver%schur,ierr)
  
  ! create SNES
  call SNESCreate(comm,solver%snes,ierr)
  !  Set runtime options (e.g., -snes_monitor -snes_rtol <rtol> -ksp_type <type>)
  call SNESSetOptionsPrefix(solver%snes, 'fsa_', ierr)
  call SNESSetFromOptions(solver%snes,ierr)
  !  Set function evaluation routine and vector
  call SNESSetFunction(solver%snes,solver%rVec2,FormFunction,solver,ierr)
  call SNESSetJacobian(solver%snes,solver%KKTmat,solver%KKTmat,FormJacobian,solver,ierr)
  call SNESSetDM(solver%snes,solver%da,ierr)

  ! set shell PC
  call SNESGetKSP(solver%snes,innerksp,ierr)
  call KSPGetPC(innerksp,spc,ierr)  
  call PCSetType(spc,PCSHELL,ierr) ! overriding opions !!!
  call PCShellSetApply(spc,PCShermanMorrisonApply,ierr)
  call PCSetApplicationContext(spc,solver,ierr)

  ! create solver::pc -- fieldsplit
  call PCCreate(solver%comm,solver%pc,ierr)
  call PCSetOptionsPrefix(solver%pc, 'fsa_', ierr)
  call PCSetFromOptions(solver%pc,ierr)
  call PCSetDM(solver%pc,solver%da,ierr)
  call PCSetOperators(solver%pc,solver%KKTmat,solver%KKTmat,ierr)
  call PCSetUp(solver%pc,ierr)
#ifndef PETSC_ED
#if (!PETSC_VERSION_LT(3,5,0))
  call PCFieldSplitSetSchurPre(solver%pc,PC_FIELDSPLIT_SCHUR_PRE_USER,solver%schur,ierr)
#endif
#endif
  call PCFieldSplitGetSubKSP(solver%pc,i,subksp,ierr)
  call KSPSetOperators(subksp(2),solver%schur,solver%schur,ierr)
  call KSPSetup(subksp(1),ierr)
  call KSPSetup(subksp(2),ierr)

  call SNESSetUp(solver%snes,ierr)
  CHKERRQ(ierr)
#if defined(PETSC_USE_LOG)
  call PetscLogStagePop(ierr)
#endif
#endif
end subroutine create_2field_solver

subroutine make_12mat(grid,Bmat,BIdentMat,FSAMass_Te,solver,bd,scale,xgc_petsc)
  use petscsnes
  use grid_class
  use psn_class
  use eq_module
  use sml_module,only:sml_mype
  use xgc_solver_module
  implicit none
  type(grid_type)::grid
  Mat::BIdentMat  ! out, but created before
  Mat,intent(in)::FSAMass_Te
  Mat,intent(out)::Bmat
  type(xgc_solver)::solver
  type(boundary2_type),intent(in)::bd
  real (8),intent(in)::scale
  PetscInt,intent(in)::xgc_petsc(grid%nnode)
  !
  PetscInt::im1,ip,low,high,idxj(1),nloc
  PetscErrorCode::ierr
  PetscInt,parameter::ione=1
  PetscScalar::value(1)
  real (8)::pn,wp
  integer::xgc_idx,ii
  PetscReal::fill

  if(sml_mype==0) print *, ' FSA: create B mat. M,N=',grid%cnv_2d_00%m,grid%cnv_2d_00%n,', First vertix in domain:',is_inside(1,bd)

  call MatGetOwnershipRange(BIdentMat,low,high,ierr)
  nloc = high-low

  ! make distributer matrix BIdentMat that make phi at vertices from <phi>
  do ii=1,nloc
     xgc_idx = solver%petscloc_xgc(ii) + 1 ! looking into global data structure
     pn = (grid%psi(xgc_idx)-grid%psi00min)/grid%dpsi00
     ip = floor(pn)+1
     if(0 < ip .and. ip < grid%npsi00 .and. is_rgn12(grid%x(1,xgc_idx),grid%x(2,xgc_idx),grid%psi(xgc_idx)) ) then
        wp=1D0 - ( pn - real(ip-1,8) )
     elseif( ip <= 0 ) then
        ip=1
        wp=1D0
     else
        ip=grid%npsi00-1
        wp=0D0
     endif

     ! BC check
     if( is_inside(xgc_idx,bd) ) then
        im1 = xgc_petsc(xgc_idx)  ! zero based global PETSc row
        ! inner
        idxj(1) = ip-1
        value(1) = wp 
        call MatSetValues(BIdentMat,ione,im1,ione,idxj,value,INSERT_VALUES,ierr)
        CHKERRQ(ierr)
        ! outer
        idxj(1) = ip
        value(1) = (1D0-wp) 
        call MatSetValues(BIdentMat,ione,im1,ione,idxj,value,INSERT_VALUES,ierr)
        CHKERRQ(ierr)
     endif
  end do
  call MatAssemblyBegin(BIdentMat,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(BIdentMat,MAT_FINAL_ASSEMBLY,ierr)

  ! B = -BMass*1/T*BIdent
  fill = 2.d0
  call MatMatMult(FSAMass_Te,BIdentMat,MAT_INITIAL_MATRIX,fill,Bmat,ierr)
  fill = -1.d0
  call MatScale(Bmat,fill,ierr)

  return 
end subroutine make_12mat

subroutine make_21mat(grid,Cmat,Bmat,CConstBCVec,bd,xgc_petsc)
  use petscsnes
  use grid_class
  use psn_class
  use sml_module, only : sml_mype, sml_plane_index
  implicit none
  type(grid_type) :: grid
  Mat :: Cmat
  Mat,intent(in) :: Bmat
  Vec :: CConstBCVec
  type(boundary2_type),intent(in)::bd
  PetscInt,intent(in)::xgc_petsc(grid%nnode)
  !
  PetscInt :: j,k,km1,im1,low,high,i,nlocRows,nGlobalCols,nLocCols
  PetscErrorCode :: ierr
  PetscInt, parameter :: ione=1
  PetscScalar :: value
  PetscInt,allocatable :: d_nnz(:),o_nnz(:)

  ! Cmat - count nnz, can't can call MatGetOwnershipRange before prealloc. Use transpose property.
  call MatGetOwnershipRangeColumn(Bmat,low,high,ierr);CHKERRQ(ierr)
  nlocRows = high-low
  allocate(d_nnz(nlocRows),o_nnz(nlocRows))
  d_nnz = 0
  o_nnz = 0
  do i=1,grid%cnv_2d_00%n ! columns of C
     do j=1,grid%cnv_2d_00%nelement(i) ! rows of C
#ifndef OPTIM_GYRO_AVG_MAT
        k = grid%cnv_2d_00%eindex(j,i)
        !y(k) = y(k) + mat%value(j,i) * x(i)
#else
        k = grid%cnv_2d_00%eindex(grid%cnv_2d_00%a(i)+j)
#endif
        km1=k-1
        if (km1 .ge.low .and. km1 .lt. high ) then ! I own this row
           d_nnz(km1-low+1) = d_nnz(km1-low+1) + 1 
           o_nnz(km1-low+1) = o_nnz(km1-low+1) + 1
        end if
     end do
  end do

  ! transpose
  call MatGetSize(Bmat,nGlobalCols,PETSC_NULL_INTEGER,ierr)
  call MatGetLocalSize(Bmat,nLocCols,PETSC_NULL_INTEGER,ierr)
  do i=1,nlocRows
     if (d_nnz(i) .gt. nLocCols) d_nnz(i) = nLocCols
     if (o_nnz(i) .gt. nGlobalCols-nLocCols) o_nnz(i) = nGlobalCols-nLocCols
  end do

  call MatSeqAIJSetPreallocation(Cmat,ione,d_nnz,ierr)
  CHKERRQ(ierr)
  call MatMPIAIJSetPreallocation(Cmat,ione,d_nnz,ione,o_nnz,ierr)
  CHKERRQ(ierr)
  deallocate(d_nnz,o_nnz)
  call MatSetOption(Cmat,MAT_IGNORE_OFF_PROC_ENTRIES,PETSC_TRUE,ierr)
  call MatSetup(Cmat,ierr)
  call VecSetOption(CConstBCVec,VEC_IGNORE_OFF_PROC_ENTRIES,PETSC_TRUE,ierr)

  ! set C, my columns only.
  do i=1,grid%cnv_2d_00%n ! column of C
     im1 = xgc_petsc(i)
     do j=1,grid%cnv_2d_00%nelement(i) ! row of C
#ifndef OPTIM_GYRO_AVG_MAT
        k = grid%cnv_2d_00%eindex(j,i)
        km1=k-1 ! zero bases row of matrix, same as XGC
        value=grid%cnv_2d_00%value(j,i)/grid%cnv_norm_1d00(k)
#else
        k = grid%cnv_2d_00%eindex(grid%cnv_2d_00%a(i)+j)
        km1=k-1 ! zero bases row of matrix, same as XGC
        value=grid%cnv_2d_00%value(grid%cnv_2d_00%a(i)+j)/grid%cnv_norm_1d00(k)
#endif
        if (im1.ge.0) then
           call MatSetValues(Cmat,ione,km1,ione,im1,value,INSERT_VALUES,ierr)
           CHKERRQ(ierr)
        else
           ! cache this for <Q> - area of BC area
           call VecSetValues(CConstBCVec,ione,km1,value,ADD_VALUES,ierr)
        end if
     enddo
  enddo
  call MatAssemblyBegin(Cmat,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(Cmat,MAT_FINAL_ASSEMBLY,ierr)
  call VecAssemblyBegin(CConstBCVec,ierr)
  call VecAssemblyEnd(  CConstBCVec,ierr)

  return  
end subroutine make_21mat

!************************************************************************
! solve_poisson - high level solve, init physics, solve
!************************************************************************
subroutine solve_poisson(grid,psn,iflag)
  use psn_class
  use grid_class
  use sml_module
  implicit none
  type(grid_type)::grid
  type(psn_type)::psn
  integer,intent(in)::iflag 

  if(sml_iter_solver) then
     call solve_poisson_iter(grid,psn,iflag)
  else
     call solve_poisson_private(grid,psn,iflag)
  endif

end subroutine solve_poisson
!
! the solver
!
subroutine solve_poisson_private(grid,psn,iflag)
  use petscsnes
  use psn_class
  use grid_class
  use sml_module
  use perf_monitor
  use smooth_module
  use omp_module, only: split_indices
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  integer,intent(in) :: iflag ! -1 for init, 1==00 mode+turb, 2==turb only, unless solver_type/=0
  !
  integer :: i,j,loop_num
  PetscErrorCode :: ierr
  integer :: ith,i_beg(sml_nthreads),i_end(sml_nthreads)
  real (kind=8),allocatable::dentmp0(:),dentmp2(:),dentmp3(:),dentmp4(:)
  real (kind=8),allocatable::sendl(:),recvr(:)
  save dentmp0, dentmp2, dentmp3, dentmp4
  save sendl,recvr
  PetscScalar, parameter :: one = 1.d0
  PetscScalar, parameter :: none = -1.d0
  PetscInt, parameter :: ione=1
  PetscInt, parameter :: itwo=2
  Vec::Fsub(2),Xsub(2),vec
  PetscScalar,pointer::phi_v(:)
  PetscInt::low,high,nnodes64
  PetscScalar::norm,norm0
  real (kind=8),allocatable::alpha(:),beta(:)
  Mat::massmat

  nnodes64 = grid%nnode
  
  if(iflag==-1) then ! this is a memory leak
     allocate(dentmp0(grid%nnode),dentmp2(grid%nnode),dentmp3(grid%nnode),dentmp4(grid%nnode))
     allocate(sendl(grid%nnode),recvr(grid%nnode))
     dentmp0=0D0 !for safety
     dentmp2=0D0
     dentmp3=0D0
     dentmp4=0D0
     if(sml_mype==0) print *, &
          '********** solve_poisson_private static variable initialized ************'
     return
  else
     !if(sml_mype==0) print *, 'solve_poisson_private: first vertex in domain:',is_inside(1,psn%pbd0_2)
  endif

  if(sml_electron_on) then
     dentmp0=psn%idensity0-psn%edensity0
  else
     dentmp0=psn%idensity0
  endif

  ! RHS density for n=0 mode
  call t_startf("POISSON_00MODE")

  psn%pot0m = 0d0 ! used as initial guess -- should use, might help some
  if(sml_poisson_solver_type==0 .and. (iflag==1 .and. sml_00_poisson) ) then
     call t_startf("POISSON_00MODE_SOLVE")
     if(sml_use_simple00) then
        call simple00(grid,psn)        
     else
        call set_boundary2_values(dentmp0,0D0,psn%pbd0_2)
        call convert_grid_2_001d(grid,dentmp0,psn%cden00_1d)
        call convert_001d_2_grid(grid,psn%cden00_1d,dentmp0)
        call set_boundary2_values(dentmp0,0D0,psn%pbd0_2)
        psn%pot0m = 0d0
        call petsc_solve(grid%nnode,dentmp0,0,psn%pot0m,psn%solver00%comm,psn%solver00,ierr)
        CHKERRQ(ierr)
     endif
     call t_stopf("POISSON_00MODE_SOLVE")

     !psn%pot0 is flux label
     if(.not.sml_use_simple00) call convert_grid_2_001d(grid,psn%pot0m,psn%pot00_1d)
     call convert_001d_2_grid(grid,psn%pot00_1d,psn%pot0)
#ifdef XGC1
  elseif(sml_poisson_solver_type/=0 .and. (sml_00_poisson .or. sml_turb_poisson) ) then
# else
  elseif(sml_poisson_solver_type/=0 .and. sml_00_poisson ) then
#endif
     call t_startf("POISSON_00MODE_SOLVE")
     call set_boundary2_values(dentmp0,0D0,psn%pbd0_2)
     !if (sml_plane_index==0) call VecView(psn%solver00%CConstBCVec,PETSC_VIEWER_STDOUT_SELF,ierr)
     psn%pot0m = 0d0
     call petsc_solve(grid%nnode,dentmp0,0,psn%pot0m,psn%solver00%comm,psn%solver00,ierr)
     CHKERRQ(ierr)
     call t_stopf("POISSON_00MODE_SOLVE")
     ! pot0m has n=0,m=aribtrary component
     ! Flux average of pot0m goes to pot00_1d (psi-space)
     ! And, flux averaged potential copied to psn%pot0
     call convert_grid_2_001d(grid,psn%pot0m,psn%pot00_1d)
     call convert_001d_2_grid(grid,psn%pot00_1d,psn%pot0)
  endif
  call t_stopf("POISSON_00MODE")

#ifdef XGC1
  ! XGC1 --->
  ! solve turbulence field********************************************************************
  call split_indices(grid%nnode, sml_nthreads, i_beg, i_end)

  if(sml_turb_poisson) then
     call t_startf("POISSON_TURB")
     ! For all poloidal plane - without n=0 mode
     !$OMP PARALLEL DO &
     !$OMP PRIVATE( ITH, I )
     do ith=1,sml_nthreads
        if(iflag==1 .and. sml_electron_on) then !use non-adiabatic electron
           do i=i_beg(ith),i_end(ith)
              dentmp3(i)= (psn%idensity(i,1) - psn%edensity(i,1) - dentmp0(i))   !????? check later convert_001d_2_grid 
              !/psn%density0_full(:)  ! ni - ne - <ni - ne>
           enddo
        else ! iflag==2 .or. no electron
           if(sml_poisson_solver_type==0) then
              call convert_001d_2_grid(grid,psn%iden00_1d,dentmp0) !dentmp0 becomes flux averaged ion 00 density
           endif
           do i=i_beg(ith),i_end(ith)
              dentmp3(i)= (psn%idensity(i,1) - dentmp0(i))
#ifdef MAX_DDEN
           if(dentmp3(i) > 2d0*dentmp0(i)) dentmp3(i)=2d0*dentmp0(i) ! factor 2 genralize later
#endif
              !  ni - <ni>
           enddo
        endif
     enddo

     ! smoothing of charge density -- RHS
     !call smooth_tr_connect(grid,dentmp3)
     !call smooth_pol(grid,dentmp3)
     call set_boundary2_values(dentmp3,0D0,psn%pbdH_2)
 
     psn%dden(:,1)=dentmp3  ! this can be optimized. redundant memory usage.
     call t_startf("POISSON_TURB_SOLVE")
     psn%dpot(:,1) = 0d0
     call petsc_solve(grid%nnode,dentmp3,0,psn%dpot(:,1),psn%solverH%comm,psn%solverH,ierr)
     CHKERRQ(ierr)
     if(sml_turb_poisson_n0_only) psn%dpot(:,1)=0D0

     call t_stopf("POISSON_TURB_SOLVE")

     if(sml_poisson_solver_type/=0) then
        ! move m/=0 component to dpot
        psn%dpot(:,1)=psn%dpot(:,1) + psn%pot0m - psn%pot0
     else
        !call smooth_tr_connect(grid,psn%dpot(:,1))

        !************************* extract 00 mode *********************************************
        call t_startf("EXTRACT_00MODE")
        call extract_00mode(grid,psn%dpot)
        call t_stopf("EXTRACT_00MODE")
     end if

     if(sml_mode_select_on==1) then
        call mode_selection(sml_mode_select_n,grid,psn)
     endif

     !apply min-max of dpot
     sml_dpot_min_factor=-1D0
     sml_dpot_max_factor=1D0
!     if(sml_dpot_min_max) then
        !$OMP PARALLEL DO &
        !$OMP PRIVATE( ITH, I )
        do ith=1,sml_nthreads
          do i=i_beg(ith),i_end(ith)
              psn%dpot(i,1)=min(max(psn%dpot(i,1),sml_dpot_min_factor*psn%tempe_ev(i)),sml_dpot_max_factor*psn%tempe_ev(i))
          enddo
        enddo
!     endif
     call t_stopf("POISSON_TURB")

  endif
  
  !********************** Phase 4 *********************************************
  !send and receive potential value to and from adjacent PE
  !***************************************************************************
  if(sml_turb_poisson) then
     call t_startf("POISSON_SR_POT")
     call send_recv_potential( psn%dpot, recvr,sendl,grid%nnode)
     call t_stopf("POISSON_SR_POT")
  endif

#else
  ! XGCa --->
  ! Smoothing has to be done here, otherwise, the weight evolution of the electrons
  ! will not be consistent with the electric field used in pushe.
  dentmp3=psn%pot0m
  if (smooth_pol_efield) then
    call t_startf("POISSON_POL_SMOOTH")
    call smooth_pol(grid,dentmp3,psn%dpot)
    call t_stopf("POISSON_POL_SMOOTH")
  else
    psn%dpot=dentmp3
  endif

  !extract 00 solution
  call convert_grid_2_001d(grid,psn%dpot,psn%pot00_1d)
  call convert_001d_2_grid(grid,psn%pot00_1d,psn%pot0)

  if(sml_sheath_mode==1) call zero_out_wall(psn%pot0)
  psn%dpot=psn%dpot-psn%pot0

  ! Ignore dpot outside some boundaries
  !call set_boundary2_values(psn%dpot,0D0,psn%pbdH_2)

#endif

  !******************** Phase 5 ***********************************************
  ! get space derivative of potential
  ! psn%E_para(node,1:nphi )  psn%E_perp(2 vec, ntriangle, 0:nphi)
  !***************************************************************************
  if(sml_add_pot0/=0) then
     if(sml_replace_pot0) then
!$OMP PARALLEL DO &
!$OMP PRIVATE( ITH, I )
        do ith=1,sml_nthreads
           do i=i_beg(ith),i_end(ith)
              psn%pot0(i)=psn%add_pot0(i)
           enddo
        enddo
     else
!$OMP PARALLEL DO &
!$OMP PRIVATE( ITH, I )
        do ith=1,sml_nthreads
           do i=i_beg(ith),i_end(ith)
              psn%pot0(i)=psn%pot0(i)+psn%add_pot0(i)
           enddo
        enddo
     endif
  endif
  
  ! debug - check maps, put partitioning in psn%dpot(:,1)
  if(.false.) then
     call VecCreateSeqWithArray(PETSC_COMM_SELF,ione,nnodes64,dentmp4,vec,ierr)
     dentmp4 = sml_mype
     call VecScatterBegin(psn%solver00%to_petsc,vec,psn%solver00%bVec,INSERT_VALUES,SCATTER_FORWARD,ierr)
     call VecScatterEnd  (psn%solver00%to_petsc,vec,psn%solver00%bVec,INSERT_VALUES,SCATTER_FORWARD,ierr) 
     call VecGetOwnershipRange(psn%solver00%xVec,low,high,ierr)
     norm=sml_mype
     do i=low,high-1
        call VecSetValue(psn%solver00%xVec,i,norm,INSERT_VALUES,ierr)
     end do
     call VecAssemblyBegin(psn%solver00%xVec,ierr)
     call VecAssemblyEnd(  psn%solver00%xVec,ierr)
     call VecAXPY(psn%solver00%xVec,-1D0,psn%solver00%bVec,ierr) 
     call VecNorm(psn%solver00%xVec,NORM_INFINITY,norm,ierr)
     if (sml_mype==0) print *,' ******* ERROR=',norm
     call VecScatterBegin(psn%solver00%from_petsc,psn%solver00%bVec,vec,INSERT_VALUES,SCATTER_FORWARD,ierr)
     call VecScatterEnd  (psn%solver00%from_petsc,psn%solver00%bVec,vec,INSERT_VALUES,SCATTER_FORWARD,ierr)
#ifdef XGC1
     if(sml_plane_index==0) psn%dpot(:,1)=dentmp4
#else
     ! What is this for?
     if (sml_plane_index==0) psn%dpot(:)=dentmp4
#endif
     call VecDestroy(vec,ierr)
  else if (.false.) then
     call DMCompositeGetAccessArray(psn%solver00%da,psn%solver00%bVec2,itwo,PETSC_NULL_INTEGER,Fsub,ierr)
     call VecCreateSeqWithArray(PETSC_COMM_SELF,ione,nnodes64,dentmp3,vec,ierr)
     call VecScatterBegin(psn%solver00%from_petsc,Fsub(1),vec,INSERT_VALUES,SCATTER_FORWARD,ierr)
     call VecScatterEnd  (psn%solver00%from_petsc,Fsub(1),vec,INSERT_VALUES,SCATTER_FORWARD,ierr)
     call VecDestroy(vec,ierr)
#ifdef XGC1
     if(sml_plane_index==0) psn%dpot(:,1) = dentmp3 ! 1: rho
#else
     ! What is this for?
     if (sml_plane_index==0) psn%dpot(:) = dentmp3 ! 1: rho
#endif
     call DMCompositeRestoreAccessArray(psn%solver00%da,psn%solver00%bVec2,itwo,PETSC_NULL_INTEGER,Fsub,ierr)
     ! matrix residaul - PETSc, we have correct X and F 
     call VecNorm(psn%solver00%bVec2,NORM_INFINITY,norm0,ierr)
     call VecScale(psn%solver00%bVec2,-1d0,ierr)
     call MatMultAdd(psn%solver00%KKTmat,psn%solver00%xVec2,psn%solver00%bVec2,psn%solver00%bVec2,ierr) ! resid
     call VecNorm(psn%solver00%bVec2,NORM_INFINITY,norm,ierr)
     if(sml_mype==0.and.norm0/=0d0) write(*,'(A,ES12.4E2)') 'Petsc |b - (A+M(I-<>))x|/|b|=',norm/norm0
     if(sml_mype==0.and.norm0==0d0) write(*,'(A,ES12.4E2)') 'Petsc |b - (A+M(I-<>))x|=',norm
     call DMCompositeGetAccessArray(psn%solver00%da,psn%solver00%bVec2,itwo,PETSC_NULL_INTEGER,Fsub,ierr)
     call VecCreateSeqWithArray(PETSC_COMM_SELF,ione,nnodes64,dentmp3,vec,ierr)
     call VecScatterBegin(psn%solver00%from_petsc,Fsub(1),vec,INSERT_VALUES,SCATTER_FORWARD,ierr)
     call VecScatterEnd  (psn%solver00%from_petsc,Fsub(1),vec,INSERT_VALUES,SCATTER_FORWARD,ierr)
     call VecDestroy(vec,ierr)
     !if(sml_plane_index==1) psn%dpot(:,1) = dentmp3 ! 2: rho - (A + M(I-<>)) phi -- debug!!!
     call DMCompositeRestoreAccessArray(psn%solver00%da,psn%solver00%bVec2,itwo,PETSC_NULL_INTEGER,Fsub,ierr)
  end if

contains
  subroutine zero_out_wall(var)
    implicit none
    real (8), intent(inout) :: var(grid%nnode)
    integer :: i
 
    do i=1, psn%nwall
      var(psn%wall_nodes(i))=0D0
    enddo
  end subroutine zero_out_wall

end subroutine solve_poisson_private


subroutine solve_poisson_iter(grid,psn,iflag)
  use psn_class
  use grid_class
  use sml_module
  use perf_monitor
  use smooth_module
  !use diag_module, only :: diag_poisson_rhs
  use omp_module, only: split_indices
  implicit none
  type(grid_type)::grid
  type(psn_type)::psn
  integer,intent(in)::iflag
  integer::i,loop_num,iter,count
  PetscErrorCode::ierr
  integer::ith,i_beg(sml_nthreads),i_end(sml_nthreads)
  real (kind=8),allocatable::dentmp(:), dentmp2(:), dentmp3(:), dentmp4(:)
  real (kind=8),allocatable::sendl(:),recvr(:)
  real (kind=8) :: t1, t2, t3
  save dentmp, dentmp2, dentmp3, dentmp4
  save sendl,recvr,count
  data count /0/

  if(iflag==-1) then
     allocate(dentmp(grid%nnode),dentmp2(grid%nnode),dentmp3(grid%nnode),dentmp4(grid%nnode))
     allocate(sendl(grid%nnode),recvr(grid%nnode))
     dentmp=0D0 !for safety
     dentmp2=0D0
     dentmp3=0D0
     dentmp4=0D0
     if(sml_mype==0) print *, &
          '**************** Poisson outer iterative solver is used ******************'
     return
  endif

  call t_startf("POISSON_00MODE")


  ! get n0*******************************************************************************
  ! RHS density for n=0 mode --> save it as dentmp
  psn%cden00_1d = psn%iden00_1d !-psn%eden00_1d -- commented out to keep the same initial condition for a given ion densiy

  if(sml_f0_grid) then   ! to include effect of f0_grid, use idensity0 and discard iden00_1d
     call convert_grid_2_001d(grid,psn%idensity0,psn%cden00_1d)
  endif

  if(iflag==1 .and. sml_electron_on) then
     dentmp=psn%idensity0-psn%edensity0
  else
     dentmp=psn%idensity0
  endif

  call set_boundary2_values(dentmp,0D0,psn%pbd0_2)
  call zero_out_wall(dentmp)
  psn%rhs1=dentmp

  ! solve n=0 field iteratively ************************************************************
  ! set initial value - from simple00
  call simple00(grid,psn)
  call convert_001d_2_grid(grid,psn%pot00_1d,psn%pot0)
  
  ! solve iterativley
  do iter=1,sml_iter_solver_niter
     ! update RHS
     call convert_grid_2_001d(grid,psn%pot0,psn%pot00_1d)
     call convert_001d_2_grid(grid,psn%pot00_1d,dentmp2)  !dentmp2 is flux averaged solution
     call zero_out_wall(dentmp2)
     ! write out 1d solutions - convergence test
     if(.false.) then
        if(sml_mype==0) then
           if(iter==1) open(unit=567,file='itersol.txt',status='replace')
           write(567,1000) psn%pot00_1d
1000       format(400(e19.12,1x))
           if(iter==sml_iter_solver_niter) close(567)
        endif
     endif
     ! solve poisson eq
     call t_startf("POISSON_00MODE_SOLVE")
     call petsc_solve(grid%nnode,dentmp,dentmp2,psn%pot0,psn%solver00%comm,psn%solver00,ierr)
     call t_stopf("POISSON_00MODE_SOLVE")
  enddo

  !extract 00 solution
  call convert_grid_2_001d(grid,psn%pot0,psn%pot00_1d)
  call convert_001d_2_grid(grid,psn%pot00_1d,dentmp2)  !dentmp2 is flux averaged solution
  call zero_out_wall(dentmp2)
  psn%pot0m=psn%pot0

  call t_stopf("POISSON_00MODE")


#ifdef XGC1
  ! solve turbulence field**********************************************************************************
  call split_indices(grid%nnode, sml_nthreads, i_beg, i_end)

  if(sml_turb_poisson) then
     call t_startf("POISSON_TURB")

     ! For all toroidal section
     ! the case when phitmp=>psn%pot0 can be eliminated. No function now.

!$OMP PARALLEL DO &
!$OMP PRIVATE( ITH, I )
     do ith=1,sml_nthreads
        if(iflag==1 .and. sml_electron_on) then
           do i=i_beg(ith),i_end(ith)
              dentmp3(i)= (psn%idensity(i,1) - psn%edensity(i,1) - dentmp(i))
              !/psn%density0_full(:)  ! ni - ne - <ni - ne>
           enddo
        else ! iflag==2 .or. no electron
!           call convert_001d_2_grid(grid,psn%iden00_1d,dentmp) !dentmp becomes flux averaged ion 00 density
           do i=i_beg(ith),i_end(ith)
              dentmp3(i)= (psn%idensity(i,1) - psn%idensity0(i))
              !  ni - <ni>
           enddo
        endif
     enddo

     ! smoothing of charge density -- RHS
     !call smooth_tr_connect(grid,dentmp3)
     !call smooth_pol(grid,dentmp3)
     call zero_out_wall(dentmp3)
     call set_boundary2_values(dentmp3,0D0,psn%pbdH_2)
     psn%rhs2=dentmp3 
     psn%dden(:,1)=dentmp3  ! this can be optimized. redundant memory usage.
     call t_startf("POISSON_TURB_SOLVE")
     psn%dpot(:,1) = 0d0
     call petsc_solve(grid%nnode,dentmp3,0,psn%dpot(:,1),psn%solverH%comm,psn%solverH,ierr)
     call t_stopf("POISSON_TURB_SOLVE")

     ! move m/=0 component to dpot
     psn%dpot(:,1)=psn%dpot(:,1)+psn%pot0 - dentmp2
     psn%pot0=dentmp2

     !apply min-max of dpot
     sml_dpot_min_factor=-1D0
     sml_dpot_max_factor=1D0
!     if(sml_dpot_min_max) then
        !$OMP PARALLEL DO &
        !$OMP PRIVATE( ITH, I )
        do ith=1,sml_nthreads
          do i=i_beg(ith),i_end(ith)
              psn%dpot(i,1)=min(max(psn%dpot(i,1),sml_dpot_min_factor*psn%tempe_ev(i)),sml_dpot_max_factor*psn%tempe_ev(i))
          enddo
        enddo
!     endif

     ! call smooth_pol0(grid,psn%dpot(:,1),smoothH)
     ! call smooth_r(psn%dpot(:,1),grid%rtmp1,smooth_r1,grid)
     ! call smooth_nearx(psn%dpot(:,1),smooth_nearx1,grid)
     !call smooth_tr_connect(grid,psn%dpot(:,1))
     !call smooth_pol(grid,psn%dpot(:,1))


     !************************* extract 00 mode **********************************************************
     !call t_startf("EXTRACT_00MODE")
     !call extract_00mode(grid,psn%dpot)
     !call t_stopf("EXTRACT_00MODE")

     !**************************** mode selection*********************************************************
     if(sml_mode_select_on==1) then
        call mode_selection(sml_mode_select_n,grid,psn)
     endif

     call t_stopf("POISSON_TURB")
  endif


  !********************** Phase 4 *********************************************
  !send and receive potential value to and from adjacent PE
  !***************************************************************************
  if(sml_turb_poisson) then
     call t_startf("POISSON_SR_POT")
     call send_recv_potential( psn%dpot, recvr,sendl,grid%nnode)
     call t_stopf("POISSON_SR_POT")
  endif

#else
  ! XGCa --->
  ! Smoothing has to be done here, otherwise, the weight evolution of the electrons
  ! will not be consistent with the electric field used in pushe.
  dentmp=psn%pot0m
  if (smooth_pol_efield) then
    call t_startf("POISSON_POL_SMOOTH")
    call smooth_pol(grid,dentmp,psn%dpot)
    call t_stopf("POISSON_POL_SMOOTH")
  else
    psn%dpot=dentmp
  endif

  !extract 00 solution
  call convert_grid_2_001d(grid,psn%dpot,psn%pot00_1d)
  call convert_001d_2_grid(grid,psn%pot00_1d,psn%pot0)

  if(sml_sheath_mode==1) call zero_out_wall(psn%pot0)
  psn%dpot=psn%dpot-psn%pot0

  ! Ignore dpot outside some boundaries
  !call set_boundary2_values(psn%dpot,0D0,psn%pbdH_2)

! on ifdef XGC1 --->
#endif

  !******************** Phase 5 ***********************************************
  ! get space derivative of potential
  ! psn%E_para(node,1:nphi )  psn%E_perp(2 vec, ntriangle, 0:nphi)
  !***************************************************************************
  if(sml_add_pot0/=0) then
     if(sml_replace_pot0) then
!$OMP PARALLEL DO &
!$OMP PRIVATE( ITH, I )
        do ith=1,sml_nthreads

           do i=i_beg(ith),i_end(ith)
              psn%pot0(i)=psn%add_pot0(i)
           enddo

        enddo
     else
!$OMP PARALLEL DO &
!$OMP PRIVATE( ITH, I )
        do ith=1,sml_nthreads

           do i=i_beg(ith),i_end(ith)
              psn%pot0(i)=psn%pot0(i)+psn%add_pot0(i)
           enddo

        enddo
     endif
  endif
contains
  subroutine zero_out_wall(var)
    implicit none
    real (8), intent(inout) :: var(grid%nnode)
    integer :: i
 
    do i=1, psn%nwall
      var(psn%wall_nodes(i))=0D0
    enddo
  end subroutine zero_out_wall
end subroutine solve_poisson_iter


