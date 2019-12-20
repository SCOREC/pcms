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
subroutine init_poisson(grid,psn,bc)
  use psn_class
  use grid_class
  use sml_module
  use eq_module
  implicit none
  type(grid_type)::grid
  type(psn_type)::psn
  type(boundary2_type),intent(in)::bc
  !
  integer::i
  real (kind=8)::psi,r,z
  PetscErrorCode::ierr

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
  endif

  call psn_init_solver(psn,grid,bc)
  
  if(sml_add_pot0>0) then
     stop 'lets not allow sml_add_pot0 for now'
     call check_point('add potential offset')
     allocate( psn%add_pot0(grid%nnode))

     ! prepare additional 0 - potential
     if(sml_add_pot0==1) then
        call read_add_pot0(grid,psn)
     elseif(sml_add_pot0==2) then
!        call neo_pot0_simple(grid,psn)
     else
        psn%add_pot0=0D0
     endif
  endif

  call check_point('prepare gyro-averaging matrix')
  if(sml_plane_mype<grid%nrho*2) then
     call init_gyro_avg_mat(grid,psn)
  endif
end subroutine init_poisson

subroutine psn_init_solver(psn,grid,bc)
  use petscsnes
  use psn_class
  use grid_class
  use sml_module
  use ptl_module,only:ptl_mass
  use eq_module
  implicit none
  type(grid_type)::grid
  type(psn_type)::psn
  type(boundary2_type),intent(in)::bc
  
 ! type(mat_type) :: matl, matr, matr2
  real (kind=8),allocatable :: alpha(:),beta(:),teev(:),den(:),b2(:)
  integer :: itr,nd(3),in_bd_node
  PetscErrorCode :: ierr
  real (kind=8) :: x_center(2),psi
  real (kind=8) :: factor,psi_in_poisson_bd,psi_out_poisson_bd
  real (kind=8), external :: gyro_radius2,psi_interpol,b_interpol
  real (kind=8), parameter :: offset=2D0
  Vec::vec
  PetscInt :: ione

  ione=1

  allocate(alpha(grid%ntriangle),beta(grid%ntriangle),teev(grid%ntriangle),den(grid%ntriangle),b2(grid%ntriangle))

  ! setup n_0,b^2,Te at centroids of triangles
  do itr=1, grid%ntriangle
     nd=grid%nd(:,itr)
     x_center=(grid%x(:,nd(1))+grid%x(:,nd(2))+grid%x(:,nd(3)))/3D0
     psi=psi_interpol(x_center(1),x_center(2),0,0)
     teev(itr)=eq_ftn(psi,x_center(1),x_center(2),eq_tempe)
     b2(itr)=b_interpol(x_center(1),x_center(2),0D0)**2
     den(itr)=eq_ftn(psi,x_center(1),x_center(2),eq_den)
  enddo

  ! 00 solver
  call check_point('Initialize 00 matrix solver')
  ! **************************************
  ! initialize petsc solver
  ! **************************************
  psn%solver00%comm = sml_plane_comm
  psn%solver00%mype = sml_plane_mype
  psn%solver00%totalpe = sml_plane_totalpe
  psn%solver00%prefix = 1
  psn%solver00%n_rhs_mat = 1 ! set RHS mode (solver_init does not need this anymore!)
  
  ! create 00 solver - this is the Poisson solver
  call check_point('- 00 matrix solver memory create')
  call solver_init(psn%solver00,grid,bc,ierr)
  CHKERRQ(ierr)
  ! LHS - alpha del^2 u
  alpha=ptl_mass(1)/sml_e_charge*den/b2 ! mn_0/eB^2
  beta=0D0
  if(sml_mype==0) print *, 'psn_init_solver: create Amat, first vertex in domain:',is_inside(1,bc)
  if (sml_sheath_mode/=0 .and. sml_mype==0) print *, 'psn_init_solvers: sml_sheath_mode/=0  --> sml_sheath_mode==0 used instead for poisson solver'
  call helm_matrix(psn%solver00%Amat,alpha,beta,grid,bc,.false.,psn%solver00%xgc_petsc,ierr)

  ! RHS mass matrix
  if(sml_mype==0) print *, 'psn_init_solver: create rhs mat, first vertex in domain:',is_inside(1,bc)
  call MatDuplicate(psn%solver00%Amat,MAT_DO_NOT_COPY_VALUES,psn%solver00%rhs_mat,ierr)
  if (sml_lumped_mass) then
     call VecCreateSeqWithArray(PETSC_COMM_SELF,ione,grid%nnode,grid%node_area,vec,ierr);CHKERRQ(ierr)
     call VecScatterBegin(psn%solver00%to_petsc, vec, psn%solver00%xVec, INSERT_VALUES,SCATTER_FORWARD,ierr)
     call VecScatterEnd(  psn%solver00%to_petsc, vec, psn%solver00%xVec, INSERT_VALUES,SCATTER_FORWARD,ierr)
     call VecDestroy(vec,ierr)

     call MatDiagonalSet(psn%solver00%rhs_mat, psn%solver00%xVec, INSERT_VALUES, ierr);CHKERRQ(ierr)
  else
     alpha=0D0  ! no del term on RHS
     beta=1d0  ! diagonal (identity)
     call helm_matrix(psn%solver00%rhs_mat,alpha,beta,grid,bc,.false.,psn%solver00%xgc_petsc,ierr)
  end if

  call create_1field_solver(grid,psn%solver00)
  
  deallocate(alpha,beta,teev,den,b2)

end subroutine psn_init_solver
!
! helper function to compute nnz per row of my rows
!
subroutine getNNZ(grid,nloc,low,high,d_nnz,o_nnz,xgc_petsc,ierr)
  use petscsnes
  use grid_class
  use sml_module
  implicit none
  type(grid_type),intent(in)::grid
  PetscInt,intent(in)::nloc
  PetscInt,intent(in)::xgc_petsc(grid%nnode)
  !
  PetscInt::d_nnz(nloc),o_nnz(nloc),low,high
  PetscErrorCode::ierr
  !
  PetscInt::ie,ind(3),itmp,idi(1),ii,jtmp,idj(1),jj

  d_nnz = 2*2
  o_nnz = 1*2
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

  o_nnz = o_nnz/2
  d_nnz = d_nnz/2
  do ii=1,nloc
     if (d_nnz(ii).gt.nloc) d_nnz(ii) = nloc
     if (o_nnz(ii).gt.grid%nnode-nloc) o_nnz(ii) = grid%nnode-nloc
  end do
end subroutine getNNZ
!
! solver_init - create Amat and its data layout, and this%petscloc_xgc, this%xgc_petsc
!
subroutine solver_init(this,grid,bc,ierr)
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
  PetscInt::nnode64,proc,iloc,nreal
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
  nnode64 = grid%nnode
  this%Amat = 0
  this%rhs_mat=0
  this%rhs2_mat=0
  ! LHS --------------
  if (allocated(this%xgc_petsc)) stop 'associated(this%xgc_petsc)'
  allocate(this%xgc_petsc(grid%nnode))  

  if(sml_mype==0) print *, 'solver_init: first vertex in domain:',is_inside(1,bc)
  
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
     if(sml_mype==0)write(*,1)'SOLVER_INIT: make partitioning with ',nreal,'/',grid%nnode,' real vertices'
1       format(A,I7,A,I7,A)
     call MatCreate(this%comm,AA,ierr);CHKERRQ(ierr)
     call MatSetSizes(AA,PETSC_DECIDE,PETSC_DECIDE,nreal,nreal,ierr);CHKERRQ(ierr)
     call MatSetType(AA,MATAIJ,ierr);CHKERRQ(ierr)
     call MatSetup(AA,ierr);CHKERRQ(ierr)
     call MatGetOwnershipRange(AA,low,high,ierr) ;CHKERRQ(ierr)
     npetscloc = high-low
     allocate(d_nnz(npetscloc),o_nnz(npetscloc))
     call getNNZ(grid,npetscloc,low,high,d_nnz,o_nnz,this%xgc_petsc,ierr);CHKERRQ(ierr)
     call MatSeqAIJSetPreallocation(AA,PETSC_NULL_INTEGER,d_nnz,ierr);CHKERRQ(ierr)
     call MatMPIAIJSetPreallocation(AA,PETSC_NULL_INTEGER,d_nnz,PETSC_NULL_INTEGER,o_nnz,ierr);CHKERRQ(ierr)
     deallocate(d_nnz,o_nnz)
     call MatSetOption(AA,MAT_IGNORE_OFF_PROC_ENTRIES,PETSC_TRUE,ierr)
     call MatSetOption(AA,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,ierr)
     call MatSetup(AA,ierr);CHKERRQ(ierr)
     beta=1d0
     alpha=1d0
     call helm_matrix(AA,alpha,beta,grid,bc,.false.,this%xgc_petsc,ierr)  
     ! partition
     call MatPartitioningCreate(this%comm,part,ierr);CHKERRQ(ierr)
     call MatPartitioningSetAdjacency(part,AA,ierr);CHKERRQ(ierr)
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
     call getNNZ(grid,npetscloc,low,high,d_nnz,o_nnz,this%xgc_petsc,ierr)
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
     if(sml_mype==0) print *, 'solver_init: warning simple partitionin gwith BCs not working, perhaps'
     ! simple chop partitioning
     allocate(petsc_xgc(grid%nnode))
     do ii=1,grid%nnode     
        this%xgc_petsc(ii) = ii - 1 ! default zero based map for assembly
        petsc_xgc(ii) = ii - 1 
     end do
     call MatCreate(this%comm,AA,ierr)
     call MatSetSizes(AA,PETSC_DECIDE,PETSC_DECIDE,nnode64,nnode64,ierr)
     call MatSetType(AA,MATAIJ,ierr)
     !call MatSeqAIJSetPreallocation(AA,maxn164,PETSC_NULL_INTEGER,ierr)  
     !call MatMPIAIJSetPreallocation(AA,maxn164,PETSC_NULL_INTEGER,maxn164/2,PETSC_NULL_INTEGER,ierr)
     stop 'solver_init: not supported'
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
     nreal = nnode64
  end if
  ! vecs
  call VecCreateMPI(this%comm,npetscloc,nreal,this%xVec,ierr)
  call VecSetFromOptions(this%xVec,ierr)
  call VecDuplicate(this%xVec,this%bVec,ierr)
  ! set scatters
  call VecCreateSeq(PETSC_COMM_SELF,nnode64,vec,ierr) ! big dumb PETSc vector
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
  ! zero out - not used for EM but zero out for destructor
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
  this%FSAMass = 0
  this%TeInv = 0
  this%snes = 0
  this%scale = 1d0
  this%schur = 0
  this%AinvB = 0
  this%ksp = 0
  this%pc = 0
end subroutine solver_init
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

  stop 'solve_poisson: not used!!!!'
  
  call solve_poisson_private(grid,psn,iflag)
  
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
  PetscScalar,pointer::phi_v(:)
  PetscInt::low,high,nnode64
  PetscScalar::norm
  Vec::vec
  real (kind=8),allocatable::alpha(:),beta(:)

  nnode64 = grid%nnode
  
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

  ! RHS density for n=0 mode
  call t_startf("POISSON_00MODE")

  psn%pot0m=0D0

  ! XGC1 --->
  ! solve turbulence field********************************************************************
  call split_indices(grid%nnode, sml_nthreads, i_beg, i_end)

  call t_startf("POISSON_TURB")
  ! For all poloidal plane - without n=0 mode
  !$OMP PARALLEL DO &
  !$OMP PRIVATE( ITH, I )
  do ith=1,sml_nthreads
     if(iflag==1) then !use non-adiabatic electron
        psn%edensity(:,1) = psn%eden_hyb
        do i=i_beg(ith),i_end(ith)
           !##              dentmp3(i)= (psn%idensity(i,1) - psn%edensity(i,1) - dentmp0(i))   !????? check later convert_001d_2_grid 
           dentmp3(i)= (psn%idensity(i,1) - psn%edensity(i,1) )   !????? check later convert_001d_2_grid 
           !/psn%density0_full(:)  ! ni - ne - <ni - ne>
        enddo
     else ! iflag==2 .or. no electron
        stop ' Error in turb_poisson'
     endif
  enddo
  
  ! smoothing of charge density -- RHS
  !call smooth_tr_connect(grid,dentmp3)
  !call smooth_pol(grid,dentmp3)
  call set_boundary2_values(dentmp3,0D0,psn%pbdH_2)
  
  psn%dden(:,1)=dentmp3  ! this can be optimized. redundant memory usage.
  call t_startf("POISSON_TURB_SOLVE")
  psn%dpot(:,1) = 0d0
  call petsc_solve(grid%nnode,dentmp3,0,psn%dpot(:,1),psn%solver00%comm,psn%solver00,ierr)
  CHKERRQ(ierr)
  call t_stopf("POISSON_TURB_SOLVE")

!!!!        call smooth_tr_connect(grid,psn%dpot(:,1))
  
  !************************* extract 00 mode *********************************************
  if(sml_mode_select_on==1) then
     call mode_selection_comb(sml_mode_select_n,grid,psn,psn%dpot(:,1))
  endif
  
  
  call convert_grid_2_001d(grid,psn%dpot(:,1),psn%pot00_1d)
  call convert_001d_2_grid(grid,psn%pot00_1d,psn%pot0)
  
  psn%dpot(:,1)=psn%dpot(:,1)-psn%pot0
  call t_stopf("POISSON_TURB")

  !********************** Phase 4 *********************************************
  !send and receive potential value to and from adjacent PE
  !***************************************************************************
  if(sml_turb_poisson) then
     call t_startf("POISSON_SR_POT")
     call send_recv_potential( psn%dpot, recvr,sendl,grid%nnode)
     call t_stopf("POISSON_SR_POT")
  endif

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

end subroutine solve_poisson_private


