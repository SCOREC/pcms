! ------------------------------------------------------------------------
! PETSc ODE time integrator (stepper) je
!
! ------------------------------------------------------------------------
!
! ts_matrix_init - create 3D Amat and its data layout:
!   a_ts: petscloc_xgc, xgc_petsc, to_petsc, from_petscscat
!   - offset global PETSc indices by (sml_plane_index-1)*grid%nnode
! ------------------------------------------------------------------------
subroutine ts_matrix_init(a_ts,grid,bc,ierr)
  use grid_class
  use xgc_ts_module
  use sml_module,only:sml_plane_totalpe,sml_mype,sml_plane_mype,sml_plane_index,sml_plane_comm,sml_comm,sml_nphi_total,sml_intpl_mype
  implicit none
#include <petscversion.h>
#if PETSC_VERSION_LT(3,6,0)
#include <finclude/petsc.h90>
#else
#include <petsc/finclude/petsc.h90>
#endif
  type(xgc_ts)::a_ts
  type(grid_type),intent(in)::grid
  type(boundary2_type),intent(in)::bc
  PetscErrorCode,intent(out)::ierr
  !
  integer::ixgc
  PetscErrorCode::ii
  PetscInt::proc,iloc,nreal2D,my0
  Mat::AA
  real (kind=8),allocatable::alpha(:),beta(:)  
  MatPartitioning::part;
  IS::is,is2;
  PetscInt,pointer::xgc_proc(:)
  PetscInt,allocatable::d_nnz(:),o_nnz(:)
  PetscInt::npetscloc,low,high,proc_eq(0:sml_plane_totalpe)
  PetscInt,parameter::ione=1
  PetscMPIInt::mype
  MPI_Comm::w_comm
  PetscViewer::viewer
  !
  call MPI_Comm_rank(sml_comm,mype,ierr)
  
  ! LHS --------------
  if (allocated(a_ts%xgc_petsc)) stop 'associated(a_ts%xgc_petsc)'
  allocate(a_ts%xgc_petsc(grid%nnode))

  if(sml_mype==0) print*,'ts_matrix_init: first vertex in domain:',is_inside(1,bc),', num planes=',sml_nphi_total

  ! get: nreal2D,npetscloc,petscloc_xgc,low,a_ts%xgc_petsc

  ! metis partitioning with one plane
  nreal2D=0 ! count locals with BC move to 0, set a_ts%xgc_petsc
  do ixgc=1,grid%nnode
     if(is_inside(ixgc,bc)) then ! put BCs on proc 0
        a_ts%xgc_petsc(ixgc) = nreal2D ! default zero based map for assembly   
        nreal2D = nreal2D + 1
     else
        a_ts%xgc_petsc(ixgc) = -1 ! petsc will skip & and fem does not look anyway
     end if
  end do 
  my0 = sml_plane_index*nreal2D ! the first index for this planes block
  allocate(alpha(grid%ntriangle),beta(grid%ntriangle))
  ! create dummy matrix for partitioning. could be done natively
  if(sml_mype==0)write(6,1)' ts_matrix_init: make partitioning with ',nreal2D,' real vertices in plane'
1 format(A,I7,A,I7,A)
  call MatCreate(sml_plane_comm,AA,ierr);CHKERRQ(ierr)
  call MatSetSizes(AA,PETSC_DECIDE,PETSC_DECIDE,nreal2D,nreal2D,ierr);CHKERRQ(ierr)
  call MatSetType(AA,MATAIJ,ierr);CHKERRQ(ierr)
  call MatSetup(AA,ierr);CHKERRQ(ierr)
  call MatGetOwnershipRange(AA,low,high,ierr);CHKERRQ(ierr)
  npetscloc = high-low
  allocate(d_nnz(npetscloc),o_nnz(npetscloc))
  ! old poisson.F90 helper
  call getNNZ(grid,npetscloc,low,high,d_nnz,o_nnz,a_ts%xgc_petsc,ierr);CHKERRQ(ierr)
  call MatSeqAIJSetPreallocation(AA,PETSC_NULL_INTEGER,d_nnz,ierr);CHKERRQ(ierr)
  call MatMPIAIJSetPreallocation(AA,PETSC_NULL_INTEGER,d_nnz,PETSC_NULL_INTEGER,o_nnz,ierr);CHKERRQ(ierr)
  deallocate(d_nnz,o_nnz)
  call MatSetup(AA,ierr);CHKERRQ(ierr)
  beta=1d0 ! not used
  alpha=0d0
  call helm_matrix(AA,alpha,beta,grid,bc,.false.,a_ts%xgc_petsc,ierr);CHKERRQ(ierr)

  ! debug
  if (sml_intpl_mype==0 .and. .false.) then
     ! write out mass
     call PetscViewerASCIIOpen(sml_plane_comm, 'M.m', viewer,ierr);CHKERRQ(ierr)
     call PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB,ierr);CHKERRQ(ierr)
     call MatView(AA,viewer,ierr);CHKERRQ(ierr)
     call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)       
  end if
  
  ! partition
  call MatPartitioningCreate(sml_plane_comm,part,ierr);CHKERRQ(ierr)
  call MatPartitioningSetAdjacency(part,AA,ierr);CHKERRQ(ierr)
  call MatPartitioningSetFromOptions(part,ierr);CHKERRQ(ierr)
  call MatPartitioningApply(part,is,ierr);CHKERRQ(ierr)

  ! clean up partitioning
  call MatPartitioningDestroy(part,ierr);CHKERRQ(ierr)
  call MatDestroy(AA,ierr);CHKERRQ(ierr)
  call ISGetLocalSize(is,npetscloc,ierr);CHKERRQ(ierr)
  if (npetscloc.ne.high-low) stop 'npetscloc.ne.high-low'
  call ISAllGather(is,is2,ierr);CHKERRQ(ierr) ! get global (2D) colors
  call ISDestroy(is,ierr);CHKERRQ(ierr)
  !if (sml_mype==0) call ISView(is2,PETSC_VIEWER_STDOUT_SELF,ierr)
  ! count proc sizes
  call ISGetIndicesF90(is2,xgc_proc,ierr);CHKERRQ(ierr)
  proc_eq = 0
  nreal2D = 0
  do ixgc=1,grid%nnode
     if(is_inside(ixgc,bc)) then
        nreal2D = nreal2D + 1
        proc = xgc_proc(nreal2D) + 1
        proc_eq(proc) = proc_eq(proc) + 1 ! increment, start at (1)
     end if
  end do
  npetscloc = proc_eq(sml_plane_mype+1) ! my new local size
  ! scan and get first 2D index for each plane proc
  proc_eq(0) = 0
  do ii=1,sml_plane_totalpe
     proc_eq(ii) = proc_eq(ii) + proc_eq(ii-1)
  end do
  ! set global PETSc equation numbers
  allocate(a_ts%petscloc_xgc(npetscloc))
  allocate(a_ts%petsc_xgc(nreal2D))
  a_ts%petsc_xgc = -2
  nreal2D = 0
  iloc = 0
  do ixgc=1,grid%nnode
     if(is_inside(ixgc,bc)) then
        nreal2D = nreal2D + 1
        proc = xgc_proc(nreal2D)
        a_ts%xgc_petsc(ixgc) = proc_eq(proc) + my0 ! zero based offset, 3D now, block index
        proc_eq(proc) = proc_eq(proc) + 1  ! increment, current one based eq on proc
        if ( a_ts%petsc_xgc( proc_eq(proc) ) /= -2 ) stop 'petsc_xgc /= -2'
        a_ts%petsc_xgc( proc_eq(proc) ) = ixgc - 1 
        if (proc==sml_plane_mype) then
           iloc = iloc + 1
           a_ts%petscloc_xgc(iloc) = ixgc - 1
        end if
     else 
        if (a_ts%xgc_petsc(ixgc) /= -1) stop 'a_ts%xgc_petsc(ixgc) /= -1'
     end if
  end do
  if (iloc.ne.npetscloc) stop 'iloc.ne.npetscloc'
  call ISRestoreIndicesF90(is2,xgc_proc,ierr);CHKERRQ(ierr)
  call ISDestroy(is2,ierr);CHKERRQ(ierr)
  deallocate(alpha,beta)

  ! create new 3D matrix for one field
  call MatCreate(sml_comm,AA,ierr);CHKERRQ(ierr) ! global 3D matrix, start with sml_comm
  call PetscObjectGetComm(AA,w_comm,ierr);CHKERRQ(ierr)
  call MatSetSizes(AA,npetscloc,npetscloc,PETSC_DECIDE,PETSC_DECIDE,ierr);CHKERRQ(ierr)
  call MatSetType(AA,MATAIJ,ierr);CHKERRQ(ierr)
  call MatSetup(AA,ierr)
  call MatGetOwnershipRange(AA,low,high,ierr);CHKERRQ(ierr)
  allocate(d_nnz(npetscloc),o_nnz(npetscloc))
  call getNNZ(grid,npetscloc,low,high,d_nnz,o_nnz,a_ts%xgc_petsc,ierr);CHKERRQ(ierr) ! low/high refer to xgc_petsc - global scalar
  call MatSeqAIJSetPreallocation(AA,PETSC_NULL_INTEGER,d_nnz,ierr)  ;CHKERRQ(ierr)
  call MatMPIAIJSetPreallocation(AA,PETSC_NULL_INTEGER,d_nnz,PETSC_NULL_INTEGER,o_nnz,ierr);CHKERRQ(ierr)
  deallocate(d_nnz,o_nnz)
  call MatSetup(AA,ierr)
  a_ts%mass= AA
  if (npetscloc.ne.high-low) stop 'npetscloc.ne.high-low'

  return
end subroutine ts_matrix_init

! ------------------------------------------------------------------------
!  subroutine addParMat
!
!  Input Parameters:
!     grid -
!     psn - 
!     bc - 
!     add_b0_terms -
!     xgc_petsc
!     scales
!     bs - block size (3)
!     idx - zero bases block row to put data in
!     jdx - zero bases block column to put data in
!     laplace_flag - if true, second parallel derivative is evaluated
!                    instead of gradient
!
!  Input/Output Parameter:
!     a_mat - useing ADD_VALUES, so can have stuff in it
!
!  Output Parameter:
!     ierr - error code
!
#undef __FUNCT__
#define __FUNCT__ "addParMat"
!> Adds the parallel gradient or parallel Laplacian to a_mat
!> only one of add_b0_terms, laplace_flag, and add_t0_n0 should be .true.
subroutine addParMat(grid,psn,a_bc,add_b0_terms,xgc_petsc,a_mat,scales,bs,idx,jdx,laplace_flag,add_t0_n0,ierr)
  use grid_class
  use psn_class
  use sml_module,only:sml_mype,sml_nphi_total,sml_Bt_sign,sml_e_charge,sml_plane_index
  use xgc_ts_module
  use boundary_class
  use eq_module
  use ptl_module,only:ptl_mass
  implicit none
#include <petscversion.h>
#if PETSC_VERSION_LT(3,6,0)
#include <finclude/petsc.h90>
#else
#include <petsc/finclude/petsc.h90>
#endif
  type(grid_type),intent(in)::grid
  type(psn_type)::psn
  type(boundary2_type),intent(in)::a_bc
  logical,intent(in)::add_b0_terms
  Mat::a_mat
  PetscErrorCode,intent(out)::ierr
  PetscInt,intent(in)::xgc_petsc(grid%nnode)
  PetscScalar,intent(in)::scales
  PetscInt,intent(in)::idx,jdx,bs
  logical, intent(in) :: laplace_flag, add_t0_n0
  !
  PetscInt::ip(1),ix,jp(1),jx,low,high,npetscplane,nd,sign1,sign2,dir,nodes(3),npetsctot
  PetscInt,parameter::ione=1
  PetscScalar::v1(1)
  PetscReal,parameter::realone=1.d0
  real (8) :: p(3), w(-1:1), dx_inv
  real (8) :: B1,B1_eB2  ! B-field at ix,  B1/(B-field times e) at corresponding field following position
  real (8) :: t0_n0
  
  call MatGetOwnershipRange(a_mat,low,high,ierr);CHKERRQ(ierr)
  call MatGetSize(a_mat,npetsctot,PETSC_NULL_INTEGER,ierr);CHKERRQ(ierr)

  if (mod(npetsctot,sml_nphi_total)/=0) stop 'addParMat: mod(npetsctot,sml_nphi_total)/=0'
  npetscplane = npetsctot/sml_nphi_total ! number of real eqs on plane

  if(sml_bt_sign>0D0) then
    sign1=1
  else
    sign1=-1
  endif

  B1=1D0
  B1_eB2=1D0
  t0_n0=1D0

  do ix=1,grid%nnode

     ip(1) = xgc_petsc(ix)*bs + idx ! can be -1*3 + 0,1,2
     if (ip(1).ge.low .and. ip(1).lt.high) then

        if(add_b0_terms) B1=grid%bfield(4,ix)
        if(add_t0_n0) t0_n0=grid%tempe(ix)/grid%dene(ix)
        if (.not. laplace_flag) then
          !rh More accurate parallel derivative with arbitrary
          !rh distances between planes -1, 0 and +1
          !rh d_parallel(X) = w(-1)*X(-1) + w(0)*X(0) + w(1)*X(1)
          !rh w(-1) = - dl(1)/(dl(-1)*(dl(-1)+dl(1)))
          !rh w(-1) =  (dl(1)-dl(-1))/(dl(-1)*dl(1))
          !rh w(1)  = + dl(-1)/(dl(1)*(dl(-1)+dl(1)))
          w(-1) = -psn%ff_1dp_dx(ix,1)/(psn%ff_1dp_dx(ix,0)*sum(psn%ff_1dp_dx(ix,:)))
          w(0)  = (psn%ff_1dp_dx(ix,1)-psn%ff_1dp_dx(ix,0))/(psn%ff_1dp_dx(ix,0)*psn%ff_1dp_dx(ix,1))
          w(1)  = psn%ff_1dp_dx(ix,0)/(psn%ff_1dp_dx(ix,1)*sum(psn%ff_1dp_dx(ix,:)))
        else
          !rh Second order parallel derivative with arbitrary
          !rh distances between planes -1, 0 and +1
          !rh d_parallel^2(X) = w(-1)*X(-1) + w(0)*X(0) + w(1)*X(1)
          !rh w(-1) = + 1/(2 * dl(-1)*(dl(-1)+dl(1)))
          !rh w(-1) = - 1/(2 * dl(-1)*dl(1))
          !rh w(1)  = + 1/(2 * dl(1)*(dl(-1)+dl(1)))
          w(-1) =  1D0/(2D0*psn%ff_1dp_dx(ix,0)*sum(psn%ff_1dp_dx(ix,:)))
          w(0)  = -1D0/(2D0*psn%ff_1dp_dx(ix,0)*psn%ff_1dp_dx(ix,1))
          w(1)  =  1D0/(2D0*psn%ff_1dp_dx(ix,1)*sum(psn%ff_1dp_dx(ix,:)))
        endif

        !rh NOTE: psn%ff_1dp is oriented along phi-direction like the petsc planes
        !rh index=0 <==> next plane in -phi-direction
        !rh index=1 <==> next plane in +phi-direction
        !rh ==> sign2=sign*dir ==> sign*sign2=dir
        do dir=-1,1
           sign2=sign1*dir
           if (dir==-1) then
             !rh plane -1: phi=phi0-delta_phi
             nodes=grid%nd(:,psn%ff_1dp_tr(ix,0))
             p=psn%ff_1dp_p(:,ix,0)
           elseif (dir==0) then
             !rh plane 0: phi=phi0
             nodes(:)=ix
             p(:)=(/ 1D0, 0D0, 0D0 /)
           elseif (dir==1) then
             !rh plane +1: phi=phi0+delta_phi
             nodes=grid%nd(:,psn%ff_1dp_tr(ix,1))
             p=psn%ff_1dp_p(:,ix,1)
           endif
           !rh *sign makes it derivative along the magnetic field
           p(:)=p(:)*w(dir)*sign1

           do nd=1,3
              if (dir==0 .and. nd .gt. 1) cycle ! mid point not used if nd>1, so 2-point (?)
              jx = nodes(nd)
              if(add_b0_terms) B1_eB2=B1/(grid%bfield(4,jx)*sml_e_charge)

              jp(1) = xgc_petsc(jx) ! can be -1
              if (jp(1).ge.0) then
                 jp(1) = jp(1)*bs + jdx 

                 jp(1) = jp(1) + npetscplane*sign1*sign2 ! forward/backward
                 jp(1) = mod(jp(1)+npetsctot,npetsctot) ! periodic boundary
                 v1(1) = scales*p(nd)*B1_eB2*t0_n0
                 if (v1(1)/=0d0) then
                    !if (ip(1).eq.6 .and. sml_plane_index==0) write(6,'(A,I3,A,I8,I8,ES12.4E2)'),'     [',sml_mype,'] addParMat: add to ',ip(1),jp(1),v1(1)
                    !print *,'addParMat:',sml_mype,v1(1),ip(1),jp(1),low,high
                    call MatSetValues(a_mat,ione,ip(1),ione,jp(1),v1(1),ADD_VALUES,ierr);CHKERRQ(ierr)
!!$                    if (sml_plane_index==0 .and. jp(1).ge.low .and. jp(1).lt.high) then 
!!$                       print *,'addParMat: diagonal entry:[',sml_mype,'] v=',v1(1),'i=',ip(1),jp(1),', lo,hi=',low,high
!!$                    end if
                 end if
              end if
           enddo
        enddo
     endif
  enddo
end subroutine addParMat

! ------------------------------------------------------------------------
!  subroutine addXGCMat
!
!  Input Parameters:
!     grid -
!     psn - 
!     bc -
!     xgc_petsc
!     scales
!     bs - block size (3)
!     idx - zero bases block row to put data in
!     jdx - zero bases block column to put data in
!
!  Input/Output Parameter:
!     a_mat - useing ADD_VALUES, so can have stuff in it
!     xgc_mat - input matrix in XGC sparse matrix representation
!
!  Output Parameter:
!     ierr - error code
!
#undef __FUNCT__
#define __FUNCT__ "addXGCMat"
!> Adds any matrix in the XGC sparse matrix format to a_mat
subroutine addXGCMat(grid,psn,a_bc,xgc_petsc,a_mat,xgc_mat,scales,bs,idx,jdx,ierr)
  use grid_class
  use mat_class
  use psn_class
  use sml_module,only:sml_mype,sml_nphi_total,sml_Bt_sign,sml_e_charge,sml_plane_index
  use xgc_ts_module
  use boundary_class
  use eq_module
  use ptl_module,only:ptl_mass
  implicit none
#include <petscversion.h>
#if PETSC_VERSION_LT(3,6,0)
#include <finclude/petsc.h90>
#else
#include <petsc/finclude/petsc.h90>
#endif
  type(grid_type),intent(in)::grid
  type(psn_type)::psn
  type(boundary2_type),intent(in)::a_bc
  type(mat_type), intent(in) :: xgc_mat
  Mat::a_mat
  PetscErrorCode,intent(out)::ierr
  PetscInt,intent(in)::xgc_petsc(grid%nnode)
  PetscScalar,intent(in)::scales
  PetscInt,intent(in):: idx, jdx, bs
  !
  PetscInt::ip(1),ix,jp(1),jx,wx,low,high,npetscplane,npetsctot
  PetscInt,parameter::ione=1
  PetscScalar::v1(1)
  PetscReal,parameter::realone=1.d0

  call MatGetOwnershipRange(a_mat,low,high,ierr);CHKERRQ(ierr)
  call MatGetSize(a_mat,npetsctot,PETSC_NULL_INTEGER,ierr);CHKERRQ(ierr)

  if (mod(npetsctot,sml_nphi_total)/=0) stop 'addParMat: mod(npetsctot,sml_nphi_total)/=0'
  npetscplane = npetsctot/sml_nphi_total ! number of real eqs on plane

  do ix=1,grid%nnode

    ip(1) = xgc_petsc(ix)*bs + idx ! can be -1*3 + 0,1,2
    if (ip(1).ge.low .and. ip(1).lt.high) then

      do wx=1,xgc_mat%nelement(ix)
        jx = xgc_mat%eindex(wx,ix)
        jp(1)=xgc_petsc(jx)
        if (jp(1) .ge. 0) then
          ! We stay on the plane --->
          jp(1) = jp(1)*bs + jdx
          v1(1) = xgc_mat%value(wx,ix)*scales
          if (jx==ix) then
            ! Add curl(nb).grad(j0/eB) to the diagonal
            v1(1) = v1(1) + grid%tearing_drive2(ix)*scales
          endif
          call MatSetValues(a_mat,ione,ip(1),ione,jp(1),v1(1),ADD_VALUES,ierr);CHKERRQ(ierr)
        endif
      enddo
    endif
  enddo

end subroutine addXGCMat

!
#undef __FUNCT__
#define __FUNCT__ "MatXPY"
subroutine MatXPY(dest,src,bs,idx,jdx,ierr)
  use sml_module,only:sml_mype,sml_plane_index
  implicit none
#include <petscversion.h>
#if PETSC_VERSION_LT(3,6,0)
#include <finclude/petsc.h90>
#else
#include <petsc/finclude/petsc.h90>
#endif
  Mat::dest,src
  PetscInt::bs,idx,jdx
  PetscErrorCode::ierr
  !
  PetscInt low,high,ncols,ii,ione,idi(1),jj
  PetscInt,allocatable::cols(:),idj(:)
  PetscScalar,allocatable::vals(:)

  ione = 1
  allocate(cols(1024),vals(1024),idj(1024))
  call MatGetOwnershipRange(src,low,high,ierr);CHKERRQ(ierr)

  do ii=low,high-1
     idi(1) = ii*bs + idx
     call MatGetRow(src,ii,ncols,PETSC_NULL_INTEGER,PETSC_NULL_SCALAR,ierr);CHKERRQ(ierr)
     if (ncols.gt.size(cols)) then
        deallocate(cols,vals,idj)
        allocate(cols(ncols),vals(ncols),idj(ncols))
        if(sml_mype==0) print *, ' MatXPY: realloc buffers to ',ncols
     end if
     call MatRestoreRow(src,ii,ncols,PETSC_NULL_INTEGER,PETSC_NULL_SCALAR,ierr);CHKERRQ(ierr)
     call MatGetRow(    src,ii,ncols,cols,vals,ierr);CHKERRQ(ierr)

     do jj=1,ncols
        idj(jj) =  cols(jj)*bs + jdx
        if (idj(jj).lt.0) stop
        !if (idi(1)==6 .and. sml_plane_index==0) write(6,'(A,I3,A,I8,I8,ES12.4E2)'),'     [',sml_mype,'] MatXPY: add to ',idi(1),idj(jj),vals(jj)
     end do
     call MatSetValues(dest,ione,idi,ncols,idj,vals,ADD_VALUES,ierr);CHKERRQ(ierr)
     call MatRestoreRow(src,ii,ncols,cols,vals,ierr);CHKERRQ(ierr)
  end do
  deallocate(cols,vals,idj)
end subroutine MatXPY
! ------------------------------------------------------------------------
!  ts_init - Setup time steper. 
!
!  Input Parameters:
!     grid - grid
!     a_bc - bc object
!
!  Output Parameter:
!     a_ts - time step object
! ------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "ts_init"
subroutine ts_init(a_grid,a_psn,a_bc,a_ts,ierr2)
  use grid_class
  use psn_class
  use sml_module,only:sml_mype,sml_plane_index,sml_intpl_mype,sml_e_charge,sml_plane_mype,  &
                      sml_plane_totalpe,sml_mu0,sml_eta,sml_dt,sml_nphi_total,sml_plane_comm, &
                      sml_comm, sml_lumped_mass, sml_use_scaling, sml_n1_diff_coef_perp, &
                      sml_n1_diff_coef_para, sml_hyb_tearing_test
  use xgc_ts_module
  use boundary_class
  use eq_module
  use sml_module, only: sml_hyb_bgraddpe_on
  use ptl_module,only:ptl_mass
  implicit none
#include <petscversion.h>
#if PETSC_VERSION_LT(3,6,0)
#include <finclude/petsc.h90>
#else
#include <petsc/finclude/petsc.h90>
#endif
  type(grid_type),intent(in)::a_grid
  type(psn_type)::a_psn
  type(boundary2_type),intent(in)::a_bc
  type(xgc_ts)::a_ts
  integer,intent(out)::ierr2
  !
  PetscErrorCode::ierr
  integer::ixgc,ix
  PetscInt::npetsc,itr,nd(3),npetscloc,low,high,ii,jj,kk,npetscplane
  PetscReal::half
  PetscScalar::val
  real (kind=8),allocatable::alpha(:),beta(:),teev(:),den(:),b2(:)
  real (kind=8) :: x_center(2),psi
  !DM::dmall,dmn1,dma,dmphi
  SNES::snes
  KSP::ksp
  PC::pc
  PetscInt,parameter::ione=1,izero=0,itwo=2,ithree=3
  PetscScalar :: eps_para
  Vec::xvec,globalvec
  PetscViewer::viewer
  external FormIJacobian,FormIFunction,test_fem, RHSFunction
  real (kind=8),external::psi_interpol,b_interpol
  IS::is,is2
  MPI_Comm::w_comm
  PetscInt,allocatable::d_nnz(:),o_nnz(:),d_nnz3(:),o_nnz3(:),d_nnzgp(:),o_nnzgp(:)
  PCCompositeType::pcctype
  Character(len=32)::strings
  PetscBool::flg
  Mat::gradparmat,lapparmat,mn1_del2_mu,A01,mat
  real (kind=8),dimension(a_grid%nnode)::n1,apar,phi
  real (kind=8)::rr,zz,r_minor, v_a, den_axis
  VecScatter::scat
  PetscRandom::rctx
  Vec::vec1,nullvec

  ! cache size for dynamic sizing of arrays
  a_ts%nnode = a_grid%nnode

  !rh This is safe -->
  pcctype = PC_COMPOSITE_MULTIPLICATIVE
  if (sml_mype==0) write(6,*) sml_mype,'TS.ts_init: make MULTIPLICATIVE FieldSplit type'

  !rh This often causes problems --->
!#if (PETSC_VERSION_GE(3,6,0) && PETSC_VERSION_RELEASE)
!  call PCFieldSplitGetType(pc,pcctype,ierr);CHKERRQ(ierr) ! does not work in 3.5.4 
!#else
!  call PetscOptionsGetString(PETSC_NULL_CHARACTER,"-adv_pc_fieldsplit_type",strings,flg,ierr)
!  CHKERRQ(ierr)
!  if (flg) then
!     if (strings == 'multiplicative') then
!        pcctype = PC_COMPOSITE_MULTIPLICATIVE
!        if (sml_mype==0) write(6,*) sml_mype,'TS.ts_init: make MULTIPLICATIVE FieldSplit type'
!     else if  (strings == 'schur') then
!        pcctype = PC_COMPOSITE_SCHUR
!        if (sml_mype==0) write(6,*) sml_mype,'TS.ts_init: make SCHUR FieldSplit type'
!     else
!        if (sml_mype==0) write(6,*) sml_mype,'TS.ts_init: warning unknown field split type:',strings
!        pcctype = PC_COMPOSITE_MULTIPLICATIVE
!     end if
!  else
!     if (sml_mype==0) write(6,*) sml_mype,'TS.ts_init: warning no -adv_pc_fieldsplit_type found'
!     pcctype = PC_COMPOSITE_MULTIPLICATIVE
!  end if
!#endif

  ! debug
  if (.false.) then
     if (sml_intpl_mype==0) call test_fem(sml_plane_comm)
     call mpi_barrier(sml_comm,ierr)
     stop 'debugging'
  end if

  !if (sml_mype==0) write(6,*) sml_mype,'TS.ts_init:',a_grid%nnode,' XGC nodes per plane'
  call PetscLogEventRegister('ADV Init',0,ADVInitEvent,ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('  ADV Mat Init',0,ADVMatrixInitEvent,ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('ADV Solve',0,ADVSolveEvent,ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('  ADV Form Jac',0,ADVFormIJacEvent,ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('  ADV Form LHS',0,ADVFormIFuncEvent,ierr);CHKERRQ(ierr)
  call PetscLogEventRegister('  ADV RHS',0,ADVFormRHSEvent,ierr);CHKERRQ(ierr)

  call PetscLogEventBegin(ADVMatrixInitEvent,ierr)
  ! create 3D mass matrix (w/o data) and its data layout
  call ts_matrix_init(a_ts,a_grid,a_bc,ierr);ierr2=ierr;CHKERRQ(ierr)
  call PetscLogEventEnd(ADVMatrixInitEvent,ierr)

  call PetscLogEventBegin(ADVInitEvent,ierr)

  ! make mass matrix - used for RHS and scaling
  call PetscObjectGetComm(a_ts%mass,w_comm,ierr);CHKERRQ(ierr)
  call MatGetLocalSize(a_ts%mass,npetscloc,PETSC_NULL_INTEGER,ierr);CHKERRQ(ierr)
  call MatGetSize(a_ts%mass,PETSC_NULL_INTEGER,npetsc,ierr);CHKERRQ(ierr)
  call MatGetOwnershipRange(a_ts%mass,low,high,ierr);CHKERRQ(ierr)
  if (npetscloc.ne.high-low) stop 'ts_init: npetscloc.ne.high-low'
#if PETSC_VERSION_LT(3,6,0)
  ! will v3.5 work for 3-field?
  call MatGetVecs(a_ts%mass,PETSC_NULL_OBJECT,a_ts%mass_inv_vec,ierr);CHKERRQ(ierr)
  stop 'you should use PETSc > v3.6. v3.5 might work ...'
#else
  call MatCreateVecs(a_ts%mass,PETSC_NULL_OBJECT,a_ts%mass_inv_vec,ierr);CHKERRQ(ierr)
#endif
  ! create mass matrix - a_ts%mass.  Diagonal matrix
  call MatZeroEntries(a_ts%mass,ierr)  
  if (sml_lumped_mass) then
     ! put a_grid%node_area into the diagonal of a_ts%mass
     do jj = 1,a_grid%nnode
        ii = a_ts%xgc_petsc(jj)
        if (ii.ge.low .and. ii.lt.high) then
           call MatSetValues(a_ts%mass,        ione,ii,ione,ii,a_grid%node_area(jj),INSERT_VALUES,ierr);CHKERRQ(ierr)
           call VecSetValues(a_ts%mass_inv_vec,ione,ii,        a_grid%node_area(jj),INSERT_VALUES,ierr);CHKERRQ(ierr)
        end if
     end do
     call MatAssemblyBegin(a_ts%mass,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
     call MatAssemblyEnd(  a_ts%mass,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
     call VecAssemblyBegin(a_ts%mass_inv_vec,ierr);CHKERRQ(ierr)
     call VecAssemblyEnd(  a_ts%mass_inv_vec,ierr);CHKERRQ(ierr)
     call VecReciprocal(a_ts%mass_inv_vec,ierr);CHKERRQ(ierr)
  else
     stop '3 field solver requires lumped mass'
     beta=1.D0
     alpha=0.D0 
     call helm_matrix(a_ts%mass,alpha,beta,a_grid,a_bc,.false.,a_ts%xgc_petsc,ierr)
     CHKERRQ(ierr)
  end if

  ! create Jacobian matrix
  allocate(d_nnz(0:npetscloc-1),o_nnz(0:npetscloc-1),d_nnzgp(0:npetscloc-1),o_nnzgp(0:npetscloc-1))
  call getNNZ(a_grid,npetscloc,low,high,d_nnz,o_nnz,a_ts%xgc_petsc,ierr);CHKERRQ(ierr) ! nnz of helm perp
  ! create Mats
  call MatCreate(w_comm,a_ts%FJacobian,ierr);CHKERRQ(ierr)
  call MatSetSizes(a_ts%FJacobian,ithree*npetscloc,ithree*npetscloc,PETSC_DECIDE,PETSC_DECIDE,ierr);CHKERRQ(ierr)
  call MatSetType(a_ts%FJacobian,MATAIJ,ierr);CHKERRQ(ierr)
  call MatSetBlockSize(a_ts%FJacobian,ithree,ierr);CHKERRQ(ierr)
  call MatCreate(w_comm,gradparmat,ierr);CHKERRQ(ierr)
  call MatSetSizes(gradparmat,ione*npetscloc,ione*npetscloc,PETSC_DECIDE,PETSC_DECIDE,ierr);CHKERRQ(ierr)
  call MatSetType(gradparmat,MATAIJ,ierr);CHKERRQ(ierr)
  call MatSetBlockSize(gradparmat,ione,ierr);CHKERRQ(ierr)
  !rh New matrix for parallel laplacian
  call MatCreate(w_comm,lapparmat,ierr);CHKERRQ(ierr)
  call MatSetSizes(lapparmat,ione*npetscloc,ione*npetscloc,PETSC_DECIDE,PETSC_DECIDE,ierr);CHKERRQ(ierr)
  call MatSetType(lapparmat,MATAIJ,ierr);CHKERRQ(ierr)
  call MatSetBlockSize(lapparmat,ione,ierr);CHKERRQ(ierr)
  call MatCreate(w_comm,mn1_del2_mu,ierr);CHKERRQ(ierr)
  call MatSetSizes(mn1_del2_mu,ione*npetscloc,ione*npetscloc,PETSC_DECIDE,PETSC_DECIDE,ierr);CHKERRQ(ierr)
  call MatSetType(mn1_del2_mu,MATAIJ,ierr);CHKERRQ(ierr)
  call MatSetBlockSize(mn1_del2_mu,ione,ierr);CHKERRQ(ierr)
  call MatCreate(w_comm,a_ts%muDel2,ierr);CHKERRQ(ierr)
  call MatSetSizes(a_ts%muDel2,npetscloc,npetscloc,PETSC_DECIDE,PETSC_DECIDE,ierr);CHKERRQ(ierr)
  call MatSetType(a_ts%muDel2,MATAIJ,ierr);CHKERRQ(ierr)
  call MatSetBlockSize(a_ts%muDel2,ione,ierr);CHKERRQ(ierr)
  ! form Jacobean for fast terms:
  !---------
  ! n_1 : diff M^-1 Del_perp^2    B_0 del_|| (1/eB_0) (mu M)^-1 Del_perp^2    0
  ! A_||: 0                                      -eta (mu M)^-1 Del_perp^2    del_||
  ! phi:  M                                                              0    del n_0 m_i/B^2 del_perp
  !---------
  ! sizes
  allocate(d_nnz3(0:ithree*npetscloc-1),o_nnz3(0:ithree*npetscloc-1))
  do ii=0,npetscloc-1,1
     !rh One additional entry for the parallel gradient
     d_nnz3(ithree*ii + 0) = d_nnz(ii) + 9 + 4*size(a_grid%tr_node,1)+2 ! M^-1 * Del_perp^2_FE + kink-term
     d_nnz3(ithree*ii + 1) = d_nnz(ii) + 1 ! nothing for daig of grad_par because hits diagonal only
     d_nnz3(ithree*ii + 2) = 2*d_nnz(ii)

     o_nnz3(ithree*ii + 0) = o_nnz(ii) + 6*(o_nnz(ii)+d_nnz(ii))+6 ! diag + del_|| * Del_perp^2 + par. diffusion
     o_nnz3(ithree*ii + 1) = o_nnz(ii) + 6 + 6 ! diag + grad_par
     o_nnz3(ithree*ii + 2) = 2*o_nnz(ii)

     d_nnzgp(ii + 0) = 1 ! we dump onto the diagonal at the edge when we blow it
     o_nnzgp(ii + 0) = 6 ! del_|| * Del_perp^2
  end do

  call MatSeqAIJSetPreallocation(a_ts%FJacobian,PETSC_NULL_INTEGER,d_nnz3,ierr);CHKERRQ(ierr)
  call MatMPIAIJSetPreallocation(a_ts%FJacobian,PETSC_NULL_INTEGER,d_nnz3,PETSC_NULL_INTEGER,o_nnz3,ierr);CHKERRQ(ierr)
  call MatSeqAIJSetPreallocation(gradparmat,PETSC_NULL_INTEGER,d_nnzgp,ierr);CHKERRQ(ierr)
  call MatMPIAIJSetPreallocation(gradparmat,PETSC_NULL_INTEGER,d_nnzgp,PETSC_NULL_INTEGER,o_nnzgp,ierr);CHKERRQ(ierr)
  !rh New matrix for parallel laplacian
  call MatSeqAIJSetPreallocation(lapparmat,PETSC_NULL_INTEGER,d_nnzgp,ierr);CHKERRQ(ierr)
  call MatMPIAIJSetPreallocation(lapparmat,PETSC_NULL_INTEGER,d_nnzgp,PETSC_NULL_INTEGER,o_nnzgp,ierr);CHKERRQ(ierr)
  call MatSeqAIJSetPreallocation(mn1_del2_mu,PETSC_NULL_INTEGER,d_nnz,ierr);CHKERRQ(ierr)
  call MatMPIAIJSetPreallocation(mn1_del2_mu,PETSC_NULL_INTEGER,d_nnz,PETSC_NULL_INTEGER,o_nnz,ierr);CHKERRQ(ierr)
  call MatSeqAIJSetPreallocation(mn1_del2_mu,PETSC_NULL_INTEGER,d_nnz,ierr);CHKERRQ(ierr)
  call MatMPIAIJSetPreallocation(mn1_del2_mu,PETSC_NULL_INTEGER,d_nnz,PETSC_NULL_INTEGER,o_nnz,ierr);CHKERRQ(ierr)
  call MatSeqAIJSetPreallocation(a_ts%muDel2,PETSC_NULL_INTEGER,d_nnz,ierr);CHKERRQ(ierr)
  call MatMPIAIJSetPreallocation(a_ts%muDel2,PETSC_NULL_INTEGER,d_nnz,PETSC_NULL_INTEGER,o_nnz,ierr);CHKERRQ(ierr)
  deallocate(d_nnz3,o_nnz3,d_nnz,o_nnz,d_nnzgp,o_nnzgp)

  ! setup n_0,b^2,Te at centroids of triangles for Poisson solver
  allocate(alpha(a_grid%ntriangle),beta(a_grid%ntriangle),teev(a_grid%ntriangle),den(a_grid%ntriangle),b2(a_grid%ntriangle))
  do itr=1, a_grid%ntriangle
     nd=a_grid%nd(:,itr)
     x_center=(a_grid%x(:,nd(1))+a_grid%x(:,nd(2))+a_grid%x(:,nd(3)))/3D0
     psi=psi_interpol(x_center(1),x_center(2),0,0)
     teev(itr)=eq_ftn(psi,x_center(1),x_center(2),eq_tempe)
     b2(itr)=b_interpol(x_center(1),x_center(2),0D0)**2
     den(itr)=eq_ftn(psi,x_center(1),x_center(2),eq_den)
  enddo
  beta = 0D0

  call MatSetOption(a_ts%FJacobian,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,ierr);CHKERRQ(ierr)
  call MatSetOption(a_ts%FJacobian,MAT_IGNORE_ZERO_ENTRIES,PETSC_FALSE,ierr);CHKERRQ(ierr) ! just make sure

  if (sml_use_scaling .and. a_ts%scales(0)==1d0) then
    den_axis=eq_ftn(0D0,eq_axis_r,eq_axis_z,eq_den)
    v_a=eq_axis_b/sqrt(sml_mu0*den_axis*ptl_mass(1))
    a_ts%scales(0)=eq_axis_r/(den_axis*v_a)
  endif

  ! K(0,0): eps * M^-1 Del^2
  call check_point('Geherate perp. Laplacian for perp. diffusion')
  alpha = sml_n1_diff_coef_perp*a_ts%scales(0) ! 1.0
  call helm_matrix_mdof(a_ts%FJacobian,alpha,beta,a_grid,a_bc,.false.,a_ts%xgc_petsc,ithree,izero,izero,a_ts%mass_inv_vec,ierr)
  CHKERRQ(ierr)
  !rh Add parallel laplacian matrix: eps_para*d^2/dl_para^2
  eps_para=-sml_n1_diff_coef_para*a_ts%scales(0)
  call addParMat(a_grid,a_psn,a_bc,.false.,a_ts%xgc_petsc,a_ts%FJacobian,eps_para,ithree,izero,izero,.true.,.false.,ierr);CHKERRQ(ierr)

  ! K(0,1) : M^-1 Del_perp^2/mu -- adding negative to compensate for built in neg. in del^2
  alpha = -a_ts%scales(0)/sml_mu0 ! scaled by diag w/o mass
  call helm_matrix_mdof(mn1_del2_mu,alpha,beta,a_grid,a_bc,.false.,a_ts%xgc_petsc,ione,izero,izero,a_ts%mass_inv_vec,ierr)
  CHKERRQ(ierr)
  call MatAssemblyBegin(mn1_del2_mu,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(  mn1_del2_mu,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

  ! B_0 grad_|| 1/eB_0:  FD -- only scale one of these two terms
  call addParMat(a_grid,a_psn,a_bc,.true.,a_ts%xgc_petsc,gradparmat,1D0,ione,izero,izero,.false.,.false.,ierr);CHKERRQ(ierr)
  call MatAssemblyBegin(gradparmat,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(  gradparmat,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

  ! K(0,1): B_0 del_|| (1/eB_0) * (mu M)^-1 Del_perp^2
  call MatMatMult(gradparmat,mn1_del2_mu,MAT_INITIAL_MATRIX,PETSC_DEFAULT_REAL,A01,ierr);CHKERRQ(ierr)
  call MatDestroy(gradparmat,ierr);CHKERRQ(ierr)
  call MatDestroy(mn1_del2_mu,ierr);CHKERRQ(ierr) ! add later
  call MatXPY(a_ts%FJacobian,A01,ithree,izero,ione,ierr);CHKERRQ(ierr) ! K(0,1)   
  call MatDestroy(A01,ierr);CHKERRQ(ierr)

  if (sml_hyb_tearing_test) then
    ! K(0,1): Add kink-drive (b x grad(j0/eB)).grad(A_par)
    call addXGCMat(a_grid,a_psn,a_bc,a_ts%xgc_petsc,a_ts%FJacobian,a_grid%kink_mat,  &
                   a_ts%scales(0),ithree,izero,ione,ierr) ; CHKERRQ(ierr)
  endif

  if (sml_hyb_bgraddpe_on) then
    !rh K(1,0) -(T_0/n_0)*b.grad(delta_n_e): Needed for energy conservation
    !rh need to multiply by T_0/n_0
    call addParMat(a_grid,a_psn,a_bc,.false.,a_ts%xgc_petsc,a_ts%FJacobian,-a_ts%scales(1),ithree,ione,izero,.false.,.true.,ierr);CHKERRQ(ierr)
  endif

  ! make muDel2: M^-1 Del^2 / mu to compute Je = M^-1 Del^2 / mu * A -- adding negative^2 to compensate for built in neg. in del^2
  alpha = 1D0/sml_mu0
  call helm_matrix_mdof(a_ts%muDel2,alpha,beta,a_grid,a_bc,.false.,a_ts%xgc_petsc,ione,izero,izero,a_ts%mass_inv_vec,ierr)
  call MatAssemblyBegin(a_ts%muDel2,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(  a_ts%muDel2,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

  ! K(1,1): FD
  if (sml_use_scaling .and. a_ts%scales(1)==1d0) then
    den_axis=eq_ftn(0D0,eq_axis_r,eq_axis_z,eq_den)
     v_A=eq_axis_b/sqrt(sml_mu0*den_axis*ptl_mass(1))
     a_ts%scales(1) = 1D0/(eq_axis_b*v_A)
  end if
  alpha = sml_eta*a_ts%scales(1)/sml_mu0 ! 1.0/mu  -- adding negative to compensate for built in neg. in del^2
  call helm_matrix_mdof(a_ts%FJacobian,alpha,beta,a_grid,a_bc,.false.,a_ts%xgc_petsc,ithree,ione,ione,a_ts%mass_inv_vec,ierr)
  CHKERRQ(ierr)

  ! grad_|| (1,2): grad_|| - FD
  call addParMat(a_grid,a_psn,a_bc,.false.,a_ts%xgc_petsc,a_ts%FJacobian,a_ts%scales(1),ithree,ione,itwo,.false.,ierr);CHKERRQ(ierr)

  ! E potential Phi matrix (2,2): FE
  alpha = ptl_mass(1)/sml_e_charge*den/b2 ! mn_0/eB^2
  if (sml_use_scaling .and. a_ts%scales(2)==1d0) then
    den_axis=eq_ftn(0D0,eq_axis_r,eq_axis_z,eq_den)
    a_ts%scales(2) = 1D0/den_axis
  endif
  alpha = alpha*a_ts%scales(2)
  call helm_matrix_mdof(a_ts%FJacobian,alpha,beta,a_grid,a_bc,.false.,a_ts%xgc_petsc, ithree,itwo,itwo,PETSC_NULL_OBJECT,ierr)
  CHKERRQ(ierr)

  ! mass (2,0): M - FE - has negative for electrons built in
  beta = 1.D0*a_ts%scales(2)
  alpha = 0.D0
  call helm_matrix_mdof(a_ts%FJacobian,alpha,beta,a_grid,a_bc,.false.,a_ts%xgc_petsc,ithree,itwo,izero,PETSC_NULL_OBJECT,ierr)
  CHKERRQ(ierr)
  beta = 0.D0

  call MatAssemblyBegin(a_ts%FJacobian,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(  a_ts%FJacobian,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

  ! setup TS
  call TSCreate(w_comm,a_ts%ts,ierr);CHKERRQ(ierr)
  call TSSetProblemType(a_ts%ts,TS_LINEAR,ierr);CHKERRQ(ierr)
  call MatDuplicate(a_ts%FJacobian,MAT_COPY_VALUES,a_ts%FJacobian2,ierr);CHKERRQ(ierr)
  call TSSetIJacobian(a_ts%ts,a_ts%FJacobian2,a_ts%FJacobian2,FormIJacobian,a_ts,ierr)
  call TSSetRHSFunction(a_ts%ts,PETSC_NULL_OBJECT,RHSFunction,a_ts,ierr);CHKERRQ(ierr)
  call TSSetIFunction(a_ts%ts,PETSC_NULL_OBJECT,FormIFunction,a_ts,ierr);CHKERRQ(ierr) 

  ! set intial to Crank-Nicolson
  !call TSSetType(a_ts%ts,TSTHETA,ierr);CHKERRQ(ierr)
  !half = 0.5d0
  !call TSThetaSetTheta(a_ts%ts,half,ierr);CHKERRQ(ierr)
  !call TSThetaSetEndpoint(a_ts%ts,PETSC_TRUE,ierr);CHKERRQ(ierr)
  ! set options
  call TSSetOptionsPrefix(a_ts%ts,'adv_',ierr);CHKERRQ(ierr)
  call TSSetFromOptions(a_ts%ts,ierr);CHKERRQ(ierr)

  ! setup field split
  call TSGetSNES(a_ts,snes,ierr);CHKERRQ(ierr)
  call SNESGetKSP(snes,ksp,ierr);CHKERRQ(ierr)
  call KSPGetPC(ksp,pc,ierr);CHKERRQ(ierr)

  call MatGetLocalSize(a_ts%mass,npetscloc,PETSC_NULL_INTEGER,ierr);CHKERRQ(ierr)
  call MatGetOwnershipRange(a_ts%FJacobian,low,PETSC_NULL_INTEGER,ierr);CHKERRQ(ierr)
  jj = low + 0
  call ISCreateStride(w_comm,npetscloc,jj,ithree,a_ts%iss(0),ierr);CHKERRQ(ierr)
  jj = low + 1
  call ISCreateStride(w_comm,npetscloc,jj,ithree,a_ts%iss(1),ierr);CHKERRQ(ierr)
  jj = low + 2
  call ISCreateStride(w_comm,npetscloc,jj,ithree,a_ts%iss(2),ierr);CHKERRQ(ierr)
  !
  if (pcctype==PC_COMPOSITE_MULTIPLICATIVE) then
     call PCFieldSplitSetIS(pc,'phi', a_ts%iss(2),ierr);CHKERRQ(ierr) ! permute to make lower triangular
     call PCFieldSplitSetIS(pc,'apar',a_ts%iss(1),ierr);CHKERRQ(ierr)
     call PCFieldSplitSetIS(pc,'n1',  a_ts%iss(0),ierr);CHKERRQ(ierr)
  else
     if (pcctype/=PC_COMPOSITE_SCHUR) stop 'pcctype/=PC_COMPOSITE_SCHUR/PC_COMPOSITE_MULTIPLICATIVE'
     stop 'PC_COMPOSITE_SCHUR not supported'
     call MatGetLocalSize(a_ts%FJacobian,npetscloc,PETSC_NULL_INTEGER,ierr);CHKERRQ(ierr)
     if (mod(npetscloc,2)/=0) stop 'ts_init: mod(npetsctot,2) /= 0'
     jj = npetscloc/2 ! two parts
     allocate(d_nnz3(jj),o_nnz3(jj))
     jj = 1
     do ii=0,npetscloc-1,ithree
        d_nnz3(jj)   = ii + 0 + low
        d_nnz3(jj+1) = ii + 1 + low
        o_nnz3(jj)   = ii + 2 + low
        jj = jj + 2
     end do
     if (mod(npetscloc,2)/=0) stop 'ts_init: 2) mod(npetsctot,2) /= 0'
     jj = npetscloc/2 ! two parts
     call ISCreateGeneral(w_comm,jj,d_nnz3,PETSC_COPY_VALUES,is,ierr);CHKERRQ(ierr)
     call PCFieldSplitSetIS(pc,'dt',is,ierr);CHKERRQ(ierr)
     call ISDestroy(is,ierr);CHKERRQ(ierr)
     
     call ISCreateGeneral(w_comm,jj,o_nnz3,PETSC_COPY_VALUES,is,ierr);CHKERRQ(ierr)
     call PCFieldSplitSetIS(pc,'aux',is,ierr);CHKERRQ(ierr)
     call ISDestroy(is,ierr);CHKERRQ(ierr)  
     deallocate(d_nnz3,o_nnz3)
  end if
  ! set scatters
#if PETSC_VERSION_LT(3,6,0)
  call MatGetVecs(a_ts%FJacobian,PETSC_NULL_OBJECT,globalvec,ierr);CHKERRQ(ierr)
#else
  call MatCreateVecs(a_ts%FJacobian,PETSC_NULL_OBJECT,globalvec,ierr);CHKERRQ(ierr)
#endif
  call VecCreateSeq(PETSC_COMM_SELF,a_ts%nnode,xvec,ierr);CHKERRQ(ierr) ! big dumb 2D PETSc vector for XGC vector

  ! forward scatter into PETSc
  call MatGetLocalSize(a_ts%mass,npetscloc,PETSC_NULL_INTEGER,ierr);CHKERRQ(ierr)
  call ISCreateGeneral(PETSC_COMM_SELF,npetscloc,a_ts%petscloc_xgc,PETSC_COPY_VALUES,is,ierr);CHKERRQ(ierr) ! is: x_petscloc(n)
  deallocate(a_ts%petscloc_xgc)
  call VecScatterCreate(xvec,is,globalvec,a_ts%iss(0),a_ts%to_petsc(0),ierr);CHKERRQ(ierr)   ! forward scatter into PETSc
  call VecScatterCreate(xvec,is,globalvec,a_ts%iss(1),a_ts%to_petsc(1),ierr);CHKERRQ(ierr)   ! forward scatter into PETSc
  call VecScatterCreate(xvec,is,globalvec,a_ts%iss(2),a_ts%to_petsc(2),ierr);CHKERRQ(ierr) ! forward scatter into PETSc
  call ISDestroy(is,ierr);CHKERRQ(ierr)

  ! scatter object from PETSc
  call MatGetSize(a_ts%FJacobian,npetscloc,PETSC_NULL_INTEGER,ierr);CHKERRQ(ierr)
  if (mod(npetscloc,sml_nphi_total)/=0) stop 'ts_init: mod(npetsctot,sml_nphi_total) /= 0'
  npetscplane = npetscloc/sml_nphi_total
  if (mod(npetscloc,ithree)/=0) stop 'ts_init: mod(npetsctot,ithree) /= 0'
  low = sml_plane_index*npetscplane ! global start my0 for plane
  npetscplane = npetscplane/ithree   ! num real vertex on one whole plane

  call ISCreateGeneral(PETSC_COMM_SELF,npetscplane,a_ts%petsc_xgc,PETSC_COPY_VALUES,is2,ierr);CHKERRQ(ierr)
  jj = low + 0
  call ISCreateStride(PETSC_COMM_SELF,npetscplane,jj,ithree,is,ierr);CHKERRQ(ierr)
  call VecScatterCreate(globalvec,is,xvec,is2,a_ts%from_petsc(0),ierr);CHKERRQ(ierr) 
  call ISDestroy(is,ierr);CHKERRQ(ierr)
  jj = low + 1
  call ISCreateStride(PETSC_COMM_SELF,npetscplane,jj,ithree,is,ierr);CHKERRQ(ierr)
  call VecScatterCreate(globalvec,is,xvec,is2,a_ts%from_petsc(1),ierr);CHKERRQ(ierr) 
  call ISDestroy(is,ierr);CHKERRQ(ierr)
  jj = low + 2
  call ISCreateStride(PETSC_COMM_SELF,npetscplane,jj,ithree,is,ierr);CHKERRQ(ierr)
  call VecScatterCreate(globalvec,is,xvec,is2,a_ts%from_petsc(2),ierr);CHKERRQ(ierr) 
  call ISDestroy(is,ierr);CHKERRQ(ierr)

  ! create single value scatter to XGC for Je
#if PETSC_VERSION_LT(3,6,0)
  call MatGetVecs(a_ts%muDel2,PETSC_NULL_OBJECT,vec1,ierr);CHKERRQ(ierr)
#else
  call MatCreateVecs(a_ts%muDel2,PETSC_NULL_OBJECT,vec1,ierr);CHKERRQ(ierr)
#endif
  jj = sml_plane_index*npetscplane ! start of my plane
  call ISCreateStride(PETSC_COMM_SELF,npetscplane,jj,ione,is,ierr);CHKERRQ(ierr)
  call VecScatterCreate(vec1,is,xvec,is2,a_ts%from_petsc_single_value,ierr);CHKERRQ(ierr) 
  call ISDestroy(is,ierr);CHKERRQ(ierr)
  call VecDestroy(vec1,ierr);CHKERRQ(ierr)
  
  call ISDestroy(is2,ierr);CHKERRQ(ierr)
  call VecDestroy(xvec,ierr);CHKERRQ(ierr)
  
  call PetscLogEventEnd(ADVInitEvent,ierr)

  if (npetscloc.le.1000) then ! total vertices
     call  PetscViewerASCIIOpen(w_comm, 'Kmat.m', viewer,ierr);CHKERRQ(ierr)
     call  PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB,ierr);CHKERRQ(ierr)
     call  MatView(a_ts%FJacobian,viewer,ierr);CHKERRQ(ierr)
     call  PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  end if
  ! debug
  if (.false.) then
     val = 0d0
     call VecSet(globalvec,val,ierr);CHKERRQ(ierr)
     if (.false. ) then
        if (sml_mype==0) write(6,*) sml_mype,'TS. test del^2 size XGC = ',a_ts%nnode
        ! test: K(0,3) B_0 grad_par/B_0  ...
        call MatGetSubMatrix(a_ts%FJacobian,a_ts%iss(0),a_ts%iss(1),MAT_INITIAL_MATRIX,mat,ierr);CHKERRQ(ierr)
        call VecDestroy(globalvec,ierr)     
#if PETSC_VERSION_LT(3,6,0)
        call MatGetVecs(mat,vec1,globalvec,ierr);CHKERRQ(ierr) 
#else
        call MatCreateVecs(mat,vec1,globalvec,ierr);CHKERRQ(ierr)
#endif
        ! get rand and A*rand in PETSc
        call PetscRandomCreate(w_comm,rctx,ierr);CHKERRQ(ierr);
        call PetscRandomSetFromOptions(rctx,ierr);CHKERRQ(ierr);
        call VecSetRandom(globalvec,rctx,ierr);CHKERRQ(ierr); ! globalvec is rand
        call PetscRandomDestroy(rctx,ierr);CHKERRQ(ierr);
        call MatMult(mat, globalvec, vec1, ierr);CHKERRQ(ierr) ! vec1 has A*rand
        ! make map to XGC
        n1 = 0d0
        call VecCreateSeqWithArray(PETSC_COMM_SELF,ione,a_ts%nnode,n1,xvec,ierr);CHKERRQ(ierr)
        call ISCreateGeneral(PETSC_COMM_SELF,npetscplane,a_ts%petsc_xgc,PETSC_COPY_VALUES,is2,ierr);CHKERRQ(ierr)
        jj = low/3
        call ISCreateStride(PETSC_COMM_SELF,npetscplane,jj,ione,is,ierr);CHKERRQ(ierr)
        call VecScatterCreate(globalvec,is,xvec,is2,scat,ierr);CHKERRQ(ierr) 
        call ISDestroy(is,ierr);CHKERRQ(ierr)
        call ISDestroy(is2,ierr);CHKERRQ(ierr)
        call VecScatterBegin(scat,globalvec,xvec,INSERT_VALUES,SCATTER_FORWARD,ierr)
        call VecScatterEnd(  scat,globalvec,xvec,INSERT_VALUES,SCATTER_FORWARD,ierr) ! rand in n1
        call VecDestroy(xvec,ierr)
        call VecScatterDestroy(scat,ierr)
        call mpi_barrier(sml_comm,ierr)
     else if (.false.) then
        !  test mass
        if (sml_mype==0) write(6,*) sml_mype,'TS.ts_init: test mass == area'
#if PETSC_VERSION_LT(3,6,0)
        call MatGetVecs(a_ts%mass,globalvec,vec1,ierr);CHKERRQ(ierr)
#else
        call MatCreateVecs(a_ts%mass,globalvec,vec1,ierr);CHKERRQ(ierr)
#endif
        n1 = 0d0
        call VecCreateSeqWithArray(PETSC_COMM_SELF,ione,a_ts%nnode,n1,xvec,ierr);CHKERRQ(ierr)
        val = 1d0
        call VecSet(globalvec,val,ierr);CHKERRQ(ierr)
        call MatMult(a_ts%mass, globalvec, vec1, ierr);CHKERRQ(ierr) ! vec has lumped mass

        call ISCreateGeneral(PETSC_COMM_SELF,npetscplane,a_ts%petsc_xgc,PETSC_COPY_VALUES,is2,ierr);CHKERRQ(ierr)
        jj = low/3
        call ISCreateStride(PETSC_COMM_SELF,npetscplane,jj,ione,is,ierr);CHKERRQ(ierr)
        call VecScatterCreate(vec1,is,xvec,is2,scat,ierr);CHKERRQ(ierr) 
        call ISDestroy(is,ierr);CHKERRQ(ierr)
        call ISDestroy(is2,ierr);CHKERRQ(ierr)
        call VecScatterBegin(scat,vec1,xvec,INSERT_VALUES,SCATTER_FORWARD,ierr)
        call VecScatterEnd(  scat,vec1,xvec,INSERT_VALUES,SCATTER_FORWARD,ierr) ! mass in n1
        if (sml_mype==0) then 
           do jj=1,a_ts%nnode
              if (n1(jj) /= 0d0 .and. abs((a_grid%node_area(jj)-n1(jj))/n1(jj)).gt.1.d-14) &
                   write (*,'(A,I4,A,I4,A,ES10.2E2,ES10.2E2,ES10.2E2)') 'ERROR A*mass != node_area', &
                   jj,'/',a_ts%nnode,') relative diff=',abs((a_grid%node_area(jj)-n1(jj))/a_grid%node_area(jj)),&
                   a_grid%node_area(jj),n1(jj)
           end do
           write(6,*) sml_mype,'TS.ts_init: test mass == area done'
        end if
        !if (sml_intpl_mype==0) call test_fem(sml_plane_comm)
        call mpi_barrier(sml_comm,ierr)
     else if (.false.) then
        call  MatGetSubMatrix(a_ts%FJacobian,a_ts%iss(0),a_ts%iss(0),MAT_INITIAL_MATRIX,mat,ierr);CHKERRQ(ierr)
        call  PetscViewerASCIIOpen(w_comm, 'N1.m', viewer,ierr);CHKERRQ(ierr)
        call  PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB,ierr);CHKERRQ(ierr)
        call  MatView(mat,viewer,ierr);CHKERRQ(ierr)
        call  PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)     
        call  MatDestroy(mat,ierr);CHKERRQ(ierr)
        
        call  MatGetSubMatrix(a_ts%FJacobian,a_ts%iss(1),a_ts%iss(1),MAT_INITIAL_MATRIX,mat,ierr);CHKERRQ(ierr)
        call  PetscViewerASCIIOpen(w_comm, 'A.m', viewer,ierr);CHKERRQ(ierr)
        call  PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB,ierr);CHKERRQ(ierr)
        call  MatView(mat,viewer,ierr);CHKERRQ(ierr)
        call  PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)     
        call  MatDestroy(mat,ierr);CHKERRQ(ierr)

        call  MatGetSubMatrix(a_ts%FJacobian,a_ts%iss(2),a_ts%iss(2),MAT_INITIAL_MATRIX,mat,ierr);CHKERRQ(ierr)
        call  PetscViewerASCIIOpen(w_comm, 'Phi.m', viewer,ierr);CHKERRQ(ierr)
        call  PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB,ierr);CHKERRQ(ierr)
        call  MatView(mat,viewer,ierr);CHKERRQ(ierr)
        call  PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)     
        call  MatDestroy(mat,ierr);CHKERRQ(ierr)
     else
        if (sml_mype==0) write(6,*) sml_mype,'TS.ts_init: write del^2 (R - r)'
        call VecDestroy(globalvec,ierr)
        ! get del^2
#if PETSC_VERSION_LT(3,6,0)
        call MatGetVecs(a_ts%mass,globalvec,vec1,ierr);CHKERRQ(ierr)
#else
        call MatCreateVecs(a_ts%mass,globalvec,vec1,ierr);CHKERRQ(ierr)
#endif
        ! get r minor
        do ixgc = 1, a_grid%nnode
           rr = a_grid%x(1,ixgc)
           zz = a_grid%x(2,ixgc)
           r_minor=dsqrt((rr-eq_axis_r)**2+(zz-eq_axis_z)**2)
           !write(6,*) ixgc,') TS.ts_init: r_minor=',r_minor,is_inside(ixgc,a_bc)
           if (.not.is_inside(ixgc,a_bc)) exit
        enddo
        ! set f(r,theta) = r_minor - r
        do ixgc = 1, a_grid%nnode
           rr = a_grid%x(1,ixgc)
           zz = a_grid%x(2,ixgc)
           rr=dsqrt((rr-eq_axis_r)**2+(zz-eq_axis_z)**2)
           if (.not.is_inside(ixgc,a_bc)) then
              n1(ixgc) = 0d0
           else
              n1(ixgc) = r_minor - rr
           end if
           if(n1(ixgc) /= n1(ixgc)) stop 'Nan'
           !write(6,*) sml_mype,'   (R - r) = ',n1(ixgc)
        enddo
        ! write function
        call VecCreateSeqWithArray(PETSC_COMM_SELF,ione,a_ts%nnode,n1,vec1,ierr);CHKERRQ(ierr)
        if (sml_mype==0) then
           ! write XGC vector
           write(strings,'("linear_",I6.6,".m  ")') a_ts%nnode
           write(6,*) 'writing file ', strings
           call PetscViewerASCIIOpen(PETSC_COMM_SELF, strings, viewer,ierr);CHKERRQ(ierr)
           call PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB,ierr)
           CHKERRQ(ierr)
           call VecView(vec1,viewer,ierr);CHKERRQ(ierr)
           call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
        end if
        call VecDestroy(vec1,ierr)        
        call VecScatterDestroy(scat,ierr)
     end if
     call mpi_barrier(sml_comm,ierr)
     stop 'debugging'
  endif

  deallocate(alpha,beta,teev,den,b2)
  call VecDestroy(globalvec,ierr)

  return
end subroutine ts_init

! ------------------------------------------------------------------------
!  subroutine ts_solve
!
!  Input Parameters:
!     a_ts - the user context
!     dt - total (ion) time step 
!
!  Output Parameter:
!     ierr - error code
!     phi - electric field potential (constraint equation, also 'advanced')
!     Je - parallel current
!
!  In/Output Parameter:
!     n1 - electron fluid density (advanced with TS)
!     apar - parallel magnetic field potential (advanced with TS)
!
#undef __FUNCT__
#define __FUNCT__ "ts_solve"
subroutine ts_solve(a_ts, a_dt, n1, apar, phi, Je, ierr)
  use sml_module,only:sml_mype
  use xgc_ts_module
  use perf_monitor
  implicit none
#include <petscversion.h>
#if PETSC_VERSION_LT(3,6,0)
#include <finclude/petsc.h90>
#else
#include <petsc/finclude/petsc.h90>
#endif
  type(xgc_ts)::a_ts
  real (kind=8),dimension(a_ts%nnode)::n1,apar,phi,Je
  real (kind=8),intent(in)::a_dt
  PetscErrorCode::ierr
  !
  PetscInt::its
  PetscReal::ftime,dt,zero
  TSConvergedReason::reason
  Vec::XX
  Mat::mat
  Vec::vec1,vec2,vec3
  PetscScalar::a
  PetscReal::norm1,norm2

  !if (sml_mype == 0) write(6,*) sml_mype,'TS.ts_solve start ',a_ts%nnode,' XGC nodes'
  call PetscLogEventBegin(ADVSolveEvent,ierr)

#if PETSC_VERSION_LT(3,6,0)
  call MatGetVecs(a_ts%FJacobian2,PETSC_NULL_OBJECT,XX,ierr);CHKERRQ(ierr)
#else
  call MatCreateVecs(a_ts%FJacobian2,PETSC_NULL_OBJECT,XX,ierr);CHKERRQ(ierr)
#endif

  call scatter_from_xgc(a_ts,XX,n1,apar,phi,ierr);CHKERRQ(ierr)
  !call VecView(XX,PETSC_VIEWER_STDOUT_WORLD,ierr)

  ! solve
  its = 1
  ftime = a_dt
  call TSSetDuration(a_ts%ts,its,ftime,ierr);CHKERRQ(ierr)
  zero = 0.d0
  dt = a_dt
  call TSSetInitialTimeStep(a_ts%ts,zero,dt,ierr);CHKERRQ(ierr)
  call TSSetSolution(a_ts%ts,XX,ierr);CHKERRQ(ierr) ! why do this?
  call TSSetUp(a_ts%ts,ierr);CHKERRQ(ierr)
  call t_startf("PETSC_TS_SOLVE")
  !call PetscFPTrapPush(PETSC_FP_TRAP_ON,ierr);CHKERRQ(ierr)
  call TSSolve(a_ts%ts,XX,ierr);CHKERRQ(ierr)
  !call PetscFPTrapPop(ierr);CHKERRQ(ierr)
  call t_stopf("PETSC_TS_SOLVE")
  call TSGetSolveTime(a_ts%ts, ftime, ierr);CHKERRQ(ierr)
  ! Get the number of steps
  call TSGetTimeStepNumber(a_ts%ts, its, ierr);CHKERRQ(ierr)
  call TSGetConvergedReason(a_ts%ts, reason, ierr);CHKERRQ(ierr)
  if (sml_mype == 0) write(6,11)  sml_mype,'ts_solve: Number of pseudo timesteps = ',its,', final time ',ftime
11 format(I2,A,I2,A,D12.4)
  call scatter_to_xgc(a_ts,XX,n1,apar,phi,ierr);CHKERRQ(ierr)

  ! compute Je for output: del^2 A / mu; can skip if not needed
#if PETSC_VERSION_LT(3,6,0)
  call MatGetVecs(a_ts%muDel2,PETSC_NULL_OBJECT,vec3,ierr);CHKERRQ(ierr)
#else
  call MatCreateVecs(a_ts%muDel2,PETSC_NULL_OBJECT,vec3,ierr);CHKERRQ(ierr)
#endif
  call VecGetSubVector(XX,a_ts%iss(1),vec1,ierr);CHKERRQ(ierr) ! A
  call MatMult(a_ts%muDel2, vec1, vec3, ierr);CHKERRQ(ierr) ! J = del^2/mu A
  call VecRestoreSubVector(XX,a_ts%iss(1),vec1,ierr);CHKERRQ(ierr)
  ! put in Je
  its = 1
  call VecCreateSeqWithArray(PETSC_COMM_SELF,its,a_ts%nnode,Je,vec1,ierr);CHKERRQ(ierr)
  call VecScatterBegin(a_ts%from_petsc_single_value,vec3,vec1,INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(  a_ts%from_petsc_single_value,vec3,vec1,INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecDestroy(vec1,ierr)
  call VecDestroy(vec3,ierr)

  ! debug - compute residuals
!!$  if (.false.) then
!!$#if PETSC_VERSION_LT(3,6,0)
!!$     call MatGetVecs(a_ts%mass,vec2,vec3,ierr);CHKERRQ(ierr)
!!$#else
!!$     call MatCreateVecs(a_ts%mass,vec2,vec3,ierr);CHKERRQ(ierr)
!!$#endif
!!$     a = 1d0 ! already has minus sign
!!$     ! residuals, del^2(3,1) A = j
!!$     call MatGetSubMatrix(a_ts%FJacobian,a_ts%iss(3),a_ts%iss(1),MAT_INITIAL_MATRIX,mat,ierr);CHKERRQ(ierr)
!!$     call VecGetSubVector(XX,a_ts%iss(1),vec1,ierr);CHKERRQ(ierr) ! A
!!$     call MatMult(mat, vec1, vec2, ierr);CHKERRQ(ierr) ! -del^2 A
!!$     call VecRestoreSubVector(XX,a_ts%iss(1),vec1,ierr);CHKERRQ(ierr) 
!!$     call MatDestroy(mat,ierr);CHKERRQ(ierr)
!!$     call VecGetSubVector(XX,a_ts%iss(3),vec1,ierr);CHKERRQ(ierr) ! j
!!$     call MatGetSubMatrix(a_ts%FJacobian,a_ts%iss(3),a_ts%iss(3),MAT_INITIAL_MATRIX,mat,ierr);CHKERRQ(ierr)
!!$     call MatMult(mat, vec1, vec3, ierr);CHKERRQ(ierr) ! RHS j
!!$     call MatDestroy(mat,ierr);CHKERRQ(ierr)
!!$     call VecRestoreSubVector(XX,a_ts%iss(3),vec1,ierr);CHKERRQ(ierr) 
!!$     call VecAXPY(vec2,a,vec3,ierr);CHKERRQ(ierr) 
!!$     call VecNorm(vec2,NORM_2,norm1,ierr);CHKERRQ(ierr)
!!$     call VecNorm(vec3,NORM_2,norm2,ierr);CHKERRQ(ierr)
!!$     if (norm2/=0d0.and.sml_mype == 0) write(*,'(A,ES9.2E2,A,ES9.2E2)'),'TS.ts_solve RESIDUAL |j - del^2 A|/|j|=',norm1/norm2,', |j|=',norm2
!!$     ! residuals, del^2(2,2) phi = n
!!$     call MatGetSubMatrix(a_ts%FJacobian,a_ts%iss(2),a_ts%iss(2),MAT_INITIAL_MATRIX,mat,ierr);CHKERRQ(ierr)
!!$     call VecGetSubVector(XX,a_ts%iss(2),vec1,ierr);CHKERRQ(ierr) ! phi
!!$     call MatMult(mat, vec1, vec2, ierr);CHKERRQ(ierr) ! -del^2 A
!!$     call VecRestoreSubVector(XX,a_ts%iss(2),vec1,ierr);CHKERRQ(ierr) 
!!$     call MatDestroy(mat,ierr);CHKERRQ(ierr)
!!$     call VecGetSubVector(XX,a_ts%iss(0),vec1,ierr);CHKERRQ(ierr) ! n
!!$     call MatGetSubMatrix(a_ts%FJacobian,a_ts%iss(2),a_ts%iss(0),MAT_INITIAL_MATRIX,mat,ierr);CHKERRQ(ierr)
!!$     call MatMult(mat, vec1, vec3, ierr);CHKERRQ(ierr) ! RHS n
!!$     call MatDestroy(mat,ierr);CHKERRQ(ierr)
!!$     call VecRestoreSubVector(XX,a_ts%iss(0),vec1,ierr);CHKERRQ(ierr) 
!!$     call VecAXPY(vec2,a,vec3,ierr);CHKERRQ(ierr) 
!!$     call VecNorm(vec2,NORM_2,norm1,ierr);CHKERRQ(ierr)
!!$     call VecNorm(vec3,NORM_2,norm2,ierr);CHKERRQ(ierr)
!!$     if (norm2/=0d0.and.sml_mype == 0) write(*,'(A,ES9.2E2,A,ES9.2E2)')'TS.ts_solve RESIDUAL |n - nm/B^2 del^2 phi|/|n|=',norm1/norm2,', |n|=',norm2
!!$     !
!!$     call VecDestroy(vec2,ierr);CHKERRQ(ierr) 
!!$     call VecDestroy(vec3,ierr);CHKERRQ(ierr) 
!!$  end if

  call VecDestroy(XX,ierr);CHKERRQ(ierr)

  call PetscLogEventEnd(ADVSolveEvent,ierr)
  return
end subroutine ts_solve
! ---------------------------------------------------------------------
!
!  FormIFunction - Compute F(U,U_dot)
!
subroutine FormIFunction(ts_dummy,t_dummy,a_XX,a_Xdot,a_FF,a_ts,ierr)
  use petscts
  use xgc_ts_module
  use sml_module,only:sml_mype
  implicit none
  TS::ts_dummy
  PetscReal::t_dummy
  Vec,intent(in)::a_XX,a_Xdot
  Vec,intent(out)::a_FF
  type(xgc_ts)::a_ts
  PetscErrorCode,intent(out)::ierr
  !
  real (kind=8),dimension(a_ts%nnode,0:2)::lhs,rhs
  Vec::Fsub(0:2),Xdotsub(0:1),Rsub(0:2),RR
  PetscReal::norm

  !if (sml_mype == 0) write(6,*)  sml_mype,' **** TS.FormIFunction: time=',t_dummy
  call PetscLogEventBegin(ADVFormIFuncEvent,ierr)

  ! compute RHS in XGC space - 
  lhs(:,:)=0D0
  !rh We need this for matrix-free components on the LHS --->
  !rh call scatter_to_xgc(a_ts, a_XX, lhs(:,0), lhs(:,1), lhs(:,2), ierr);CHKERRQ(ierr) ! stay in native ordering
  !rh call FormRHSFunctionLocal(a_ts, lhs(:,0), lhs(:,1), lhs(:,2), rhs(:,0), rhs(:,1), rhs(:,2), a_ts%nnode)
  !rh Otherwise, we set rhs(:,:) to 0 --->
  rhs(:,:)=0D0
  rhs(:,0) = -rhs(:,0)*a_ts%scales(0) ! in xgc space so no permute, needs mass
  rhs(:,1) = -rhs(:,1)*a_ts%scales(1)
  rhs(:,2) = -rhs(:,2)*a_ts%scales(2)
  call VecDuplicate(a_FF,RR,ierr);CHKERRQ(ierr) ! permute to put in RHS (row)
  call scatter_from_xgc(a_ts,RR,rhs(:,0),rhs(:,1),rhs(:,2),ierr);
  CHKERRQ(ierr)

  !  mult by mass, put in F
  call VecGetSubVector(a_FF,a_ts%iss(0), Fsub(0), ierr);CHKERRQ(ierr)
  call VecGetSubVector(a_FF,a_ts%iss(1), Fsub(1), ierr);CHKERRQ(ierr)
  call VecGetSubVector(a_FF,a_ts%iss(2), Fsub(2), ierr);CHKERRQ(ierr)
  call VecGetSubVector(RR,  a_ts%iss(0), Rsub(0), ierr);CHKERRQ(ierr)
  call VecGetSubVector(RR,  a_ts%iss(1), Rsub(1), ierr);CHKERRQ(ierr)
  call VecGetSubVector(RR,  a_ts%iss(2), Rsub(2), ierr);CHKERRQ(ierr)
  call VecGetSubVector(a_Xdot,a_ts%iss(0),Xdotsub(0),ierr);CHKERRQ(ierr)
  call VecGetSubVector(a_Xdot,a_ts%iss(1),Xdotsub(1),ierr);CHKERRQ(ierr)

  ! add xdot to first two FD eqs
  ! y = alpha x + y.
  !PetscErrorCode  VecAXPY(Vec y,PetscScalar alpha,Vec x)
  call VecAXPY(Rsub(0),a_ts%scales(0),Xdotsub(0),ierr);CHKERRQ(ierr) 
  call VecAXPY(Rsub(1),a_ts%scales(1),Xdotsub(1),ierr);CHKERRQ(ierr)
  ! Copies a vector. y <- x
  ! PetscErrorCode  VecCopy(Vec x,Vec y)
  call VecCopy(           Rsub(0), Fsub(0), ierr);CHKERRQ(ierr)
  call VecCopy(           Rsub(1), Fsub(1), ierr);CHKERRQ(ierr)
  ! mult by mass for last two FE equations
  call MatMult(a_ts%mass, Rsub(2), Fsub(2), ierr);CHKERRQ(ierr)

  call VecRestoreSubVector(a_FF,a_ts%iss(0), Fsub(0), ierr);CHKERRQ(ierr)
  call VecRestoreSubVector(a_FF,a_ts%iss(1), Fsub(1), ierr);CHKERRQ(ierr)
  call VecRestoreSubVector(a_FF,a_ts%iss(2), Fsub(2), ierr);CHKERRQ(ierr)
  call VecRestoreSubVector(RR,  a_ts%iss(0), Rsub(0), ierr);CHKERRQ(ierr)
  call VecRestoreSubVector(RR,  a_ts%iss(1), Rsub(1), ierr);CHKERRQ(ierr)
  call VecRestoreSubVector(RR,  a_ts%iss(2), Rsub(2), ierr);CHKERRQ(ierr)
  call VecRestoreSubVector(a_Xdot,a_ts%iss(0),Xdotsub(0),ierr);CHKERRQ(ierr)
  call VecRestoreSubVector(a_Xdot,a_ts%iss(1),Xdotsub(1),ierr);CHKERRQ(ierr)

  call VecDestroy(RR,ierr);CHKERRQ(ierr)

  !Computes v3 = v2 + A * v1.
  !PetscErrorCode  MatMultAdd(Mat mat,Vec v1,Vec v2,Vec v3)
  ! add linear Jacobian matrix w/o time terms
  call MatMultAdd(a_ts%FJacobian,a_XX,a_FF,a_FF,ierr);CHKERRQ(ierr)

  !call VecNorm(a_FF,NORM_2,norm,ierr);CHKERRQ(ierr)
  !if (sml_mype==0.and.abs(norm).gt.1d30) write(6,*)  sml_mype,' ++++++++ FormIFunction |f|=',norm

  if (.false.) then ! debug
     call VecNorm(a_FF,NORM_2,norm,ierr);CHKERRQ(ierr)
     if (sml_mype==0.and.norm/=norm) write(6,*) sml_mype,' *** TS.FormIFunction: |F|=',norm
     if (sml_mype==0.and.abs(norm).gt.1d30) write(6,*) sml_mype,' *** TS.FormIFunction: |F|=',norm
  end if
  call PetscLogEventEnd(ADVFormIFuncEvent,ierr) 
  return
end subroutine FormIFunction
! ---------------------------------------------------------------------
!
!  FormIJacobian - Compute IJacobian = dF/dU + shift*dF/dUdot
!  f(TS ts,PetscReal t,Vec U,Vec U_t,PetscReal a,Mat Amat,Mat Pmat,void *ctx);
!
subroutine FormIJacobian(ts,t_dummy,a_XX,a_Xdot_dummy,shift,J,Jpre,a_ts,ierr)
  use petscts
  use xgc_ts_module
  use sml_module,only:sml_mype
  implicit none
  TS::ts
  PetscReal::t_dummy,shift
  Vec::a_XX,a_Xdot_dummy
  Mat::J,Jpre
  type(xgc_ts)::a_ts
  PetscErrorCode,intent(out)::ierr
  !
  PetscInt,parameter::ione=1,izero=0,bs=3
  Vec::vec
  PetscInt::ii,ix,jj,low,high
  PetscScalar::val,zero
  
  !if (sml_mype == 0) write(6,*)  sml_mype,' +++++++++ TS.FormIJacobian: shift=',shift
  call PetscLogEventBegin(ADVFormIJacEvent,ierr)
  ! add (lumped) mass * shift
#if PETSC_VERSION_LT(3,6,0)
  call MatGetVecs(Jpre,PETSC_NULL_OBJECT,vec,ierr);CHKERRQ(ierr)
#else
  call MatCreateVecs(Jpre,PETSC_NULL_OBJECT,vec,ierr);CHKERRQ(ierr)
#endif
  zero = 0.d0
  call MatGetOwnershipRange(a_ts%mass,low,high,ierr);CHKERRQ(ierr)
  do jj = 1,a_ts%nnode
     ii = a_ts%xgc_petsc(jj)
     if (ii.ge.low .and. ii.lt.high) then
        ii = ii*bs
        val = shift*a_ts%scales(0)
        call VecSetValue(vec,ii + 0, val,  INSERT_VALUES, ierr);CHKERRQ(ierr)
        val = shift*a_ts%scales(1)
        call VecSetValue(vec,ii + 1, val,  INSERT_VALUES, ierr);CHKERRQ(ierr)
        call VecSetValue(vec,ii + 2, zero, INSERT_VALUES, ierr);CHKERRQ(ierr)
     end if
  end do
  call VecAssemblyBegin(vec,ierr);CHKERRQ(ierr)
  call VecAssemblyEnd(vec,ierr);CHKERRQ(ierr)
  call MatCopy(a_ts%FJacobian,Jpre,SAME_NONZERO_PATTERN,ierr);CHKERRQ(ierr) 
  call MatDiagonalSet(Jpre,vec,ADD_VALUES,ierr);CHKERRQ(ierr)
  call VecDestroy(vec,ierr);CHKERRQ(ierr)
  
  !if (sml_mype==0) write(6,*)  sml_mype,'TS.FormIJacobian: shift = ', shift

  call MatAssemblyBegin(Jpre,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(Jpre,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

  if (J /= Jpre) then
     stop 'J /= Jpre ?'
  end if

  call PetscLogEventEnd(ADVFormIJacEvent,ierr) 

  return
end subroutine FormIJacobian

! ---------------------------------------------------------------------
!
!  FormLHSFunctionLocal - Compute F(t=0,U)
!
!   Input Parameters:
!      n1, apar, phi - U
!      nn - size XGC array (grid%nnode)
!
!  Output parameters:
!      lhsn1, lhsapar, lhsphi - value of function G(t=0,U)
!
!
subroutine FormLHSFunctionLocal(n1,apar,phi,lhsn1,lhsapar,lhsphi,nn)
  use sml_module,only:sml_mype
  implicit none
  integer::nn
  real (kind=8),dimension(nn),intent(in)::n1,apar,phi
  real (kind=8),dimension(nn),intent(out)::lhsn1,lhsapar,lhsphi
  
  !if (sml_mype==0) write(6,*)  sml_mype,'FormLHSFunctionLocal start N=',nn

  lhsn1 = 0d0 
  lhsapar = 0d0
  lhsphi = 0d0
  stop 'FormLHSFunctionLocal not used'

  return
end subroutine FormLHSFunctionLocal

! ---------------------------------------------------------------------
!
!  FormRHSFunctionLocal - Compute G(t=0,U)
!
!   Input Parameters:
!      n1, apar, phi - U
!      ui
!      nn - size XGC array (grid%nnode)
!
!  Output parameters:
!      rhsn1, rhsapar_dummy, rhsphi_dummy - value of function G(t=0,U)
!
subroutine FormRHSFunctionLocal(a_ts,n1,apar,phi,rhsn1,rhsapar,rhsphi,nn)
  use sml_module, only: sml_mype, sml_hyb_tearing_test, sml_e_charge
  use xgc_ts_module
  implicit none
  integer::nn
  type(xgc_ts)::a_ts
  real (kind=8),dimension(nn),intent(out)::rhsn1,rhsapar,rhsphi
  real (kind=8),dimension(nn),intent(in)::n1,apar,phi
  real (kind=8) :: deltape(nn), delta_B_perp(2,nn), dpar_deltape(nn), dB_grad_u0(nn)

  !if (sml_mype==0) write(6,*) sml_mype,'TS.FormRHSFunctionLocal start N=',nn

  rhsn1  = 0d0 ! -B_0 del_|| n_0 u_i / B_0
  rhsapar = 0d0
  rhsphi = 0d0 ! -dn_i

  !if (sml_hyb_tearing_test) then
  if (.false.) then
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!! Evaluate kink term dB_perp.grad(j_0/(e*B))
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !delta_B_perp(:,:)=0D0
    !call get_dbperp(a_ts%grid,apar,delta_B_perp)
    !call kink_drive(a_ts%grid,delta_B_perp,dB_grad_u0)
    call kink_drive2(grid,apar,dB_grad_u0)

    !rh "-" !!!
    rhsn1 = -dB_grad_u0
    !rh 2016/01/13: sign is questionable, try both versions
    ! rhsn1 = dB_grad_u0

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!! Evaluate b.grad(dP_e)/(e n_0)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !rh is needed for energy conservation of Shear-Alfven wave system
    !rh therefore, this term should be on the LHS

    !deltape(:)=0D0
    !deltape=sml_e_charge*n1*a_ts%grid%tempe
    !dpar_deltape(:)=0D0
    !call grad_par_delta_pe(a_ts%grid,a_ts%tr,a_ts%p,a_ts%dx,a_ts%bd,deltape,dpar_deltape)

    !!rh "-" !!!
    !rhsapar = -dpar_deltape
  endif


  return
end subroutine FormRHSFunctionLocal

! ------------------------------------------------------------------------
!
!   RHSFunction - User-provided routine that evaluates the RHS function
!   in the ODE.  This routine is set TSSetRHSFunction().
!     We compute: a_FF = G(t=0,U)
!
!   Input Parameters:
!      ts     - timestep context
!      time   - current time
!      a_XX   - input vector to function (n1,Apar,phi)
!      a_ts - user-provided context for function evaluation
!
!   Output Parameter:
!      a_FF - value of function
!      ierr - error code
!
subroutine RHSFunction(ts,time,a_XX,a_FF,a_ts,ierr)  ! NOT USED!!!
  use petscts
  use xgc_ts_module
  use sml_module,only:sml_mype
  implicit none
  !  Input/output parameters:
  TS::ts
  PetscReal::time
  Vec::a_XX,a_FF
  type(xgc_ts)::a_ts
  PetscErrorCode,intent(out)::ierr
  ! Local variables
  real (kind=8),dimension(a_ts%nnode,0:2)::lhs,rhs !n1,apar,phi,rhsn1,rhsapar,rhsphi
  Vec::Fsub(0:2),Rsub(0:2),RR
  integer::nn
  PetscReal::norm
  
  !stop 'RHSFunction not used'

  !if (sml_mype==0) write(6,*) sml_mype,' **** TS.RHSFunction: ',a_ts%nnode,' XGC vertices on plane'
  call PetscLogEventBegin(ADVFormRHSEvent,ierr) 

  !call VecNorm(a_XX,NORM_2,norm,ierr);CHKERRQ(ierr)
  !if (sml_mype==0.and.abs(norm).gt.1d30) write(6,*)  sml_mype,'RHSFunction |x|=',norm

  ! compute RHS in XGC space - 
  lhs(:,:)=0D0
  call scatter_to_xgc(a_ts, a_XX, lhs(:,0), lhs(:,1), lhs(:,2), ierr);CHKERRQ(ierr) ! stay in native ordering
  nn = a_ts%nnode
  call FormRHSFunctionLocal(a_ts, lhs(:,0), lhs(:,1), lhs(:,2), rhs(:,0), rhs(:,1), rhs(:,2), nn)
  rhs(:,0) = rhs(:,0)*a_ts%scales(0) ! in xgc space so no permute
  rhs(:,1) = rhs(:,1)*a_ts%scales(1)
  rhs(:,2) = rhs(:,2)*a_ts%scales(2)

  call VecDuplicate(a_FF,RR,ierr);CHKERRQ(ierr) ! permute to put in RHS (row)
  call scatter_from_xgc(a_ts,RR,rhs(:,0),rhs(:,1),rhs(:,2),ierr);
  CHKERRQ(ierr)

  !call VecNorm(RR,NORM_2,norm,ierr);CHKERRQ(ierr)
  !if (sml_mype==0.and.abs(norm).gt.1d30) write(6,*)  sml_mype,'RHSFunction 1) |f|=',norm

  !  mult by mass, put in F
  call VecGetSubVector(a_FF,a_ts%iss(0), Fsub(0), ierr);CHKERRQ(ierr)
  call VecGetSubVector(a_FF,a_ts%iss(1), Fsub(1), ierr);CHKERRQ(ierr)
  call VecGetSubVector(a_FF,a_ts%iss(2), Fsub(2), ierr);CHKERRQ(ierr)
  call VecGetSubVector(RR,  a_ts%iss(0), Rsub(0), ierr);CHKERRQ(ierr)
  call VecGetSubVector(RR,  a_ts%iss(1), Rsub(1), ierr);CHKERRQ(ierr)
  call VecGetSubVector(RR,  a_ts%iss(2), Rsub(2), ierr);CHKERRQ(ierr)

  call VecCopy(           Rsub(0), Fsub(0), ierr);CHKERRQ(ierr)
  call VecCopy(           Rsub(1), Fsub(1), ierr);CHKERRQ(ierr)
  call MatMult(a_ts%mass, Rsub(2), Fsub(2), ierr);CHKERRQ(ierr)

  call VecRestoreSubVector(a_FF,a_ts%iss(0), Fsub(0), ierr);CHKERRQ(ierr)
  call VecRestoreSubVector(a_FF,a_ts%iss(1), Fsub(1), ierr);CHKERRQ(ierr)
  call VecRestoreSubVector(a_FF,a_ts%iss(2), Fsub(2), ierr);CHKERRQ(ierr)
  call VecRestoreSubVector(RR,  a_ts%iss(0), Rsub(0), ierr);CHKERRQ(ierr)
  call VecRestoreSubVector(RR,  a_ts%iss(1), Rsub(1), ierr);CHKERRQ(ierr)
  call VecRestoreSubVector(RR,  a_ts%iss(2), Rsub(2), ierr);CHKERRQ(ierr)

  call VecDestroy(RR,ierr);CHKERRQ(ierr)

  !  debugging information if desired
  if (.false.) then ! debug
     call VecNorm(a_FF,NORM_2,norm,ierr);CHKERRQ(ierr)
     if (sml_mype==0.and.norm/=norm) write(6,*)  sml_mype,'RHSFunction done, RHS function vector: |f|=',norm
     if (sml_mype==0.and.abs(norm).gt.1d30) write(6,*)  sml_mype,'RHSFunction done, RHS function vector: |f|=',norm
  endif

  call PetscLogEventEnd(ADVFormRHSEvent,ierr) 

  return
end subroutine RHSFunction

! ------------------------------------------------------------------------
!  subroutine scatter_to_xgc
!
!  Input Parameters:
!     a_ts - the user context
!     XX - PETSc fields
!     
!  Output Parameter:
!     n1,apar,phi - XGC fields
!     ierr - error code
!
#undef __FUNCT__
#define __FUNCT__ "scatter_to_xgc"
subroutine scatter_to_xgc(a_ts,a_XX,a_n1,a_apar,a_phi,ierr)
  use petscts
  use sml_module,only:sml_mype
  use xgc_ts_module
  implicit none
  type(xgc_ts),intent(in)::a_ts
  Vec,intent(in)::a_XX
  real (kind=8),dimension(a_ts%nnode)::a_n1,a_apar,a_phi
  PetscErrorCode,intent(out)::ierr
  ! locals
  PetscInt,parameter::ione=1
  PetscScalar,dimension(a_ts%nnode)::n1,apar,phi
  Vec::xxvec(0:2)

  ! scatter solution back - n1
  n1 = a_n1
  call VecCreateSeqWithArray(PETSC_COMM_SELF,ione,a_ts%nnode,n1,xxvec(0),ierr);CHKERRQ(ierr)
  call VecScatterBegin(a_ts%from_petsc(0),a_XX,xxvec(0),INSERT_VALUES,SCATTER_FORWARD,ierr)

  ! scatter solution back - apar
  apar = a_apar
  call VecCreateSeqWithArray(PETSC_COMM_SELF,ione,a_ts%nnode,apar,xxvec(1),ierr);CHKERRQ(ierr)
  call VecScatterBegin(a_ts%from_petsc(1),a_XX,xxvec(1),INSERT_VALUES,SCATTER_FORWARD,ierr)

  ! scatter solution back - phi
  phi = a_phi
  call VecCreateSeqWithArray(PETSC_COMM_SELF,ione,a_ts%nnode,phi,xxvec(2),ierr);CHKERRQ(ierr)
  call VecScatterBegin(a_ts%from_petsc(2),a_XX,xxvec(2),INSERT_VALUES,SCATTER_FORWARD,ierr)

  ! end
  call VecScatterEnd(  a_ts%from_petsc(0),a_XX,xxvec(0),INSERT_VALUES,SCATTER_FORWARD,ierr)
  a_n1 = n1
  call VecDestroy(xxvec(0),ierr)

  call VecScatterEnd(  a_ts%from_petsc(1),a_XX,xxvec(1),INSERT_VALUES,SCATTER_FORWARD,ierr)
  a_apar = apar
  call VecDestroy(xxvec(1),ierr)

  call VecScatterEnd(  a_ts%from_petsc(2),a_XX,xxvec(2),INSERT_VALUES,SCATTER_FORWARD,ierr)
  a_phi = phi
  call VecDestroy(xxvec(2),ierr)

  return
end subroutine scatter_to_xgc

! ------------------------------------------------------------------------
!  subroutine scatter_from_xgc
!
!  Input Parameters:
!     a_ts - the user context
!     n1,Aparpot,epot - XGC fields
!     
!  Output Parameter:
!     ierr - error code
!     XX - PETSc fields
!
#undef __FUNCT__
#define __FUNCT__ "scatter_from_xgc"
subroutine scatter_from_xgc(a_ts,a_XX,a_n1,a_apar,a_phi,ierr)
  use petscts
  use sml_module,only:sml_mype
  use xgc_ts_module
  implicit none
  type(xgc_ts)::a_ts
  Vec::a_XX
  real (kind=8),dimension(a_ts%nnode)::a_n1,a_apar,a_phi
  PetscErrorCode,intent(out)::ierr
  ! locals
  PetscInt,parameter::ione=1
  PetscScalar,dimension(a_ts%nnode)::n1,apar,phi
  Vec::xxvec(0:2)

  ! scatter current solution - n1
  n1 = a_n1
  call VecCreateSeqWithArray(PETSC_COMM_SELF,ione,a_ts%nnode,n1,xxvec(0),ierr);CHKERRQ(ierr)
  call VecScatterBegin(a_ts%to_petsc(0),xxvec(0),a_XX,INSERT_VALUES,SCATTER_FORWARD,ierr)

  ! scatter current solution - a_par
  apar = a_apar
  call VecCreateSeqWithArray(PETSC_COMM_SELF,ione,a_ts%nnode,apar,xxvec(1),ierr);CHKERRQ(ierr)
  call VecScatterBegin(a_ts%to_petsc(1),xxvec(1),a_XX,INSERT_VALUES,SCATTER_FORWARD,ierr)

  ! scatter current solution - phi
  phi = a_phi
  call VecCreateSeqWithArray(PETSC_COMM_SELF,ione,a_ts%nnode,phi,xxvec(2),ierr);CHKERRQ(ierr)
  call VecScatterBegin(a_ts%to_petsc(2),xxvec(2),a_XX,INSERT_VALUES,SCATTER_FORWARD,ierr)

  ! end
  call VecScatterEnd(  a_ts%to_petsc(0),xxvec(0),a_XX,INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecDestroy(xxvec(0),ierr)

  call VecScatterEnd(  a_ts%to_petsc(1),xxvec(1),a_XX,INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecDestroy(xxvec(1),ierr)

  call VecScatterEnd(  a_ts%to_petsc(2),xxvec(2),a_XX,INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecDestroy(xxvec(2),ierr)

  return
end subroutine scatter_from_xgc
