! Interface methods for FEM operators
! Mark Adams October 2012
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     helm_matrix
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
subroutine helm_matrix( &
     Amat, &                   ! matrix (out)
     alpha, &                  ! operator: 
     beta, &                   !   alpha del^2 u + beta u (in)
     grid, &                   ! grid (in)
     bd, &                     ! BCs (in)
     bcflag, &                 ! BC flag (in)
     xgc_petsc, &              ! map to PETSc equation number
     ierr  &                   ! error code (out)
     )
  use grid_class
  use boundary_class
#include <petscversion.h>
#if PETSC_VERSION_LT(3,6,0)
#include <finclude/petscdef.h>
#else
#include <petsc/finclude/petscdef.h>
#endif
  use petscmat
  implicit  none
  Mat::Amat
  type(grid_type),intent(in)::grid
  real (kind=8),intent(in)::alpha(grid%ntriangle), beta(grid%ntriangle)
  logical,intent(in)::bcflag
  type(boundary2_type),intent(in)::bd
  PetscInt,intent(in)::xgc_petsc(grid%nnode)
  PetscErrorCode,intent(out)::ierr
  !
  PetscInt::bs,idx,jdx
  
  bs = 1; idx = 0; jdx = 0
  call helm_matrix_private(Amat,alpha,beta,grid%ntriangle,grid%nd,grid%nnode,grid%x,bd,bcflag,xgc_petsc,bs,idx,jdx,PETSC_NULL_OBJECT,ierr)

  call MatAssemblyBegin(Amat,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(Amat,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

end subroutine helm_matrix
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     helm_matrix
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
subroutine helm_matrix_mdof( &
     Amat, &                   ! matrix (out)
     alpha, &                  ! operator: 
     beta, &                   !   alpha del^2 u + beta u (in)
     grid, &                   ! grid (in)
     bd, &                     ! BCs (in)
     bcflag, &                 ! BC flag (in)
     xgc_petsc, &              ! map to PETSc equation number
     bs, idx, jdx, &           ! blocks size (>1) and row/col equation index (zero bases, <bs)
     scale_vec, &              ! vector to use to scale
     ierr  &                   ! error code (out)
     )
  use grid_class
  use boundary_class
#include <petscversion.h>
#if PETSC_VERSION_LT(3,6,0)
#include <finclude/petscdef.h>
#else
#include <petsc/finclude/petscdef.h>
#endif
  use petscmat
  implicit  none  
  Mat::Amat
  Vec::scale_vec
  type(grid_type),intent(in)::grid
  real (kind=8),intent(in)::alpha(grid%ntriangle), beta(grid%ntriangle)
  logical,intent(in)::bcflag
  type(boundary2_type),intent(in)::bd
  PetscInt,intent(in)::xgc_petsc(grid%nnode), bs, idx, jdx
  PetscErrorCode,intent(out)::ierr
  !
  if (idx.ge.bs.or.jdx.ge.bs) stop 'helm_matrix_mdof bad eq index'
  !
  call helm_matrix_private(Amat,alpha,beta,grid%ntriangle,grid%nd,grid%nnode,grid%x,bd,bcflag,xgc_petsc,bs,idx,jdx,scale_vec,ierr)
  
end subroutine helm_matrix_mdof

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     fem_del2: -del^2 u [need to change to fem_neg_del2()]
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
subroutine fem_del2( &
     fu, &                     ! f(u) (out)
     uu, &                     ! u (in)
     grid, &                   ! grid (in)
     bc &
     )
  use grid_class
  implicit  none
  
  type(grid_type), intent(in):: grid
  real (kind=8), intent(out), dimension(grid%nnode) :: fu
  real (kind=8), intent(in) :: uu(grid%nnode)
  type(boundary2_type),intent(in)::bc
  
  call fem_op_private(fu,uu,grid%ntriangle,grid%nd,grid%nnode,grid%x,bc,2)
  
end subroutine fem_del2

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     fem_mass
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
subroutine fem_mass( &
     fu, &                     ! f(u) (out)
     uu, &                     ! u (in)
     grid, &                   ! grid (in)
     bc &
     )
  use grid_class
  implicit  none
  type(grid_type), intent(in):: grid
  real (kind=8), intent(out), dimension(grid%nnode) :: fu
  real (kind=8), intent(in) :: uu(grid%nnode)
  type(boundary2_type),intent(in)::bc
  
  call fem_op_private(fu,uu,grid%ntriangle,grid%nd,grid%nnode,grid%x,bc,1)
  !stop 'fem_mass not used'
end subroutine fem_mass

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     b_dot_grad_u
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
subroutine b_dot_grad_u( &
     bdotgradu, &              ! f(u) (out)
     b1, &                     ! b_x (in)
     b2, &                     ! b_y (in)
     uu, &                     ! u (in)
     grid &                    ! grid (in)
     )
  use grid_class
  implicit  none
  !
  type(grid_type),intent(in)::grid
  real (kind=8),intent(in),dimension(grid%nnode)::b1,b2,uu
  real (kind=8),intent(out),dimension(grid%nnode)::bdotgradu
  !
  call b_dot_grad_u_private(bdotgradu,b1,b2,uu,grid%ntriangle, grid%nd, grid%nnode, grid%x)
  
end subroutine b_dot_grad_u

! private methods
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     helm_matrix_private - make Helmoltz (alpha del^2 u + beta u) matrix
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
subroutine helm_matrix_private( &
     Amat, &                   ! matrix (out)
     alpha, &                  ! operator: 
     beta, &                   !   alpha del^2 u + beta u (in)
     iemax, &                  ! number of triangles (in)
     trind, &                  ! indedices for triangles (in)
     npts,  &                  ! number of vertices (in)
     points, &                 ! coordinates of vertices (in)
     bd, &                     ! BCs (in)
     bcflag, &                 ! BC flag (in)
     xgc_petsc, &              ! map to PETSc equation number
     bs, idx, jdx, &           ! blocks size (>1) and row/col equation index (<bs)
     scale_vec, &              ! vector to use to scale
     ierr  &                   ! error code (out)
     )      
  use boundary_class
!#include <finclude/petscdef.h>
  use petscmat
  use sml_module,only:sml_mype,sml_plane_index
  implicit  none
  !
  Mat::Amat
  Vec::scale_vec
  integer,intent(in)::iemax,npts,trind(3,iemax)
  real (kind=8),intent(in)::points(2,npts),alpha(iemax),beta(iemax)
  logical, intent(in)::bcflag
  type(boundary2_type),intent(in)::bd
  PetscInt,intent(in)::xgc_petsc(npts), bs, idx, jdx
  PetscErrorCode,intent(out)::ierr
  !
  real (kind=8)::c1,c2,xlt(2,3),ul(3),p(3)
  integer::ii,jj,kk
  PetscInt::ie,ind(3),itmp,jtmp,idi(1),idj(1)
  PetscInt,parameter::ione=1
  PetscScalar::arr(3,3),val
  PetscScalar,pointer::scale_v(:)
  PetscInt::low1,lowbs,highbs

  call MatGetOwnershipRange(Amat,lowbs,highbs,ierr)
  if (scale_vec/=0) then
     call VecGetOwnershipRange(scale_vec,low1,PETSC_NULL_INTEGER,ierr)
     call VecGetArrayReadF90(scale_vec,scale_v,ierr)
  end if

  ul  = 0D0 ! displacement not used for matrix
  do ie=1,iemax
     ind(:) = trind(:,ie)
     if(ind(1) .le. 0) stop 'ind(1) == 0'
     
     xlt(1,1) = points(1,ind(1))
     xlt(1,2) = points(1,ind(2))
     xlt(1,3) = points(1,ind(3))
     xlt(2,1) = points(2,ind(1))
     xlt(2,2) = points(2,ind(2))
     xlt(2,3) = points(2,ind(3))

     c1 = alpha(ie) ! del^2
     c2 = beta(ie)  ! mass

     call helm2dElem(c1,c2,ul,xlt,arr,p,1)

     do itmp=1,3
        ii = ind(itmp)
        if (is_inside(ii,bd)) then
           do jtmp=1,3
              jj = ind(jtmp)
              if (is_inside(jj,bd) .or. bcflag) then
                 kk = xgc_petsc(ii)
                 if (kk*bs.ge.lowbs .and. kk*bs.lt.highbs) then
                    idi(1) =            kk*bs + idx
                    idj(1) = xgc_petsc(jj)*bs + jdx
                    if (scale_vec/=0) then
                       val = arr(itmp,jtmp)*scale_v(kk-low1+1) ! scale, one based
                    else 
                       val = arr(itmp,jtmp)
                    end if
                    !if (sml_plane_index==0 .and. idi(1).eq.8)write(6,'(A,I3,A,I8,I8,ES12.4E2)'),'     [',sml_mype,'] helm_matrix_private: add to ',idi(1),idj(1),val
                    call MatSetValues(Amat,ione,idi,ione,idj,val,ADD_VALUES,ierr)
                    CHKERRQ(ierr)
                    if (abs(arr(itmp,jtmp)-arr(jtmp,itmp))/arr(1,1).gt.1.e-8 .or. ierr/=0) then
                       write(*,'(I6,A,I12,I12,ES12.4E2)'),sml_mype,') I,J,diff=',idi(1),idj(1),(arr(itmp,jtmp)-arr(jtmp,itmp))/arr(1,1)
                       write(*,'(ES12.4E2,ES12.4E2,ES12.4E2)'),arr(1,1),arr(1,2),arr(1,3)
                       write(*,'(ES12.4E2,ES12.4E2,ES12.4E2)'),arr(2,1),arr(2,2),arr(2,3)
                       stop 'helm_matrix_private: not symmetric'
                    end if
                 end if
              endif
           enddo
        end if
     enddo
  enddo
  if (scale_vec/=0) then
     call VecRestoreArrayReadF90(scale_vec,scale_v,ierr) 
  end if
  !    xs BCs  
!!$  do ii=1,npts
!!$     if( .not.is_inside(ii,bd) ) then
!!$        idi(1) = xgc_petsc(ii)*bs + idx
!!$        call MatSetValues(Amat,ione,idi,ione,idi,diagonal,ADD_VALUES,ierr)
!!$        CHKERRQ(ierr)
!!$        if (sml_mype==0.and.idi(1).ge.0) write(*,'(A,I6,A,I6)'),'helm_matrix_private: zeroing out xgc row=',ii,', PETSc row=',idi(1),'.  should not happen'
!!$     endif
!!$  enddo
  
end subroutine helm_matrix_private


! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     fem_op_private - apply Helmoltz (alpha del^2 u + beta u) operator
!       just apply mass or del^2 term w/o coeficient (isw switch)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
subroutine fem_op_private( &
     fu, &                     ! f(u) (out)
     uu, &                     ! u (in)
     iemax, &                  ! number of triangles (in)
     trind, &                  ! indedices for triangles (in)
     npts,  &                  ! number of vertices (in)
     points, &                 ! coordinates of vertices (in)
     bc, & 
     isw &                     ! 1:u, 2:del^2(u)
     )
  use boundary_class
  implicit  none
  !     
  integer,intent(in)::iemax,npts,trind(3,iemax),isw
  real (kind=8),intent(out)::fu(npts)
  real (kind=8),intent(in)::uu(npts),points(2,npts)
  type(boundary2_type),intent(in)::bc
  !     
  real (kind=8) ::  c1,c2,xlt(2,3),arr(3,3),ul(3),pp(3)
  integer :: ie,ii,jj,ind(3),itmp,jtmp

  arr = 1d300 ! not used

  fu = 0.d0
  do ie=1,iemax
     ind(:) = trind(:,ie)
     if(ind(1) .le. 0) stop 'ind(1) == 0'
     
     xlt(1,1) = points(1,ind(1))
     xlt(1,2) = points(1,ind(2))
     xlt(1,3) = points(1,ind(3))
     xlt(2,1) = points(2,ind(1))
     xlt(2,2) = points(2,ind(2))
     xlt(2,3) = points(2,ind(3))
     
     do itmp=1,3
        ii = ind(itmp)
        ul(itmp) = uu(ii) ! grab dispacements
     enddo
     
     if( isw == 1 ) then
        c1 = 0.d0
        c2 = 1.d0
     else ! if( isw .eq. 2 ) then
        c1 = 1.d0 
        c2 = 0.d0
     endif

     call helm2dElem(c1,c2,ul,xlt,arr,pp,2)
     
     do itmp=1,3
        ii = ind(itmp)
        if (is_inside(ii,bc)) then
           fu(ii) = fu(ii) + pp(itmp)
        end if
     enddo
  enddo
  
end subroutine fem_op_private

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     b_dot_grad_u_private - apply Helmoltz (alpha del^2 u + beta u) operator
!       just apply mass or del^2 term w/o coeficient (isw switch)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
subroutine b_dot_grad_u_private( &
     bdotgradu, &              ! f(u) (out)
     b1, &                     ! b_x (in)
     b2, &                     ! b_y (in)
     uu, &                     ! u (in)
     iemax, &                  ! number of triangles (in)
     trind, &                  ! indedices for triangles (in)
     npts,  &                  ! number of vertices (in)
     points &                  ! coordinates of vertices (in)
     )
  implicit  none
  !
  integer,intent(in)::iemax,npts,trind(3,iemax)
  real (kind=8),intent(in)::b1(npts),b2(npts),uu(npts),points(2,npts) 
  real (kind=8),intent(out)::bdotgradu(npts)
  !
  real (kind=8)::c1,c2,xl(2,3),bb(3,3),ul(3),pp(3),tarea(3),area(npts)
  !real (kind=8),allocatable::area(:)
  integer::ie,ii,jj,ind(3),itmp,jtmp
  !
  c1 = 1.d300                 ! dummy args 
  c2 = 1.d300
  !allocate(area(npts))
  bdotgradu = 0.d0
  area = 0.d0
  do ie=1,iemax
     ind(:) = trind(:,ie)
     if(ind(1) .le. 0) stop 'ind(1) == 0'
     xl(1,1) = points(1,ind(1))
     xl(1,2) = points(1,ind(2))
     xl(1,3) = points(1,ind(3))
     xl(2,1) = points(2,ind(1))
     xl(2,2) = points(2,ind(2))
     xl(2,3) = points(2,ind(3))
     do itmp=1,3
        ii = ind(itmp)
        ul(itmp) = uu(ii)     ! grab dispacements
        bb(itmp,1) = b1(ii)   ! grab B
        bb(itmp,2) = b2(ii)
     enddo
     !call helm2dElem(c1,c2,ul,xl,bb,pp,3)
     call bdotgradelem(ul,xl,bb,pp,tarea)
     do itmp=1,3
        ii = ind(itmp)
        bdotgradu(ii) = bdotgradu(ii) + pp(itmp)
        !area(ii) = area(ii) + bb(itmp,1)
        area(ii) = area(ii) + tarea(itmp)
     enddo

  enddo
  bdotgradu = bdotgradu / area
  !deallocate(area)
  return
end subroutine b_dot_grad_u_private

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     main
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
subroutine test_fem(comm)
  ! program main
  use boundary_class
  use petscksp
  use grid_class
  implicit none
  MPI_Comm::comm
  !
  PetscErrorCode::ierr
  PetscMPIInt::npe,mype
  Mat::Amat,MassMat
  KSP::ksp
  Vec::xx,ff
  PetscInt,parameter::ntri=8,nnodes=9
  integer :: trind(3,ntri),ii,jj,ix,jx,kx
  real (kind=8)::crds(2,nnodes), uu(nnodes), bb1(nnodes), bb2(nnodes),pp(nnodes)
  real (kind=8)::alpha(ntri), beta(ntri), dom_len
  type(boundary2_type)::bd
  PetscInt::idx(1),ione
  PetscScalar::val(1)
  PetscInt::xgc_petsc(nnodes),idj
  PetscViewer::viewer
  MatNullSpace::nullsp
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  call MPI_Comm_size(comm,npe,ierr)
  call MPI_Comm_rank(comm,mype,ierr)
  if (mype==0) print *,'test_fem: npe=',npe
 
  ! setup grid -- coords
  dom_len = 2.d0
  crds(1,1) = 0.d0*dom_len;       crds(2,1) = 0.d0*dom_len
  crds(1,2) = .5d0*dom_len;       crds(2,2) = 0.d0*dom_len
  crds(1,3) = 1.d0*dom_len;       crds(2,3) = 0.d0*dom_len
  crds(1,4) = 0.d0*dom_len;       crds(2,4) = .5d0*dom_len
  !crds(1,5) = .6d0*dom_len;       crds(2,5) = .7d0*dom_len ! patch test
  crds(1,5) = .5d0*dom_len;       crds(2,5) = .5d0*dom_len  ! regular grid
  crds(1,6) = 1.d0*dom_len;       crds(2,6) = .5d0*dom_len
  crds(1,7) = 0.d0*dom_len;       crds(2,7) = 1.d0*dom_len
  crds(1,8) = .5d0*dom_len;       crds(2,8) = 1.d0*dom_len
  crds(1,9) = 1.d0*dom_len;       crds(2,9) = 1.d0*dom_len
  !     mesh
  trind(1,1) = 1;  trind(2,1) = 2; trind(3,1) = 4
  trind(1,2) = 2;  trind(2,2) = 5; trind(3,2) = 4
  trind(1,3) = 2;  trind(2,3) = 3; trind(3,3) = 5
  trind(1,4) = 3;  trind(2,4) = 6; trind(3,4) = 5
  trind(1,5) = 4;  trind(2,5) = 5; trind(3,5) = 7
  trind(1,6) = 5;  trind(2,6) = 8; trind(3,6) = 7
  trind(1,7) = 5;  trind(2,7) = 6; trind(3,7) = 8
  trind(1,8) = 6;  trind(2,8) = 9; trind(3,8) = 8
  !     BCs
  !      bd%in%start = 1; bd%in%end = 1;
  !      bd%out1%start = 4; bd%out1%end = 4
  !      bd%out2%start = 7; bd%out2%end = 7
  
  ! matrix ops
  if (.false.) then
     ! create matrix
     call MatCreate(comm,Amat,ierr)
     call MatSetSizes(Amat,PETSC_DECIDE,PETSC_DECIDE,nnodes,nnodes,ierr)
     call MatSetType( Amat, MATMPIAIJ, ierr )
     call MatCreate(comm,MassMat,ierr)
     call MatSetSizes(MassMat,PETSC_DECIDE,PETSC_DECIDE,nnodes,nnodes,ierr)
     call MatSetType( MassMat, MATMPIAIJ, ierr )
     
     call MatSeqAIJSetPreallocation(Amat,9,PETSC_NULL_INTEGER,ierr)
     call MatMPIAIJSetPreallocation(Amat,9,PETSC_NULL_INTEGER,6,PETSC_NULL_INTEGER, ierr)
     call MatSeqAIJSetPreallocation(MassMat,9,PETSC_NULL_INTEGER,ierr)
     call MatMPIAIJSetPreallocation(MassMat,9,PETSC_NULL_INTEGER,6,PETSC_NULL_INTEGER, ierr)
     
     call MatSetOption(Amat,MAT_IGNORE_OFF_PROC_ENTRIES,PETSC_TRUE,ierr)
     call MatSetOption(MassMat,MAT_IGNORE_OFF_PROC_ENTRIES,PETSC_TRUE,ierr)
     !     get Amat
     idx(1) = 1
     idj = 0
     alpha = 0.0d0 ! lapacian
     beta = 1.d0   ! mass
     call helm_matrix_private(MassMat,alpha,beta,ntri,trind,nnodes,crds,bd,.false.,xgc_petsc,idx(1),idj,idj,PETSC_NULL_OBJECT,ierr)
     call MatAssemblyBegin(MassMat,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
     call MatAssemblyEnd(MassMat,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
     alpha = 1.0d0 ! lapacian
     beta = 0.d0   ! mass
     call helm_matrix_private(Amat,alpha,beta,ntri,trind,nnodes,crds,bd,.false.,xgc_petsc,idx(1),idj,idj,PETSC_NULL_OBJECT,ierr)
     call MatAssemblyBegin(Amat,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
     call MatAssemblyEnd(Amat,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
     
     call PetscViewerASCIIOpen(comm, 'testA.m', viewer,ierr);CHKERRQ(ierr)
     call PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB,ierr);CHKERRQ(ierr)
     call MatView(Amat,viewer,ierr);CHKERRQ(ierr)
     call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
     
     call PetscViewerASCIIOpen(comm, 'testM.m', viewer,ierr);CHKERRQ(ierr)
     call PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB,ierr);CHKERRQ(ierr)
     call MatView(MassMat,viewer,ierr);CHKERRQ(ierr)
     call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
     
     call KSPCreate(comm,ksp,ierr);CHKERRQ(ierr)
     call KSPSetOperators(ksp,Amat,Amat,ierr);CHKERRQ(ierr)
     call MatNullSpaceCreate(comm, PETSC_TRUE, 0, PETSC_NULL_OBJECT, nullsp,ierr);CHKERRQ(ierr)
     call MatSetNullSpace(Amat,nullsp,ierr);CHKERRQ(ierr)
     call MatNullSpaceDestroy(nullsp,ierr);CHKERRQ(ierr)
     call KSPSetFromOptions(ksp,ierr);CHKERRQ(ierr)
     
     !     solve
#if PETSC_VERSION_LT(3,6,0)
     call MatGetVecs(Amat,PETSC_NULL_OBJECT,xx,ierr);CHKERRQ(ierr)
#else
     call MatCreateVecs(Amat,PETSC_NULL_OBJECT,xx,ierr);CHKERRQ(ierr)
#endif
     call VecDuplicate(xx,ff,ierr);CHKERRQ(ierr)
     !     f
     ione = 1;
     if (.true.) then
        idx(1) = 3-1; val = 1d0
        call VecSetValues(ff,ione,idx(1),val,ADD_VALUES,ierr)
        idx(1) = 9-1; val = 1d0
        call VecSetValues(ff,ione,idx(1),val,ADD_VALUES,ierr)
        idx(1) = 6-1; val = 4d0
        call VecSetValues(ff,ione,idx(1),val,ADD_VALUES,ierr)
        idx(1) = 1-1; val = -1d0
        call VecSetValues(ff,ione,idx(1),val,ADD_VALUES,ierr)
        idx(1) = 7-1; val = -1d0
        call VecSetValues(ff,ione,idx(1),val,ADD_VALUES,ierr)
        idx(1) = 4-1; val = -4d0
        call VecSetValues(ff,ione,idx(1),val,ADD_VALUES,ierr)
        call VecAssemblyBegin(ff,ierr);
        call VecAssemblyEnd(ff,ierr);
     else ! use a mass matrix - very very bad
        idx(1) = 3-1; val = 1d0
        call VecSetValues(xx,ione,idx(1),val,ADD_VALUES,ierr)
        idx(1) = 9-1; 
        call VecSetValues(xx,ione,idx(1),val,ADD_VALUES,ierr)
        idx(1) = 6-1; 
        call VecSetValues(xx,ione,idx(1),val,ADD_VALUES,ierr)
        idx(1) = 1-1; val = -1d0
        call VecSetValues(xx,ione,idx(1),val,ADD_VALUES,ierr)
        idx(1) = 7-1; 
        call VecSetValues(xx,ione,idx(1),val,ADD_VALUES,ierr)
        idx(1) = 4-1; 
        call VecSetValues(xx,ione,idx(1),val,ADD_VALUES,ierr)
        call VecAssemblyBegin(xx,ierr);
        call VecAssemblyEnd(xx,ierr);
        call MatMult(MassMat,xx,ff,ierr);CHKERRQ(ierr)
        call VecZeroEntries(xx,ierr)
     end if
     call KSPSolve(ksp,ff,xx,ierr);CHKERRQ(ierr)
     call PetscViewerASCIIOpen(comm, 'testX.m', viewer,ierr);CHKERRQ(ierr)
     call PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB,ierr);CHKERRQ(ierr)
     call VecView(xx,viewer,ierr);CHKERRQ(ierr)
     call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
     call PetscViewerASCIIOpen(comm, 'testF.m', viewer,ierr);CHKERRQ(ierr)
     call PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB,ierr);CHKERRQ(ierr)
     call VecView(ff,viewer,ierr);CHKERRQ(ierr)
     call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
     call MatDestroy(Amat,ierr)
     call MatDestroy(MassMat,ierr)
     call KSPDestroy(ksp,ierr)
     call VecDestroy(xx,ierr)
     call VecDestroy(ff,ierr)
  end if
  ! matrix free ops
  if (.true. .and. mype==0) then
     print *, 'u='
     jj = 1
     do ii=1,3
        print *, uu(jj:(jj+2))
        jj = jj + 3
     end do
     call fem_op_private( bb1, uu, ntri, trind, nnodes, crds, bd, 1 )
     print *, 'M*u='
     jj = 1
     do ii=1,3
        print *, bb1(jj:(jj+2))
        jj = jj + 3
     end do
 
     call fem_op_private( bb1, uu, ntri, trind, nnodes, crds, bd, 2 )
     print *, 'del^2(u) = '
     jj = 1
     do ii=1,3
        print *, bb1(jj:(jj+2))
        jj = jj + 3
     end do
     bb1 = 5.d0
     bb2 = 1.d30
     call b_dot_grad_u_private( pp, bb1, bb2, uu, ntri, trind, nnodes, crds )
     print *, 'b_dot_grad_u='
     jj = 1
     do ii=1,3
        print *, pp(jj:(jj+2))
        jj = jj + 3
     end do
     ! mass with bases
     kx = 0
     do ix=1,3
        do jx=1,3
           kx = kx + 1
           uu = 0d0
           uu(kx) = 1d0
           call fem_op_private( bb1, uu, ntri, trind, nnodes, crds, bd, 1 )           
           write(*,'(A,I2,A,I2,A,ES9.2E2,ES9.2E2,ES9.2E2,ES9.2E2,ES9.2E2,ES9.2E2,ES9.2E2,ES9.2E2,ES9.2E2,A,ES9.2E2)'),'M*u(',ix,',',jx,')=',bb1(1:9),', sum=',sum(bb1(1:9))
        end do
     end do
  end if

end subroutine test_fem
