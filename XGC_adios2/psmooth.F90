

!!$module psmth_module
!!$#include "finclude/petscdef.h"
!!$#if !defined(PETSC_31)
!!$#include "finclude/petscvecdef.h" 
!!$#include "finclude/petscmatdef.h" 
!!$#include "finclude/petsckspdef.h" 
!!$#endif
!!$  type psmth_type
!!$     Mat smoother, sm_projector
!!$     Vec tmp1
!!$  end type psmth_type
!!$  type(psmth_type) :: psmthH, psmthr, psmth00
!!$
!!$contains
!!$  subroutine psmth_setup(psmth, grid, mode, istart, pnloc,smooth, smooth_r)
!!$    use smooth_module    
!!$    use eq_module
!!$    use grid_class
!!$    use sml_module
!!$    implicit none
!!$#include "finclude/petsc.h"
!!$#if !defined(PETSC_31)
!!$#include "finclude/petscvec.h" 
!!$#include "finclude/petscmat.h" 
!!$#endif
!!$    type(psmth_type) :: psmth
!!$    type(grid_type) :: grid
!!$    type(smooth_type) :: smooth
!!$    type(smooth_r_type) :: smooth_r
!!$!    type(smooth_nearx_type) :: smooth_nearx
!!$    integer, intent(in) :: mode, istart, pnloc
!!$    !
!!$    integer :: i,j,j2,nfs,n1,nnodes,nsmnodes,nqloc,nqtot,ierr,nd1,nd2,nnz,onnz
!!$    integer :: nlocfsi,nlocalsmooth, iqeq,icolid
!!$    integer :: icolid2, cnd, rl
!!$    real (8) :: asum,ar, wsum(2)
!!$    real (8), parameter :: one=1D0
!!$    !
!!$    nnodes = grid%nnode
!!$    nsmnodes =  grid%itheta0(grid%npsi) - 1 
!!$    ! get number of flux surfaces that we are doing flux surface averaging
!!$    nnz = 180; onnz = 90 
!!$    if(mode == 0) then
!!$       nfs = grid%npsi - 1
!!$       n1 = nsmnodes
!!$       nnz=1000
!!$       onnz=10
!!$    elseif(mode>= 2) then
!!$       nfs=0
!!$       n1 = 0
!!$    else ! mode == 1
!!$       do i=1,grid%npsi-1
!!$          nd1 = grid%itheta0(i) 
!!$          if( grid%psi(nd1) < 0.999D0*eq_x_psi ) then
!!$          else
!!$             nfs = i - 1
!!$             n1 = nd1 - 1 
!!$             exit
!!$          endif
!!$       enddo
!!$       if(sml_plane_mype==0) nnz = 1000
!!$    endif
!!$
!!$call MPI_Barrier(sml_comm,ierr)
!!$if(sml_mype==0) print *,sml_mype, ') n1=',n1, ' nsmnodes=',nsmnodes, ' nnodes=',nnodes
!!$
!!$    ! get number of smooth nodes local
!!$    nqtot = nfs + (nnodes - n1)
!!$    if(istart >= n1) then 
!!$       nqloc = pnloc
!!$    elseif( istart + pnloc -1 >= n1 ) then ! iend = istart + pnloc -1
!!$       nqloc = istart + pnloc - n1 
!!$    else
!!$       nqloc=0
!!$    endif
!!$    nlocalsmooth = nqloc
!!$    ! add in number of fux surfeace equations local
!!$    nlocfsi = 0
!!$    if(sml_plane_mype==0) then
!!$       nqloc = nqloc + nfs
!!$       nlocfsi = nfs
!!$    endif
!!$call MPI_Barrier(sml_comm,ierr)
!!$if(sml_mype==0) print *,sml_mype, ') nqloc=',nqloc,nqtot
!!$    ! allocate data structures
!!$      if(nnz.gt.pnloc) nnz = pnloc
!!$      if(onnz.gt.(nnodes-pnloc)) onnz = nnodes-pnloc
!!$    call MatCreateMPIAIJ( sml_plane_comm, nqloc, pnloc, nqtot, nnodes, &
!!$      nnz, PETSC_NULL_INTEGER, onnz, PETSC_NULL_INTEGER, psmth%smoother, ierr )
!!$    call MatCreateMPIAIJ( sml_plane_comm, nqloc, pnloc, nqtot, nnodes, &
!!$         nnz, PETSC_NULL_INTEGER, onnz, PETSC_NULL_INTEGER, psmth%sm_projector, ierr )
!!$    call VecCreateMPI(sml_plane_comm, nqloc, nqtot, psmth%tmp1, ierr )
!!$    ! create matrices -- flux surface average part
!!$    icolid = 0 ! all one proc 0
!!$
!!$call MPI_Barrier(sml_comm,ierr)
!!$if(sml_mype==0) print *,sml_mype, ') MatSetValues ...'
!!$
!!$    do i=1,nlocfsi
!!$       iqeq = i - 1
!!$       nd1=grid%itheta0(i) 
!!$       nd2=grid%itheta0(i)+grid%ntheta(i)-1
!!$       asum = sum( grid%node_area(nd1:nd2) )
!!$if(sml_mype==0) print *,sml_mype, ') setting row ',i
!!$       do j = nd1,nd2
!!$          ar = grid%node_area(j)/asum
!!$          call MatSetValues(psmth%smoother, 1, iqeq, 1, icolid, ar, ADD_VALUES,ierr)
!!$          call MatSetValues(psmth%sm_projector,1,iqeq,1,icolid, one,ADD_VALUES,ierr)
!!$          icolid = icolid + 1
!!$       enddo
!!$    enddo
!!$print *,sml_mype, ') mat set values done'
!!$call MPI_Barrier(sml_comm,ierr)
!!$if(sml_mype==0) print *,sml_mype, ') ***************************** nlocalsmooth=',nlocalsmooth
!!$    ! "smoothing" part
!!$    do j=1,nlocalsmooth
!!$       icolid = (istart+pnloc-nlocalsmooth) + j - 1
!!$       iqeq = nfs + (icolid - n1)
!!$       call MatSetValues(psmth%sm_projector,1,iqeq,1,icolid,one,ADD_VALUES,ierr)
!!$       ! smoother
!!$       if( icolid < grid%itheta0(grid%npsi)-1 ) then
!!$          !if(sml_mype==0)print*,'smoothing eq:',icolid
!!$          if(mode==0) then
!!$             print *, 'real smoothing mode==0, icolid', icolid, grid%itheta0(grid%npsi)-1
!!$             stop
!!$          elseif(mode==2.or.mode==1) then
!!$             wsum=0D0
!!$             do rl=1,2
!!$                cnd=icolid+1  ! make one based
!!$                NODES : do i=1,smooth%n
!!$                   wsum(rl)=wsum(rl)+smooth%weight(i)*grid%node_area(cnd)
!!$                   !find next node
!!$                   cnd=grid%nn(rl,cnd)
!!$                   if(cnd <=0 .or. cnd > grid%nnode) exit NODES
!!$                enddo NODES
!!$             enddo
!!$             do rl=1,2
!!$                cnd=icolid+1  ! make one based
!!$                NODES2 : do i=1,smooth%n             
!!$                   icolid2=cnd-1 ! make zero based
!!$                   ar=smooth%weight(i)*grid%node_area(cnd)/(2D0*wsum(rl))
!!$                   call MatSetValues(psmth%smoother,1,iqeq,1,icolid2,ar,ADD_VALUES,ierr)
!!$                   !find next node
!!$                   cnd=grid%nn(rl,cnd)
!!$                   if(cnd <=0 .or. cnd > grid%nnode) exit NODES2
!!$                enddo NODES2
!!$             enddo
!!$          elseif(mode==3) then
!!$             cnd=icolid+1  ! make one based
!!$             do j2=1,smooth_r%mat%nelement(cnd)
!!$                icolid2=smooth_r%mat%eindex(j2,cnd)-1 ! make zero based
!!$                ar=smooth_r%mat%value(j2,cnd)
!!$                call MatSetValues(psmth%smoother,1,iqeq,1,icolid2,ar,ADD_VALUES,ierr) 
!!$             enddo
!!$          else
!!$             print *, 'real smoothing mode==',mode
!!$             stop 
!!$          endif
!!$       else
!!$          ! Wall 
!!$          call MatSetValues(psmth%smoother,1,iqeq,1,icolid,one,ADD_VALUES,ierr)
!!$       endif
!!$    enddo
!!$call MPI_Barrier(sml_comm,ierr)
!!$if(sml_mype==0) print *,sml_mype, ') mat set values DONE'
!!$    call MatAssemblyBegin(psmth%smoother, MAT_FINAL_ASSEMBLY, ierr )
!!$    call MatAssemblyBegin(psmth%sm_projector, MAT_FINAL_ASSEMBLY, ierr )
!!$    call MatAssemblyEnd(psmth%smoother, MAT_FINAL_ASSEMBLY, ierr )
!!$    call MatAssemblyEnd(psmth%sm_projector, MAT_FINAL_ASSEMBLY, ierr )
!!$call MPI_Barrier(sml_comm,ierr)
!!$if(sml_mype==0)print *, 'assembly DONED'
!!$  end subroutine psmth_setup
!!$  !
!!$  !
!!$  !
!!$  subroutine psmth_apply(psmth,den)
!!$    use sml_module
!!$    implicit none
!!$#include "finclude/petsc.h"
!!$#include "finclude/petscviewer.h"
!!$    type(psmth_type) :: psmth    
!!$    Vec den
!!$    !
!!$    integer ierr,nn,mm
!!$    PetscViewer  viewer
!!$
!!$    if(sml_plane_index == -1)then
!!$       print *, 'printing *******************************', sml_mype
!!$       call PetscViewerASCIIOpen(sml_plane_comm,"mat.output0",viewer,ierr);
!!$       !call PetscViewerMatLabOpen(sml_plane_comm,"mat.m",FILE_MODE_WRITE,viewer,ierr);
!!$       call MatView(psmth%smoother, viewer,ierr)
!!$       call PetscViewerDestroy(viewer,ierr)
!!$       !
!!$       call PetscViewerASCIIOpen(sml_plane_comm,'vec.output0',viewer,ierr);
!!$       call VecView(den,viewer,ierr)
!!$       call PetscViewerDestroy(viewer,ierr)
!!$       
!!$       call MatGetSize(psmth%smoother,mm,nn,ierr)
!!$       print *, sml_mype,') mat global N=',mm,nn
!!$       call MatGetLocalSize(psmth%smoother,mm,nn,ierr)
!!$       print *, sml_mype,') mat local n=',mm,nn
!!$       call VecGetSize(den,mm,ierr)
!!$       print *, sml_mype,') den vec global N=',mm
!!$       call VecGetLocalSize(den,mm,ierr)
!!$       print *, sml_mype,') den vec local n=',mm
!!$       call VecGetSize(psmth%tmp1,mm,ierr)
!!$       print *, sml_mype,') tmp vec global N=',mm
!!$       call VecGetLocalSize(psmth%tmp1,mm,ierr)
!!$       print *, sml_mype,') tmp vec local n=',mm
!!$    endif
!!$    !call check_point2('---')
!!$
!!$    call MatMult( psmth%smoother, den, psmth%tmp1, ierr )
!!$    call MatMultTranspose( psmth%sm_projector, psmth%tmp1, den, ierr )
!!$
!!$    if(sml_plane_index == -1)then
!!$       print *, 'printing *******************************', sml_mype
!!$       call PetscViewerASCIIOpen(sml_plane_comm,"mat.output",viewer,ierr);
!!$       !call PetscViewerMatLabOpen(sml_plane_comm,"mat.m",FILE_MODE_WRITE,viewer,ierr);
!!$       call MatView(psmth%smoother, viewer,ierr)
!!$       call PetscViewerDestroy(viewer,ierr)
!!$       !
!!$       call PetscViewerASCIIOpen(sml_plane_comm,'vec.output',viewer,ierr);
!!$       call VecView(den,viewer,ierr)
!!$       call PetscViewerDestroy(viewer,ierr)
!!$    endif
!!$
!!$  end subroutine psmth_apply
!!$  !
!!$  !
!!$  !
!!$  subroutine psmth_destroy(psmth)
!!$    implicit none
!!$    type(psmth_type) :: psmth    
!!$    !
!!$    integer ierr
!!$    call MatDestroy(psmth%smoother,ierr)
!!$    call MatDestroy(psmth%sm_projector,ierr)
!!$    call VecDestroy(psmth%tmp1,ierr)
!!$  end subroutine psmth_destroy
!!$end module psmth_module
!!$
