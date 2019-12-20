       subroutine gen_perm_gpu( nx, ny, xmin,ymin, inv_dx,inv_dy,        &
     &                          n, gid, xy, xydim, iperm, xstart )
       use cudafor
       use sml_module
       use precision_mod_gpu
       implicit none

       integer, intent(in) :: nx, ny
       real(kind=work_p), intent(in) :: xmin,ymin, inv_dx,inv_dy
       integer, intent(in) :: n, xydim
       integer*8, intent(in) :: gid(n)
       real(kind=work_p) :: xy(xydim, 2)
       integer, intent(inout) :: iperm(n)
       integer, intent(inout) :: xstart(nx*ny+1)


       integer, parameter :: nblock = 32
       integer, parameter :: ncopy = nblock
       integer, parameter :: nthread = 32*6

       integer(kind=4), allocatable, dimension(:,:) :: lcount

#ifdef USE_GPU
       attributes(device) :: gid, xy, iperm, xstart
       attributes(device) :: lcount
#endif

       integer, parameter :: block_limit = 2**15-1

       integer :: mm,nn,lda
       integer :: istat, istart_val
       
       integer, parameter :: idebug = 0
       integer(kind=4), allocatable, dimension(:,:) :: h_lcount
       integer :: i,j,k,iblock
 

!      ---------------------------------------------
!      1st pass to count number of particles in cell
!      ---------------------------------------------
       allocate( lcount( nx*ny,ncopy ) , stat=istat)
       if (istat.ne.0) then
         write(*,*) sml_mype, 'gen_perm_gpu: allocate lcount return istat=',istat
         write(*,*) 'nx, ny, nblock ', nx, ny, nblock
         stop '** error ** '
       endif

       call isetval_gpu( size(lcount), lcount, 0)


#ifdef USE_GPU
       call gen_perm_gpu_pass1<<<nblock,nthread>>>(                       &
     &               nx,ny, xmin,ymin,                                    &
     &               inv_dx,inv_dy, n,gid,xy,xydim, ncopy, lcount)
#else
       call gen_perm_gpu_pass1(                                           &
     &               nx,ny,xmin,ymin,                                     &
     &               inv_dx,inv_dy, n,gid,xy,xydim, ncopy, lcount)
#endif

       if (idebug >= 2) then
          allocate( h_lcount(nx*ny,ncopy) )
          h_lcount = lcount

!          print*,'after gen_perm_gpu_pass1 '
          do iblock=1,ncopy
            do j=1,ny
            do i=1,nx
              k = i + (j-1)*nx
!              print 9010,i,j,iblock,h_lcount(k,iblock)
 9010         format(1x,'h_lcount(',i8,',',i8,',',i8,') = ',i9 )
            enddo
            enddo
          enddo
!          print*,'sum(h_lcount) ', sum(h_lcount)
          deallocate( h_lcount )
         endif

!      ------------------------
!      setup pointers in xstart
!      ------------------------

!      ---------------------------------------
!      compute xstart(k)  = sum( lcount(k,:) )
!      ---------------------------------------

       mm = size(lcount,1)
       nn = size(lcount,2)
       lda = size(lcount,1)

#ifdef USE_GPU
       call isum_col_gpu<<<nblock,nthread>>>(mm,nn,lda, lcount, xstart)
#else
       call isum_col_gpu(mm,nn,lda, lcount, xstart)
#endif
!      --------------------------------------------
!      compute xstart(:) as prefix sum on icount(:)
!      --------------------------------------------
       istart_val = 1
       mm = size(xstart)-1
       call iprefix_sum_gpu( mm, xstart, istart_val)
     

!      ------------------------------------------
!      setup lcount(:,:) as pointers into iperm
!      ------------------------------------------
#ifdef USE_GPU
       call setup_lcount_gpu<<<nblock,nthread>>>(nx,ny,                 &
     &                  ncopy,lcount,xstart)
#else
       call setup_lcount_gpu(nx,ny,ncopy,lcount,xstart)
#endif

!      --------------------------
!      2nd pass to populate iperm
!      --------------------------

#ifdef USE_GPU
       call gen_perm_gpu_pass2<<<nblock,nthread>>>( nx,ny,xmin,ymin,      &
     &       inv_dx,inv_dy,n,gid,xy,xydim,ncopy,lcount,iperm,xstart )
#else
       call gen_perm_gpu_pass2( nx,ny,xmin,ymin,                          &
     &       inv_dx,inv_dy,n,gid,xy,xydim,ncopy,lcount,iperm,xstart )
#endif


       deallocate( lcount, stat=istat )
       if (istat.ne.0) then
          write(*,*) sml_mype, 'gen_perm_gpu: dealloc(lcount),istat=',istat
          stop '** error ** '
       endif
       return
       end subroutine gen_perm_gpu
