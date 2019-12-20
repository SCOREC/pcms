#ifdef USE_GPU
attributes(global) &
#endif
subroutine gen_perm_tri_gpu_pass1(n,ntriangle,ncopy,gpu_ibegin,gpu_iend,lcount,tr_save_gpu)
  use cudafor
  !use dimensions_mod_gpu !, only : nthreads_dim
  use ptl_module_gpu, only : ptl_nphase, ptl_gid_gpu, ptl_ph_gpu
  use precision_mod_gpu
  use sml_module_gpu, only : sml_2pi
  implicit none
  integer, value, intent(in) :: n,ntriangle,gpu_ibegin,gpu_iend,ncopy
  integer, intent(in) :: tr_save_gpu(n)
  integer, intent(inout) :: lcount(ntriangle,ncopy)
  real(kind=work_p) :: p(3), xff(2)
  integer :: ith, nthread, iblock, nblock
  integer :: gith, gnthread, pith
  integer :: i, k, idummy,icopy, itr
#ifdef USE_GPU
  attributes(device) :: lcount
#endif
!      ------------------------------------
!      assume lcount(nx*ny,:) to be zero
!      ------------------------------------

  ith = threadIdx%x + (threadIdx%y - 1)*blockDim%x
  nthread = blockDim%x * blockDim%y
  iblock = blockIdx%x + (blockIdx%y-1)*gridDim%x

  nblock = gridDim%x * gridDim%y
  gnthread = nblock * nthread
  gith = ith + (iblock-1)*nthread

  icopy = 1 + mod(iblock-1,ncopy)

  do k=gith,n,gnthread
    pith = k + gpu_ibegin - 1
    if (ptl_gid_gpu(pith) <= 0) then
      itr = ntriangle
    else
! Check guess won't actually work here, since we don't do
! the field following. (and don't want to, because it's expensive).
! Instead, we'll check the saved triangle. If it doesn't exist, we'll
! Run the search and use that, but that triangle probably won't be in
! the correct place.
#ifdef USE_TR_CHECK
      itr = tr_save_gpu(k)
      if (itr <= 0) then
        call search_tr2_gpu(xff,itr,p)
      endif
#else
      call search_tr2_gpu(xff,itr,p)
#endif
      if (itr <= 0) then
        itr = ntriangle
      endif
    endif

    idummy = atomicadd( lcount(itr,icopy), 1 )
  enddo

return
end subroutine

#ifdef USE_GPU
attributes(global) &
#endif
subroutine gen_perm_tri_gpu_pass2(n,ntriangle,ncopy,gpu_ibegin,gpu_iend,lcount,iperm,tr_save_gpu)
  use cudafor
  !use dimensions_mod_gpu !, only : nthreads_dim
  use ptl_module_gpu, only : ptl_nphase, ptl_gid_gpu, ptl_ph_gpu
  use precision_mod_gpu
  implicit none
  integer, value, intent(in) :: n,ntriangle,gpu_ibegin,gpu_iend,ncopy
  integer, intent(inout) :: lcount(ntriangle,ncopy)
  integer, intent(in) :: tr_save_gpu(n)
  integer, intent(out) :: iperm(n)
  real(kind=work_p) :: p(3), xff(2)
  integer :: ith, nthread, iblock, nblock
  integer :: gith, gnthread, pith
  integer :: i, k, ip, idummy,icopy, itr

#ifdef USE_GPU
  attributes(device) :: lcount, iperm
#endif

!      ------------------------------------
!      assume lcount(nx*ny,:) to be zero
!      ------------------------------------

  ith = threadIdx%x + (threadIdx%y - 1)*blockDim%x
  nthread = blockDim%x * blockDim%y
  iblock = blockIdx%x + (blockIdx%y-1)*gridDim%x

  nblock = gridDim%x * gridDim%y
  gnthread = nblock * nthread
  gith = ith + (iblock-1)*nthread

  icopy = 1 + mod(iblock-1,ncopy)

  do k=gith,n,gnthread
    pith = k + gpu_ibegin - 1
    if (ptl_gid_gpu(pith) <= 0) then
      itr = ntriangle
    else
#ifdef USE_TR_CHECK
      itr = tr_save_gpu(k)
      if (itr <= 0) then
        call search_tr2_gpu(xff,itr,p)
      endif
#else
      call search_tr2_gpu(xff,itr,p)
#endif
      if (itr <= 0) then
        itr = ntriangle
      endif
      ! Figure out the phi plane
    endif

    ip = atomicadd( lcount(itr,icopy), 1 )
    iperm(ip) = pith
  enddo

return
end subroutine

subroutine gen_perm_tri_gpu(n,ntriangle,gpu_ibegin,gpu_iend,iperm,xstart)
       use cudafor
       use sml_module
       use ptl_module_gpu
       use precision_mod_gpu
       use gen_perm_gpu_mod
       implicit none

       integer, intent(in) :: gpu_ibegin, gpu_iend
       integer, intent(in) :: n, ntriangle
       integer, intent(inout) :: iperm(n)
       integer, intent(inout) :: xstart(ntriangle+1)


       integer, parameter :: nblock = 32
       integer(kind=cuda_count_kind) :: freemem, totalmem
       integer :: ncopy
       integer, parameter :: nthread = 32*6

       integer(kind=4), allocatable, dimension(:,:) :: lcount

#ifdef USE_GPU
       attributes(device) :: iperm, xstart
       attributes(device) :: lcount, current_loc
#endif

       integer, parameter :: block_limit = 2**15-1

       integer :: mm,nn,lda
       integer :: istat, istart_val
       
       integer(kind=4), allocatable, dimension(:,:) :: h_lcount
       integer :: i,j,k,iblock
 

!      ---------------------------------------------
!      1st pass to count number of particles in cell
!      ---------------------------------------------
       istat = cudaMemGetInfo(freemem,totalmem)
       if (istat.ne.0) then
         write(*,*) sml_mype, 'gen_perm_tri_gpu: get free mem return istat=',istat
         stop '** error ** '
       endif

       ncopy = min(32,freemem/4/ntriangle - 1)
       if (ncopy.eq.0) then
         write(*,*) sml_mype, 'gen_perm_tri_gpu: not enough free memory for lcount! ncopy:', ncopy
         write(*,*) 'ntriangle, ncopy ', ntriangle, ncopy
         stop '** error ** '
       endif

       allocate( lcount( ntriangle,ncopy ) , stat=istat)
       if (istat.ne.0) then
         write(*,*) sml_mype, 'gen_perm_tri_gpu: allocate lcount return istat=',istat
         write(*,*) 'ntriangle,  nblock ', ntriangle, nblock
         stop '** error ** '
       endif

       call isetval_gpu( size(lcount), lcount, 0)


#ifdef USE_GPU
       call gen_perm_tri_gpu_pass1<<<nblock,nthread>>>(                       &
     &               n,ntriangle,ncopy, &
     &               gpu_ibegin,gpu_iend,lcount,tr_save_gpu)
#else
       call gen_perm_tri_gpu_pass1(                       &
     &               n,ntriangle,ncopy, &
     &               gpu_ibegin,gpu_iend,lcount,tr_save_gpu)
#endif
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
       call setup_lcount_gpu<<<nblock,nthread>>>(ntriangle,1,&
     &                  ncopy,lcount,xstart)
#else
       call setup_lcount_gpu(ntriangle,1,ncopy,lcount,xstart)
#endif

!      --------------------------
!      2nd pass to populate iperm
!      --------------------------

#ifdef USE_GPU
       call gen_perm_tri_gpu_pass2<<<nblock,nthread>>>(                       &
     &               n,ntriangle,ncopy, &
     &               gpu_ibegin,gpu_iend,lcount,iperm,tr_save_gpu)
#else
       call gen_perm_tri_gpu_pass2(                       &
     &               n,ntriangle,ncopy, &
     &               gpu_ibegin,gpu_iend,lcount,iperm,tr_save_gpu)
#endif

       deallocate( lcount, stat=istat )
       if (istat.ne.0) then
          write(*,*) sml_mype, 'gen_perm_tri_gpu: dealloc(lcount),istat=',istat
          stop '** error ** '
       endif

       return
       end subroutine gen_perm_tri_gpu
