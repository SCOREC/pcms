       subroutine gen_perm( ilo,ihi,jlo,jhi,guess_min,inv_guess_d,        &
     &                n, gid, x,y, iperm )

       use sml_module, only : sml_nthreads
       use omp_module, only : split_indices

       implicit none
       integer, intent(in) :: ilo,ihi, jlo,jhi
       real(8), intent(in) :: guess_min(2)
       real(8), intent(in) :: inv_guess_d(2)

       integer, intent(in) :: n
       integer*8, intent(in) :: gid(n)
       real(8), intent(in) :: x(n), y(n)
       integer, intent(inout) :: iperm(n)

       integer :: ith,isize, i,j,k,  ifree,ipos
       real(8) :: guess_min_1, guess_min_2
       real(8) :: inv_guess_d_1, inv_guess_d_2

       integer :: i_beg(sml_nthreads),i_end(sml_nthreads)

       integer, allocatable, dimension(:,:,:) :: icount

       integer :: nerr
       integer, parameter :: idebug = 0
       integer :: icheck(n)


       guess_min_1 = guess_min(1)
       guess_min_2 = guess_min(2)
       inv_guess_d_1 = inv_guess_d(1)
       inv_guess_d_2 = inv_guess_d(2)

       call split_indices( n, sml_nthreads, i_beg, i_end)

       allocate( icount(ilo:ihi,jlo:jhi,sml_nthreads ) )

!$omp  parallel do private(ith)
       do ith=1,sml_nthreads
          icount(:,:,ith) = 0
       enddo


!      ----------------------------------------
!      first pass to count number of particles
!      in each cell
!      ----------------------------------------

!$omp  parallel do private(ith,k,i,j)
       do ith=1,sml_nthreads
        do k=i_beg(ith),i_end(ith)
         if (gid(k) > 0) then
           i = max(ilo,min(ihi,                                          &
     &            int((x(k) - guess_min_1) * inv_guess_d_1 + 1)))

           j = max(jlo,min(jhi,                                          &
     &             int((y(k) - guess_min_2) * inv_guess_d_2 + 1)))
         else
           i = ihi
           j = jhi
         endif

         icount(i,j,ith) = icount(i,j,ith) + 1
        enddo
       enddo

!      ---------------------
!      setup pointer 
!      similar to prefix sum
!      ---------------------

       ifree = 1
       do j=jlo,jhi
       do i=ilo,ihi
         do ith=1,sml_nthreads
           isize = icount(i,j,ith)
           icount(i,j,ith) = ifree
           ifree = ifree + isize
         enddo
       enddo
       enddo

!      ----------------------
!      second pass to setup perm
!      ----------------------

!$omp  parallel do private(ith,k,i,j,ipos)
       do ith=1,sml_nthreads
        do k=i_beg(ith),i_end(ith)
         if (gid(k) > 0) then
           i = max(ilo,min(ihi,                                          &
     &             int((x(k) - guess_min_1) * inv_guess_d_1 + 1 )))
           j = max(jlo,min(jhi,                                          &
     &             int((y(k) - guess_min_2) * inv_guess_d_2 + 1 )))
         else
!          ----------------------------------------
!          try to map invalid particles at the end
!          ----------------------------------------
           i = ihi
           j = jhi
         endif

         ipos = icount(i,j,ith)
         iperm(ipos) = k
         icount(i,j,ith) = ipos + 1
        enddo
       enddo

       deallocate( icount )

!      ------------
!      double check
!      ------------

       if (idebug >= 1) then

       icheck(1:n) = 0
       nerr = 0
       do i=1,n
         if ((1 <= iperm(i)).and.(iperm(i) <= n)) then
             icheck(iperm(i)) = icheck(iperm(i)) + 1
         else
             write(*,*) 'iperm out of range '
             write(*,*) 'i, iperm(i) ', i, iperm(i)
             nerr = nerr + 1
         endif
       enddo

       do i=1,n
         if (icheck(i) .ne. 1) then
           write(*,*) 'invalid iperm value '
           write(*,*) 'i, icheck(i) ', i,icheck(i)
           nerr = nerr + 1
         endif
       enddo

       if (nerr.ne.0) then
         stop '** error ** '
       endif

       endif


       return
       end subroutine gen_perm
