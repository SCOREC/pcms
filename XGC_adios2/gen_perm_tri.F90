! gen_perm_tri : genarate the permutation vector based on triangle number
       subroutine gen_perm_tri( ptl, grid, iperm, n )

       use sml_module, only : sml_nthreads
       use omp_module, only : split_indices
       use grid_class
       use ptl_module
       implicit none

       type(grid_type), intent(in) :: grid
       type(ptl_type), intent(in) :: ptl(n)
       integer, intent(out) :: iperm(n)
       integer, intent(in) :: n
       integer :: itr, ith,isize, i,j,k,  ifree,ipos
       real(8) :: p(3)
       integer :: i_beg(sml_nthreads),i_end(sml_nthreads)

       ! Since tr_save might be used on the cpu, 
       ! we don't want to overwrite it
       integer, allocatable, dimension(:) :: current_tr
       integer, allocatable, dimension(:,:) :: icount

       integer(8) :: gauss, checksum
       integer :: nerr
       integer, parameter :: idebug = 0
       integer :: icheck(n)

       call split_indices( n, sml_nthreads, i_beg, i_end)

       allocate( icount(grid%ntriangle,sml_nthreads) , current_tr(n))

!$omp  parallel do private(ith)
       do ith=1,sml_nthreads
          icount(:,ith) = 0
       enddo


!      ----------------------------------------
!      First pass
!      Get the current triangle, and count # ptl in that triangle
!      (could be sped up it we read check tr_save first)
!      ----------------------------------------

!$omp  parallel do private(ith,k,itr)
       do ith=1,sml_nthreads
        do k=i_beg(ith),i_end(ith)
          call search_tr2(grid, ptl(k)%ph(1:2), itr, p)
          if (itr > 0 .and. ptl(k)%gid > 0) then
            current_tr(k) = itr
          else
            current_tr(k) = grid%ntriangle
          endif
          icount(current_tr(k), ith) = icount(current_tr(k), ith) + 1
        enddo
      enddo

!      ---------------------
!      setup pointer 
!      similar to prefix sum
!      ---------------------

       ifree = 1
       do i=1,grid%ntriangle
         do ith=1,sml_nthreads
           isize = icount(i,ith)
           icount(i,ith) = ifree
           ifree = ifree + isize
         enddo
       enddo

!      ----------------------
!      second pass to setup perm
!      ----------------------

!$omp  parallel do private(ith,k,ipos)
       do ith=1,sml_nthreads
        do k=i_beg(ith),i_end(ith)
         ipos = icount(current_tr(k),ith)
         iperm(ipos) = k
         icount(current_tr(k),ith) = ipos + 1
        enddo
       enddo

       deallocate( icount, current_tr)

!      ------------
!      double check
!      ------------

       if (idebug > 0) then
         checksum = sum(int(iperm,8)) 
         gauss = (int(n,8) * int((n+1),8)) / 2
         print *, "Iperm checksum:", checksum
         print *, "Gauss check: ", gauss

         if ((idebug > 1) .or. (checksum /= gauss)) then
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
      endif

     return
     end subroutine gen_perm_tri
