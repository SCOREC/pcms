

        module perm_module
        implicit none
! -------------------------------------------------------
! module to perform copy with permutation
! operation similar to
!
! Atmp(:,:) = A(:,:)
! A(:,:) = Atmp(:,iperm(:))
!
! other routines to generate cycles from permutation and
! perform permuted copy "in place" without another copy
!
! -------------------------------------------------------
        public

        interface permute_inplace
        module procedure &
     & permute_d1, permute_d2, permute_d3, &
     & permute_i1, permute_i2, permute_i3, &
     & permute_j1, permute_j2, permute_j3, &
     & perm_inplace_d1, perm_inplace_d2, perm_inplace_d3, &
     & perm_inplace_i1, perm_inplace_i2, perm_inplace_i3, &
     & perm_inplace_j1, perm_inplace_j2, perm_inplace_j3
        end interface

! --------------------------------------
! use temporary storage in heap
! call in place version if memory is low
! --------------------------------------
        interface permcopy
        module procedure &
     & permcopy_d1, permcopy_d2, permcopy_d3, &
     & permcopy_i1, permcopy_i2, permcopy_i3, &
     & permcopy_j1, permcopy_j2, permcopy_j3
        end interface

! ------------------------------
! use temporary storage on stack
! ------------------------------
        interface permscopy
        module procedure &
     & permscopy_d1, permscopy_d2, permscopy_d3, &
     & permscopy_i1, permscopy_i2, permscopy_i3, &
     & permscopy_j1, permscopy_j2, permscopy_j3
        end interface


        interface copy
        module procedure copyi, copyj, copys, copyd, copyc, copyz
        end interface copy


        contains
        subroutine assert(isok,mesg,ivalue)
        implicit none
        logical isok
        character*(*) mesg
        integer ivalue
        if (.not.isok) then
           write(*,*) mesg,ivalue
           stop '**** assertion failed **** '
        endif
        return
        end subroutine
        subroutine copyi(n,x,y)
        implicit none
        integer, intent(in) :: n
        integer, parameter :: i4 = selected_int_kind(8)
        integer(kind=i4), dimension(*), intent(in) :: x
        integer(kind=i4), dimension(*), intent(inout) :: y
        integer :: i
        do i=1,n
          y(i) = x(i)
        enddo
        end subroutine copyi
        subroutine copyj(n,x,y)
        implicit none
        integer, intent(in) :: n
        integer, parameter :: i8 = selected_int_kind(12)
        integer(kind=i8), dimension(*), intent(in) :: x
        integer(kind=i8), dimension(*), intent(inout) :: y
        integer :: i
        do i=1,n
          y(i) = x(i)
        enddo
        end subroutine copyj
        subroutine copys(n,x,y)
        implicit none
        integer, intent(in) :: n
        integer, parameter :: sp = kind(1.0)
        real(kind=sp), dimension(*), intent(in) :: x
        real(kind=sp), dimension(*), intent(inout) :: y
        integer :: i
        do i=1,n
          y(i) = x(i)
        enddo
        end subroutine copys
        subroutine copyd(n,x,y)
        implicit none
        integer, intent(in) :: n
        integer, parameter :: dp = kind(1.0d0)
        real(kind=dp), dimension(*), intent(in) :: x
        real(kind=dp), dimension(*), intent(inout) :: y
        integer :: i
        do i=1,n
          y(i) = x(i)
        enddo
        end subroutine copyd
        subroutine copyc(n,x,y)
        implicit none
        integer, intent(in) :: n
        integer, parameter :: sp = kind(1.0)
        complex(kind=sp), dimension(*), intent(in) :: x
        complex(kind=sp), dimension(*), intent(inout) :: y
        integer :: i
        do i=1,n
          y(i) = x(i)
        enddo
        end subroutine copyc
        subroutine copyz(n,x,y)
        implicit none
        integer, intent(in) :: n
        integer, parameter :: dp = kind(1.0d0)
        complex(kind=dp), dimension(*), intent(in) :: x
        complex(kind=dp), dimension(*), intent(inout) :: y
        integer :: i
        do i=1,n
          y(i) = x(i)
        enddo
        end subroutine copyz
        subroutine gencycle(n, iperm, ncycle, leaders,ldim )
        implicit none
        integer, intent(in) :: n
        integer, dimension(n), intent(inout) :: iperm
        integer, intent(inout) :: ncycle
        integer, intent(inout) :: ldim
        integer, dimension(*), intent(inout) :: leaders
        integer, parameter :: idebug = 0
        logical :: isok, has_touched
        integer :: istart, i, j
        integer :: total_swap, max_swap, nswap
        ncycle = 0
        total_swap = 0
        max_swap = 0
! ---------------------------------
! use sign bit of iperm has marker
! in finding leaders of cycles
!
! initially set sign bit
! ---------------------------------
        do i=1,n
          iperm(i) = -iperm(i)
        enddo
        isok = .true.
        do i=1,n
           has_touched = (iperm(i) .gt. 0)
           if (.not. has_touched) then
                  iperm(i) = abs(iperm(i))
                  istart = i
                  j = iperm(i)
                  if (j .ne. i) then
                    ncycle = ncycle + 1
                    isok = (ncycle.le.ldim)
                    if (isok) then
                      leaders(ncycle) = i
                    endif
                   nswap = 0
                   do while (j .gt. istart)
                     nswap = nswap + 1
                     iperm(j) = abs(iperm(j))
                     j = iperm(j)
                   enddo
                   max_swap = max( max_swap, nswap)
                   total_swap = total_swap + nswap
                 endif
           endif
        enddo
        if (.not.isok) then
                ldim = -ncycle
        endif
        if (idebug.ge.1) then
                write(*,*) 'gencycle: total_swap, max_swap ', &
     & total_swap,max_swap
        endif
        return
        end subroutine
        subroutine permute_d1(n, iperm, ncycle, leaders, A )
        implicit none
        integer, intent(in) :: n
        real(kind=8), dimension(n), intent(inout) :: A
        integer, dimension(n), intent(in) :: iperm
        integer, intent(in) :: ncycle
        integer, dimension(n), intent(in) :: leaders
        integer :: k, i,j,istart
        real(kind=8) :: Ai
        do k=1,ncycle
          istart = leaders(k)
          i = istart
          j = iperm(istart)
          Ai = A(i)
          do while (j .ne. istart)
             A(i) = A(j)
             i = j
             j = iperm(i)
          enddo
          A(i) = Ai
         enddo
         return
         end subroutine
        subroutine permute_d2(n, iperm, ncycle, leaders, A )
        implicit none
        integer, intent(in) :: n
        real(kind=8), dimension(:,:), intent(inout) :: A
        integer, dimension(n), intent(in) :: iperm
        integer, intent(in) :: ncycle
        integer, dimension(n), intent(in) :: leaders
        integer :: k, i,j,istart
        real(kind=8), dimension( size(A,1)) :: Ai
        do k=1,ncycle
          istart = leaders(k)
          i = istart
          j = iperm(istart)
          Ai = A(:,i)
          do while (j .ne. istart)
             A(:,i) = A(:,j)
             i = j
             j = iperm(i)
          enddo
          A(:,i) = Ai
         enddo
         return
         end subroutine
        subroutine permute_d3(n, iperm, ncycle, leaders, A )
        implicit none
        integer, intent(in) :: n
        real(kind=8), dimension(:,:, :), intent(inout) :: A
        integer, dimension(n), intent(in) :: iperm
        integer, intent(in) :: ncycle
        integer, dimension(n), intent(in) :: leaders
        integer :: k, i,j,istart
        real(kind=8), dimension(size(A,1),size(A,2)) :: Ai
        do k=1,ncycle
          istart = leaders(k)
          i = istart
          j = iperm(istart)
          Ai = A(:,:,i)
          do while (j .ne. istart)
             A(:,:,i) = A(:,:,j)
             i = j
             j = iperm(i)
          enddo
          A(:,:,i) = Ai
         enddo
         return
         end subroutine
        subroutine permute_i1(n, iperm, ncycle, leaders, A )
        implicit none
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: iperm
        integer, intent(in) :: ncycle
        integer, dimension(n), intent(in) :: leaders
        integer, parameter :: i4 = selected_int_kind(8)
        integer(kind=i4), dimension(n), intent(inout) :: A
        integer :: k, i,j,istart
        integer(kind=i4) :: Ai
        do k=1,ncycle
          istart = leaders(k)
          i = istart
          j = iperm(istart)
          Ai = A(i)
          do while (j .ne. istart)
             A(i) = A(j)
             i = j
             j = iperm(i)
          enddo
          A(i) = Ai
         enddo
         return
         end subroutine
        subroutine permute_i2(n, iperm, ncycle, leaders, A )
        implicit none
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: iperm
        integer, intent(in) :: ncycle
        integer, dimension(n), intent(in) :: leaders
        integer, parameter :: i4 = selected_int_kind(8)
        integer(kind=i4), dimension(:,:), intent(inout) :: A
        integer :: k, i,j,istart
        integer(kind=i4), dimension(size(A,1)) :: Ai
        do k=1,ncycle
          istart = leaders(k)
          i = istart
          j = iperm(istart)
          Ai(:) = A(:,i)
          do while (j .ne. istart)
             A(:,i) = A(:,j)
             i = j
             j = iperm(i)
          enddo
          A(:,i) = Ai(:)
         enddo
         return
         end subroutine
        subroutine permute_i3(n, iperm, ncycle, leaders, A )
        implicit none
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: iperm
        integer, intent(in) :: ncycle
        integer, dimension(n), intent(in) :: leaders
        integer, parameter :: i4 = selected_int_kind(8)
        integer(kind=i4), dimension(:,:,:), intent(inout) :: A
        integer :: k, i,j,istart
        integer(kind=i4), dimension(size(A,1),size(A,2)) :: Ai
        do k=1,ncycle
          istart = leaders(k)
          i = istart
          j = iperm(istart)
          Ai = A(:,:,i)
          do while (j .ne. istart)
             A(:,:,i) = A(:,:,j)
             i = j
             j = iperm(i)
          enddo
          A(:,:,i) = Ai
         enddo
         return
         end subroutine
        subroutine permute_j1(n, iperm, ncycle, leaders, A )
        implicit none
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: iperm
        integer, intent(in) :: ncycle
        integer, dimension(n), intent(in) :: leaders
        integer, parameter :: i8 = selected_int_kind(12)
        integer(kind=i8), dimension(n), intent(inout) :: A
        integer :: k, i,j,istart
        integer(kind=i8) :: Ai
        do k=1,ncycle
          istart = leaders(k)
          i = istart
          j = iperm(istart)
          Ai = A(i)
          do while (j .ne. istart)
             A(i) = A(j)
             i = j
             j = iperm(i)
          enddo
          A(i) = Ai
         enddo
         return
         end subroutine
        subroutine permute_j2(n, iperm, ncycle, leaders, A )
        implicit none
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: iperm
        integer, intent(in) :: ncycle
        integer, dimension(n), intent(in) :: leaders
        integer, parameter :: i8 = selected_int_kind(12)
        integer(kind=i8), dimension(:,:), intent(inout) :: A
        integer :: k, i,j,istart
        integer(kind=i8), dimension(size(A,1)) :: Ai
        do k=1,ncycle
          istart = leaders(k)
          i = istart
          j = iperm(istart)
          Ai(:) = A(:,i)
          do while (j .ne. istart)
             A(:,i) = A(:,j)
             i = j
             j = iperm(i)
          enddo
          A(:,i) = Ai(:)
         enddo
         return
         end subroutine
        subroutine permute_j3(n, iperm, ncycle, leaders, A )
        implicit none
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: iperm
        integer, intent(in) :: ncycle
        integer, dimension(n), intent(in) :: leaders
        integer, parameter :: i8 = selected_int_kind(12)
        integer(kind=i8), dimension(:,:,:), intent(inout) :: A
        integer :: k, i,j,istart
        integer(kind=i8), dimension(size(A,1),size(A,2)) :: Ai
        do k=1,ncycle
          istart = leaders(k)
          i = istart
          j = iperm(istart)
          Ai(:,:) = A(:,:,i)
          do while (j .ne. istart)
             A(:,:,i) = A(:,:,j)
             i = j
             j = iperm(i)
          enddo
          A(:,:,i) = Ai(:,:)
         enddo
         return
         end subroutine
        subroutine perm_inplace_d1(n, iperm, A)
        implicit none
        integer, intent(in) :: n
        integer, dimension(n), intent(inout) :: iperm
        real(kind=8), dimension(n), intent(inout) :: A
        integer :: i, j, k, istart
        logical :: has_touched_k
        real(kind=8) :: Ai
! ----------------------------------
! use the sign bit in iperm(:) array
! to help find leaders of cycles
! ----------------------------------
        iperm = -iperm
        do k=1,n
          has_touched_k = (iperm(k) .gt. 0)
          if (.not. has_touched_k) then
                  iperm(k) = abs(iperm(k))
                  istart = k
                  i = istart
                  j = iperm(i)
                  Ai = A(i)
                  do while (j .gt. istart)
                     A(i) = A(j)
                     iperm(j) = abs(iperm(j))
                     i = j
                     j = iperm(j)
                   enddo
                   A(i) = Ai
           endif
         enddo
         return
         end subroutine
        subroutine perm_inplace_d2(n, iperm, A)
        implicit none
        integer, intent(in) :: n
        integer, dimension(n), intent(inout) :: iperm
        real(kind=8), dimension(:,:), intent(inout) :: A
        integer :: i, j, k, istart
        logical :: has_touched_k
        real(kind=8), dimension(size(A,1)) :: Ai
! ----------------------------------
! use the sign bit in iperm(:) array
! to help find leaders of cycles
! ----------------------------------
        iperm = -iperm
        do k=1,n
          has_touched_k = (iperm(k) .gt. 0)
          if (.not. has_touched_k) then
                  iperm(k) = abs(iperm(k))
                  istart = k
                  i = istart
                  j = iperm(i)
                  Ai = A(:,i)
                  do while (j .gt. istart)
                     A(:,i) = A(:,j)
                     iperm(j) = abs(iperm(j))
                     i = j
                     j = iperm(j)
                   enddo
                   A(:,i) = Ai
           endif
         enddo
         return
         end subroutine
        subroutine perm_inplace_d3(n, iperm, A)
        implicit none
        integer, intent(in) :: n
        integer, dimension(n), intent(inout) :: iperm
        real(kind=8), dimension(:,:,:), intent(inout) :: A
        integer :: i, j, k, istart
        logical :: has_touched_k
        real(kind=8), dimension(size(A,1),size(A,2)) :: Ai
! ----------------------------------
! use the sign bit in iperm(:) array
! to help find leaders of cycles
! ----------------------------------
        iperm = -iperm
        do k=1,n
          has_touched_k = (iperm(k) .gt. 0)
          if (.not. has_touched_k) then
                  iperm(k) = abs(iperm(k))
                  istart = k
                  i = istart
                  j = iperm(i)
                  Ai = A(:,:,i)
                  do while (j .gt. istart)
                     A(:,:,i) = A(:,:,j)
                     iperm(j) = abs(iperm(j))
                     i = j
                     j = iperm(j)
                   enddo
                   A(:,:,i) = Ai
           endif
         enddo
         return
         end subroutine
        subroutine perm_inplace_i1(n, iperm, A)
        implicit none
        integer, intent(in) :: n
        integer, dimension(n), intent(inout) :: iperm
        integer, parameter :: i4 = selected_int_kind(8)
        integer(kind=i4), dimension(n), intent(inout) :: A
        integer :: i, j, k, istart
        logical :: has_touched_k
        integer(kind=i4) :: Ai
! ----------------------------------
! use the sign bit in iperm(:) array
! to help find leaders of cycles
! ----------------------------------
        iperm = -iperm
        do k=1,n
          has_touched_k = (iperm(k) .gt. 0)
          if (.not. has_touched_k) then
                  iperm(k) = abs(iperm(k))
                  istart = k
                  i = istart
                  j = iperm(i)
                  Ai = A(i)
                  do while (j .gt. istart)
                     A(i) = A(j)
                     iperm(j) = abs(iperm(j))
                     i = j
                     j = iperm(j)
                   enddo
                   A(i) = Ai
           endif
         enddo
         return
         end subroutine
        subroutine perm_inplace_i2(n, iperm, A)
        implicit none
        integer, intent(in) :: n
        integer, dimension(n), intent(inout) :: iperm
        integer, parameter :: i4 = selected_int_kind(8)
        integer(kind=i4), dimension(:,:), intent(inout) :: A
        integer :: i, j, k, istart
        logical :: has_touched_k
        integer(kind=i4), dimension(size(A,1)) :: Ai
! ----------------------------------
! use the sign bit in iperm(:) array
! to help find leaders of cycles
! ----------------------------------
        iperm = -iperm
        do k=1,n
          has_touched_k = (iperm(k) .gt. 0)
          if (.not. has_touched_k) then
                  iperm(k) = abs(iperm(k))
                  istart = k
                  i = istart
                  j = iperm(i)
                  Ai(:) = A(:,i)
                  do while (j .gt. istart)
                     A(:,i) = A(:,j)
                     iperm(j) = abs(iperm(j))
                     i = j
                     j = iperm(j)
                   enddo
                   A(:,i) = Ai(:)
           endif
         enddo
         return
         end subroutine
        subroutine perm_inplace_i3(n, iperm, A)
        implicit none
        integer, intent(in) :: n
        integer, dimension(n), intent(inout) :: iperm
        integer, parameter :: i4 = selected_int_kind(8)
        integer(kind=i4), dimension(:,:,:), intent(inout) :: A
        integer :: i, j, k, istart
        logical :: has_touched_k
        integer(kind=i4), dimension(size(A,1),size(A,2)) :: Ai
! ----------------------------------
! use the sign bit in iperm(:) array
! to help find leaders of cycles
! ----------------------------------
        iperm = -iperm
        do k=1,n
          has_touched_k = (iperm(k) .gt. 0)
          if (.not. has_touched_k) then
                  iperm(k) = abs(iperm(k))
                  istart = k
                  i = istart
                  j = iperm(i)
                  Ai(:,:) = A(:,:,i)
                  do while (j .gt. istart)
                     A(:,:,i) = A(:,:,j)
                     iperm(j) = abs(iperm(j))
                     i = j
                     j = iperm(j)
                   enddo
                   A(:,:,i) = Ai(:,:)
           endif
         enddo
         return
         end subroutine
        subroutine perm_inplace_j1(n, iperm, A)
        implicit none
        integer, intent(in) :: n
        integer, dimension(n), intent(inout) :: iperm
        integer, parameter :: i8 = selected_int_kind(12)
        integer(kind=i8), dimension(n), intent(inout) :: A
        integer :: i, j, k, istart
        logical :: has_touched_k
        integer(kind=i8) :: Ai
! ----------------------------------
! use the sign bit in iperm(:) array
! to help find leaders of cycles
! ----------------------------------
        iperm = -iperm
        do k=1,n
          has_touched_k = (iperm(k) .gt. 0)
          if (.not. has_touched_k) then
                  iperm(k) = abs(iperm(k))
                  istart = k
                  i = istart
                  j = iperm(i)
                  Ai = A(i)
                  do while (j .gt. istart)
                     A(i) = A(j)
                     iperm(j) = abs(iperm(j))
                     i = j
                     j = iperm(j)
                   enddo
                   A(i) = Ai
           endif
         enddo
         return
         end subroutine
        subroutine perm_inplace_j2(n, iperm, A)
        implicit none
        integer, intent(in) :: n
        integer, dimension(n), intent(inout) :: iperm
        integer, parameter :: i8 = selected_int_kind(12)
        integer(kind=i8), dimension(:,:), intent(inout) :: A
        integer :: i, j, k, istart
        logical :: has_touched_k
        integer(kind=i8), dimension(size(A,1)) :: Ai
! ----------------------------------
! use the sign bit in iperm(:) array
! to help find leaders of cycles
! ----------------------------------
        iperm = -iperm
        do k=1,n
          has_touched_k = (iperm(k) .gt. 0)
          if (.not. has_touched_k) then
                  iperm(k) = abs(iperm(k))
                  istart = k
                  i = istart
                  j = iperm(i)
                  Ai(:) = A(:,i)
                  do while (j .gt. istart)
                     A(:,i) = A(:,j)
                     iperm(j) = abs(iperm(j))
                     i = j
                     j = iperm(j)
                   enddo
                   A(:,i) = Ai(:)
           endif
         enddo
         return
         end subroutine
        subroutine perm_inplace_j3(n, iperm, A)
        implicit none
        integer, intent(in) :: n
        integer, dimension(n), intent(inout) :: iperm
        integer, parameter :: i8 = selected_int_kind(12)
        integer(kind=i8), dimension(:,:,:), intent(inout) :: A
        integer :: i, j, k, istart
        logical :: has_touched_k
        integer(kind=i8), dimension(size(A,1),size(A,2)) :: Ai
! ----------------------------------
! use the sign bit in iperm(:) array
! to help find leaders of cycles
! ----------------------------------
        iperm = -iperm
        do k=1,n
          has_touched_k = (iperm(k) .gt. 0)
          if (.not. has_touched_k) then
                  iperm(k) = abs(iperm(k))
                  istart = k
                  i = istart
                  j = iperm(i)
                  Ai(:,:) = A(:,:,i)
                  do while (j .gt. istart)
                     A(:,:,i) = A(:,:,j)
                     iperm(j) = abs(iperm(j))
                     i = j
                     j = iperm(j)
                   enddo
                   A(:,:,i) = Ai(:,:)
           endif
         enddo
         return
         end subroutine
        subroutine permcopy_d1(n, iperm, ncycle, leader, A )
        implicit none
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: iperm
        integer, intent(in) :: ncycle
        integer, dimension(ncycle), intent(in) :: leader
        real(kind=8), dimension(:), intent(inout) :: A
        logical :: isok
        real(kind=8), dimension(:), allocatable :: Atmp
        integer :: i, ierr
        allocate( Atmp( size(A) ), stat=ierr)
        isok = (ierr.eq.0)
        if (isok) then
           do i=1,n
             Atmp(i) = A(i)
           enddo
           do i=1,n
             A(i) = Atmp( iperm(i) )
           enddo
           deallocate( Atmp, stat=ierr)
        else
! -------------------------------------------
! insufficient memory, use in place algorithm
! -------------------------------------------
           call permute_inplace(n,iperm,ncycle,leader,A)
        endif
        return
        end subroutine
        subroutine permcopy_d2(n, iperm, ncycle, leader, A )
        implicit none
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: iperm
        integer, intent(in) :: ncycle
        integer, dimension(ncycle), intent(in) :: leader
        real(kind=8), dimension(:,:), intent(inout) :: A
        logical :: isok
        real(kind=8), dimension(:,:), allocatable :: Atmp
        integer :: i, ierr
        allocate( Atmp( size(A,1), size(A,2) ), stat=ierr)
        isok = (ierr.eq.0)
        if (isok) then
           Atmp(:,1:n) = A(:,1:n)
           do i=1,n
             A(:,i) = Atmp( :, iperm(i) )
           enddo
           deallocate( Atmp, stat=ierr)
        else
! -------------------------------------------
! insufficient memory, use in place algorithm
! -------------------------------------------
           call permute_inplace(n,iperm,ncycle,leader,A)
        endif
        return
        end subroutine
        subroutine permcopy_d3(n, iperm, ncycle, leader, A )
        implicit none
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: iperm
        integer, intent(in) :: ncycle
        integer, dimension(ncycle), intent(in) :: leader
        real(kind=8), dimension(:,:,:), intent(inout) :: A
        logical :: isok
        real(kind=8), dimension(:,:,:), allocatable :: Atmp
        integer :: i, ierr
        allocate( Atmp( size(A,1), size(A,2), n ), stat=ierr)
        isok = (ierr.eq.0)
        if (isok) then
           Atmp(:,:,1:n) = A(:,:,1:n)
           do i=1,n
             A(:,:,i) = Atmp(:,:,iperm(i) )
           enddo
           deallocate( Atmp, stat=ierr)
        else
! -------------------------------------------
! insufficient memory, use in place algorithm
! -------------------------------------------
           call permute_inplace(n,iperm,ncycle,leader,A)
        endif
        return
        end subroutine
        subroutine permcopy_i1(n, iperm, ncycle, leader, A )
        implicit none
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: iperm
        integer, intent(in) :: ncycle
        integer, dimension(ncycle), intent(in) :: leader
        integer, parameter :: i4 = selected_int_kind(8)
        integer(kind=i4), dimension(:), intent(inout) :: A
        logical :: isok
        integer(kind=i4), dimension(:), allocatable :: Atmp
        integer :: i, ierr
        allocate( Atmp( size(A) ), stat=ierr)
        isok = (ierr.eq.0)
        if (isok) then
           do i=1,n
             Atmp(i) = A(i)
           enddo
           do i=1,n
             A(i) = Atmp( iperm(i) )
           enddo
           deallocate( Atmp, stat=ierr)
        else
! -------------------------------------------
! insufficient memory, use in place algorithm
! -------------------------------------------
           call permute_inplace(n,iperm,ncycle,leader,A)
        endif
        return
        end subroutine
        subroutine permcopy_i2(n, iperm, ncycle, leader, A )
        implicit none
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: iperm
        integer, intent(in) :: ncycle
        integer, dimension(ncycle), intent(in) :: leader
        integer, parameter :: i4 = selected_int_kind(8)
        integer(kind=i4), dimension(:,:), intent(inout) :: A
        logical :: isok
        integer(kind=i4), dimension(:,:), allocatable :: Atmp
        integer :: i, ierr
        allocate( Atmp(size(A,1), 1:n ), stat=ierr)
        isok = (ierr.eq.0)
        if (isok) then
           Atmp(:,1:n) = A(:,1:n)
           do i=1,n
             A(:,i) = Atmp(:,iperm(i) )
           enddo
           deallocate( Atmp, stat=ierr)
        else
! -------------------------------------------
! insufficient memory, use in place algorithm
! -------------------------------------------
           call permute_inplace(n,iperm,ncycle,leader,A)
        endif
        return
        end subroutine
        subroutine permcopy_i3(n, iperm, ncycle, leader, A )
        implicit none
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: iperm
        integer, intent(in) :: ncycle
        integer, dimension(ncycle), intent(in) :: leader
        integer, parameter :: i4 = selected_int_kind(8)
        integer(kind=i4), dimension(:,:,:), intent(inout) :: A
        logical :: isok
        integer(kind=i4), dimension(:,:,:), allocatable :: Atmp
        integer :: i, ierr
        allocate( Atmp(size(A,1),size(A,2), 1:n ), stat=ierr)
        isok = (ierr.eq.0)
        if (isok) then
           Atmp(:,:,1:n) = A(:,:,1:n)
           do i=1,n
             A(:,:,i) = Atmp( :,:,iperm(i) )
           enddo
           deallocate( Atmp, stat=ierr)
        else
! -------------------------------------------
! insufficient memory, use in place algorithm
! -------------------------------------------
           call permute_inplace(n,iperm,ncycle,leader,A)
        endif
        return
        end subroutine
        subroutine permcopy_j1(n, iperm, ncycle, leader, A )
        implicit none
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: iperm
        integer, intent(in) :: ncycle
        integer, dimension(ncycle), intent(in) :: leader
        integer, parameter :: i8 = selected_int_kind(12)
        integer(kind=i8), dimension(:), intent(inout) :: A
        logical :: isok
        integer(kind=i8), dimension(:), allocatable :: Atmp
        integer :: i, ierr
        allocate( Atmp( size(A) ), stat=ierr)
        isok = (ierr.eq.0)
        if (isok) then
           do i=1,n
             Atmp(i) = A(i)
           enddo
           do i=1,n
             A(i) = Atmp( iperm(i) )
           enddo
           deallocate( Atmp, stat=ierr)
        else
! -------------------------------------------
! insufficient memory, use in place algorithm
! -------------------------------------------
           call permute_inplace(n,iperm,ncycle,leader,A)
        endif
        return
        end subroutine
        subroutine permcopy_j2(n, iperm, ncycle, leader, A )
        implicit none
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: iperm
        integer, intent(in) :: ncycle
        integer, dimension(ncycle), intent(in) :: leader
        integer, parameter :: i8 = selected_int_kind(12)
        integer(kind=i8), dimension(:,:), intent(inout) :: A
        logical :: isok
        integer(kind=i8), dimension(:,:), allocatable :: Atmp
        integer :: i, ierr
        allocate( Atmp(size(A,1), 1:n ), stat=ierr)
        isok = (ierr.eq.0)
        if (isok) then
           Atmp(:,1:n) = A(:,1:n)
           do i=1,n
             A(:,i) = Atmp(:,iperm(i) )
           enddo
           deallocate( Atmp, stat=ierr)
        else
! -------------------------------------------
! insufficient memory, use in place algorithm
! -------------------------------------------
           call permute_inplace(n,iperm,ncycle,leader,A)
        endif
        return
        end subroutine
        subroutine permcopy_j3(n, iperm, ncycle, leader, A )
        implicit none
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: iperm
        integer, intent(in) :: ncycle
        integer, dimension(ncycle), intent(in) :: leader
        integer, parameter :: i8 = selected_int_kind(12)
        integer(kind=i8), dimension(:,:,:), intent(inout) :: A
        logical :: isok
        integer(kind=i8), dimension(:,:,:), allocatable :: Atmp
        integer :: i, ierr
        allocate( Atmp(size(A,1),size(A,2), 1:n ), stat=ierr)
        isok = (ierr.eq.0)
        if (isok) then
           Atmp(:,:,1:n) = A(:,:,1:n)
           do i=1,n
             A(:,:,i) = Atmp( :,:,iperm(i) )
           enddo
           deallocate( Atmp, stat=ierr)
        else
! -------------------------------------------
! insufficient memory, use in place algorithm
! -------------------------------------------
           call permute_inplace(n,iperm,ncycle,leader,A)
        endif
        return
        end subroutine
        subroutine permscopy_d1(n,iperm,A)
        implicit none
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: iperm
        real(kind=8), dimension(:), intent(inout) :: A
        real(kind=8), dimension(n) :: Atmp
        Atmp(1:n) = A(1:n)
        A(1:n) = Atmp( iperm(1:n) )
        return
        end subroutine
        subroutine permscopy_d2(n,iperm,A)
        implicit none
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: iperm
        real(kind=8), dimension(:,:), intent(inout) :: A
        real(kind=8), dimension(size(A,1),n) :: Atmp
        Atmp(:,1:n) = A(:,1:n)
        A(:,1:n) = Atmp(:, iperm(1:n) )
        return
        end subroutine
        subroutine permscopy_d3(n,iperm,A)
        implicit none
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: iperm
        real(kind=8), dimension(:,:,:), intent(inout) :: A
        real(kind=8), dimension(size(A,1),size(A,2),n) :: Atmp
        Atmp(:,:,1:n) = A(:,:,1:n)
        A(:,:,1:n) = Atmp(:,:, iperm(1:n) )
        return
        end subroutine
        subroutine permscopy_i1(n,iperm,A)
        implicit none
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: iperm
        integer, parameter :: i4 = selected_int_kind(8)
        integer(kind=i4), dimension(:), intent(inout) :: A
        integer(kind=i4), dimension(n) :: Atmp
        Atmp(1:n) = A(1:n)
        A(1:n) = Atmp( iperm(1:n) )
        return
        end subroutine
        subroutine permscopy_i2(n,iperm,A)
        implicit none
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: iperm
        integer, parameter :: i4 = selected_int_kind(8)
        integer(kind=i4), dimension(:,:), intent(inout) :: A
        integer(kind=i4), dimension(size(A,1),n) :: Atmp
        Atmp(:,1:n) = A(:,1:n)
        A(:,1:n) = Atmp(:, iperm(1:n) )
        return
        end subroutine
        subroutine permscopy_i3(n,iperm,A)
        implicit none
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: iperm
        integer, parameter :: i4 = selected_int_kind(8)
        integer(kind=i4), dimension(:,:,:), intent(inout) :: A
        integer(kind=i4), dimension(size(A,1),size(A,2),n) :: Atmp
        Atmp(:,:,1:n) = A(:,:,1:n)
        A(:,:,1:n) = Atmp(:,:, iperm(1:n) )
        return
        end subroutine
        subroutine permscopy_j1(n,iperm,A)
        implicit none
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: iperm
        integer, parameter :: i8 = selected_int_kind(12)
        integer(kind=i8), dimension(:), intent(inout) :: A
        integer(kind=i8), dimension(n) :: Atmp
        Atmp(1:n) = A(1:n)
        A(1:n) = Atmp( iperm(1:n) )
        return
        end subroutine
        subroutine permscopy_j2(n,iperm,A)
        implicit none
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: iperm
        integer, parameter :: i8 = selected_int_kind(12)
        integer(kind=i8), dimension(:,:), intent(inout) :: A
        integer(kind=i8), dimension(size(A,1),n) :: Atmp
        Atmp(:,1:n) = A(:,1:n)
        A(:,1:n) = Atmp(:, iperm(1:n) )
        return
        end subroutine
        subroutine permscopy_j3(n,iperm,A)
        implicit none
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: iperm
        integer, parameter :: i8 = selected_int_kind(12)
        integer(kind=i8), dimension(:,:,:), intent(inout) :: A
        integer(kind=i8), dimension(size(A,1),size(A,2),n) :: Atmp
        Atmp(:,:,1:n) = A(:,:,1:n)
        A(:,:,1:n) = Atmp(:,:, iperm(1:n) )
        return
        end subroutine
        subroutine test_perm_module(n)
        implicit none
        integer, intent(in) :: n
        integer, parameter :: i4 = selected_int_kind(8)
        integer, parameter :: i8 = selected_int_kind(12)
        real(kind=8), dimension(n) :: A, Asort, Atmp
        integer(kind=i4), dimension(n) :: iA,iAsort,iAtmp
        integer(kind=i8), dimension(n) :: jA,jAsort,jAtmp
        integer, parameter :: k1 = 4, k2 = 4
        real(kind=8), dimension(k1,k2,n) :: A3, A3sort, A3tmp
        real(kind=8), dimension(k1,n) :: A2, A2sort, A2tmp
        integer(kind=i4), dimension(k1,k2,n) :: iA3, iA3sort, iA3tmp
        integer(kind=i4), dimension(k1,n) :: iA2, iA2sort, iA2tmp
        integer(kind=i8), dimension(k1,k2,n) :: jA3, jA3sort, jA3tmp
        integer(kind=i8), dimension(k1,n) :: jA2, jA2sort, jA2tmp
        integer, dimension(n) :: iperm
        integer :: ldim
        integer :: ncycle
        integer, dimension(n) :: leader
        logical :: isok
        real(kind=8) :: time_A, time_A2, time_A3
        real(kind=8) :: time_iA, time_iA2, time_iA3
        real(kind=8) :: time_jA, time_jA2, time_jA3
        real(kind=8) :: t1, t2
        integer :: it, niter, i,j
        real(kind=8) :: ij(2)
        integer :: itemp
        integer, parameter :: ncase = 4
        integer :: icase
        write(*,*) 'n = ',n
        call random_number( A )
        Atmp = A
        call random_number( A2 )
        A2tmp = A2
        call random_number( A3 )
        A3tmp = A3
        iA = int(A * (2.0d0**30))
        iA2 = int(A2 * (2.0d0**30))
        iA3 = int(A3 * (2.0d0**30))
        iAtmp = iA
        iA2tmp = iA2
        iA3tmp = iA3
        jA = int(A * (2.0d0**30))
        jA2 = int(A2 * (2.0d0**30))
        jA3 = int(A3 * (2.0d0**30))
        jAtmp = jA
        jA2tmp = jA2
        jA3tmp = jA3
        do icase=1,ncase
        write(*,*) 'icase  == ', icase
        if (icase.eq.1) then
          niter = 4*n
          write(*,*) 'random exchange niter = ',niter
          do i=1,n
            iperm(i) = i
          enddo
          do it=1,niter
            call random_number(ij)
            i = abs(int( ij(1)*(2.0d0**20) ))
            j = abs(int( ij(2)*(2.0d0**20) ))
            i = mod( i, n ) + 1
            j = mod( j, n ) + 1
            itemp = iperm(i)
            iperm(i) = iperm(j)
            iperm(j) = itemp
          enddo
        else if (icase.eq.2) then
          write(*,*) '2 cycles '
          do i=1,n,2
             itemp = iperm(i)
             iperm(i) = iperm(i+1)
             iperm(i+1) = itemp
          enddo
        else if (icase.eq.3) then
            write(*,*) 'reversal '
            do i=1,n
              iperm(i) = n-i+1
            enddo
        else if (icase .eq. 4) then
           write(*,*) 'identity permutation'
           do i=1,n
             iperm(i) = i
           enddo
        else
           call assert(.false.,'unknown case ',icase )
        endif
! --------------------------------
! straight forward way for A,A2,A3
! --------------------------------
        call cpu_time(t1)
        Atmp(1:n) = A(1:n)
        Asort(1:n) = Atmp( iperm(1:n) )
        call cpu_time(t2)
        time_A = t2-t1
        write(*,*) ' Atmp(iperm) took ', time_A
        call cpu_time(t1)
        A2tmp = A2
        A2sort = A2tmp(:,iperm(1:n) )
        call cpu_time(t2)
        time_A2 = t2-t1
        write(*,*) ' A2tmp(iperm) took ', time_A2
        call cpu_time(t1)
        A3tmp = A3
        A3sort = A3tmp(:,:,iperm(1:n) )
        call cpu_time(t2)
        time_A3 = t2-t1
        write(*,*) ' A3tmp(iperm) took ', time_A3
! --------------------------------
! straight forward way for iA,iA2,iA3
! --------------------------------
        call cpu_time(t1)
        iAtmp(1:n) = iA(1:n)
        iAsort(1:n) = iAtmp( iperm(1:n) )
        call cpu_time(t2)
        time_iA = t2-t1
        write(*,*) ' iAtmp(iperm) took ', time_iA
        call cpu_time(t1)
        iA2tmp = iA2
        iA2sort = iA2tmp(:,iperm(1:n) )
        call cpu_time(t2)
        time_iA2 = t2-t1
        write(*,*) ' iA2tmp(iperm) took ', time_iA2
        call cpu_time(t1)
        iA3tmp = iA3
        iA3sort = iA3tmp(:,:,iperm(1:n) )
        call cpu_time(t2)
        time_iA3 = t2-t1
        write(*,*) ' iA3tmp(iperm) took ', time_iA3
! --------------------------------
! straight forward way for jA,jA2,jA3
! --------------------------------
        call cpu_time(t1)
        jAtmp(1:n) = jA(1:n)
        jAsort(1:n) = jAtmp( iperm(1:n) )
        call cpu_time(t2)
        time_jA = t2-t1
        write(*,*) ' jAtmp(iperm) took ', time_jA
        call cpu_time(t1)
        jA2tmp = jA2
        jA2sort = jA2tmp(:,iperm(1:n) )
        call cpu_time(t2)
        time_jA2 = t2-t1
        write(*,*) ' jA2tmp(iperm) took ', time_jA2
        call cpu_time(t1)
        jA3tmp = jA3
        jA3sort = jA3tmp(:,:,iperm(1:n) )
        call cpu_time(t2)
        time_jA3 = t2-t1
        write(*,*) ' jA3tmp(iperm) took ', time_jA3
        call cpu_time(t1)
        ldim = n
        call gencycle( n, iperm, ncycle, leader, ldim )
        call cpu_time(t2)
        write(*,*) 'gencycle took ', t2-t1
        write(*,*) 'ncycle ', ncycle
! ----------------------
! test permcopy A,A2,A3
! ----------------------
        A = Atmp
        call cpu_time(t1)
        call permcopy(n,iperm,ncycle,leader,A)
        call cpu_time(t2)
        write(*,*) 'permcopy A took ', t2-t1
        isok = all( A .eq. Asort )
        call assert(isok,' permcopy A failed ',n)
        A2 = A2tmp
        call cpu_time(t1)
        call permcopy(n,iperm,ncycle,leader,A2)
        call cpu_time(t2)
        write(*,*) 'permcopy A2 took ', t2-t1
        isok = all( A2 .eq. A2sort )
        call assert(isok,' permcopy A2 failed ',n)
        A3 = A3tmp
        call cpu_time(t1)
        call permcopy(n,iperm,ncycle,leader,A3)
        call cpu_time(t2)
        write(*,*) 'permcopy A3 took ', t2-t1
        isok = all( A3 .eq. A3sort )
        call assert(isok, 'permcopy A3 failed ',n)
! --------------------------
! test permcopy iA,iA2,iA3
! --------------------------
        iA = iAtmp
        call cpu_time(t1)
        call permcopy(n,iperm,ncycle,leader,iA)
        call cpu_time(t2)
        write(*,*) 'permcopy iA took ', t2-t1
        isok = all( iA .eq. iAsort )
        call assert(isok,' permcopy iA failed ',n)
        iA2 = iA2tmp
        call cpu_time(t1)
        call permcopy(n,iperm,ncycle,leader,iA2)
        call cpu_time(t2)
        write(*,*) 'permcopy iA2 took ', t2-t1
        isok = all( iA2 .eq. iA2sort )
        call assert(isok,' permcopy iA2 failed ',n)
        iA3 = iA3tmp
        call cpu_time(t1)
        call permcopy(n,iperm,ncycle,leader,iA3)
        call cpu_time(t2)
        write(*,*) 'permcopy iA3 took ', t2-t1
        isok = all( iA3 .eq. iA3sort )
        call assert(isok, 'permcopy iA3 failed ',n)
! --------------------------
! test permcopy jA,jA2,jA3
! --------------------------
        jA = jAtmp
        call cpu_time(t1)
        call permcopy(n,iperm,ncycle,leader,jA)
        call cpu_time(t2)
        write(*,*) 'permcopy jA took ', t2-t1
        isok = all( jA .eq. jAsort )
        call assert(isok,' permcopy jA failed ',n)
        jA2 = jA2tmp
        call cpu_time(t1)
        call permcopy(n,iperm,ncycle,leader,jA2)
        call cpu_time(t2)
        write(*,*) 'permcopy jA2 took ', t2-t1
        isok = all( jA2 .eq. jA2sort )
        call assert(isok,' permcopy jA2 failed ',n)
        jA3 = jA3tmp
        call cpu_time(t1)
        call permcopy(n,iperm,ncycle,leader,jA3)
        call cpu_time(t2)
        write(*,*) 'permcopy jA3 took ', t2-t1
        isok = all( jA3 .eq. jA3sort )
        call assert(isok, 'permcopy jA3 failed ',n)
! =====================================================
! ----------------------
! test permscopy A,A2,A3
! ----------------------
        A = Atmp
        call cpu_time(t1)
        call permscopy(n,iperm,A)
        call cpu_time(t2)
        write(*,*) 'permscopy A took ', t2-t1
        isok = all( A .eq. Asort )
        call assert(isok,' permscopy A failed ',n)
        A2 = A2tmp
        call cpu_time(t1)
        call permscopy(n,iperm,A2)
        call cpu_time(t2)
        write(*,*) 'permscopy A2 took ', t2-t1
        isok = all( A2 .eq. A2sort )
        call assert(isok,' permscopy A2 failed ',n)
        A3 = A3tmp
        call cpu_time(t1)
        call permscopy(n,iperm,A3)
        call cpu_time(t2)
        write(*,*) 'permscopy A3 took ', t2-t1
        isok = all( A3 .eq. A3sort )
        call assert(isok, 'permscopy A3 failed ',n)
! --------------------------
! test permscopy iA,iA2,iA3
! --------------------------
        iA = iAtmp
        call cpu_time(t1)
        call permscopy(n,iperm,iA)
        call cpu_time(t2)
        write(*,*) 'permscopy iA took ', t2-t1
        isok = all( iA .eq. iAsort )
        call assert(isok,' permscopy iA failed ',n)
        iA2 = iA2tmp
        call cpu_time(t1)
        call permscopy(n,iperm,iA2)
        call cpu_time(t2)
        write(*,*) 'permscopy iA2 took ', t2-t1
        isok = all( iA2 .eq. iA2sort )
        call assert(isok,' permscopy iA2 failed ',n)
        iA3 = iA3tmp
        call cpu_time(t1)
        call permscopy(n,iperm,iA3)
        call cpu_time(t2)
        write(*,*) 'permscopy iA3 took ', t2-t1
        isok = all( iA3 .eq. iA3sort )
        call assert(isok, 'permscopy iA3 failed ',n)
! --------------------------
! test permscopy jA,jA2,jA3
! --------------------------
        jA = jAtmp
        call cpu_time(t1)
        call permscopy(n,iperm,jA)
        call cpu_time(t2)
        write(*,*) 'permscopy jA took ', t2-t1
        isok = all( jA .eq. jAsort )
        call assert(isok,' permscopy jA failed ',n)
        jA2 = jA2tmp
        call cpu_time(t1)
        call permscopy(n,iperm,jA2)
        call cpu_time(t2)
        write(*,*) 'permscopy jA2 took ', t2-t1
        isok = all( jA2 .eq. jA2sort )
        call assert(isok,' permscopy jA2 failed ',n)
        jA3 = jA3tmp
        call cpu_time(t1)
        call permscopy(n,iperm,jA3)
        call cpu_time(t2)
        write(*,*) 'permscopy jA3 took ', t2-t1
        isok = all( jA3 .eq. jA3sort )
        call assert(isok, 'permscopy jA3 failed ',n)
! =====================================================
! ----------------------
! test permute_inplace A,A2,A3
! ----------------------
        A = Atmp
        call cpu_time(t1)
        call permute_inplace(n,iperm,ncycle,leader,A)
        call cpu_time(t2)
        write(*,*) 'permute_inplace(ncyle) A took ', t2-t1
        isok = all( A .eq. Asort )
        call assert(isok,' permute_inplace A failed ',n)
        A2 = A2tmp
        call cpu_time(t1)
        call permute_inplace(n,iperm,ncycle,leader,A2)
        call cpu_time(t2)
        write(*,*) 'permute_inplace(ncyle) A2 took ', t2-t1
        isok = all( A2 .eq. A2sort )
        call assert(isok,' permute_inplace A2 failed ',n)
        A3 = A3tmp
        call cpu_time(t1)
        call permute_inplace(n,iperm,ncycle,leader,A3)
        call cpu_time(t2)
        write(*,*) 'permute_inplace(ncyle) A3 took ', t2-t1
        isok = all( A3 .eq. A3sort )
        call assert(isok, 'permute_inplace A3 failed ',n)
! --------------------------
! test permute_inplace iA,iA2,iA3
! --------------------------
        iA = iAtmp
        call cpu_time(t1)
        call permute_inplace(n,iperm,ncycle,leader,iA)
        call cpu_time(t2)
        write(*,*) 'permute_inplace(ncyle) iA took ', t2-t1
        isok = all( iA .eq. iAsort )
        call assert(isok,' permute_inplace iA failed ',n)
        iA2 = iA2tmp
        call cpu_time(t1)
        call permute_inplace(n,iperm,ncycle,leader,iA2)
        call cpu_time(t2)
        write(*,*) 'permute_inplace(ncyle) iA2 took ', t2-t1
        isok = all( iA2 .eq. iA2sort )
        call assert(isok,' permute_inplace iA2 failed ',n)
        iA3 = iA3tmp
        call cpu_time(t1)
        call permute_inplace(n,iperm,ncycle,leader,iA3)
        call cpu_time(t2)
        write(*,*) 'permute_inplace(ncyle) iA3 took ', t2-t1
        isok = all( iA3 .eq. iA3sort )
        call assert(isok, 'permute_inplace iA3 failed ',n)
! --------------------------
! test permute_inplace jA,jA2,jA3
! --------------------------
        jA = jAtmp
        call cpu_time(t1)
        call permute_inplace(n,iperm,ncycle,leader,jA)
        call cpu_time(t2)
        write(*,*) 'permute_inplace(ncyle) jA took ', t2-t1
        isok = all( jA .eq. jAsort )
        call assert(isok,' permute_inplace jA failed ',n)
        jA2 = jA2tmp
        call cpu_time(t1)
        call permute_inplace(n,iperm,ncycle,leader,jA2)
        call cpu_time(t2)
        write(*,*) 'permute_inplace(ncyle) jA2 took ', t2-t1
        isok = all( jA2 .eq. jA2sort )
        call assert(isok,' permute_inplace jA2 failed ',n)
        jA3 = jA3tmp
        call cpu_time(t1)
        call permute_inplace(n,iperm,ncycle,leader,jA3)
        call cpu_time(t2)
        write(*,*) 'permute_inplace(ncyle) jA3 took ', t2-t1
        isok = all( jA3 .eq. jA3sort )
        call assert(isok, 'permute_inplace jA3 failed ',n)
!======================================================
! ----------------------
! test permute_inplace A,A2,A3
! ----------------------
        A = Atmp
        call cpu_time(t1)
        call permute_inplace(n,iperm,A)
        call cpu_time(t2)
        write(*,*) 'permute_inplace A took ', t2-t1
        isok = all( A .eq. Asort )
        call assert(isok,' permute_inplace A failed ',n)
        A2 = A2tmp
        call cpu_time(t1)
        call permute_inplace(n,iperm,A2)
        call cpu_time(t2)
        write(*,*) 'permute_inplace A2 took ', t2-t1
        isok = all( A2 .eq. A2sort )
        call assert(isok,' permute_inplace A2 failed ',n)
        A3 = A3tmp
        call cpu_time(t1)
        call permute_inplace(n,iperm,A3)
        call cpu_time(t2)
        write(*,*) 'permute_inplace A3 took ', t2-t1
        isok = all( A3 .eq. A3sort )
        call assert(isok, 'permute_inplace A3 failed ',n)
! --------------------------
! test permute_inplace iA,iA2,iA3
! --------------------------
        iA = iAtmp
        call cpu_time(t1)
        call permute_inplace(n,iperm,iA)
        call cpu_time(t2)
        write(*,*) 'permute_inplace iA took ', t2-t1
        isok = all( iA .eq. iAsort )
        call assert(isok,' permute_inplace iA failed ',n)
        iA2 = iA2tmp
        call cpu_time(t1)
        call permute_inplace(n,iperm,iA2)
        call cpu_time(t2)
        write(*,*) 'permute_inplace iA2 took ', t2-t1
        isok = all( iA2 .eq. iA2sort )
        call assert(isok,' permute_inplace iA2 failed ',n)
        iA3 = iA3tmp
        call cpu_time(t1)
        call permute_inplace(n,iperm,iA3)
        call cpu_time(t2)
        write(*,*) 'permute_inplace iA3 took ', t2-t1
        isok = all( iA3 .eq. iA3sort )
        call assert(isok, 'permute_inplace iA3 failed ',n)
! --------------------------
! test permute_inplace jA,jA2,jA3
! --------------------------
        jA = jAtmp
        call cpu_time(t1)
        call permute_inplace(n,iperm,jA)
        call cpu_time(t2)
        write(*,*) 'permute_inplace jA took ', t2-t1
        isok = all( jA .eq. jAsort )
        call assert(isok,' permute_inplace jA failed ',n)
        jA2 = jA2tmp
        call cpu_time(t1)
        call permute_inplace(n,iperm,jA2)
        call cpu_time(t2)
        write(*,*) 'permute_inplace jA2 took ', t2-t1
        isok = all( jA2 .eq. jA2sort )
        call assert(isok,' permute_inplace jA2 failed ',n)
        jA3 = jA3tmp
        call cpu_time(t1)
        call permute_inplace(n,iperm,jA3)
        call cpu_time(t2)
        write(*,*) 'permute_inplace jA3 took ', t2-t1
        isok = all( jA3 .eq. jA3sort )
        call assert(isok, 'permute_inplace jA3 failed ',n)
        enddo
        return
        end subroutine
        end module perm_module
