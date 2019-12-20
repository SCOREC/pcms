

      module comm_mod
      use assert_mod
      implicit none
      include 'mpif.h'

      private
      integer, parameter :: idebug = 0

    
      public :: allgatherv_2d,allgatherv_2i
      public :: allgatherv_2s,allgatherv_2l

      public :: allgatherv_2d_ring,allgatherv_2i_ring
      public :: allgatherv_2s_ring,allgatherv_2l_ring

      contains
      subroutine allgatherv_2d( array, n1,n2, segsize, comm )
!     ----------------------------
!     Perform Allgatherv on 2D array
!     layer to call MPI_Allgatherv
!     ----------------------------
      implicit none
      integer, intent(in) :: n1, n2, segsize
      integer, intent(in) :: comm

      real(8), intent(inout) :: array(n1*n2)
      integer, parameter :: mpi_type = MPI_REAL8
      real(8), allocatable, dimension(:) :: sendbuf, recvbuf

      integer, allocatable, dimension(:) :: displs, recvcounts
      integer :: mype, npes
      integer :: sendcount, sendtype, recvtype
      integer :: ierr, istat, sendbuf_size, recvbuf_size 
      integer :: istart,iend,isize, pe, j,i1,i2, j1,j2
   

      ierr = MPI_SUCCESS
      call MPI_Comm_rank( comm, mype, ierr )
      if (ierr.ne.MPI_SUCCESS) then
         write(*,*) 'allgatherv_2d: mpi_comm_rank return ',ierr
         return
      endif

      call MPI_Comm_size( comm, npes, ierr )
      if (ierr.ne.MPI_SUCCESS) then
         write(*,*) 'allgatherv_2d: mpi_comm_size return ',ierr
         return
      endif

      sendbuf_size = max(1, n1*segsize)
 
      isize = min( segsize*npes, n2 )
      recvbuf_size = max(1, n1*isize)
      allocate( sendbuf(sendbuf_size),                                   &
     &          recvbuf(recvbuf_size),                                   &
     &          displs(0:(npes-1)),                                      &
     &          recvcounts(0:(npes-1)),                                  &
     &          stat=istat )
      if (istat.ne.0) then
        write(*,*) 'allgatherv_2d: alloc recvbuf, n1,n2 ',n1,n2
        stop '** allocation error ** '
      endif

      do pe=0,npes-1
         istart = 1 + pe*segsize
         iend = min( n2, istart + segsize-1)
         isize = max(0,iend - istart + 1)

         displs(pe) =  pe*(n1*segsize)
         recvcounts(pe) = (n1*isize)
      enddo

!     ----------------
!     copy into send buffer
!     ----------------
      pe = mype
      istart = 1 + pe*segsize
      iend = min(n2, istart + segsize-1)
      isize = max(0,iend-istart+1)

      sendcount = isize*n1
      
      do j=istart,iend
        i1 = (j-istart)*n1 + 1
        i2 = i1 + n1-1

        
        j1 = 1 + (j-1)*n1
        j2 = j1 + n1 - 1
        sendbuf( i1:i2 ) = array( j1:j2 )
      enddo


      if (idebug >= 1) then
        if (mype .eq. (npes-1)) then
          do pe=0,npes-1
           write(*,*) 'pe,displs,recvcount ',                           &
     &                 pe,displs(pe),recvcounts(pe)
          enddo
        endif
      endif

      sendtype = mpi_type
      recvtype = sendtype
      call MPI_Allgatherv( sendbuf, sendcount, sendtype,                 &
     &                     recvbuf, recvcounts, displs, recvtype,        &
     &                     comm, ierr)
      if (ierr.ne.MPI_SUCCESS) then
         write(*,*) 'allgatherv_2d: mpi_allgatherv return ',ierr
         return
      endif

!    ------------------------------
!    copy from recv buffer to array
!    ------------------------------

      do j=1,min(n2, npes*segsize)
        i1 = 1 + (j-1)*n1
        i2 = i1 + n1 - 1

        j1 = 1 + (j-1)*n1
        j2 = j1 + n1 - 1
        array( j1:j2 ) = recvbuf( i1:i2 )
      enddo

!     --------
!     clean up
!     --------
      deallocate( sendbuf,                                               &
     &            recvbuf,                                               &
     &            displs,                                                &
     &            recvcounts,                                            &
     &            stat=istat )

      return
      end subroutine allgatherv_2d
      subroutine allgatherv_2i( array, n1,n2, segsize, comm )
!     ----------------------------
!     Perform Allgatherv on 2D array
!     layer to call MPI_Allgatherv
!     ----------------------------
      implicit none
      integer, intent(in) :: n1, n2, segsize
      integer, intent(in) :: comm

      integer, intent(inout) :: array(n1*n2)
      integer, parameter :: mpi_type = MPI_INTEGER
      integer, allocatable, dimension(:) :: sendbuf, recvbuf

      integer, allocatable, dimension(:) :: displs, recvcounts
      integer :: mype, npes
      integer :: sendcount, sendtype, recvtype
      integer :: ierr, istat, sendbuf_size, recvbuf_size 
      integer :: istart,iend,isize, pe, j,i1,i2, j1,j2
   

      ierr = MPI_SUCCESS
      call MPI_Comm_rank( comm, mype, ierr )
      if (ierr.ne.MPI_SUCCESS) then
         write(*,*) 'allgatherv_2i: mpi_comm_rank return ',ierr
         return
      endif

      call MPI_Comm_size( comm, npes, ierr )
      if (ierr.ne.MPI_SUCCESS) then
         write(*,*) 'allgatherv_2i: mpi_comm_size return ',ierr
         return
      endif

      sendbuf_size = max(1, n1*segsize)
 
      isize = min( segsize*npes, n2 )
      recvbuf_size = max(1, n1*isize)
      allocate( sendbuf(sendbuf_size),                                   &
     &          recvbuf(recvbuf_size),                                   &
     &          displs(0:(npes-1)),                                      &
     &          recvcounts(0:(npes-1)),                                  &
     &          stat=istat )
      if (istat.ne.0) then
        write(*,*) 'allgatherv_2i: alloc recvbuf, n1,n2 ',n1,n2
        stop '** allocation error ** '
      endif

      do pe=0,npes-1
         istart = 1 + pe*segsize
         iend = min( n2, istart + segsize-1)
         isize = max(0,iend - istart + 1)

         displs(pe) =  pe*(n1*segsize)
         recvcounts(pe) = (n1*isize)
      enddo

!     ----------------
!     copy into send buffer
!     ----------------
      pe = mype
      istart = 1 + pe*segsize
      iend = min(n2, istart + segsize-1)
      isize = max(0,iend-istart+1)

      sendcount = isize*n1
      
      do j=istart,iend
        i1 = (j-istart)*n1 + 1
        i2 = i1 + n1-1

        
        j1 = 1 + (j-1)*n1
        j2 = j1 + n1 - 1
        sendbuf( i1:i2 ) = array( j1:j2 )
      enddo


      if (idebug >= 1) then
        if (mype .eq. (npes-1)) then
          do pe=0,npes-1
           write(*,*) 'pe,displs,recvcount ',                           &
     &                 pe,displs(pe),recvcounts(pe)
          enddo
        endif
      endif

      sendtype = mpi_type
      recvtype = sendtype
      call MPI_Allgatherv( sendbuf, sendcount, sendtype,                 &
     &                     recvbuf, recvcounts, displs, recvtype,        &
     &                     comm, ierr)
      if (ierr.ne.MPI_SUCCESS) then
         write(*,*) 'allgatherv_2i: mpi_allgatherv return ',ierr
         return
      endif

!    ------------------------------
!    copy from recv buffer to array
!    ------------------------------

      do j=1,min(n2, npes*segsize)
        i1 = 1 + (j-1)*n1
        i2 = i1 + n1 - 1

        j1 = 1 + (j-1)*n1
        j2 = j1 + n1 - 1
        array( j1:j2 ) = recvbuf( i1:i2 )
      enddo

!     --------
!     clean up
!     --------
      deallocate( sendbuf,                                               &
     &            recvbuf,                                               &
     &            displs,                                                &
     &            recvcounts,                                            &
     &            stat=istat )

      return
      end subroutine allgatherv_2i
      subroutine allgatherv_2l( array, n1,n2, segsize, comm )
!     ----------------------------
!     Perform Allgatherv on 2D array
!     layer to call MPI_Allgatherv
!     ----------------------------
      implicit none
      integer, intent(in) :: n1, n2, segsize
      integer, intent(in) :: comm

      logical, intent(inout) :: array(n1*n2)
      integer, parameter :: mpi_type = MPI_LOGICAL
      logical, allocatable, dimension(:) :: sendbuf, recvbuf

      integer, allocatable, dimension(:) :: displs, recvcounts
      integer :: mype, npes
      integer :: sendcount, sendtype, recvtype
      integer :: ierr, istat, sendbuf_size, recvbuf_size 
      integer :: istart,iend,isize, pe, j,i1,i2, j1,j2
   

      ierr = MPI_SUCCESS
      call MPI_Comm_rank( comm, mype, ierr )
      if (ierr.ne.MPI_SUCCESS) then
         write(*,*) 'allgatherv_2l: mpi_comm_rank return ',ierr
         return
      endif

      call MPI_Comm_size( comm, npes, ierr )
      if (ierr.ne.MPI_SUCCESS) then
         write(*,*) 'allgatherv_2l: mpi_comm_size return ',ierr
         return
      endif

      sendbuf_size = max(1, n1*segsize)
 
      isize = min( segsize*npes, n2 )
      recvbuf_size = max(1, n1*isize)
      allocate( sendbuf(sendbuf_size),                                   &
     &          recvbuf(recvbuf_size),                                   &
     &          displs(0:(npes-1)),                                      &
     &          recvcounts(0:(npes-1)),                                  &
     &          stat=istat )
      if (istat.ne.0) then
        write(*,*) 'allgatherv_2l: alloc recvbuf, n1,n2 ',n1,n2
        stop '** allocation error ** '
      endif

      do pe=0,npes-1
         istart = 1 + pe*segsize
         iend = min( n2, istart + segsize-1)
         isize = max(0,iend - istart + 1)

         displs(pe) =  pe*(n1*segsize)
         recvcounts(pe) = (n1*isize)
      enddo

!     ----------------
!     copy into send buffer
!     ----------------
      pe = mype
      istart = 1 + pe*segsize
      iend = min(n2, istart + segsize-1)
      isize = max(0,iend-istart+1)

      sendcount = isize*n1
      
      do j=istart,iend
        i1 = (j-istart)*n1 + 1
        i2 = i1 + n1-1

        
        j1 = 1 + (j-1)*n1
        j2 = j1 + n1 - 1
        sendbuf( i1:i2 ) = array( j1:j2 )
      enddo


      if (idebug >= 1) then
        if (mype .eq. (npes-1)) then
          do pe=0,npes-1
           write(*,*) 'pe,displs,recvcount ',                           &
     &                 pe,displs(pe),recvcounts(pe)
          enddo
        endif
      endif

      sendtype = mpi_type
      recvtype = sendtype
      call MPI_Allgatherv( sendbuf, sendcount, sendtype,                 &
     &                     recvbuf, recvcounts, displs, recvtype,        &
     &                     comm, ierr)
      if (ierr.ne.MPI_SUCCESS) then
         write(*,*) 'allgatherv_2l: mpi_allgatherv return ',ierr
         return
      endif

!    ------------------------------
!    copy from recv buffer to array
!    ------------------------------

      do j=1,min(n2, npes*segsize)
        i1 = 1 + (j-1)*n1
        i2 = i1 + n1 - 1

        j1 = 1 + (j-1)*n1
        j2 = j1 + n1 - 1
        array( j1:j2 ) = recvbuf( i1:i2 )
      enddo

!     --------
!     clean up
!     --------
      deallocate( sendbuf,                                               &
     &            recvbuf,                                               &
     &            displs,                                                &
     &            recvcounts,                                            &
     &            stat=istat )

      return
      end subroutine allgatherv_2l
      subroutine allgatherv_2s( array, n1,n2, segsize, comm )
!     ----------------------------
!     Perform Allgatherv on 2D array
!     layer to call MPI_Allgatherv
!     ----------------------------
      implicit none
      integer, intent(in) :: n1, n2, segsize
      integer, intent(in) :: comm

      real, intent(inout) :: array(n1*n2)
      integer, parameter :: mpi_type = MPI_REAL
      real, allocatable, dimension(:) :: sendbuf, recvbuf

      integer, allocatable, dimension(:) :: displs, recvcounts
      integer :: mype, npes
      integer :: sendcount, sendtype, recvtype
      integer :: ierr, istat, sendbuf_size, recvbuf_size 
      integer :: istart,iend,isize, pe, j,i1,i2, j1,j2
   

      ierr = MPI_SUCCESS
      call MPI_Comm_rank( comm, mype, ierr )
      if (ierr.ne.MPI_SUCCESS) then
         write(*,*) 'allgatherv_2s: mpi_comm_rank return ',ierr
         return
      endif

      call MPI_Comm_size( comm, npes, ierr )
      if (ierr.ne.MPI_SUCCESS) then
         write(*,*) 'allgatherv_2s: mpi_comm_size return ',ierr
         return
      endif

      sendbuf_size = max(1, n1*segsize)
 
      isize = min( segsize*npes, n2 )
      recvbuf_size = max(1, n1*isize)
      allocate( sendbuf(sendbuf_size),                                   &
     &          recvbuf(recvbuf_size),                                   &
     &          displs(0:(npes-1)),                                      &
     &          recvcounts(0:(npes-1)),                                  &
     &          stat=istat )
      if (istat.ne.0) then
        write(*,*) 'allgatherv_2s: alloc recvbuf, n1,n2 ',n1,n2
        stop '** allocation error ** '
      endif

      do pe=0,npes-1
         istart = 1 + pe*segsize
         iend = min( n2, istart + segsize-1)
         isize = max(0,iend - istart + 1)

         displs(pe) =  pe*(n1*segsize)
         recvcounts(pe) = (n1*isize)
      enddo

!     ----------------
!     copy into send buffer
!     ----------------
      pe = mype
      istart = 1 + pe*segsize
      iend = min(n2, istart + segsize-1)
      isize = max(0,iend-istart+1)

      sendcount = isize*n1
      
      do j=istart,iend
        i1 = (j-istart)*n1 + 1
        i2 = i1 + n1-1

        
        j1 = 1 + (j-1)*n1
        j2 = j1 + n1 - 1
        sendbuf( i1:i2 ) = array( j1:j2 )
      enddo


      if (idebug >= 1) then
        if (mype .eq. (npes-1)) then
          do pe=0,npes-1
           write(*,*) 'pe,displs,recvcount ',                           &
     &                 pe,displs(pe),recvcounts(pe)
          enddo
        endif
      endif

      sendtype = mpi_type
      recvtype = sendtype
      call MPI_Allgatherv( sendbuf, sendcount, sendtype,                 &
     &                     recvbuf, recvcounts, displs, recvtype,        &
     &                     comm, ierr)
      if (ierr.ne.MPI_SUCCESS) then
         write(*,*) 'allgatherv_2s: mpi_allgatherv return ',ierr
         return
      endif

!    ------------------------------
!    copy from recv buffer to array
!    ------------------------------

      do j=1,min(n2, npes*segsize)
        i1 = 1 + (j-1)*n1
        i2 = i1 + n1 - 1

        j1 = 1 + (j-1)*n1
        j2 = j1 + n1 - 1
        array( j1:j2 ) = recvbuf( i1:i2 )
      enddo

!     --------
!     clean up
!     --------
      deallocate( sendbuf,                                               &
     &            recvbuf,                                               &
     &            displs,                                                &
     &            recvcounts,                                            &
     &            stat=istat )

      return
      end subroutine allgatherv_2s


      subroutine allgatherv_2d_ring( array, n1,n2, segsize, comm )
!     ----------------------------
!     Perform Allgatherv on 2D array
!     layer to use ring algorithm 
!     ----------------------------
      implicit none
      integer, intent(in) :: n1, n2, segsize
      integer, intent(in) :: comm

      real(8), intent(inout) :: array(n1*n2)
      integer, parameter :: mpi_type = MPI_REAL8
      real(8) :: sendbuf(n1*segsize),recvbuf(n1*segsize)

      integer :: mype, npes
      integer :: sendcount, recvcount 
      integer :: ierr 
      integer :: istart,iend,isize,  j,i1,i2, j1,j2, istep
      integer :: recv_request, send_request
      integer :: pe , source_pe, dest_pe
      integer :: tag
   

      ierr = MPI_SUCCESS
      call MPI_Comm_rank( comm, mype, ierr )
      if (ierr.ne.MPI_SUCCESS) then
         write(*,*) 'allgatherv_2d_ring: mpi_comm_rank return ',ierr
         return
      endif

      call MPI_Comm_size( comm, npes, ierr )
      if (ierr.ne.MPI_SUCCESS) then
         write(*,*) 'allgatherv_2d_ring: mpi_comm_size return ',ierr
         return
      endif


      
      do istep=1,npes-1

        pe = mod( (mype-1+npes)  + (istep-1), npes)
        istart = 1 + pe * segsize
        iend = min(n2,  istart + segsize-1)
        isize = max(0, iend - istart + 1)
        source_pe = mod( (mype+1+npes), npes)
        tag = source_pe + istep

        recvcount = n1*segsize
        call MPI_Irecv( recvbuf, recvcount, mpi_type, source_pe, tag,     &
     &                  comm, recv_request, ierr )
        call assert( ierr.eq.MPI_SUCCESS,                                 &
     &       'allgather_2d_ring:mpi_irecv',ierr)


        pe = mod( mype + (istep-1), npes)
        istart = 1+ pe * segsize
        iend = min(n2, istart + segsize-1)
        isize = max(0, iend - istart + 1)
        dest_pe = mod( (mype-1+npes), npes )
        tag = mype + istep


!       ----------------------        
!       copy data into sendbuf
!       ----------------------        
        do j=istart,iend
          i1 = 1 + (j-istart)*n1
          i2 = i1 + n1-1

          j1 = 1 + (j-1)*n1
          j2 = j1 + n1-1
          sendbuf( i1:i2 ) = array( j1:j2 )
        enddo

        sendcount = n1*segsize
        call MPI_Isend(sendbuf, sendcount, mpi_type, dest_pe, tag,       &
     &                 comm, send_request,ierr)
        call assert(ierr.eq.MPI_SUCCESS,                                 &
     &              'allgather_2d_ring:mpi_isend',ierr)


!       -------------------------         
!       wait for incoming message
!       -------------------------         
        call MPI_Wait( recv_request, MPI_STATUS_IGNORE, ierr )
        call assert(ierr.eq.MPI_SUCCESS,                                 &
     &              'allgather_2d_ring:mpi_wait recv',ierr)
         
!       ----------------------
!       copy data from recvbuf
!       ----------------------
        pe = mod( (mype+1+npes)  + (istep-1), npes)
        istart = 1 + pe * segsize
        iend = min(n2,  istart + segsize-1)
        do  j=istart,iend
           i1 = 1 + (j-istart)*n1
           i2 = i1 + n1-1

           j1 = 1 + (j-1)*n1
           j2 = j1 + n1 - 1

           array( j1:j2 ) = recvbuf( i1:i2 )
        enddo



!       -------------------------
!       wait for outgoing message
!       -------------------------
        call MPI_Wait( send_request, MPI_STATUS_IGNORE, ierr )
        call assert(ierr.eq.MPI_SUCCESS,                                 &
     &              'allgather_2d_ring:mpi_wait send',ierr)

      enddo
        
      return
      end subroutine allgatherv_2d_ring
      subroutine allgatherv_2i_ring( array, n1,n2, segsize, comm )
!     ----------------------------
!     Perform Allgatherv on 2D array
!     layer to use ring algorithm
!     ----------------------------
      implicit none
      integer, intent(in) :: n1, n2, segsize
      integer, intent(in) :: comm

      integer, intent(inout) :: array(n1*n2)
      integer, parameter :: mpi_type = MPI_INTEGER
      integer :: sendbuf(n1*segsize),recvbuf(n1*segsize)

      integer :: mype, npes
      integer :: sendcount, recvcount 
      integer :: ierr 
      integer :: istart,iend,isize,  j,i1,i2, j1,j2, istep
      integer :: recv_request, send_request
      integer :: pe , source_pe, dest_pe
      integer :: tag
   

      ierr = MPI_SUCCESS
      call MPI_Comm_rank( comm, mype, ierr )
      if (ierr.ne.MPI_SUCCESS) then
         write(*,*) 'allgatherv_2i_ring: mpi_comm_rank return ',ierr
         return
      endif

      call MPI_Comm_size( comm, npes, ierr )
      if (ierr.ne.MPI_SUCCESS) then
         write(*,*) 'allgatherv_2i_ring: mpi_comm_size return ',ierr
         return
      endif


      
      do istep=1,npes-1

        pe = mod( (mype-1+npes)  + (istep-1), npes)
        istart = 1 + pe * segsize
        iend = min(n2,  istart + segsize-1)
        isize = max(0, iend - istart + 1)
        source_pe = mod( (mype+1+npes), npes)
        tag = source_pe + istep

        recvcount = n1*segsize
        call MPI_Irecv( recvbuf, recvcount, mpi_type, source_pe, tag,     &
     &                  comm, recv_request, ierr )
        call assert( ierr.eq.MPI_SUCCESS,                                 &
     &       'allgather_2i_ring:mpi_irecv',ierr)


        pe = mod( mype + (istep-1), npes)
        istart = 1+ pe * segsize
        iend = min(n2, istart + segsize-1)
        isize = max(0, iend - istart + 1)
        dest_pe = mod( (mype-1+npes), npes )
        tag = mype + istep


!       ----------------------        
!       copy data into sendbuf
!       ----------------------        
        do j=istart,iend
          i1 = 1 + (j-istart)*n1
          i2 = i1 + n1-1

          j1 = 1 + (j-1)*n1
          j2 = j1 + n1-1
          sendbuf( i1:i2 ) = array( j1:j2 )
        enddo

        sendcount = n1*segsize
        call MPI_Isend(sendbuf, sendcount, mpi_type, dest_pe, tag,       &
     &                 comm, send_request,ierr)
        call assert(ierr.eq.MPI_SUCCESS,                                 &
     &              'allgather_2i_ring:mpi_isend',ierr)


!       -------------------------         
!       wait for incoming message
!       -------------------------         
        call MPI_Wait( recv_request, MPI_STATUS_IGNORE, ierr )
        call assert(ierr.eq.MPI_SUCCESS,                                 &
     &              'allgather_2i_ring:mpi_wait recv',ierr)
         
!       ----------------------
!       copy data from recvbuf
!       ----------------------
        pe = mod( (mype+1+npes)  + (istep-1), npes)
        istart = 1 + pe * segsize
        iend = min(n2,  istart + segsize-1)
        do  j=istart,iend
           i1 = 1 + (j-istart)*n1
           i2 = i1 + n1-1

           j1 = 1 + (j-1)*n1
           j2 = j1 + n1 - 1

           array( j1:j2 ) = recvbuf( i1:i2 )
        enddo



!       -------------------------
!       wait for outgoing message
!       -------------------------
        call MPI_Wait( send_request, MPI_STATUS_IGNORE, ierr )
        call assert(ierr.eq.MPI_SUCCESS,                                 &
     &              'allgather_2i_ring:mpi_wait send',ierr)

      enddo
        
      return
      end subroutine allgatherv_2i_ring
      subroutine allgatherv_2l_ring( array, n1,n2, segsize, comm )
!     ----------------------------
!     Perform Allgatherv on 2D array
!     layer to use ring algorithm
!     ----------------------------
      implicit none
      integer, intent(in) :: n1, n2, segsize
      integer, intent(in) :: comm

      logical, intent(inout) :: array(n1*n2)
      integer, parameter :: mpi_type = MPI_LOGICAL
      logical :: sendbuf(n1*segsize),recvbuf(n1*segsize)

      integer :: mype, npes
      integer :: sendcount, recvcount 
      integer :: ierr 
      integer :: istart,iend,isize,  j,i1,i2, j1,j2, istep
      integer :: recv_request, send_request
      integer :: pe , source_pe, dest_pe
      integer :: tag
   

      ierr = MPI_SUCCESS
      call MPI_Comm_rank( comm, mype, ierr )
      if (ierr.ne.MPI_SUCCESS) then
         write(*,*) 'allgatherv_2l_ring: mpi_comm_rank return ',ierr
         return
      endif

      call MPI_Comm_size( comm, npes, ierr )
      if (ierr.ne.MPI_SUCCESS) then
         write(*,*) 'allgatherv_2l_ring: mpi_comm_size return ',ierr
         return
      endif


      
      do istep=1,npes-1

        pe = mod( (mype-1+npes)  + (istep-1), npes)
        istart = 1 + pe * segsize
        iend = min(n2,  istart + segsize-1)
        isize = max(0, iend - istart + 1)
        source_pe = mod( (mype+1+npes), npes)
        tag = source_pe + istep

        recvcount = n1*segsize
        call MPI_Irecv( recvbuf, recvcount, mpi_type, source_pe, tag,     &
     &                  comm, recv_request, ierr )
        call assert( ierr.eq.MPI_SUCCESS,                                 &
     &       'allgather_2d_ring:mpi_irecv',ierr)


        pe = mod( mype + (istep-1), npes)
        istart = 1+ pe * segsize
        iend = min(n2, istart + segsize-1)
        isize = max(0, iend - istart + 1)
        dest_pe = mod( (mype-1+npes), npes )
        tag = mype + istep


!       ----------------------        
!       copy data into sendbuf
!       ----------------------        
        do j=istart,iend
          i1 = 1 + (j-istart)*n1
          i2 = i1 + n1-1

          j1 = 1 + (j-1)*n1
          j2 = j1 + n1-1
          sendbuf( i1:i2 ) = array( j1:j2 )
        enddo

        sendcount = n1*segsize
        call MPI_Isend(sendbuf, sendcount, mpi_type, dest_pe, tag,       &
     &                 comm, send_request,ierr)
        call assert(ierr.eq.MPI_SUCCESS,                                 &
     &              'allgather_2d_ring:mpi_isend',ierr)


!       -------------------------         
!       wait for incoming message
!       -------------------------         
        call MPI_Wait( recv_request, MPI_STATUS_IGNORE, ierr )
        call assert(ierr.eq.MPI_SUCCESS,                                 &
     &              'allgather_2d_ring:mpi_wait recv',ierr)
         
!       ----------------------
!       copy data from recvbuf
!       ----------------------
        pe = mod( (mype+1+npes)  + (istep-1), npes)
        istart = 1 + pe * segsize
        iend = min(n2,  istart + segsize-1)
        do  j=istart,iend
           i1 = 1 + (j-istart)*n1
           i2 = i1 + n1-1

           j1 = 1 + (j-1)*n1
           j2 = j1 + n1 - 1

           array( j1:j2 ) = recvbuf( i1:i2 )
        enddo



!       -------------------------
!       wait for outgoing message
!       -------------------------
        call MPI_Wait( send_request, MPI_STATUS_IGNORE, ierr )
        call assert(ierr.eq.MPI_SUCCESS,                                 &
     &              'allgather_2d_ring:mpi_wait send',ierr)

      enddo
        
      return
      end subroutine allgatherv_2l_ring
      subroutine allgatherv_2s_ring( array, n1,n2, segsize, comm )
!     ----------------------------
!     Perform Allgatherv on 2D array
!     layer to use ring algorithm
!     ----------------------------
      implicit none
      integer, intent(in) :: n1, n2, segsize
      integer, intent(in) :: comm

      real, intent(inout) :: array(n1*n2)
      integer, parameter :: mpi_type = MPI_REAL
      real :: sendbuf(n1*segsize),recvbuf(n1*segsize)

      integer :: mype, npes
      integer :: sendcount, recvcount 
      integer :: ierr 
      integer :: istart,iend,isize,  j,i1,i2, j1,j2, istep
      integer :: recv_request, send_request
      integer :: pe , source_pe, dest_pe
      integer :: tag
   

      ierr = MPI_SUCCESS
      call MPI_Comm_rank( comm, mype, ierr )
      if (ierr.ne.MPI_SUCCESS) then
         write(*,*) 'allgatherv_2s_ring: mpi_comm_rank return ',ierr
         return
      endif

      call MPI_Comm_size( comm, npes, ierr )
      if (ierr.ne.MPI_SUCCESS) then
         write(*,*) 'allgatherv_2s_ring: mpi_comm_size return ',ierr
         return
      endif


      
      do istep=1,npes-1

        pe = mod( (mype-1+npes)  + (istep-1), npes)
        istart = 1 + pe * segsize
        iend = min(n2,  istart + segsize-1)
        isize = max(0, iend - istart + 1)
        source_pe = mod( (mype+1+npes), npes)
        tag = source_pe + istep

        recvcount = n1*segsize
        call MPI_Irecv( recvbuf, recvcount, mpi_type, source_pe, tag,     &
     &                  comm, recv_request, ierr )
        call assert( ierr.eq.MPI_SUCCESS,                                 &
     &       'allgather_2s_ring:mpi_irecv',ierr)


        pe = mod( mype + (istep-1), npes)
        istart = 1+ pe * segsize
        iend = min(n2, istart + segsize-1)
        isize = max(0, iend - istart + 1)
        dest_pe = mod( (mype-1+npes), npes )
        tag = mype + istep


!       ----------------------        
!       copy data into sendbuf
!       ----------------------        
        do j=istart,iend
          i1 = 1 + (j-istart)*n1
          i2 = i1 + n1-1

          j1 = 1 + (j-1)*n1
          j2 = j1 + n1-1
          sendbuf( i1:i2 ) = array( j1:j2 )
        enddo

        sendcount = n1*segsize
        call MPI_Isend(sendbuf, sendcount, mpi_type, dest_pe, tag,       &
     &                 comm, send_request,ierr)
        call assert(ierr.eq.MPI_SUCCESS,                                 &
     &              'allgather_2s_ring:mpi_isend',ierr)


!       -------------------------         
!       wait for incoming message
!       -------------------------         
        call MPI_Wait( recv_request, MPI_STATUS_IGNORE, ierr )
        call assert(ierr.eq.MPI_SUCCESS,                                 &
     &              'allgather_2s_ring:mpi_wait recv',ierr)
         
!       ----------------------
!       copy data from recvbuf
!       ----------------------
        pe = mod( (mype+1+npes)  + (istep-1), npes)
        istart = 1 + pe * segsize
        iend = min(n2,  istart + segsize-1)
        do  j=istart,iend
           i1 = 1 + (j-istart)*n1
           i2 = i1 + n1-1

           j1 = 1 + (j-1)*n1
           j2 = j1 + n1 - 1

           array( j1:j2 ) = recvbuf( i1:i2 )
        enddo



!       -------------------------
!       wait for outgoing message
!       -------------------------
        call MPI_Wait( send_request, MPI_STATUS_IGNORE, ierr )
        call assert(ierr.eq.MPI_SUCCESS,                                 &
     &              'allgather_2s_ring:mpi_wait send',ierr)

      enddo
        
      return
      end subroutine allgatherv_2s_ring
      end module comm_mod
