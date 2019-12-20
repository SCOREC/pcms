      integer function get_gpu_streamid()
      use cudafor

      integer, allocatable, dimension(:) :: streamid_array
      integer :: ierr
      save

#ifdef USE_GPU_ASYNC
      if (.not.allocated(streamid_array)) then
         allocate( streamid_array(1), stat=ierr )
         call assert(ierr.eq.0,'alloc(streamid_array)',ierr)

         ierr = cudaStreamCreate( streamid_array(1) )
         call assert(ierr.eq.0,'cudaStreamCreate(streamid_array)',ierr)
      endif

      get_gpu_streamid = streamid_array(1)
#else
!     ----------------------
!     use the default stream
!     ----------------------
      get_gpu_streamid = 0
#endif
      return
      end function get_gpu_streamid

       

