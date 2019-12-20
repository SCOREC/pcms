!MPI related function 
! for absoft compiler call external name is UPPERCASE and 
! all mpi name is lowercase


! mympif.h is the same as mpif.h except the all text is lowercase.

subroutine MY_MPI_INIT
  use sml_module
  use rem_module
  implicit none
  include 'mpif.h'
  integer :: ierror
  logical :: exist

  call mpi_init(ierror)
  call mpi_comm_rank(mpi_comm_world,sml_mype,ierror)
  if (sml_mype == 0) then
    rem_starttime = mpi_wtime()

    ! xgc_mpi_init used only in conjunction with Walltime.Remaining 
    ! file, which should have been recreated before execution began
    inquire(FILE="Walltime.Remaining",EXIST=exist)
    if (exist) then
      open(unit=14,file='xgc_mpi_init',status='replace')
      write(14,*) rem_starttime
      close(14)
    endif

  endif
  call mpi_comm_size(mpi_comm_world,sml_totalpe,ierror)
  call mpi_comm_dup(mpi_comm_world,sml_comm,ierror)
  sml_comm_null = mpi_comm_null

end subroutine MY_MPI_INIT


subroutine MY_MPI_FINALIZE
  use sml_module
  implicit none
  include 'mpif.h'
  integer :: ierror
  
  call mpi_finalize(ierror)

end subroutine MY_MPI_FINALIZE

subroutine MY_MPI_REDUCE(a,b,n)
  implicit none
  include 'mpif.h'
  integer ,intent(in):: n
  real (kind=8), intent(inout) :: a(n), b(n)
  integer :: ierror
  
  call mpi_reduce(a,b,n,mpi_real8,mpi_sum,0,mpi_comm_world,ierror)
  if (ierror.ne.mpi_success) then
    write(*,*) 'my_mpi_reduce return ierror=',ierror
  endif

end subroutine MY_MPI_REDUCE


subroutine MY_MPI_ALLREDUCE(a,b,n)
  implicit none
  include 'mpif.h'
  integer ,intent(in):: n
  real (kind=8), intent(inout) :: a(n), b(n)
  integer :: ierror
  
  call mpi_allreduce(a,b,n,mpi_real8,mpi_sum,mpi_comm_world,ierror)
  if (ierror.ne.mpi_success) then
    write(*,*) 'my_mpi_allreduce return ierror = ',ierror
  endif

end subroutine MY_MPI_ALLREDUCE
