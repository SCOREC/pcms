
program main
  IMPLICIT NONE
  include 'mpif.h'
  INTEGER :: ierr
  CALL MPI_INIT(ierr)
  CALL test_htable(.false.)
  CALL test_htable(.true.)
  CALL MPI_FINALIZE(ierr)
end program

subroutine test_htable(restart)
  use futils
  use hashtable
  IMPLICIT NONE
  include 'mpif.h'
  TYPE(BUFFER_TYPE) :: hbuf
  LOGICAL, INTENT(IN) :: restart
  INTEGER :: fresid
  TYPE(lel), POINTER :: foundel
  CHARACTER(len=32) :: text
  CHARACTER(len=256) :: filename='htable.h5'   ! Result file name
  CHARACTER(len=256) :: groupname='/zero_dim'
  INTEGER :: i,j,k
  INTEGER :: nrun, ierr, me_world
  nrun = 20

  CALL mpi_comm_rank(MPI_COMM_WORLD, me_world, ierr)

  IF(me_world==0) THEN
     IF(.NOT.restart) THEN
        CALL creatf(filename, fresid, 'Test of 0D buffer module htable')
        PRINT *, 'Create HDF5 file  ', TRIM(filename), ' :fresid: ',fresid
        CALL creatg(fresid, groupname, 'Simulation results')
     ELSE
        CALL openf(filename, fresid)
        PRINT *, 'Reopen HDF5 file  ', TRIM(filename), ' :fresid: ',fresid
     END IF
  END IF

  ! BUFFER module does not need to know about restart.
  CALL htable_init(hbuf,4)
  CALL set_htable_fileid(hbuf,fresid,groupname)

  DO i=1,nrun
     CALL add_record(hbuf,"time",   "Simulation time",      (1.0d0*i) )
     CALL add_record(hbuf,"time_add","Simulation time add", (1.0d0*i), MPI_COMM_WORLD,MPI_SUM   )
     CALL add_record(hbuf,"nnodes", "Number of processors", 1.0d0,          MPI_COMM_WORLD         )
     CALL add_record(hbuf,"max",    "MAX(3.0*me_world)",    3.0d0*me_world, MPI_COMM_WORLD,MPI_SUM )
     CALL htable_endstep(hbuf)
  END DO

  CALL htable_hdf5_flush(hbuf)

  IF(me_world==0) THEN
     CALL closef(fresid)
  END IF

end subroutine test_htable
