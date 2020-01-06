PROGRAM main
!
!   Particle array using fixed dims
!
  USE futils
  IMPLICIT NONE
  CHARACTER(len=32) :: file='part.h5'
  CHARACTER(len=32) :: name
  INTEGER, PARAMETER :: npart=8192
  INTEGER :: ierr, n, fid, istep, nrun=10
  DOUBLE PRECISION :: r0=100.0, a=20.0, time, pi, temp1(npart), temp2(npart)
  DOUBLE PRECISION, TARGET  :: part(3, npart)   ! (r, z, phi) coordinates
  DOUBLE PRECISION, DIMENSION(:), POINTER :: r, z, phi
!===========================================================================
!                   1. Prologue
!
!
  IF( command_argument_count() > 0 ) THEN
     CALL get_command_argument(1, file, n, ierr)
  END IF
!
  CALL creatf(file, fid, 'A simple simulation')
  WRITE(*,'(a)') file(1:LEN_TRIM(file))//' open'
!
  CALL creatg(fid, "/part", "1D profiles")        ! Group
!
  r   => part(1,:)
  z   => part(2,:)
  phi => part(3,:)
!===========================================================================
!                   2. Time loop
!
  pi = 4.0d0*ATAN(1.0d0)
  DO istep=1,nrun
     time = istep
!
     CALL RANDOM_NUMBER(temp1); temp1 = a*temp1
     CALL RANDOM_NUMBER(temp2); temp2 = 2.*pi*temp2
     r = r0 + temp1*COS(temp2)
     z = temp1*SIN(temp2)
     CALL RANDOM_NUMBER(temp1(1:npart))
     phi = 2.0*pi*temp1
!
     WRITE(name,'(a,i3.3)') "/part/", istep
     CALL putarr(fid, name, part, "(r, z, phi) coordinates")  ! dataset
     CALL attach(fid, name, "time", time)                    ! Attr on dataset
     CALL attach(fid, name, "step", istep)
     WRITE(*,'(a,i4,a)') 'Step', istep, ' done'
  END DO
  CALL attach(fid, "/part", "nsteps", nrun)                  ! Attr on group
!===========================================================================
!                   9. Epilogue
!
  CALL closef(fid)
END PROGRAM main
