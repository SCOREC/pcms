PROGRAM main
  IMPLICIT NONE
  REAL*8, DIMENSION(:,:), ALLOCATABLE :: a, b
  INTEGER :: i, mb=1024*1024
  DOUBLE PRECISION :: mem
!
  WRITE(*,'(a,f10.3)')'Mem at start (MB): ', mem()
  ALLOCATE(b(10,mb), a(10,mb))
  WRITE(*,'(a,f10.3)') 'Size of A (MB)', REAL(SIZE(a))*8/mb
  WRITE(*,'(a,f10.3)') 'Size of B (MB)', REAL(SIZE(b))*8/mb
!
  CALL RANDOM_NUMBER(a)
  b=2*a-1
!
  WRITE(*,'(a,f10.3)') 'Mem after allocating A & B (MB): ', mem()
  DEALLOCATE(a)
  WRITE(*,'(a,f10.3)') 'Mem after deallocating A (MB): ', mem()
!
  ALLOCATE(a(50,mb))
  WRITE(*,'(a,f10.3)') 'Size of A (MB)', REAL(SIZE(a))*8/mb
  a=0.0
  WRITE(*,'(a,f10.3)') 'Mem after allocating larger A (MB): ', mem()
END PROGRAM main
