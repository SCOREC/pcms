#include "redef.h"
MODULE boundary_exchange_XGC
  Use par_mod
  USE BoundaryDescriptionModule
  USE boundary_exchange_general
  USE discretization, ONLY: pi1gl, pi2gl, pi1,pi2, li1,li2, lj1,lj2, pj0,pj1,pj2, &
       & nzb,lbz,lk1,lk2,ubz,&
       & hky, ky0_ind, my_pez, n_procs_z, pmi0, nky0
  USE communications, ONLY: communicators, COMM_X, COMM_Z, COMM_V, COMM_W,&
       & N_PE_DIMS, comm_cart, mpi_comm_z, mpi_comm_w, mpi_comm_v, mpi_comm_xy
  use par_poloidal_planes,only: mylk0, mylk1, mylk2
  use mpi
  implicit none

  PUBLIC :: initialize_boundary_exchange_XGC, finalize_boundary_exchange_XGC
  PUBLIC :: z_boundary_real_XGC, exchange_z_1D_real_XGC, exchange_zbase_1D_XGC
  PUBLIC :: z_boundary_real_GENE, exchange_z_1D_real_GENE

  TYPE(BoundaryDescription),dimension(:), allocatable :: z_boundary_real_XGC
  TYPE(BoundaryDescription) :: z_boundary_real_GENE

  PRIVATE

CONTAINS


  SUBROUTINE initialize_boundary_exchange_XGC
    ! given the unstructured grid I cannot put this with the other boundaries
    ! can call only after the grids have been set. Actually I can put it together, but
    ! I don't want

    ! Local variables
    INTEGER, DIMENSION(N_PE_DIMS) :: pe_dims, my_pe_coords
    LOGICAL, dimension(N_PE_DIMS) :: periods
    integer :: i,ierr

    CALL mpi_cart_get(comm_cart,N_PE_DIMS, pe_dims, periods, my_pe_coords, ierr)

    n_procs_w = pe_dims(2)
    n_procs_v = pe_dims(3)
    n_procs_z = pe_dims(4)
    n_procs_x = pe_dims(6)
    my_pew    = my_pe_coords(2)
    my_pev    = my_pe_coords(3)
    my_pez    = my_pe_coords(4)
    my_pex    = my_pe_coords(6)


    !z boundary
    allocate(z_boundary_real_XGC(li1:li2))
    DO i=li1,li2
       CALL initialize_type(z_boundary_real_XGC(i),1,mylk0(i)+2*nzb,nzb,nzb,1)
       IF (n_procs_z.GT.1) THEN
          CALL set_mpi_type_real(z_boundary_real_XGC(i))
       END IF
    END DO

    call initialize_type(z_boundary_real_GENE,1,lz0,nzb,nzb,1)
    IF (n_procs_z.GT.1) THEN
       CALL set_mpi_type_real(z_boundary_real_GENE)
    END IF

  END SUBROUTINE initialize_boundary_exchange_XGC


  SUBROUTINE finalize_boundary_exchange_XGC
    integer ::i

    DO i=li1,li2
       CALL finalize_type(z_boundary_real_XGC(i))
    ENDDO

    deallocate(z_boundary_real_XGC)

    call finalize_type(z_boundary_real_GENE)

  END SUBROUTINE finalize_boundary_exchange_XGC


  !Exchanges on a plane to get the -pi point
  SUBROUTINE exchange_z_1D_real_XGC(bdesc,u)
    TYPE(BoundaryDescription) :: bdesc
    real, DIMENSION(:) :: u

    IF ((bdesc%n_points.NE.SIZE(u))&
         &.OR.bdesc%count.ne.1) THEN
       PRINT*,"Wrong BoundaryDescription in call to exchange_z_1D_XGC!"
       PRINT*,"bdesc%n_points = ",bdesc%n_points,", bdesc%count = ", bdesc%count
       PRINT*,"size = ",SIZE(u)
       STOP
    END IF

    CALL exchange_z_general_nopb_real(bdesc, u,1)

  END SUBROUTINE exchange_z_1D_real_XGC

  !Exchanges on a plane to get the -pi point
  SUBROUTINE exchange_z_1D_real_GENE(bdesc,u)
    TYPE(BoundaryDescription) :: bdesc
    real, DIMENSION(:) :: u

    IF ((bdesc%n_points.NE.SIZE(u))&
         &.OR.bdesc%count.ne.1) THEN
       PRINT*,"Wrong BoundaryDescription in call to exchange_z_1D_XGC!"
       PRINT*,"bdesc%n_points = ",bdesc%n_points,", bdesc%count = ", bdesc%count
       PRINT*,"size = ",SIZE(u)
       STOP
    END IF

    CALL exchange_z_general_nopb_real(bdesc, u,1)

  END SUBROUTINE exchange_z_1D_real_GENE

  !exchangeof the basis, if XGC is truly equidistant we can make it simpler
  SUBROUTINE exchange_zbase_1D_XGC(bdesc,u)
    TYPE(BoundaryDescription) :: bdesc
    real, DIMENSION(:) :: u

    IF ((bdesc%n_points.NE.SIZE(u))&
         &.OR.bdesc%count.ne.1) THEN
       PRINT*,"Wrong BoundaryDescription in call to exchange_z_1D_XGC!"
       PRINT*,"bdesc%n_points = ",bdesc%n_points,", bdesc%count = ", bdesc%count
       PRINT*,"size = ",SIZE(u)
       STOP
    END IF

    CALL exchange_z_general_base_real(bdesc, u,1)

  END SUBROUTINE exchange_zbase_1D_XGC


  ! =======================================================
  !Exchange in z with no ppbc, useless
  SUBROUTINE exchange_z_general_nopb_real(bdesc,u, n_dim_2)
    TYPE(BoundaryDescription) :: bdesc
    INTEGER :: n_dim_2
    REAL, DIMENSION(bdesc%n_points,1:n_dim_2) :: u

    ! Local variables
    INTEGER :: n_exchanges, i_exchange

    PERFON('exz_nopb')
    IF (MOD(n_dim_2,bdesc%count).NE.0) THEN
       WRITE(*,"(2(A,I6),A)") "Wrong BoundaryDescription in call to exchange_z_general! We cannot exchange ",&
            & n_dim_2," times in z direction in blocks of ",bdesc%count,". Aborting!"
       PRINT*, "bdesc = ",bdesc
       STOP
    END IF

    bdesc%exchange_direction = 3 ! set to z direction
    n_exchanges = n_dim_2/bdesc%count
    DO i_exchange=0,n_exchanges-1
       CALL exchange_general_real(bdesc,u(:,1+i_exchange*bdesc%count:(i_exchange+1)*bdesc%count))
    END DO
    PERFOFF
  END SUBROUTINE exchange_z_general_nopb_real

  ! =======================================================
  !Exchange in z of poloidal points
  SUBROUTINE exchange_z_general_base_real(bdesc,u, n_dim_2)
    TYPE(BoundaryDescription) :: bdesc
    INTEGER :: n_dim_2
    REAL, DIMENSION(bdesc%n_points,1:n_dim_2) :: u

    ! Local variables
    INTEGER :: n_exchanges, i_exchange

    PERFON('exz_XGC')
    IF (MOD(n_dim_2,bdesc%count).NE.0) THEN
       WRITE(*,"(2(A,I6),A)") "Wrong BoundaryDescription in call to exchange_z_general! We cannot exchange ",&
            & n_dim_2," times in z direction in blocks of ",bdesc%count,". Aborting!"
       PRINT*, "bdesc = ",bdesc
       STOP
    END IF

    bdesc%exchange_direction = 3 ! set to z direction
    n_exchanges = n_dim_2/bdesc%count
    DO i_exchange=0,n_exchanges-1
       CALL exchange_general_real(bdesc,u(:,1+i_exchange*bdesc%count:(i_exchange+1)*bdesc%count))
    END DO

    IF (my_pez.EQ.0) THEN
       CALL post_apply_boundary_condition_zbase(bdesc,u,2)
    END IF
    IF (my_pez.EQ.n_procs_z-1) THEN
       CALL post_apply_boundary_condition_zbase(bdesc,u,1)
    END IF

    PERFOFF
  END SUBROUTINE exchange_z_general_base_real


  SUBROUTINE post_apply_boundary_condition_zbase(bdesc,u,side)
    TYPE(BoundaryDescription) :: bdesc
    REAL, DIMENSION(0:,0:) :: u
    INTEGER :: side

    ! Local variables
    INTEGER :: i_index, i_subarray, start_index, end_index

    IF (side.EQ.1) THEN
       DO i_index=0,SIZE(u,2)-1
          DO i_subarray = 1,bdesc%n_upper_subarrays
             start_index = bdesc%innerlast+(i_subarray-1)*bdesc%subarray_size+1
             end_index   = bdesc%innerlast+i_subarray*bdesc%subarray_size
             u(start_index:end_index,i_index)= u(start_index:end_index,i_index)+2*pi
          END DO
       END DO
    ELSE
       DO i_index=0,SIZE(u,2)-1
          DO i_subarray = 1,bdesc%n_lower_subarrays
             start_index = (i_subarray-1)*bdesc%subarray_size
             end_index   = i_subarray*bdesc%subarray_size-1
             u(start_index:end_index,i_index)= u(start_index:end_index,i_index)-2*pi
          END DO
       END DO
    END IF

  END SUBROUTINE post_apply_boundary_condition_zbase


END MODULE boundary_exchange_XGC
