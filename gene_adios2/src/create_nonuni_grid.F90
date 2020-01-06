#include "redef.h"
#include "intrinsic_sizes.h"

module create_nonuni_grid
  use discretization, only: li1,li2, n_procs_x, n_procs_z, lk0,mype, &
                          & my_pex, my_pey, my_pez, my_pev, my_pew, my_pespec
  use coordinates, only: deli
  use par_other, only: pi
  use par_geom, only: minor_r, rhostar
  use communications
  use coordinates, only: xval_a
  use par_poloidal_planes, only: nz_out_nonuni,z_out_nonuni, R, Z, nz_out, &
          & counts_nonuni, displ_nonuni, pts_nonuni, tot_pts_nonuni, pts_nonuni, &
          & mylk1, mylk2, mylk0, counts_nonuni_z, displ_nonuni_z
  use geometry, only: geomdir,geomfile, x_def
  use par_in, only: file_extension, diagdir
  use lagrange_interpolation, only: lag3interp_2d
#ifdef WITHFUTILS
  use futils
#endif

  implicit none

  public :: generate_grid_nonuni, finalize_zinds_nonuni

  private

contains

  subroutine generate_grid_nonuni

     call initialize_zinds_nonuni

     call generate_RZ_values

     call write_generated_grid

  end subroutine generate_grid_nonuni


  ! intiilize z indexes, not distributed for now
  subroutine initialize_zinds_nonuni
    real :: dz, unbalance
    integer :: i,k,ierr, extra
    integer, dimension(0:n_procs_z-1) :: tmparr
    integer :: n_procs_z_act


    !points per flux surface
    allocate(nz_out_nonuni(li1:li2))

    !indexes of the points per flux surface
    allocate(mylk1(li1:li2))
    allocate(mylk2(li1:li2))
    allocate(mylk0(li1:li2))

    !mpi IO related
    allocate(counts_nonuni_z(0:n_procs_z-1,li1:li2))
    allocate(displ_nonuni_z(0:n_procs_z-1,li1:li2))

    pts_nonuni=0

    do i=li1,li2
       if (nz_out.lt.0) then
          nz_out_nonuni(i)=ceiling(2*pi*(xval_a(i)*minor_r)/&
                          & (deli*minor_r*rhostar)/abs(nz_out))
       else
          nz_out_nonuni(i)=nz_out
       endif

       !distribute in z
       mylk0(i)=nz_out_nonuni(i)/n_procs_z

       n_procs_z_act=n_procs_z
       do while (mylk0(i).lt.lk0)
          mylk0(i)=nz_out_nonuni(i)/n_procs_z_act
          n_procs_z_act=n_procs_z_act-1
       end do

       !this works only if external grid is uniform
       if (n_procs_z_act.lt.n_procs_z) then
          if (mype.eq.0) write(*,"(A,I4,A)")"flux surface ",i," is coarser than GENE"
          mylk1(i)=my_pez*mylk0(i)
          mylk2(i)=min(mylk1(i) + mylk0(i) - 1,nz_out_nonuni(i))
          if (my_pez.gt. (n_procs_z_act-1)) then
             mylk0(i)=0
             mylk1(i)=0
             mylk2(i)=0
          elseif (my_pez.eq. (n_procs_z_act-1)) then
             extra=mod(nz_out_nonuni(i),n_procs_z_act)
             unbalance=real(extra)/real(mylk2(i))
             mylk0(i)=mylk0(i)+extra
             mylk2(i)=mylk1(i)+extra
          endif
       else
          mylk1(i)=my_pez*mylk0(i)
          mylk2(i)=mylk1(i) + mylk0(i) - 1
          extra=mod(nz_out_nonuni(i),n_procs_z)
          unbalance=real(extra)/real(mylk2(i))
          !last takes the remaining
          if ((my_pez.eq.(n_procs_z-1)).and.(extra.gt.0)) then
             mylk0(i)=mylk0(i)+extra
             mylk2(i)=mylk1(i)+extra
!          if ((my_pey+my_pez+my_pev+my_pew+my_pespec).eq.(n_procs_z-1)) then
!             write(*,'(G12.4,A,I4)')unbalance, &
!                       &"% unbalance in the parallel distribution for fs ",i
!          endif
          endif
       endif

       !for gatherv
       tmparr(my_pez)=mylk0(i)

       call mpi_allgather(MPI_IN_PLACE,1,MPI_integer,tmparr,1,&
                        & MPI_integer,mpi_comm_z,ierr)
       counts_nonuni_z(:,i)=tmparr

       tmparr=0
       do k=1,n_procs_z-1
          tmparr(k)=sum(counts_nonuni_z(0:k-1,i),1)
       end do

       call mpi_allgather(MPI_IN_PLACE,1,MPI_integer,tmparr,1,&
                       &MPI_integer,mpi_comm_z,ierr)
       displ_nonuni_z(:,i)=tmparr


       !total points on current x block
       pts_nonuni=pts_nonuni+nz_out_nonuni(i)

!   if ((my_pey+my_pez+my_pev+my_pew+my_pespec).eq.(0)) &
!           &  write(*,*)i,'pts_nonuni',pts_nonuni,&
!           & 'counts',counts_nonuni_z(:,i),'displ',displ_nonuni_z(:,i), &
!           & mylk0(i),mylk1(i),mylk2(i)
    end do

    !mpi IO related
    allocate(counts_nonuni(0:n_procs_x-1))
    allocate(displ_nonuni(0:n_procs_x-1))

    counts_nonuni(my_pex)=pts_nonuni

    call mpi_allgather(MPI_IN_PLACE,1,MPI_integer,counts_nonuni,1,&
            & MPI_integer,mpi_comm_x,ierr)

    tot_pts_nonuni=sum(counts_nonuni)

    displ_nonuni=0
    do i=1,n_procs_x-1
        displ_nonuni(i)=sum(counts_nonuni(0:i-1))
    end do
    call mpi_allgather(MPI_IN_PLACE,1,MPI_integer,displ_nonuni,1,&
            &MPI_integer,mpi_comm_x,ierr)


    !if (mype.eq.0) write(*,*)i,'pts_nonuni',tot_pts_nonuni,counts_nonuni,displ_nonuni


    !polidal grid output
    allocate(z_out_nonuni(li1:li2))
    do i=li1,li2
       allocate(z_out_nonuni(i)%data(mylk0(i)))
       dz = (2.0 * pi ) / nz_out_nonuni(i)
       do k=1,mylk0(i)
          z_out_nonuni(i)%data(k) = -pi + (k-1+mylk1(i))*dz
       enddo

    end do


  end subroutine initialize_zinds_nonuni


  subroutine finalize_zinds_nonuni
    integer :: i

    deallocate(counts_nonuni,displ_nonuni)

    do i=li1,li2
       deallocate(z_out_nonuni(i)%data)
    end do
    deallocate(z_out_nonuni,nz_out_nonuni)

    deallocate(mylk1,mylk2,mylk0)
    deallocate(counts_nonuni_z,displ_nonuni_z)

  end subroutine finalize_zinds_nonuni


  subroutine generate_RZ_values
    integer :: hdf5_ioutgyro
    integer :: NPSI, NCHI
    real, dimension(:), allocatable:: x_in, z_in
    real, dimension(:,:), allocatable:: R_chease, Z_chease, tmparr
    integer :: i,k,kpi
#ifdef WITHFUTILS
    call openf(trim(geomdir)//'/'//trim(geomfile), hdf5_ioutgyro)

    call getatt(hdf5_ioutgyro, "/data", "NPSI", npsi)
    call getatt(hdf5_ioutgyro, "/data", "NCHI", nchi)

    allocate(x_in(1:npsi))
    allocate(R_chease(1:npsi,1:nchi),Z_chease(1:npsi,1:nchi))

    call getarr(hdf5_ioutgyro, "/data/var2d/R",   R_chease)
    call getarr(hdf5_ioutgyro, "/data/var2d/Z",   Z_chease)

    select case (x_def)
    case('arho_t')
       call getarr(hdf5_ioutgyro, "/data/var1d/rho_tor"  , x_in)
       x_in=x_in(:)/maxval(x_in)
    case default
       stop 'Use rhotoror code something else'
    end select

    call closef(hdf5_ioutgyro)
#else
    stop*,"compile with futils to read a grid from CHEASE"
#endif

    ! add one point before interpolating
    allocate(z_in(1:nchi+1))

    do k=1,nchi+1
       z_in(k)=(k-1)*2*pi/nchi-pi
    enddo

    allocate(tmparr(1:npsi,1:nchi+1))

    ! generate R values
    do k=0,nchi-1
       kpi=modulo(nchi/2+k,nchi)+1
       tmparr(:,k+1)=R_chease(:,kpi)
    end do
    tmparr(:,nchi+1)=tmparr(:,1)

    allocate(R(li1:li2))

    do i=li1,li2
       allocate(R(i)%data(mylk0(i)))

       call lag3interp_2d(tmparr,x_in,npsi,z_in,nchi+1,&
                        & R(i)%data,xval_a(i),1, &
                        & z_out_nonuni(i)%data,mylk0(i))
    end do

    deallocate(R_chease)

    ! generate Z values
    do k=0,nchi-1
       kpi=modulo(nchi/2+k,nchi)+1
       tmparr(:,k+1)=Z_chease(:,kpi)
    end do
    tmparr(:,nchi+1)=tmparr(:,1)

    allocate(Z(li1:li2))

    do i=li1,li2
       allocate(Z(i)%data(mylk0(i)))
       call lag3interp_2d(tmparr,x_in,npsi,z_in,nchi+1,&
                        & Z(i)%data,xval_a(i),1, &
                        & z_out_nonuni(i)%data,mylk0(i))
    end do

    deallocate(Z_chease)
    deallocate(tmparr, x_in)

  end subroutine generate_RZ_values


  subroutine  write_generated_grid
#ifdef WITHFUTILS
     character(len=FILENAME_MAX) :: dset_name
     integer :: i,ierr,current_data
     real, dimension(0:tot_pts_nonuni-1) :: data_layer
     real, dimension(:), allocatable:: data_fs
     real,dimension(0:pts_nonuni-1) :: data_block

     integer :: fidout_h5

     !R
     current_data=0
     do i=li1,li2

        allocate(data_fs(0:nz_out_nonuni(i)-1))
        !gather to proc 0 the distributed along z
        call mpi_gatherv(R(i)%data(1),mylk0(i), MPI_REAL_TYPE,&
                & data_fs,counts_nonuni_z(:,i),displ_nonuni_z(:,i),&
                & MPI_REAL_TYPE,0, mpi_comm_z, ierr)
        !reshape into 1D
        data_block(i-li1+current_data:i-li1+current_data+nz_out_nonuni(i)-1)=data_fs
        current_data=current_data+nz_out_nonuni(i)-1
        deallocate(data_fs)
     enddo

     call mpi_gatherv(data_block(0),pts_nonuni, MPI_REAL_TYPE,&
            & data_layer(0),counts_nonuni,displ_nonuni,&
            & MPI_REAL_TYPE,0, mpi_comm_x, ierr)

     if (mype.eq.0) then
        call creatf(trim(diagdir)//'/dummy_grid'//trim(file_extension)//'.h5', &
               fidout_h5, "Poloidal_plane", 'd')!,mpi_comm_x)
        !R points
        write(dset_name, "(A)") "/R"
        call putarr(fidout_h5, dset_name, data_layer)
        call flushh5(fidout_h5)
     endif

     !Z
     current_data=0
     do i=li1,li2
        allocate(data_fs(0:nz_out_nonuni(i)-1))
        !gather to proc 0 the distributed along z
        call mpi_gatherv(Z(i)%data,mylk0(i), MPI_REAL_TYPE,&
                & data_fs,counts_nonuni_z(:,i),displ_nonuni_z(:,i),&
                & MPI_REAL_TYPE,0, mpi_comm_z, ierr)
        !reshape into 1D
        data_block(i-li1+current_data:i-li1+current_data+nz_out_nonuni(i)-1)=data_fs
        current_data=current_data+nz_out_nonuni(i)-1
        deallocate(data_fs)
     enddo


     call mpi_gatherv(data_block(0),pts_nonuni, MPI_REAL_TYPE,&
                & data_layer(0),counts_nonuni,displ_nonuni,&
                & MPI_REAL_TYPE,0, mpi_comm_x, ierr)

     if (mype.eq.0) then
        !Z points
        write(dset_name, "(A)") "/Z"
        call putarr(fidout_h5, dset_name,data_layer)

        call flushh5(fidout_h5)

     end if
#endif
  end subroutine write_generated_grid


end module create_nonuni_grid
