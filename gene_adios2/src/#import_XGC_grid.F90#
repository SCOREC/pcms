#include "intrinsic_sizes.h"
#include "redef.h"

module import_XGC_grid
  use communications
  use par_other, only: pi
  use discretization, only: mype, my_pex, my_pey, my_pez, my_pev, my_pew,&
                         &  n_procs_x, n_procs_z,li0,li1,li2,nx0,lk1,lk2,nz0,nzb
  use geometry, only: geomdir,geomfile, x_def, q_prof
  use coordinates, only:  xval_a,zval
  use par_in, only: file_extension, diagdir
  use file_io, only: get_unit_nr
  use par_poloidal_planes, only: XGC_file, data_fs_int, z_out_nonuni, pre_proc,&
                               & mylk0, mylk1, mylk2, &
                               & R, Z, q_XGC, x_XGC, x_GENE_wb, z_GENE_wb,&
                               & counts_nonuni, displ_nonuni, &
                               & counts_nonuni_z, displ_nonuni_z, extra_zero,&
                               & block_start, block_end, block_count, &
                               & active_nodes, XGC_inboard
  use lagrange_interpolation, only: lag3interp, lag3interp_2D
  use boundary_exchange_XGC, only: initialize_boundary_exchange_XGC, finalize_boundary_exchange_XGC, &
                                 & exchange_zbase_1D_XGC, z_boundary_real_XGC
!#ifdef COUPLE_XGC
  use coupling_core_gene, only: check_coupler
!#endif
  use futils

  implicit none

  public :: read_XGC_grid, finalize_import_XGC_grid
  public :: reshuffle_nodes, ReshuffleNodes

  private

  integer :: nnodes
  integer :: nsurfs
  real, dimension(:), allocatable ::R_XGC, Z_XGC
  integer, dimension(:), allocatable :: nodes_per_fs
  type(data_fs_int),dimension(:), allocatable :: sorted_nodes
contains

  subroutine read_XGC_grid
    integer :: shift_x

    if (mype.eq.0) then
       WRITE(*,'(A)') "**********************************************  "
       WRITE(*,'(A)') "*** Importing XGC grid "//trim(XGC_file)//"***  "
       WRITE(*,'(A)') "**********************************************  "
    endif


    if (pre_proc) then
       call get_radial_grid_fast
    else
       call read_XGC_nodes

       call read_XGC_aif

       call sort_nodes

       call convert_radial_grid
    endif

    if (extra_zero) then
       if (mype.eq.0) write(*,'(A)')'  Adding an extra point to XGC '
       allocate(XGC_inboard(4,li1:li2))
    endif

    call distribute_flux_surfaces(shift_x)

    call distribute_poloidal_points(shift_x)

    call initialize_boundary_exchange_XGC

    call exchange_boundaries_for_XGC
#ifdef COUPLE_XGC
    call check_coupler(shift_x, nx0, active_nodes)
#endif
  endsubroutine read_XGC_grid


  subroutine read_XGC_nodes
    integer :: i
    integer :: XGC_file_id, iostat
    integer :: dummy(3), dummy_ind
    real :: dummy_phi

    CALL get_unit_nr(XGC_file_id)
    open(XGC_file_id, file=trim(geomdir)//trim(XGC_file)//'.node',&
        access='sequential',Status='old', Action='read')

    read (XGC_file_id,'(I6, I1, I1, I1)',iostat=iostat)nnodes, dummy

    if (mype.eq.0) write(*,'(A,I6,A)')'   reading in ',nnodes,' nodes'

    allocate(R_XGC(nnodes))
    allocate(Z_XGC(nnodes))

    do i=1,nnodes
       read (XGC_file_id,*,iostat=iostat) &
           & dummy(1),R_XGC(i),Z_XGC(i),dummy_phi
    end do

  end subroutine read_XGC_nodes


  subroutine read_XGC_aif
    integer :: XGC_file_id, iostat

    CALL get_unit_nr(XGC_file_id)
    open(XGC_file_id, file=trim(geomdir)//trim(XGC_file)//'.flx.aif',&
         access='sequential',Status='old', Action='read')

    read(XGC_file_id,'(I6)',iostat=iostat)nsurfs

    if (mype.eq.0) write(*,'(A,I4,A)')'  Found ',nsurfs,' flux-surfaces'

    allocate(nodes_per_fs(nsurfs))

    read(XGC_file_id,*,iostat=iostat)nodes_per_fs

    close(XGC_file_id)

!    if (mype.eq.0) write(*,*)nodes_per_fs

  end subroutine read_XGC_aif


  subroutine convert_radial_grid
    integer :: hdf5_ioutgyro
    integer :: NRBOX, NZBOX, NPSI
    real, dimension(:), allocatable:: x_in, psi_in, psi_XGC, RBOX, ZBOX
    real, dimension(:,:), allocatable:: psiRZ
    integer :: surf

    call openf(trim(geomdir)//'/'//trim(geomfile), hdf5_ioutgyro)

    !psi
    call getatt(hdf5_ioutgyro, "/data", "NPSI", NPSI)

    allocate(x_in(1:NPSI), psi_in(1:NPSI))

    !radial coordinate
    call getarr(hdf5_ioutgyro, "/data/grid/PSI"  , psi_in)
    select case (x_def)
    case('arho_t')
       call getarr(hdf5_ioutgyro, "/data/var1d/rho_tor"  , x_in)
       x_in=x_in(:)/maxval(x_in)
    case default
       if (mype.eq.0) print*, 'I work only with arho_t'
       stop
    end select

    !inverse mappings
    call getatt(hdf5_ioutgyro, "/data", "NRBOX", NRBOX)
    call getatt(hdf5_ioutgyro, "/data", "NZBOX", NZBOX)

    allocate(RBOX(1:NRBOX))
    call getarr(hdf5_ioutgyro, "/data/var1d/rmesh",   RBOX)

    allocate(ZBOX(1:NZBOX))
    call getarr(hdf5_ioutgyro, "/data/var1d/zmesh",   ZBOX)

    allocate(psiRZ(1:NRBOX,1:NZBOX))
    call getarr(hdf5_ioutgyro, "/data/var2d/psiRZ",   psiRZ)

    call closef(hdf5_ioutgyro)

    !now interpolate
    allocate(x_XGC(nsurfs),psi_XGC(nsurfs))

    do surf=1,nsurfs
       call lag3interp_2d(psiRZ,RBOX,NRBOX,ZBOX,NZBOX,&
          & psi_XGC(surf), R_XGC(sorted_nodes(surf)%data(1)),1, &
          & Z_XGC(sorted_nodes(surf)%data(1)),1)
!       if (mype.eq.0) write(*,*)'Rb',R_XGC(ind_fs(surf)),'Zb',Z_XGC(ind_fs(surf)),'psi',psi_XGC(surf)
    end do

    call lag3interp(x_in,psi_in,NPSI,x_XGC,psi_XGC,nsurfs)

    deallocate(x_in, psi_in)
    deallocate(psiRZ)

    if (mype.eq.0)write(*,'(A,F4.2,A,F4.2)')"   XGC extends from ",x_XGC(1)," to ",x_XGC(nsurfs)

  end subroutine convert_radial_grid


  subroutine get_radial_grid_fast
    integer :: hdf5_ioutgyro

    real, dimension(:), allocatable:: x_in, psi_in, psi_XGC, RBOX, ZBOX
    real, dimension(:,:), allocatable:: psiRZ
!    real,dimension(nx0) :: q_in, x_tmp
    integer :: surf

    call openf(trim(geomdir)//trim(XGC_file), hdf5_ioutgyro)

    call getatt(hdf5_ioutgyro, "/data", 'NSURF', nsurfs)

    allocate(x_XGC(0:nsurfs-1))

    !radial coordinate
    call getarr(hdf5_ioutgyro, "/data/x"  , x_XGC)

    call closef(hdf5_ioutgyro)

    if (mype.eq.0) write(*,'(A,I4,A)')'  Found ',nsurfs,' flux_surfaces'

  end subroutine get_radial_grid_fast


  subroutine distribute_flux_surfaces(shift_x)
    integer, intent(out):: shift_x

    integer :: i,ind_s, ierr
    logical :: inside
    integer :: myli0, myli1,myli2

    call distribute_points(x_XGC,1,nsurfs,nx0, &
                   & myli1,myli2,myli0)

    ind_s=minloc(abs(x_XGC-xval_a(0)),1)
    shift_x=ind_s-1

!if ((my_pez+my_pev+my_pew).eq.0) then
!write(*,*)'surf',my_pex,myli1,myli2,myli0,shift_x
!endif
    if (mype.eq.0) write(*,"(A,I4,A)")'  Starting at XGC flux-surface', shift_x,' (0 is axis)'

  end subroutine distribute_flux_surfaces


  subroutine sort_nodes
    integer :: i,ip,i1,i2,i_surface=1,i0

    allocate(sorted_nodes(nsurfs))
    ip=1
    do i_surface=1,nsurfs
       i1=ip
       i2=i1+nodes_per_fs(i_surface)-1
       i0=i2-i1+1
       sorted_nodes(i_surface)%num=i0

       allocate(sorted_nodes(i_surface)%data(1:i0))
       sorted_nodes(i_surface)%data= (/(i,i=i1,i2)/) !
       ip=i2+1

    enddo

  end subroutine sort_nodes


  subroutine distribute_poloidal_points(shift_x)
    integer, intent(in) :: shift_x
    integer,dimension(0:n_procs_x-1) :: inds
    integer :: i, ierr

    !allocate the indexes of poloidal points
    allocate(mylk0(li1:li2))
    allocate(mylk1(li1:li2))
    allocate(mylk2(li1:li2))

    call define_zinds(shift_x)

    block_count=0
    if (extra_zero) then
       do i=li1,li2
          block_count=block_count+z_out_nonuni(i)%num-1
       end do
    else
       do i=li1,li2
          block_count=block_count+z_out_nonuni(i)%num
       end do
    endif

    !total number of nodes
    call MPI_allreduce(block_count, active_nodes, 1, MPI_INTEGER, MPI_SUM, mpi_comm_x, ierr);

    inds=0
    inds(my_pex)=block_count

    !Totla number of nodes per processor
    call MPI_allgather(MPI_IN_PLACE,1,MPI_integer,inds,1,&
                     & MPI_integer,mpi_comm_x,ierr)


    block_start=0
    do i=0,my_pex-1
       block_start=block_start+inds(i)
    end do
    block_end = block_start+block_count-1
    if ((my_pey+my_pez+my_pev+my_pew).eq.0) write(*,'(A,I2,A,I7,A,I7,A)') '  x-proc',i,': starts at node ', &
                                        & block_start,' and reads ', inds(my_pex),' nodes'
    if (mype.eq.0) write(*,'(A,I6,A)')'  Using ',active_nodes,' nodes'

    if ((my_pey+my_pez+my_pev+my_pew).eq.0) then
          write(*,'(A,4I6, I8)')"    x-proc. ",my_pex,li0,li1,li2,block_count
          write(*,*)'blocks',my_pex,block_start,block_end
          write(*,"(A,I6,4G18.4)")"    x-proc. ",my_pex,xval_a(li1),xval_a(li2), &
                                        & x_XGC(shift_x+li1),x_XGC(shift_x+li2)
    end if
#if 0
    if ((my_pex+my_pey+my_pev+my_pew).eq.0) then
          write(*,'(A,I2,3I6)')"    z-proc. ",my_pez,mylk0(li1), mylk1(li1),mylk2(li1)
          write(*,'(A,4G18.4)')"    z-proc: ",zval(lk1),zval(lk2),z_out_nonuni(li1)%data(mylk1(li1)), &
                                            & z_out_nonuni(li1)%data(mylk2(li1))
    end if
#endif
  end subroutine distribute_poloidal_points


  subroutine define_zinds(shift_x)
    integer, intent(in) :: shift_x

    integer :: hdf5_ioutgyro
    character(len=FILENAME_MAX) :: dset_name
    integer :: i,k
    integer :: NRBOX, NZBOX
    real, dimension(:), allocatable :: RBOX, ZBOX, data
    real, dimension(:,:), allocatable :: chiRZ
    integer :: surf

    real,dimension(:), allocatable :: pts, tmp_pts

    if (.not.pre_proc) then
       call openf(trim(geomdir)//'/'//trim(geomfile), hdf5_ioutgyro)
       !inverse mappings
       call getatt(hdf5_ioutgyro, "/data", "NRBOX", NRBOX)
       call getatt(hdf5_ioutgyro, "/data", "NZBOX", NZBOX)

       allocate(RBOX(1:NRBOX))
       call getarr(hdf5_ioutgyro, "/data/var1d/rmesh",   RBOX)

       allocate(ZBOX(1:NZBOX))
       call getarr(hdf5_ioutgyro, "/data/var1d/zmesh",   ZBOX)

       allocate(chiRZ(1:NRBOX,1:NZBOX))
       call getarr(hdf5_ioutgyro, "/data/var2d/chiRZ",   chiRZ)

       call closef(hdf5_ioutgyro)

    else
        call openf(trim(geomdir)//trim(XGC_file), hdf5_ioutgyro)
    endif

    !loop on local flux surfaces
    allocate(z_out_nonuni(li1:li2))

    allocate(counts_nonuni_z(0:n_procs_z-1,li1:li2))
    allocate(displ_nonuni_z(0:n_procs_z-1,li1:li2))

    do i=li1,li2

       if (pre_proc) then
          write(dset_name,"(A, i10.10)") '/data/',i+shift_x

          call getatt(hdf5_ioutgyro, dset_name, 'NCHI', z_out_nonuni(i)%num)

          if (extra_zero) then
             allocate(pts(z_out_nonuni(i)%num+1))

             call getarr(hdf5_ioutgyro, dset_name, pts(2:z_out_nonuni(i)%num+1))
          else
             allocate(pts(z_out_nonuni(i)%num))

             call getarr(hdf5_ioutgyro, dset_name, pts)
          endif

       else
          z_out_nonuni(i)%num = sorted_nodes(i+1)%num

          !get chiz
          allocate(pts(z_out_nonuni(i)%num))

          do k=1,z_out_nonuni(i)%num
             call lag3interp_2d(chiRZ,RBOX,NRBOX,ZBOX,NZBOX,&
                  & pts(k), &
                  & R_XGC(sorted_nodes(i+1)%data(k)), 1,&
                  & Z_XGC(sorted_nodes(i+1)%data(k)), 1)
          end do

       end if
       !XGC starts from chi=0, GENE from -pi. Need to sort before distributing
       ! no fancy sorting since the points are already sorted in their domains
       if (extra_zero) then
          z_out_nonuni(i)%GENE_starts=minloc(pts(2:z_out_nonuni(i)%num),1)
          pts(2:z_out_nonuni(i)%num+1)=reshuffle_nodes(pts(2:z_out_nonuni(i)%num+1), &
                           & z_out_nonuni(i)%num, z_out_nonuni(i)%GENE_starts)

          !build here a base array for inteprolating at -pi
          XGC_inboard(:,i)=(/pts(z_out_nonuni(i)%num)-2*pi,  &
                           pts(z_out_nonuni(i)%num+1)-2*pi,&
                           pts(2), &
                           pts(3)/)

          pts(1)=-pi
          z_out_nonuni(i)%num =z_out_nonuni(i)%num+1
       else
          z_out_nonuni(i)%GENE_starts=minloc(pts,1)
          pts=reshuffle_nodes(pts, z_out_nonuni(i)%num, z_out_nonuni(i)%GENE_starts)
       endif


       call distribute_points(pts,3,z_out_nonuni(i)%num, nz0,&
                      & mylk1(i),mylk2(i),mylk0(i))

       call set_MPI_counts_z(i)

       allocate(z_out_nonuni(i)%data(mylk1(i)-nzb:mylk2(i)+nzb))

       z_out_nonuni(i)%data(mylk1(i):mylk2(i)) = pts(mylk1(i)+1:mylk2(i)+1)

       deallocate(pts)

    end do

  end subroutine define_zinds

  !for now this is kept separate, but is the same thign as distribute radially
  subroutine distribute_points(external_grid,direction,n_ext,n_int,myi1,myi2,myi0)
    real, dimension(:), intent(in) :: external_grid
    integer, intent(in) :: direction
    integer, intent(out) ::myi1,myi2,myi0
    integer,intent(in) ::n_ext,n_int

    real, dimension(0:n_int-1) :: internal_grid
    real :: internal_ub
    integer :: i1,i2,ind_s
    logical :: inside


    call set_distribute_indexes(direction,i1,i2,internal_grid,internal_ub)

    ind_s=minloc(abs(external_grid-internal_grid(i1)),1)

    !ind_s must be in my domain or will duplicate
    do while(external_grid(ind_s).lt.internal_grid(i1))
       ind_s=ind_s+1
    end do

    myi1=ind_s
    myi2=ind_s

    inside=.true.
    do while (inside)
       !upper limit
       if (myi2+1.gt.n_ext) exit
       !chunk
       if (external_grid(myi2+1).lt.internal_ub) then
          myi2=myi2+1
       else
          inside=.false.
       endif
    end do

    myi2=myi2-1
    myi1=myi1-1

    myi0=myi2-myi1 +1
!!if ((my_pey+my_pez+my_pew+my_pev).eq.0) write(*,*)'proc x',my_pex, myi1,myi2
  end subroutine distribute_points


  subroutine set_distribute_indexes(direction,i1,i2,internal_grid, internal_ub)
    integer,intent(in) :: direction
    integer, intent(out) :: i1,i2
    real,dimension(:), intent(out) :: internal_grid
    real, intent(out) :: internal_ub

    select case (direction)
      case(1)
        i1=li1
        i2=li2
        internal_grid=xval_a
        if (i2.eq.nx0-1) then
           internal_ub = xval_a(li2)
        else
           internal_ub = xval_a(li2+1)
        endif
      case(3)
        i1=lk1
        i2=lk2
        internal_grid=zval
        if (i2.eq.nz0-1) then
           internal_ub = pi
        else
           internal_ub = zval(i2+1)
        endif
      case default
           print *, "Wrong direction in distribute indexes"
        stop
    end select

  end subroutine set_distribute_indexes


  function ReshuffleNodes(pts, how_many, start, maxplane) result(pts_out)
    integer, intent(in) :: how_many
    integer, intent(in) :: start
    integer, intent(in) :: maxplane
    real, dimension(0:maxplane, how_many) :: pts, pts_out

    pts_out(0:maxplane, 1:how_many-start+1) = pts(0:maxplane, start:how_many)
    pts_out(0:maxplane, how_many-start+2:how_many) = pts(0:maxplane, 1:start-1)

  end function ReshuffleNodes


  function reshuffle_nodes(pts,how_many,start) result(pts_out)
    integer, intent(in) :: how_many
    integer, intent(in) :: start
    real, dimension(how_many) :: pts, pts_out

    pts_out(1:how_many-start+1)= pts(start:how_many)

    pts_out(how_many-start+2:how_many)= pts(1:start-1)

  end function reshuffle_nodes


  subroutine set_MPI_counts_z(i_fs)
    integer,intent(in) ::i_fs

    integer :: k,ierr
    integer, dimension(0:n_procs_z-1) :: tmparr

    !for gatherv
    tmparr(my_pez)=mylk0(i_fs)

    call mpi_allgather(MPI_IN_PLACE,1,MPI_integer,tmparr,1,&
                     & MPI_integer,mpi_comm_z,ierr)
    counts_nonuni_z(:,i_fs)=tmparr

    displ_nonuni_z(:,i_fs)=0
    do k=1,n_procs_z-1
       displ_nonuni_z(k,i_fs)= displ_nonuni_z(k-1,i_fs)+ &
            & counts_nonuni_z(k-1,i_fs)
    end do


!    tmparr=0
!    do k=1,n_procs_z-1
!       tmparr(k)=sum(counts_nonuni_z(0:k-1,i_fs))
!if (mype.eq.50) write(*,*)'inds',k,sum(counts_nonuni_z(0:k-1,i_fs)),counts_nonuni_z(0:k-1,i_fs)
!    end do
!if (mype.eq.50) write(*,*)'i_fs',i_fs
!if (mype.eq.50) write(*,*)'counts_nonuni_z',counts_nonuni_z(:,i_fs)
!if (mype.eq.50) write(*,*)'tmparr',tmparr
!if (mype.eq.50) write(*,*)'displ',displ_nonuni_z(:,i_fs)
!if (mype.eq.50) write(*,*)'counts_nonuni_z',counts_nonuni_z
!    call mpi_allgather(MPI_IN_PLACE,1,MPI_integer,tmparr,1,&
!                     & MPI_integer,mpi_comm_z,ierr)
!    displ_nonuni_z(:,i_fs)=tmparr

  end subroutine set_MPI_counts_z


  subroutine exchange_boundaries_for_XGC
    integer ::i

    !boudary exchange in parallel
    do i=li1,li2
!if (mype.eq.0) then
!       if (i.eq.167) print*,'167 before exc',z_out_nonuni(i)%data
!       if (i.eq.168) print*,'168 before exc',z_out_nonuni(i)%data
!       if (i.eq.169) print*,'169 before exc',z_out_nonuni(i)%data
!endif
       call exchange_zbase_1D_XGC(z_boundary_real_XGC(i), z_out_nonuni(i)%data)
!if (mype.eq.0) then
!       if (i.eq.168) then
!         print*,'after exchange',z_out_nonuni(i)%data
!       endif
!endif
    end do

  end subroutine exchange_boundaries_for_XGC


  subroutine finalize_import_XGC_grid
    integer ::i

    call finalize_boundary_exchange_XGC

    do i=li1,li2
       deallocate(z_out_nonuni(i)%data)
    end do
    deallocate(z_out_nonuni)

    deallocate(mylk0,mylk1,mylk2)
    deallocate(counts_nonuni_z,displ_nonuni_z)

    if (allocated(XGC_inboard)) deallocate(XGC_inboard)

  end subroutine finalize_import_XGC_grid

end module import_XGC_grid
