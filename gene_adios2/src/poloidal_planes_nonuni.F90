#include "redef.h"
#include "intrinsic_sizes.h"

#undef CHECKNAN
#undef WRITE_ALL
#define H5_INIT
!TODO memory estimate
module poloidal_planes_nonuni
  use par_mod, only: time, itime
  use par_in, only: file_extension, diagdir, n0_global,write_hac,write_h5,spec
  use par_geom, only: minor_r, rhostar
  use communications
  use coordinates, only: xval_a
  use discretization, only: li1,li2,li0,lj1,lj2,lj0,lk1,lk2,lk0,pmi1gl,pmi2gl,pmx0,&
      & ky0_ind, nky0,nj0,nz0,lbz,ubz,lz0, nzb, lbi, ubi, nib, lx0, pj1,&
      & n_procs_x,n_procs_y,n_procs_z,&
      & my_pex,my_pey,my_pez,my_pev,my_pew,my_pespec,mype
  use coordinates, only: deli,delj,zval
  use lagrange_interpolation, only: lag3interp,lag3interp_complex
  use geometry, only: geom, C_y, q_prof
  use par_other, only: pi, imag
  use fourier, only: initialize_fourier_arb, finalize_fourier_arb
  use par_poloidal_planes, only: n_cuts, y_BC, phi_cut, BC_in_fourier,&
                            & z_out_nonuni,nz_out_nonuni,z_GENE_wb, &
                            & mylj0, mylj1, mylj2, &
                            & mylk0, mylk1, mylk2, &
                            & mat_from_plane, mat_to_plane, &
                            & data_unst, data_3d_nonuni, res_fact,&
                            & norm_fact_field, norm_fact_mom, istep_field_in,XGC_inboard
  use poloidal_planes_aux, only: compute_pseudo_inv, &
                            & gto_real_n, gto_fourier_n
  use par_poloidal_planes, only: pts_nonuni, tot_pts_nonuni, &
                            & counts_nonuni, displ_nonuni, &
                            & block_start, block_end, block_count, &
                            & counts_nonuni_z, displ_nonuni_z, extra_zero
  USE boundaries, only: exchange_z
  USE boundary_exchange_XGC, only: z_boundary_real_XGC, exchange_z_1D_real_XGC, &
                                 & z_boundary_real_GENE, exchange_z_1D_real_GENE

  use profiles_mod, only: nref
  use import_XGC_grid, only: reshuffle_nodes, ReshuffleNodes
  use perf_monitor

#ifdef WITHFUTILS
  use futils
  use par_poloidal_planes, only: isnap_planes, fidplanes_h5, &
       isnap_planes_from_XGC, fidplanes_h5_from_XGC, &
       fidfield_in_h5, isnap_field_in
#endif

#ifdef COUPLE_XGC
  use coupling_core_gene, only: &
      send_density, &
      receive_field, &
      cce_process_field, &
      cce_process_density, &
      write_check_file, &
      receive_gene_density
#endif

#ifdef INIT_XGC
  use coupling_core_gene, only: receive_density, receive_density_3d_nogc
#endif


  implicit none

  public :: prepare_matrixes_nonuni, xyz_to_planes_nonuni, planes_to_xyz_nonuni, &
          & write_planes, dump_original_field_wb

  public :: map_XGC_init_cond

  private

  integer :: istep_interp=0
  integer :: istep_orig, istep_orig_wb, istep_orig_real_wb, istep_interp_real
  integer :: fidorig_h5, fidorig_wb_h5, fidorig_real_wb_h5, fidinterp_real_h5

contains


  subroutine prepare_matrixes_nonuni
    real :: dy
    integer :: y_res
    integer :: i

    ! y direction in real space no dealiasing
    mylj0=2*nj0*res_fact/n_procs_y
    mylj1=my_pey*mylj0
    mylj2=mylj1+mylj0-1

    call define_y_inds_nonuni(y_res,dy)

    call initialize_fourier_arb

    !data to be exchanged.
    allocate(data_unst(li1:li2))

    !GENE data in real space
    allocate(data_3d_nonuni(li1:li2))

    do i=li1,li2
       allocate(data_unst(i)%data(0:n_cuts-1,mylk1(i):mylk2(i)))
       allocate(data_3d_nonuni(i)%data(0:2*nj0*res_fact-1,mylk1(i):mylk2(i)))
    enddo

#ifdef WRITEALL
    call initialize_diag_debug
#endif
  end subroutine prepare_matrixes_nonuni


  subroutine finalize_planes_nonuni
    integer :: i

    call finalize_fourier_arb

    do i=li1,li2
       deallocate(data_unst(i)%data,data_3d_nonuni(i)%data)
    end do

    deallocate(data_unst,data_3d_nonuni)

#ifdef WRITEALL
    call finalize_diag_debug
#endif
  end subroutine finalize_planes_nonuni


  subroutine define_y_inds_nonuni(y_res,dy)
    real, intent(OUT) :: dy
    integer, intent(OUT) :: y_res

    real :: y_cut
    real :: plane_start,plane_end
    integer :: i,k,o,tmp_ind,ind_l_tmp,ind_h_tmp

    y_res = 2 * nky0 * res_fact
    dy = delj * rhostar * minor_r / res_fact

    Call get_systime(plane_start)

    allocate(mat_to_plane(li1:li2))
    allocate(mat_from_plane(li1:li2))

    do i=li1,li2
       allocate(mat_to_plane(i)%data(0:n_cuts-1,0:mylj0-1,mylk1(i):mylk2(i)))
       allocate(mat_from_plane(i)%data(0:mylj0-1,0:n_cuts-1,mylk1(i):mylk2(i)))
    end do

    do i=li1,li2
       mat_to_plane(i)%data=0.0
       do k=mylk1(i),mylk2(i)
          do o=0,n_cuts-1
             y_cut=C_y(0)*(q_prof(i)*z_out_nonuni(i)%data(k)-phi_cut(o))/dy
             !uncomment following for shift in y
             !y_cut=y_cut+real(y_res/2)
             y_cut=mod(mod(y_cut,real(y_res))+real(y_res), real(y_res))

             tmp_ind=int(y_cut)
             ind_l_tmp=mod(mod(tmp_ind,y_res)+y_res,y_res)
             ind_h_tmp=mod(mod(tmp_ind+1,y_res)+y_res,y_res)

             mat_to_plane(i)%data(o,ind_h_tmp,k)=y_cut-tmp_ind
             mat_to_plane(i)%data(o,ind_l_tmp,k)=1.0-(y_cut-tmp_ind)

          end do

          call compute_pseudo_inv(mat_to_plane(i)%data(:,:,k),mat_from_plane(i)%data(:,:,k))

           !distribute the matrix to y_procs
!! y parll NOT SUPPORTED
!           do o=0,n_cuts-1
!              call mpi_scatter(A_inv(:,o), mylj0,MPI_REAL_TYPE,inv_mat(:,o,k,i),&
!                 & mylj0,MPI_REAL_TYPE,0,mpi_comm_y,ierr)
!           end do

       end do
    end do

    Call get_systime(plane_end)
    If (mype.eq.0) write(*,'(A,F10.3,a)')&
       &'Time for matrix inversion:',plane_end-plane_start,' sec'

  end  subroutine define_y_inds_nonuni


  subroutine xyz_to_planes_nonuni(data_to_dump)
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2), intent(IN) ::data_to_dump
    call start_timer("MAP_TO_PLANES")
    call invert_extend_nonuni(data_to_dump)
    call stop_timer("MAP_TO_PLANES")

    call start_timer("WRITE_PLANES")
    call write_planes
    call stop_timer("WRITE_PLANES")
  end subroutine xyz_to_planes_nonuni


  subroutine invert_extend_nonuni(data_to_dump_cpx)
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2), intent(IN) ::data_to_dump_cpx   !complex data

    complex, dimension(li1:li2,lj1:lj2,lbz:ubz) :: data_to_dump_cpx_wb           !with z-boundaries
    real, dimension(li1:li2,mylj1:mylj2,lbz:ubz) :: data_to_dump         !with z-boundaries
    real, dimension(lbz:ubz) :: data_extended   !local chunk in z to avoid temporary copy in laginterp
    real, dimension(:), allocatable :: loc_data   !data on XGC poloidal grid, different ech surface
    integer :: i,j,k,ierr

#ifdef WRITE_ALL
    call dump_original_field(data_to_dump_cpx)
#endif
    !extend GENE dain Fourier
    data_to_dump_cpx_wb(li1:li2,lj1:lj2,lk1:lk2)=data_to_dump_cpx
    call exchange_z(data_to_dump_cpx_wb)

#ifdef WRITE_ALL
    call dump_original_field_wb(data_to_dump_cpx_wb)
#endif
    !real space
    do k=lbz,ubz
       call gto_real_n(data_to_dump_cpx_wb(:,:,k),data_to_dump(:,:,k),li0)
    end do

#ifdef WRITE_ALL
    call dump_original_field_real(data_to_dump)
#endif

    !project
    do i=li1,li2
       allocate(loc_data(mylk1(i):mylk2(i)))

       do j=0,mylj0-1
          data_extended=data_to_dump(i,j,:)
          !interpolate on finer poloidal grid
          call lag3interp(data_extended,z_GENE_wb,lz0, &
                        & loc_data,z_out_nonuni(i)%data(mylk1(i):mylk2(i)),mylk0(i))
!          call linear_interp(z_GENE_wb,data_extended,lz0, &
!                        & z_out_nonuni(i)%data(mylk1(i):mylk2(i)),loc_data,mylk0(i))
          data_3d_nonuni(i)%data(j,:)=loc_data
       end do

       !build output
       do k=mylk1(i),mylk2(i)
          data_unst(i)%data(:,k)=matmul(mat_to_plane(i)%data(:,:,k),data_3d_nonuni(i)%data(:,k))
       enddo

       deallocate(loc_data)

    end do

  end subroutine invert_extend_nonuni


  subroutine planes_to_xyz_nonuni(data_imported)
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2), intent(out) :: data_imported

    real, dimension(:,:), allocatable :: data_fs
    real, dimension(li1:li2,mylj1:mylj2,lk1:lk2) :: data_xyz
    real, dimension(:), allocatable :: lag_in
     real, dimension(lk1:lk2) :: lag_out
    integer :: i,j,k
    real :: dz

    integer :: lk=0, ierr
#ifdef COUPLE_XGC
    call start_timer("READ_PLANES")
    call import_planes_nonuni
    call stop_timer("READ_PLANES")
#endif
#ifdef PROJ_DENS
    call import_GENE_density
#endif

    call start_timer("MAP_FROM_PLANES")
    !poloidal grid, local chunk on XGC grid
    do i=li1,li2
       allocate(lag_in(mylk1(i)-nzb:mylk2(i)+nzb))
       allocate(data_fs(mylj1:mylj2,mylk1(i):mylk2(i)))

       do k=mylk1(i),mylk2(i)
           data_fs(:,k)=matmul(mat_from_plane(i)%data(:,:,k),data_unst(i)%data(:,k))
       end do

#ifdef CHECKNAN
       do j=mylj1,mylj2
          do k=mylk1(i),mylk2(i)
             if (isnan(real(data_fs(j,k)))) then
             print*,'nan after projection',i,j,k,real(data_fs(j,k))
             stop
             end if
          end do
       end do
#endif

       do j=mylj1,mylj2
          lag_in(mylk1(i):mylk2(i))=data_fs(j,:)
          !exchange boundary for inner interpolation
          call exchange_z_1D_real_XGC(z_boundary_real_XGC(i),lag_in)

          !downsample in z
!          call linear_interp(z_out_nonuni(i)%data(mylk1(i):mylk2(i)),lag_in,mylk0(i), &
!                        & zval(lk1:lk2),lag_out,lk0)
          call lag3interp(lag_in ,z_out_nonuni(i)%data,mylk0(i)+2*nzb, &
                        & lag_out,zval(lk1:lk2),lk0)

!          call lag3interp(lag_in ,z_out_nonuni(i)%data(mylk1(i):mylk2(i)),mylk0(i), &
!                        & lag_out,zval(lk1:lk2),lk0)

          data_xyz(i,j,lk1:lk2)= lag_out

#ifdef CHECKNAN
          do k=lk1,lk2
             if (isnan(real(data_xyz(i,j,k)))) then
             print*,'nan after interp',i,j,k,real(data_xyz(i,j,k))
             print*,'lag_in',lag_in
             print*,'lag_out', lag_out
             print*,'z_in',z_out_nonuni(i)%data
             print*,'z_out',zval(lk1:lk2)
             stop
             end if
          end do
#endif
       end do

       deallocate(lag_in,data_fs)

    end do

#ifdef WRITE_ALL
    call dump_interpolated_field_real(data_xyz)
#endif
    ! ifft in y
    do k=lk1,lk2
       call gto_fourier_n(data_xyz(:,:,k),data_imported(:,:,k),li0)
    end do

#ifdef CHECKNAN
    do i=li1,li2
       do j=lj1,lj2
          do k=lk1,lk2
          if (isnan(real(data_imported(i,j,k))).or.(isnan(aimag(data_imported(i,j,k))))) then
          print*,'nan after fft',i,j,k,real(data_imported(i,j,k)),aimag(data_imported(i,j,k))
          stop
          end if
        end do
      end do
    end do
#endif
    call stop_timer("MAP_FROM_PLANES")

  end subroutine planes_to_xyz_nonuni


  subroutine import_planes_nonuni
    real, dimension(:, :), allocatable :: data_block, data_fs, tmp
    real, dimension(:), allocatable :: first_plane
    real, dimension(1:4) :: tmp_inboard
    integer, dimension(0:n_procs_z-1) :: my_count, my_displ
    integer :: i, o, current_data, ierr, maxplane
    real :: lag_in(1), point(1)


    maxplane = n_cuts - 1
#ifdef COUPLE_XGC
    call cce_process_field
#endif
    allocate(data_block(0:maxplane, 0:block_count-1))
    allocate(first_plane(0:block_count-1))

#ifdef COUPLE_XGC
    call start_timer("RECEIVE_FIELD")
    ! There is a phase shift in the field of XGC and GENE, which is why we need the index trick
    call receive_field(data_block, block_start, block_end, block_count, n_cuts, my_mpi_comm_world)
    call stop_timer("RECEIVE_FIELD")
!   first_plane = data_block(0, :)
!   data_block(0:maxplane-1, :) = data_block(1:maxplane, :)
!   data_block(maxplane, :) = first_plane(:)
    data_block = data_block / norm_fact_field

#endif

#ifdef CHECKNAN
    do i=0, block_count-1
       do o=0, maxplane
          if (isnan(data_block(i, 1))) then
             print *, 'NaN from XGC', i, 'plane', o
             stop
          end if
       end do
    end do
#endif

    current_data = 0
    do i=li1, li2
         data_unst(i)%data(0:maxplane, :) = 0.0

        allocate(data_fs(0:maxplane, 0:z_out_nonuni(i)%num-1))
        allocate(tmp(0:maxplane, 0:mylk0(i)-1))

        my_count = counts_nonuni_z(:, i)
        my_displ = displ_nonuni_z(:, i)

        if (extra_zero) then
           ! get only the local fs from whole bock
           data_fs(0:maxplane, 1:z_out_nonuni(i)%num-1) = data_block(0:maxplane, &
                      i-li1+current_data:i-li1+current_data+z_out_nonuni(i)%num-2)

           ! reshuffle nodes
           data_fs(0:maxplane, 1:z_out_nonuni(i)%num-1) = ReshuffleNodes(&
           data_fs(0:maxplane, 1:z_out_nonuni(i)%num-1), &
           z_out_nonuni(i)%num-1, &
           z_out_nonuni(i)%GENE_starts, maxplane)

           !ridiculous workaround
           point = -pi
           do o=0, maxplane
              ! interpolate -pi
              tmp_inboard = (/&
                          data_fs(o, z_out_nonuni(i)%num-2), &
                          data_fs(o, z_out_nonuni(i)%num-1), &
                          data_fs(o, 2), data_fs(o, 3) /)

              ! I didn't know how to remove this from the loop, if it can be
              call lag3interp(tmp_inboard, XGC_inboard(:, i), 4, lag_in, point, 1)
              data_fs(o, 0) = lag_in(1)
              !if ((i.eq.80).and.(mype.eq.0)) write(*,*) tmp_inboard, XGC_inboard(:,i), lag_in, point
           end do

           ! scatter to all processors in z their
           call mpi_scatterv(data_fs(0:maxplane, 0), my_count*n_cuts, my_displ*n_cuts, MPI_REAL_TYPE, &
                             tmp(0, 0), mylk0(i)*n_cuts, MPI_REAL_TYPE, &
                             0, mpi_comm_z, ierr)

           data_unst(i)%data(0:maxplane, mylk1(i):mylk2(i)) = tmp
           current_data = current_data + z_out_nonuni(i)%num-2

        else
           ! get my part from the whole block
           data_fs = data_block(0:maxplane, i-li1+current_data:i-li1+current_data+z_out_nonuni(i)%num-1)

           ! reshuffle nodes
           data_fs = ReshuffleNodes(data_fs, z_out_nonuni(i)%num, z_out_nonuni(i)%GENE_starts, maxplane)

           ! passing data strucures
           call mpi_scatterv(data_fs(0, 0), my_count*n_cuts, my_displ*n_cuts, MPI_REAL_TYPE, &
                            tmp(0, 0), mylk0(i)*n_cuts, MPI_REAL_TYPE, 0, &
                            mpi_comm_z, ierr)

           data_unst(i)%data(0:maxplane, mylk1(i):mylk2(i)) = tmp
           current_data = current_data + z_out_nonuni(i)%num-1
        endif
        deallocate(data_fs,tmp)
     enddo


     deallocate(data_block, first_plane)
     isnap_planes_from_XGC = isnap_planes_from_XGC + 1

  end subroutine import_planes_nonuni


  subroutine import_GENE_density
 ! testing purpose only, read what written in existing file
    real,dimension(:,:), allocatable ::data_block
    real,dimension(:), allocatable :: data_fs, tmp
    real,dimension(1:4) :: tmp_inboard
    integer, dimension(0:n_procs_z-1) :: my_count,my_displ
    real :: lag_in(1), point(1)
    integer:: i,o, current_data, ierr

    allocate(data_block(0:block_count-1,1))

    do o=0,n_cuts-1
#ifdef COUPLE_XGC
       call receive_GENE_density(data_block,block_start,block_end,block_count,o,mpi_comm_x)
       data_block=data_block/norm_fact_mom

#endif
       current_data=0

       do i=li1, li2

          allocate(data_fs(0:z_out_nonuni(i)%num-1))
          allocate (tmp(0:mylk0(i)-1))

          my_count=counts_nonuni_z(:,i)
          my_displ=displ_nonuni_z(:,i)

          !reshuffle data
          if (extra_zero) then
              !get only the local fs from whole bock
              data_fs(1:z_out_nonuni(i)%num-1)=data_block(i-li1+current_data: &
                       & i-li1+current_data+z_out_nonuni(i)%num-2,1)

              !reshuffle nodes
              data_fs(1:z_out_nonuni(i)%num-1)=reshuffle_nodes(data_fs(1:z_out_nonuni(i)%num-1), &
                                           & z_out_nonuni(i)%num-1, &
                                           & z_out_nonuni(i)%GENE_starts)
              !interpolate -pi
              tmp_inboard=(/data_fs(z_out_nonuni(i)%num-2),  &
                            data_fs(z_out_nonuni(i)%num-1),  &
                            data_fs(2),                      &
                            data_fs(3)/)
              !ridiculous workaround
              point=-pi
              call lag3interp(tmp_inboard,XGC_inboard(:,i), 4, lag_in, point,1)
              data_fs(0)=lag_in(1)

              !scatter to all processors in z their data
              call mpi_scatterv(data_fs(0),my_count,my_displ, MPI_REAL_TYPE,&
                       & tmp(0), mylk0(i),&
                       & MPI_REAL_TYPE,0, mpi_comm_z, ierr)

              data_unst(i)%data(o,mylk1(i):mylk2(i))=tmp
              current_data=current_data+z_out_nonuni(i)%num-2

          else
              !get my part from the whole block
              data_fs=data_block(i-li1+current_data:i-li1+current_data+z_out_nonuni(i)%num-1,1)
              !reshuffle nodes
              data_fs=reshuffle_nodes(data_fs, z_out_nonuni(i)%num, &
                      & z_out_nonuni(i)%GENE_starts)

              !assing data strucures
              call mpi_scatterv(data_fs(0),my_count,my_displ, MPI_REAL_TYPE,&
                       & tmp(0), mylk0(i),&
                       & MPI_REAL_TYPE,0, mpi_comm_z, ierr)

              data_unst(i)%data(o,mylk1(i):mylk2(i))=tmp
              current_data=current_data+z_out_nonuni(i)%num-1
          endif

          deallocate(data_fs,tmp)
       enddo
    enddo

    deallocate(data_block)

  end subroutine import_GENE_density


  subroutine import_init_cond
    real, dimension(:, :), allocatable :: data_block, data_fs, tmp
    real, dimension(:), allocatable :: first_plane
    real, dimension(1:4) :: tmp_inboard
    integer, dimension(0:n_procs_z-1) :: my_count, my_displ
    integer :: i, o, current_data, ierr, maxplane
    real :: lag_in(1), point(1)
    real, dimension(32,322909)::data_in

    maxplane = n_cuts - 1

    allocate(data_block(0:maxplane, 0:block_count-1))
    allocate(first_plane(0:block_count-1))

#ifdef INIT_XGC
     call receive_density_3D_noGC(data_block, block_start, block_end, block_count, n_cuts, my_mpi_comm_world)
#endif
#ifdef H5_init
    CALL openf('init_cond_XGC.h5', hdf5_init)
    CALL getarr(hdf5_init, "/iden" ,   data_in)
    call closef(hdf5_init) 
    data_block =data_in(:,block_start:block_end)
#endif
    current_data = 0.0

    do i=li1, li2
       data_unst(i)%data(0:maxplane, :) = 0.0

       allocate(data_fs(0:maxplane, 0:z_out_nonuni(i)%num-1))
       allocate(tmp(0:maxplane, 0:mylk0(i)-1))

       my_count = counts_nonuni_z(:, i)*n_cuts
       my_displ = displ_nonuni_z(:, i)*n_cuts

       if (extra_zero) then
          ! get only the local fs from whole bock
          data_fs(0:maxplane, 1:z_out_nonuni(i)%num-1) = data_block(0:maxplane, &
               i-li1+current_data:i-li1+current_data+z_out_nonuni(i)%num-2)

          ! reshuffle nodes
          data_fs(0:maxplane, 1:z_out_nonuni(i)%num-1) = ReshuffleNodes(&
                           data_fs(0:maxplane, 1:z_out_nonuni(i)%num-1), &
                           z_out_nonuni(i)%num-1, &
                           z_out_nonuni(i)%GENE_starts, maxplane)

          ! ridiculous workaround
          point = -pi
          do o=0, maxplane

             ! interpolate -pi
             tmp_inboard = (/&
             data_fs(o, z_out_nonuni(i)%num-2), &
             data_fs(o, z_out_nonuni(i)%num-1), &
             data_fs(o, 2), data_fs(o, 3) /)

             ! I didn't know how to remove this from the loop, if it can be
             call lag3interp(tmp_inboard, XGC_inboard(:, i), 4, lag_in, point, 1)
             data_fs(o, 0) = lag_in(1)
             !if ((i.eq.91).and.(mype.eq.0)) write(*,*) tmp_inboard, XGC_inboard(:,i), lag_in, point

          end do

          ! scatter to all processors in z their
          call mpi_scatterv(data_fs(0:maxplane, 0), my_count, my_displ, MPI_REAL_TYPE, &
                            tmp(0, 0), mylk0(i)*n_cuts, MPI_REAL_TYPE, &
                            0, mpi_comm_z, ierr)

          data_unst(i)%data(0:maxplane, mylk1(i):mylk2(i)) = tmp / &
                      & (rhostar * minor_r) / &
                      & (spec(0)%dens_prof(i) * nref * 1.0e19)

          current_data = current_data + z_out_nonuni(i)%num-2

       else

          ! get my part from the whole block
          data_fs = data_block(0:maxplane, i-li1+current_data:i-li1+current_data+z_out_nonuni(i)%num-1)

          ! reshuffle nodes
          data_fs = ReshuffleNodes(data_fs, z_out_nonuni(i)%num, z_out_nonuni(i)%GENE_starts, maxplane)

          ! passing data strucures
          call mpi_scatterv(data_fs(0, 0), my_count, my_displ, MPI_REAL_TYPE, &
                            tmp(0, 0), mylk0(i)*n_cuts, MPI_REAL_TYPE, 0, &
                            mpi_comm_z, ierr)

          data_unst(i)%data(0:maxplane, mylk1(i):mylk2(i)) = tmp / &
                      & (rhostar * minor_r) / &
                      & (spec(0)%dens_prof(i) * nref * 1.0e19)

          current_data = current_data + z_out_nonuni(i)%num-1
       endif

       deallocate(data_fs,tmp)
    enddo


    deallocate(data_block, first_plane)

  end subroutine import_init_cond


  subroutine map_XGC_init_cond(data_imported)
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2), intent(out) :: data_imported

    real, dimension(:,:), allocatable :: data_fs
    real, dimension(li1:li2,mylj1:mylj2,lk1:lk2) :: data_xyz
    real, dimension(:), allocatable :: lag_in
     real, dimension(lk1:lk2) :: lag_out
    integer :: i,j,k
    real :: dz

    integer :: lk=0, ierr
#ifdef WITHFUTILS
    integer :: rank
    character(len=FILENAME_MAX) :: dset_name,fname
    integer :: fidinterp_h5
#endif
    call import_init_cond

    !poloidal grid, local chunk on XGC grid
    do i=li1,li2
       allocate(lag_in(mylk1(i)-nzb:mylk2(i)+nzb))
       allocate(data_fs(mylj1:mylj2,mylk1(i):mylk2(i)))

       do k=mylk1(i),mylk2(i)
           data_fs(:,k)=matmul(mat_from_plane(i)%data(:,:,k),data_unst(i)%data(:,k))
       end do

       do j=mylj1,mylj2
          lag_in(mylk1(i):mylk2(i))=data_fs(j,:)
          !exchange boundary for inner interpolation
          call exchange_z_1D_real_XGC(z_boundary_real_XGC(i),lag_in)

          !downsample in z
!          call linear_interp(z_out_nonuni(i)%data(mylk1(i):mylk2(i)),lag_in,mylk0(i), &
!                        & zval(lk1:lk2),lag_out,lk0)
          call lag3interp(lag_in ,z_out_nonuni(i)%data,mylk0(i)+2*nzb, &
                        & lag_out,zval(lk1:lk2),lk0)

!          call lag3interp(lag_in ,z_out_nonuni(i)%data(mylk1(i):mylk2(i)),mylk0(i), &
!                        & lag_out,zval(lk1:lk2),lk0)

          data_xyz(i,j,lk1:lk2)= lag_out
if ((j.eq.1).and.(i.eq.91).and.(mype.eq.0)) then
write(*,*) 'data',data_unst(i)%data(1,:)
write(*,*) 'lag in',lag_in
write(*,*) 'lag out',lag_out
write(*,*) 'z XGC',z_out_nonuni(i)%data
write(*,*) 'z',zval(lk1:lk2)
endif

       end do

       deallocate(lag_in,data_fs)

!if ((i.eq.90).and.((my_pey+my_pez+my_pew+my_pev).eq.0)) print *,'surface 91',data_xyz(i,1,:), my_pex!data_unst(i)%data(0,:),my_pex
    end do


    ! ifft in y
    do k=lk1,lk2
       call gto_fourier_n(data_xyz(:,:,k),data_imported(:,:,k),li0)
 !      do j=lj1,lj2
 !         data_imported(:,j,k)=data_imported(:,j,k)*geom%Bfield(:,pj1,k)
 !      end do
    end do

!    call xyz_to_planes_nonuni(data_imported)
if (((li1.le.90).and.(li2.ge.90)).and.((my_pey+my_pez+my_pew+my_pev).eq.0)) print *,'surface 91',data_imported(90,1,:), my_pex

#ifdef WITHFUTILS
    if ((my_pev+my_pew+my_pespec).eq.0) then
       rank = 0

       write(fname, "(A,A)") trim(diagdir)//'/init_cond', trim(file_extension)//'.h5'

       call creatf(fname, fidinterp_h5, "Initial condition", 'd', mpi_comm_xyz)
!       call creatg(fidinterp_h5, '/init_cond')
       call putarrnd(fidinterp_h5, '/init_cond', data_imported, (/3, 2, 1/))

       call closef(fidinterp_h5)
    end if
#endif

  end subroutine map_XGC_init_cond


  subroutine dump_interpolated_field(interpolated_field)
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2) :: interpolated_field
    ! testing purpose only
#ifdef WITHFUTILS
    integer :: rank
    character(len=FILENAME_MAX) :: dset_name,fname
    integer :: fidinterp_h5

    if ((my_pev+my_pew+my_pespec).eq.0) then
       write(fname, "(A, i10.10,A)") trim(diagdir)//'/field_interp_',istep_interp, trim(file_extension)//'.h5'

       call creatf(fname, fidinterp_h5, "field diagnostics", 'd', mpi_comm_xyz)
       call creatg(fidinterp_h5, '/field')
       rank = 0

       write(dset_name, "(A, '/', i10.10)") "/field/", 0

       call putarrnd(fidinterp_h5, dset_name, interpolated_field, (/3, 2, 1/))

       call closef(fidinterp_h5)

    end if
#endif

    istep_interp=istep_interp+1

  end subroutine dump_interpolated_field


	! Originally, each plane was written as a different file, with a loop over
	! the planes. This loop has been removed, and all planes write to a single
	! ADIOS file.
	!
	! GENE is working in a 3D system:
	!
	! (1)	mpi_comm_x has the domain split such that each rank has a
	!		different subset of the flux surfaces. i iterates over these flux
	!		surfaces
	!
	! (2) 	Going around each flux surface is the z-direction (each
	!		point going around is a different node from XGC).
	!
	! (3)	The different poloidal planes (o)

  subroutine write_planes(diagnostic,dset_name)
     logical, intent(in),optional :: diagnostic
     character(len=FILENAME_MAX), intent(in), optional :: dset_name

     real, dimension(:, :), allocatable :: data_fs, data_block, tmp
     integer, dimension(0:n_procs_z-1) :: my_count, my_displ
     integer :: o, i, l, ierr, current_data, maxplane

     real, dimension(:,:), allocatable :: tmp_h5


#ifdef COUPLE_XGC
     if (.not.present(diagnostic)) call cce_process_density
#endif

     maxplane = n_cuts - 1
     allocate(data_block(0:maxplane, block_start:block_end))
     current_data = block_start

     do i=li1, li2
       ! How the nodes at each flux surface need to be arranged
       my_count = counts_nonuni_z(:, i)*n_cuts
       my_displ = displ_nonuni_z(:, i)*n_cuts

       ! Collect over planes and over z simulataneously (z is mpi_gatherv, planes are just from array slice)
       allocate(tmp(0:maxplane, 0:mylk0(i)-1))
       allocate(data_fs(0:maxplane, z_out_nonuni(i)%num))
       tmp(0:maxplane, 0:mylk0(i)-1) = data_unst(i)%data(0:maxplane, mylk1(i):mylk2(i))

       call mpi_gatherv(&
                        tmp(0,0), mylk0(i)*n_cuts, MPI_REAL_TYPE, &
                        data_fs, my_count, my_displ, MPI_REAL_TYPE, &
                        0, mpi_comm_z, ierr)

       if (extra_zero) then
          ! Reshuffle data -- special handling when cross 0
          data_fs(0:maxplane, 2:z_out_nonuni(i)%num) = ReshuffleNodes(&
                   data_fs(0:maxplane, 2:z_out_nonuni(i)%num), &
                   z_out_nonuni(i)%num-1, &
                   z_out_nonuni(i)%num - z_out_nonuni(i)%GENE_starts+1, &
                   maxplane)

          data_block(0:maxplane, i-li1+current_data:i-li1+current_data+z_out_nonuni(i)%num-2) = &
                             data_fs(0:maxplane, 2:z_out_nonuni(i)%num)
          current_data = current_data + z_out_nonuni(i)%num-2
       else
          ! Reshuffle data -- normal scenarios
          data_fs(0:maxplane, :) = ReshuffleNodes(&
                              data_fs(0:maxplane, :), &
                              z_out_nonuni(i)%num, &
                              z_out_nonuni(i)%num - z_out_nonuni(i)%GENE_starts +2, &
                              maxplane)

          data_block(0:maxplane, i-li1+current_data:i-li1+current_data+z_out_nonuni(i)%num-1) = data_fs(0:maxplane, :)
          current_data = current_data+z_out_nonuni(i)%num-1
       endif

       deallocate(data_fs, tmp)

    end do

#ifdef COUPLE_XGC
    if ((my_pez+my_pev+my_pew+my_pespec).eq.0) then
       data_block = data_block * norm_fact_mom

!       if (my_pex .eq. 0) then
       do i = 1,10
       print *,  "The first 10 data_block content at ", i," is: ", data_block(i, 0)
       end do

       do i = 1,10
       print *,  "The last 10 data_block content at ", n_cuts-1 - 10 + i, " is :", data_block((n_cuts-1)-10+i, block_end)
       !print *,  "The last 10 data_block content at ", maxplane-10+i, " is :", data_block(maxplane-10+i, block_count-1)! last 10rowz
       end do
 !      endif

       if (.not.present(diagnostic)) call send_density(data_block, o, n_cuts, block_count, block_start, block_end, mpi_comm_x)
    endif
#endif

     if (present(diagnostic)) then
#ifdef WITHFUTILS
        if ((my_pey+my_pez+my_pev+my_pew+my_pespec).eq.0) then
           call putarrnd(fidplanes_h5, dset_name, data_block, (/2/))
           call attach(fidplanes_h5, dset_name, "time", time)
           call flushh5(fidplanes_h5)
        endif
#endif
     endif

     deallocate(data_block)

#ifdef COUPLE_XGC
     ! This creates the unloack file
     if (.not.present(diagnostic)) call write_check_file(mpi_comm_x)
#endif

  end subroutine write_planes


  subroutine linear_interp (x_in,y_in,n_in,x_out,y_out,n_out)
    integer, intent(IN) :: n_in
    integer, intent(in) :: n_out
    real, dimension(n_in), intent(in) :: x_in,y_in
    real, dimension(n_out), intent(in) :: x_out
    real, dimension(n_out), intent(out) :: y_out

    integer :: i,left,right
    real :: x_act

    do i=1,n_out
       x_act=x_out(i)

       call place_point(n_in, x_in,x_act,left,right)

       y_out(i)= ((x_in(right)-x_act)*y_in(left) + &
               &  (x_act-x_in(left))*y_in(right))/ &
               &  (x_in(right)-x_in(left))

     end do

  end subroutine linear_interp


  subroutine place_point(n_in,x_in,x,left,right)
    integer,intent(in) :: n_in
    real,dimension(n_in), intent(in) :: x_in
    real,intent(in) :: x
    integer,intent(out) :: left,right

    integer :: i

    do i=2,n_in-1
       if (x < x_in(i) ) then
          left=i-1
          right=i
          return
       end if
    end do

    left=n_in-1
    right=n_in

  end subroutine place_point


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine initialize_diag_debug
#ifdef WITHFUTILS
    integer :: rank
    character(len=FILENAME_MAX) :: fname

    if ((my_pev+my_pew+my_pespec).eq.0) then
       !original GENE field
       write(fname, "(A,A)") trim(diagdir)//'/field_original',trim(file_extension)//'.h5'
       call creatf(fname, fidorig_h5, "field diagnostics", 'd', mpi_comm_xyz)
       call creatg(fidorig_h5, '/field')
       call creatg(fidorig_h5, '/field/phi')
       rank = 0
       istep_orig=0
       !original GENE field wb
       write(fname, "(A,A)") trim(diagdir)//'/field_original_wb',trim(file_extension)//'.h5'
       call creatf(fname, fidorig_wb_h5, "field diagnostics", 'd', mpi_comm_xyz)
       call creatg(fidorig_wb_h5, '/field')
       call creatg(fidorig_wb_h5, '/field/phi')
       rank = 0
       istep_orig_wb=0
       !original GENE real field wb
       write(fname, "(A,A)") trim(diagdir)//'/field_original_real_wb',trim(file_extension)//'.h5'
       call creatf(fname, fidorig_real_wb_h5, "field diagnostics", 'd', mpi_comm_xyz)
       call creatg(fidorig_real_wb_h5, '/field')
       call creatg(fidorig_real_wb_h5, '/field/phi')
       rank = 0
       istep_orig_real_wb=0
       !interpolated field real
       write(fname, "(A,A)") trim(diagdir)//'/field_interp_real', trim(file_extension)//'.h5'
       call creatf(fname, fidinterp_real_h5, "field diagnostics", 'd', mpi_comm_xyz)
       call creatg(fidinterp_real_h5, '/field')
       call creatg(fidinterp_real_h5, '/field/phi')
       rank = 0
       istep_interp_real=0

    endif

#endif
  end subroutine initialize_diag_debug


  subroutine dump_original_field(field)
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2) :: field
    ! testing purpose only
#ifdef WITHFUTILS
    character(len=FILENAME_MAX) :: dset_name
    if ((my_pev+my_pew+my_pespec).eq.0) then

      write(dset_name, "(A, '/', i10.10)") "/field/phi", istep_orig

      call putarrnd(fidorig_h5, dset_name, field, (/3, 2, 1/))

      istep_orig=istep_orig+1

      call flushh5(fidorig_h5)

    endif

#endif
  end subroutine dump_original_field


  subroutine dump_original_field_wb(field)
    complex, dimension(li1:li2,lj1:lj2,lbz:ubz) :: field
    ! testing purpose only
#ifdef WITHFUTILS
    character(len=FILENAME_MAX) :: dset_name
    if ((my_pev+my_pew+my_pespec).eq.0) then

      write(dset_name, "(A, '/', i10.10)") "/field/phi", istep_orig_wb

      call putarrnd(fidorig_wb_h5, dset_name, field, (/3, 2, 1/))

      istep_orig_wb=istep_orig_wb+1

      call flushh5(fidorig_wb_h5)

    endif

#endif

  end subroutine dump_original_field_wb


  subroutine dump_original_field_real(field)
    real, dimension(li1:li2,mylj1:mylj2,lbz:ubz) :: field
    ! testing purpose only
#ifdef WITHFUTILS
    character(len=FILENAME_MAX) :: dset_name

    if ((my_pev+my_pew+my_pespec).eq.0) then

       write(dset_name, "(A, '/', i10.10)") "/field/phi", istep_orig_real_wb

       call putarrnd(fidorig_real_wb_h5, dset_name, field, (/3, 2, 1/))

       istep_orig_real_wb=istep_orig_real_wb+1

       call flushh5(fidorig_real_wb_h5)
    endif
#endif

  end subroutine dump_original_field_real


  subroutine dump_interpolated_field_real(field)
    real, dimension(li1:li2,mylj1:mylj2,lk1:lk2) ::field
    ! testing purpose only
#ifdef WITHFUTILS
    character(len=FILENAME_MAX) :: dset_name

    if ((my_pev+my_pew+my_pespec).eq.0) then

       write(dset_name, "(A, '/', i10.10)") "/field/phi", istep_interp_real

       call putarrnd(fidinterp_real_h5, dset_name, field, (/3, 2, 1/))

       istep_interp_real=istep_interp_real+1

       call flushh5(fidinterp_real_h5)
    end if
#endif
  end subroutine dump_interpolated_field_real


  subroutine finalize_diag_debug

    if ((my_pev+my_pew+my_pespec).eq.0) then

       call closef(fidorig_h5)
       call closef(fidorig_wb_h5)
       call closef(fidorig_real_wb_h5)
       call closef(fidinterp_real_h5)

    end if

  end subroutine finalize_diag_debug


end module poloidal_planes_nonuni
