#include "intrinsic_sizes.h"
#include "redef.h"
#undef WRITE_ALL
module poloidal_planes
  use mpi
  use par_other, only: pi
  use geometry, only: rhostar, minor_r, Lref
  use communications, only: get_systime, mpi_comm_x, mpi_comm_z
  use discretization, only: mype, lz0, nz0, nx0,lx0, nzb, nib, lbi,ubi,lbz,ubz, &
                          & li0,li1,li2,lj1,lj2,lk1,lk2, ln1,ln2, &
                          & pi1gl, pi2gl,&
                          & my_pex,my_pey,my_pez,my_pev,my_pew,my_pespec
  use coordinates, only: zval, xval_a, deli, dz
  use par_poloidal_planes, only: nz_out,n_cuts,istep_planes,phi_cut, res_fact,XGC_file, &
                         & x_GENE_wb, z_GENE_wb, weight,norm_fact_field, norm_fact_mom,active_nodes, &
                         & sign_phi, istep_field_in, cref, norm_fact_init,&
                         & block_start, block_count
  use create_nonuni_grid, only: generate_grid_nonuni, finalize_zinds_nonuni
  use poloidal_planes_aux, only: define_cuts
  use poloidal_planes_nonuni, only: prepare_matrixes_nonuni, &
                             & xyz_to_planes_nonuni, write_planes
  use par_in, only: file_extension, diagdir, n_pol, dt_max, spec
  use import_XGC_grid, only: read_XGC_grid, finalize_import_XGC_grid
  use profiles_mod, only: Tref,mref,nref
  use file_io, only: get_unit_nr
  use par_mod, only: charge_electron, m_proton
#ifdef WITHFUTILS
  use futils
  use par_poloidal_planes, only: isnap_planes, fidplanes_h5, isnap_planes_from_XGC, fidplanes_h5_from_XGC
#endif
#ifdef COUPLE_XGC
  use coupling_base_setup, only: cce_first_surface,cce_dt,cce_node_number 
  use coupling_core_gene, only: cce_varpi_grid, initialize_coupling_engines
#endif

  implicit none

  public :: initialize_planes, finalize_planes, set_poloidal_planes_defaults,&
         &  diag_planes

  private

contains

  subroutine set_poloidal_planes_defaults
    nz_out=400
    n_cuts=1
    istep_planes=0
    res_fact=1
    sign_phi=-1.0
    istep_field_in=0
  end subroutine set_poloidal_planes_defaults


  subroutine initialize_planes
    real :: plane_start, plane_end
    integer :: k,n

    call get_systime(plane_start)

    norm_fact_field = Tref * 1.0e3* rhostar * minor_r
    norm_fact_mom = nref * 1.0e19 * rhostar* minor_r

    allocate(norm_fact_init(pi1gl:pi2gl,ln1:ln2))

    do n=ln1,ln2
       norm_fact_init(:,n) = spec(n)%dens*spec(n)%dens_prof * nref*1e19 /sqrt(pi * &
               & 2 * spec(n)%temp*spec(n)%temp_prof * Tref * charge_electron * 1e3 / &
               & ( mref * m_proton ))
    end do

    cref=sqrt(Tref*1e3*charge_electron/(mref*m_proton))
#ifdef COUPLE_XGC
    dt_max=cce_dt/(Lref/cref)
    if (mype.eq.0) write(*,"(A,G17.8)")"Override timestep to ", dt_max
#endif
    call define_cuts

    call set_GENE_aux_grids

    if ((XGC_file == " " ).or. (nz_out.gt.0)) then
       If (mype.eq.0) write(*,'(A)')"Creating a poloidal grid"
       call generate_grid_nonuni
    else
       call read_XGC_grid
    endif

    If (mype.eq.0) write(*,'(A)')"Preparing weights..."
    call define_weight

    If (mype.eq.0) write(*,'(A)')"Preparing interpolation matrixes..."
    call prepare_matrixes_nonuni

    call initialize_diag_planes
#ifdef WRITE_ALL
    call initialize_rewrite_XGC_planes
#endif
    call get_systime(plane_end)

#ifdef ADIOS2
    If (mype.eq.0) write(*,'(A)')"Preparing ADIOS engines..."
    call initialize_coupling_engines(cce_node_number, n_cuts, block_start, block_count)
#endif

    If (mype.eq.0) write(*,'(A,F10.3,a)')&
       &'Time for initializing poloidal planes:',plane_end-plane_start,' sec'

  end subroutine initialize_planes


  subroutine finalize_planes

    deallocate(phi_cut)
    deallocate(z_GENE_wb,x_GENE_wb)
    if (allocated(weight)) deallocate(weight)

    call finalize_import_XGC_grid

    call finalize_diag_planes
#ifdef WRITE_ALL
    call finalize_rewrite_XGC_planes
#endif
  end subroutine finalize_planes


  subroutine set_GENE_aux_grids
    integer :: i,k

    !parallel coordinate wb
    allocate(z_GENE_wb(lbz:ubz))
    z_GENE_wb(lk1:lk2)=zval(lk1:lk2)
    do k=1,nzb
       z_GENE_wb(lk2+k)=zval(lk2)+dz*k
       z_GENE_wb(lk1-k)=zval(lk1)-dz*k
    end do

    !radial GENE grid with boundaries
    allocate(x_GENE_wb(lbi:ubi))
    x_GENE_wb(li1:li2)=xval_a(li1:li2)
    do i=1,nib
       x_GENE_wb(li2+i)=xval_a(li2)+deli* rhostar * i
       x_GENE_wb(li1-i)=xval_a(li1)-deli* rhostar * i
    end do

  end subroutine set_GENE_aux_grids


  subroutine define_weight
    integer :: i
#ifdef COUPLE_XGC
    allocate(weight(li1:li2))

    do i=li1,li2
       weight(i)=cce_varpi_grid(i+cce_first_surface)
    end do

    call write_weight
#endif
  end subroutine define_weight


  subroutine  write_weight
    CHARACTER(len=512) :: tmpstr
    character(len=128) :: str_fmt
    integer :: ierr, par_coupling
    real, dimension(0:nx0-1) :: tmp
#ifdef WITHFUTILS
    integer, dimension(2) :: dims
    integer :: rank=0
    integer :: fid_h5
    real, dimension(li1:li2,1) :: tmp_h5

    if ((my_pez+my_pev+my_pew+my_pespec).eq.0) then

       call creatf(trim(diagdir)//'/coupling_weight'//trim(file_extension)//'.h5', &
               fid_h5, "Coupling weights", 'd', mpi_comm_x)

       call creatg(fid_h5, '/coupling')

!       call creatd(fid_h5, rank, dims, '/coupling/weight')
       tmp_h5(li1:li2,1)=weight(li1:li2)
       call putarrnd(fid_h5, '/coupling/weight',tmp_h5, (/1/))

!       call creatd(fid_h5, rank, dims,'/coupling/x')
       tmp_h5(li1:li2,1)=xval_a(li1:li2)
       call putarrnd(fid_h5, '/coupling/x_o_a',tmp_h5, (/1/))

!       call putarrnd(fid_h5, '/coupling/norm_init',norm_fact_init, (/2, 1/))
#ifdef COUPLE_XGC
       call attach(fid_h5, "/coupling", "first_fs", cce_first_surface )
#endif
       call attach(fid_h5, "/coupling", "nodes", active_nodes)
       call attach(fid_h5, "/coupling", "norm_fact_mom", norm_fact_mom )
       call attach(fid_h5, "/coupling", "norm_fact_field", norm_fact_field)

       call closef(fid_h5)

    endif
#endif

  end subroutine write_weight


  subroutine arr2str(arr, str)
    real,dimension(0:), intent(in) :: arr
    character(len=*), intent(out) :: str
    integer :: i

    Write(str,"(G16.8)") arr(0)
    do i=1,size(arr)
       Write(str,"(2A,G16.8)") TRIM(str),',',arr(i)
    enddo

  end subroutine arr2str


  !-----------------------------------------------------------------------------------
  ! the data can be dumped on  a file. H5 only for the moment
  ! can introduce a xz communicator to do parallel IO
  subroutine  initialize_diag_planes
#ifdef WITHFUTILS
     integer :: rank,o,ierr
     integer, dimension(2) :: dims
     character(len=FILENAME_MAX) :: dset_name

     isnap_planes=0

     if ((my_pez+my_pev+my_pew+my_pespec).eq.0) then

        call creatf(trim(diagdir)//'/planes'//trim(file_extension)//'.h5', &
               fidplanes_h5, "Poloidal_plane", 'd', mpi_comm_x)
        call creatg(fidplanes_h5, '/planes')
             rank = 0
        call creatd(fidplanes_h5, rank, dims, "/planes/time", "time")

        call creatg(fidplanes_h5, '/planes/n')
        call creatg(fidplanes_h5, '/planes/n_GC')
        call creatg(fidplanes_h5, '/planes/T')
        call creatg(fidplanes_h5, '/planes/phi')

        call flushh5(fidplanes_h5)

     endif
#endif
   end subroutine initialize_diag_planes


   subroutine diag_planes(data_in,dset_name)
#ifdef WITHFUTILS
     complex,dimension(li1:li2,lj1:lj2,lk1:lk2), intent(in) :: data_in
     character(len=FILENAME_MAX) :: dset_name

     call xyz_to_planes_nonuni(data_in)

     call write_planes(.true.,dset_name)
#endif
   end subroutine diag_planes


   subroutine finalize_diag_planes
#ifdef WITHFUTILS
     if ((my_pez+my_pev+my_pew).eq.0) then
        call closef(fidplanes_h5)
     endif
#endif
   end subroutine finalize_diag_planes

  subroutine  initialize_rewrite_XGC_planes
#ifdef WITHFUTILS
     integer :: rank,o,ierr
     integer, dimension(2) :: dims
     character(len=FILENAME_MAX) :: dset_name

     isnap_planes_from_XGC=0

     if ((my_pez+my_pev+my_pew+my_pespec).eq.0) then

        call creatf(trim(diagdir)//'/planes_from_XGC'//trim(file_extension)//'.h5', &
               fidplanes_h5_from_XGC, "Poloidal_plane", 'd', mpi_comm_x)
        call creatg(fidplanes_h5_from_XGC, '/planes')
             rank = 0
        call creatd(fidplanes_h5_from_XGC, rank, dims, "/planes/time", "time")

        do o = 0, n_cuts-1
           write(dset_name,"(A, '/', i10.10)") "/planes",o+1
           call creatg(fidplanes_h5_from_XGC, trim(dset_name))
           call attach(fidplanes_h5_from_XGC, trim(dset_name), "position", phi_cut(o))
        end do

        call flushh5(fidplanes_h5_from_XGC)

     endif
#endif
   end subroutine initialize_rewrite_XGC_planes


   subroutine finalize_rewrite_XGC_planes
#ifdef WITHFUTILS
     if ((my_pez+my_pev+my_pew).eq.0) then
        call closef(fidplanes_h5_from_XGC)
     endif
#endif
   end subroutine finalize_rewrite_XGC_planes


end module poloidal_planes
