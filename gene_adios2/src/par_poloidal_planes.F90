#include "intrinsic_sizes.h"

module par_poloidal_planes
  use discretization, only: lbi, ubi, lbz, ubz

  !data on a flux surface (size(data)=ntheta[i])
  type data_fs
    real, allocatable :: data(:)
    integer :: num
    integer :: GENE_starts
  end type data_fs

  !data on a flux surface (size(data)=ncuts,2nky0,ntheta[i])
  type data_plane!
    real, allocatable :: data(:,:,:)
  end type data_plane

  type data_2d
    real, allocatable :: data(:,:)
  end type data_2d

  type data_fs_int
    integer, allocatable :: data(:)
    integer :: num
    integer :: GENE_starts
  end type data_fs_int

  integer, dimension(:), allocatable :: nz_out_nonuni !number of points per flux-surface
  type(data_fs), allocatable :: R(:)
  type(data_fs), allocatable :: Z(:)
  type(data_fs), allocatable :: z_out_nonuni(:)

  type(data_2d), allocatable :: data_3d_nonuni(:)
  type(data_2d), allocatable :: data_unst(:)
  type(data_plane), allocatable :: mat_to_plane(:)
  type(data_plane), allocatable :: mat_from_plane(:)

  integer, dimension(:,:),allocatable :: y_BC
  real, dimension(:), allocatable :: phi_cut
  real, dimension(:), allocatable :: q_XGC
  real, dimension(:), allocatable :: x_XGC
  real, dimension(:), allocatable :: x_GENE_wb
  real, dimension(:), allocatable :: z_GENE_wb
  real, dimension(:), allocatable :: weight


  integer :: mylj1, mylj2, mylj0 !local indexes y real space, no dealiasing
  integer, dimension(:), allocatable :: mylk1, mylk2, mylk0 !local indexes poloidal points
!  integer :: myli1, myli2, myli0 !local indexes radial points
  integer :: block_start,block_end,block_count
  real :: norm_fact_field, norm_fact_mom, cref
  real, dimension(:,:), allocatable :: norm_fact_init
  integer :: active_nodes
  real,dimension(:,:), allocatable :: XGC_inboard

  !namelist
  integer :: nz_out
  integer :: n_cuts
  integer :: istep_planes
  integer :: istep_field_in
  integer :: istep_mom_in
  logical :: BC_in_fourier=.true.
  logical :: pre_proc
  integer :: res_fact
  real    :: sign_phi
  character(len=FILENAME_MAX) :: XGC_file=''
  logical :: extra_zero=.false.


!!this shit is for using hdf5 futils in the dumbest way
  integer, dimension(:), allocatable :: counts_nonuni!how many x planes are to be sent
  integer, dimension(:), allocatable :: displ_nonuni  !wherefinalize_zin to place them
  integer, dimension(:,:), allocatable :: counts_nonuni_z!how many x planes are to be sent
  integer, dimension(:,:), allocatable :: displ_nonuni_z  !wherefinalize_zin to place them
  integer :: pts_nonuni        !number of points per x processor
  integer :: tot_pts_nonuni    !total number of points
  integer :: fidplanes_h5, isnap_planes
  integer :: fidplanes_h5_from_XGC, isnap_planes_from_XGC
  integer :: fidfield_in_h5,isnap_field_in=0, isnap_mom_in=0

contains


end module par_poloidal_planes
