module bicub_mod_gpu
  use precision_mod_gpu
  use cudafor
  implicit none

  integer, parameter :: ndeg = 3
  integer, parameter :: max_nr = eq_mr-1
  integer, parameter :: max_nz = eq_mz-1

  integer :: nr_gpu, nz_gpu

  real (kind=work_p) :: dr_inv_gpu, dz_inv_gpu, rmin_gpu, zmin_gpu

  ! Using texture memory here provides such a speedup, and they
  ! are small enough to always fit, that the only reason to turn
  ! them off would be debugging.
#ifndef NO_BICUB_TEXTURE
  real (kind=work_p),texture,pointer :: rc_cub_gpu(:)
  real (kind=work_p),texture,pointer :: zc_cub_gpu(:)
  real (kind=work_p),texture,pointer :: acoef_cub_gpu(:,:,:,:)

  real (kind=work_p),target,allocatable :: rc_cub_gpu_(:)
  real (kind=work_p),target,allocatable :: zc_cub_gpu_(:)
  real (kind=work_p),target,allocatable :: acoef_cub_gpu_(:,:,:,:)

#define RC_CUB_GPU rc_cub_gpu_
#define ZC_CUB_GPU zc_cub_gpu_
#define ACOEF_CUB_GPU acoef_cub_gpu_
  
#else
  real (kind=work_p),allocatable :: rc_cub_gpu(:)
  real (kind=work_p),allocatable :: zc_cub_gpu(:)
  real (kind=work_p),allocatable :: acoef_cub_gpu(:,:,:,:)

#define RC_CUB_GPU rc_cub_gpu
#define ZC_CUB_GPU zc_cub_gpu
#define ACOEF_CUB_GPU acoef_cub_gpu
#endif

!  real (kind=work_p), dimension(max_nr) :: rc_cub_gpu
!  real (kind=work_p), dimension(max_nz) :: zc_cub_gpu
!  real (kind=work_p), dimension(0:ndeg,0:ndeg,max_nr,max_nz) :: acoef_cub_gpu

#ifndef NO_BICUB_TEXTURE
!  ------------------------------------------------
!  note  attributes(texture) may not work correctly
!  ------------------------------------------------
  attributes(device) :: rc_cub_gpu_, zc_cub_gpu_, acoef_cub_gpu_
#else
  attributes(device) :: rc_cub_gpu, zc_cub_gpu, acoef_cub_gpu
#endif
  ! GPU constant memory
  attributes(constant) :: nr_gpu, nz_gpu
  attributes(constant) :: dr_inv_gpu, dz_inv_gpu, rmin_gpu, zmin_gpu

  public :: ndeg
  public :: update_device_bicub
  public :: acoef_cub_gpu
  public :: dr_inv_gpu, dz_inv_gpu, rmin_gpu, zmin_gpu
! Texture pointers only behave as textures within the same module,
! outside they'll be regular pointers. Private helps avoid confusion.
  private :: rc_cub_gpu, zc_cub_gpu, nr_gpu, nz_gpu

contains

  attributes(host) &
  subroutine update_device_bicub()
  use sml_module, only : sml_mype
  use precision_mod_gpu
  use bicub_mod, only : psi_bicub
#ifdef USE_C_BICUB
  use ex_bicub_mod, only : init_external_bicub
#endif
  use cudafor
  
implicit none

!  type(bicub_type), intent(in) :: bicub
  integer :: icount
  integer :: ierr  
  integer :: lb1,ub1, lb2,ub2, lb3,ub3, lb4,ub4

  integer, parameter :: idebug = 0
#if USE_GPU_EMU
  logical, parameter :: use_cudaMemcpy = .false.
#else
  logical, parameter :: use_cudaMemcpy = .true.
#endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!! copy data from CPU to GPU !!!
#ifdef USE_C_BICUB
  ! The C routine uses 2D cudaArrays and textures, which are slightly (5% or so) faster. 
  call init_external_bicub(psi_bicub%rc_cub, psi_bicub%zc_cub, psi_bicub%acoef_cub, &
        size(psi_bicub%acoef_cub,1), size(psi_bicub%rc_cub,1), size(psi_bicub%zc_cub,1), &
        psi_bicub%rmin, psi_bicub%zmin, psi_bicub%dr_inv, psi_bicub%dz_inv)
#else

  if (.not.allocated(RC_CUB_GPU)) then
    lb1 = lbound(psi_bicub%rc_cub,1)
    ub1 = ubound(psi_bicub%rc_cub,1)
    allocate( RC_CUB_GPU(lb1:ub1), stat=ierr )
    call assert(ierr.eq.0,'alloc(rc_cub_gpu)',ierr)

    RC_CUB_GPU = 0
  endif

  if (.not.allocated(ZC_CUB_GPU)) then
    lb1 = lbound(psi_bicub%zc_cub,1)
    ub1 = ubound(psi_bicub%zc_cub,1)
    allocate( ZC_CUB_GPU(lb1:ub1),stat=ierr )
    call assert(ierr.eq.0,'alloc(zc_cub_gpu)',ierr)

    ZC_CUB_GPU = 0
  endif

  if (.not.allocated(ACOEF_CUB_GPU)) then
    lb1 = lbound(psi_bicub%acoef_cub,1)
    lb2 = lbound(psi_bicub%acoef_cub,2)
    lb3 = lbound(psi_bicub%acoef_cub,3)
    lb4 = lbound(psi_bicub%acoef_cub,4)

    ub1 = ubound(psi_bicub%acoef_cub,1)
    ub2 = ubound(psi_bicub%acoef_cub,2)
    ub3 = ubound(psi_bicub%acoef_cub,3)
    ub4 = ubound(psi_bicub%acoef_cub,4)

    allocate( ACOEF_CUB_GPU(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4),stat=ierr )
    call assert(ierr.eq.0,'alloc(acoef_cub_gpu)',ierr)

    ACOEF_CUB_GPU = 0
  endif
#endif

  
  dr_inv_gpu = psi_bicub%dr_inv
  dz_inv_gpu = psi_bicub%dz_inv
  rmin_gpu = psi_bicub%rmin
  zmin_gpu = psi_bicub%zmin   
  nr_gpu=psi_bicub%nr
  nz_gpu=psi_bicub%nz

#ifndef USE_C_BICUB
 if(use_cudaMemcpy) then


     icount = size(RC_CUB_GPU)
     ierr = cudaMemcpy( RC_CUB_GPU, psi_bicub%rc_cub,                    &
                             icount, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(rc_cub_gpu)',ierr)

     icount = size(ZC_CUB_GPU)
     ierr = cudaMemcpy( ZC_CUB_GPU, psi_bicub%zc_cub,                    &
                             icount, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(zc_cub_gpu)',ierr)


     icount = size(ACOEF_CUB_GPU)
     ierr = cudaMemcpy( ACOEF_CUB_GPU, psi_bicub%acoef_cub,               &
                             icount, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(acoef_cub_gpu)',ierr)

 else

     RC_CUB_GPU = psi_bicub%rc_cub

     ZC_CUB_GPU = psi_bicub%zc_cub

     ACOEF_CUB_GPU = psi_bicub%acoef_cub
     
 endif

#ifndef NO_BICUB_TEXTURE
    rc_cub_gpu => rc_cub_gpu_
    zc_cub_gpu => zc_cub_gpu_
    acoef_cub_gpu => acoef_cub_gpu_
#endif
#endif
    return
  end subroutine update_device_bicub

!cdir$r unroll 
#ifdef USE_GPU
  attributes(device) &
#endif
  subroutine bicub_interpol2(r,z, f00,f10,f01,f11,f20,f02)
    use precision_mod_gpu
#ifdef USE_C_BICUB
    use ex_bicub_mod, only : eval_2_external
#endif
    implicit none
    real (kind=work_p), intent(in) :: r,z
    real (kind=work_p), intent(out) :: f00,f10,f01,f11,f20,f02
#ifdef USE_C_BICUB
    call eval_2_external(r,z,f00,f10,f01,f11,f20,f02)
#else
    real(kind=work_p) :: fx_i,dfx_i,dfx2_i, A(0:ndeg,0:ndeg)
    real (kind=work_p) :: dx,dy
    integer :: nr,nz, i,j
    !real (8) :: psi_acoef_local(0:ndeg,0:ndeg)
    ! ---------------
    ! perform hashing
    ! ---------------
    nr = nr_gpu
    nz = nz_gpu
    i = max(1,min(nr, 1 + int( (r-rmin_gpu)*dr_inv_gpu ) ))
    j = max(1,min(nz, 1 + int( (z-zmin_gpu)*dz_inv_gpu ) ))
    A(:,:) = acoef_cub_gpu(:,:,i,j)
!     real (kind=work_p), dimension(0:ndeg) :: xv,yv
!     real (kind=work_p), dimension(0:ndeg) :: fx 
!     real (kind=work_p), dimension(0:ndeg) :: dfx
!     real (kind=work_p), dimension(0:ndeg) :: dfx2




    dx = (r-rc_cub_gpu(i))
    dy = (z-zc_cub_gpu(j))


    f00 = 0
    f01 = 0
    f02 = 0

!        fx_0 = A(0,0) + dx*A(1,0) + (dx*dx)*A(2,0) + ((dx*dx)*dx)*A(3,0)
 
    fx_i = ((A(3,0)*dx + A(2,0))*dx + A(1,0))*dx + A(0,0)
    f00 = f00 + fx_i


!        fx_1 = A(0,1) + dx*A(1,1) + (dx*dx)*A(2,1) + ((dx*dx)*dx)*A(3,1)

    fx_i = ((A(3,1)*dx + A(2,1))*dx + A(1,1))*dx + A(0,1)
    f00 = f00 + dy*fx_i
    f01 = f01 +    fx_i

!        fx_2 = A(0,2) + dx*A(1,2) + (dx*dx)*A(2,2) + ((dx*dx)*dx)*A(3,2)

    fx_i = ((A(3,2)*dx + A(2,2))*dx + A(1,2))*dx + A(0,2)
    f00 = f00 + dy*(dy*fx_i)
    f01 = f01 + 2.0d0*(dy*fx_i)
    f02 = f02 + 2.0d0*fx_i
        

!        fx_3 = A(0,3) + dx*A(1,3) + (dx*dx)*A(2,3) + ((dx*dx)*dx)*A(3,3)
    fx_i = ((A(3,3)*dx + A(2,3))*dx + A(1,3))*dx + A(0,3)
    f00 = f00 + dy*(dy*(dy*fx_i))
    f01 = f01 + 3.0d0*(dy*(dy*fx_i))
    f02 = f02 + 6.0d0*(dy*fx_i)


    f10 = 0
    f11 = 0


!        dfx_0 = A(1,0) + 2*dx*A(2,0) + 3*(dx*dx)*A(3,0)

    dfx_i = (A(3,0)*3.0d0*dx + A(2,0)*2.0d0)*dx + A(1,0)
    f10 = f10 + dfx_i

!        dfx_1 = A(1,1) + 2*dx*A(2,1) + 3*(dx*dx)*A(3,1)

    dfx_i = (A(3,1)*3.0d0*dx + A(2,1)*2.0d0)*dx + A(1,1)
    f10 = f10 + dy*dfx_i
    f11 = f11 +    dfx_i


!        dfx_2 = A(1,2) + 2*dx*A(2,2) + 3*(dx*dx)*A(3,2)

    dfx_i = (A(3,2)*3.0d0*dx + A(2,2)*2.0d0)*dx + A(1,2)
    f10 = f10 + dy*(dy*dfx_i)
    f11 = f11 + 2.0d0*(dy*dfx_i)

!        dfx_3 = A(1,3) + 2*dx*A(2,3) + 3*(dx*dx)*A(3,3)

    dfx_i = (A(3,3)*3.0d0*dx + A(2,3)*2.0d0)*dx + A(1,3)
    f10 = f10 + dy*(dy*dy*dfx_i)
    f11 = f11 + 3.0d0*(dy*dy*dfx_i)


    f20 = 0

    dfx2_i = (3*2*A(3,0)*dx + 2*1*A(2,0))
    f20 = f20 + dfx2_i

    dfx2_i = (3*2*A(3,1)*dx + 2*1*A(2,1))
    f20 = f20 + dy*dfx2_i

    dfx2_i = (3*2*A(3,2)*dx + 2*1*A(2,2))
    f20 = f20 + (dy*dy)*dfx2_i

    dfx2_i = (3*2*A(3,3)*dx + 2*1*A(2,3))
    f20 = f20 + (dy*dy)*dy*dfx2_i
#endif
    return
  end subroutine bicub_interpol2


!cdir$r unroll 
#ifdef USE_GPU
   attributes(device) &
#endif
  subroutine bicub_interpol1(r,z, f00,f10,f01)
    use precision_mod_gpu
#ifdef USE_C_BICUB
    use ex_bicub_mod, only : eval_1_external
#endif
    implicit none
    real (kind=work_p), intent(in) :: r,z
    real (kind=work_p), intent(out) :: f00,f10,f01
#ifdef USE_C_BICUB
    call eval_1_external(r, z, f00, f10, f01)
#else
    real (kind=work_p) :: dx,dy
    real (kind=work_p) :: dfx_i
    real (kind=work_p) :: fx_i, A(0:ndeg,0:ndeg)
    integer :: nr,nz, i,j
    !real (8) :: psi_acoef_local(0:ndeg,0:ndeg)
    ! ---------------
    ! perform hashing
    ! ---------------
    nr = nr_gpu
    nz = nz_gpu
    i = max(1,min(nr, 1 + int( (r-rmin_gpu)*dr_inv_gpu ) ))
    j = max(1,min(nz, 1 + int( (z-zmin_gpu)*dz_inv_gpu ) ))
    A(:,:) = acoef_cub_gpu(:,:,i,j)
 
!    f00 = 0.0d0
!    f01 = 0.0d0
!    f10 = 0.0d0
    dx = (r-rc_cub_gpu(i))
    dy = (z-zc_cub_gpu(j))


    f00 = 0
    f01 = 0
!       fx_0 = A(0,0) + dx*A(1,0) + (dx*dx)*A(2,0) + ((dx*dx)*dx)*A(3,0)
    fx_i = ((A(3,0)*dx + A(2,0))*dx + A(1,0))*dx + A(0,0)
    f00 = f00 + fx_i

!       fx_1 = A(0,1) + dx*A(1,1) + (dx*dx)*A(2,1) + ((dx*dx)*dx)*A(3,1)
    fx_i = ((A(3,1)*dx + A(2,1))*dx + A(1,1))*dx + A(0,1)
    f00 = f00 + dy*fx_i
    f01 = f01 + fx_i


!       fx_2 = A(0,2) + dx*A(1,2) + (dx*dx)*A(2,2) + ((dx*dx)*dx)*A(3,2)
    fx_i = ((A(3,2)*dx + A(2,2))*dx + A(1,2))*dx + A(0,2)
    f00 = f00 + (dy*dy)*fx_i
    f01 = f01 + 2.0d0*dy*fx_i


!       fx_3 = A(0,3) + dx*A(1,3) + (dx*dx)*A(2,3) + ((dx*dx)*dx)*A(3,3)
    fx_i = ((A(3,3)*dx + A(2,3))*dx + A(1,3))*dx + A(0,3)
    f00 = f00 + dy*((dy*dy)*fx_i)
    f01 = f01 + 3.0d0*((dy*dy)*fx_i)

!       f00 = fx_0 + dy*fx_1 + (dy*dy)*fx_2 + ((dy*dy)*dy)*fx_3
!       f01 =            fx_1 + 2*dy*fx_2    + 3*(dy*dy)*fx_3

    f10 = 0

!       dfx_0 =         A(1,0) + 2*dx*A(2,0)   + 3*(dx*dx)*A(3,0)

    dfx_i = (A(3,0)*3.0d0*dx + A(2,0)*2.0d0)*dx + A(1,0)
    f10 = f10 + dfx_i

!       dfx_1 =         A(1,1) + 2*dx*A(2,1)   + 3*(dx*dx)*A(3,1)
    dfx_i = (A(3,1)*3.0d0*dx + A(2,1)*2.0d0)*dx + A(1,1)
    f10 = f10 + dy*dfx_i

!       dfx_2 =         A(1,2) + 2*dx*A(2,2)   + 3*(dx*dx)*A(3,2)
    dfx_i = (A(3,2)*3.0d0*dx + A(2,2)*2.0d0)*dx + A(1,2)
    f10 = f10 + (dy*dy)*dfx_i

!       dfx_3 =         A(1,3) + 2*dx*A(2,3)   + 3*(dx*dx)*A(3,3)
    dfx_i = (A(3,3)*3.0d0*dx + A(2,3)*2.0d0)*dx + A(1,3)
    f10 = f10 + (dy*dy)*dy*dfx_i
#endif
    return
  end subroutine bicub_interpol1

!cdir$r unroll 
#ifdef USE_GPU
  attributes(device) &
#endif
  real (8) function psi_interpol_00(r,z )
    use precision_mod_gpu
#ifdef USE_C_BICUB
    use ex_bicub_mod, only : eval_0_external
#endif
    implicit none
    real (kind=work_p), intent(in) :: r,z
    logical :: isok
    real (kind=work_p) :: f00
    integer :: nr,nz, i,j
    real (kind=work_p) :: dx, dy, A(0:ndeg,0:ndeg)
    real (kind=work_p) :: fy_i

#ifdef USE_C_BICUB
    call eval_0_external(r,z,f00)
#else
    !real(8) :: psi_acoef_local(0:ndeg,0:ndeg)
    ! ---------------
    ! perform hashing
    ! ---------------
    nr = nr_gpu
    nz = nz_gpu
    i = max(1,min(nr, 1+int( (r-rmin_gpu)*dr_inv_gpu ) ))
    j = max(1,min(nz, 1+int( (z-zmin_gpu)*dz_inv_gpu ) ))
    A(:,:) = acoef_cub_gpu(:,:,i,j)
    dx = (r-rc_cub_gpu(i))
    dy = (z-zc_cub_gpu(j))

    f00 = 0

!        fy_0 = A(0,0) + A(0,1)*dy + A(0,2)*(dy*dy) + A(0,3)*((dy*dy)*dy)

    fy_i = ((A(0,3)*dy + A(0,2))*dy + A(0,1))*dy + A(0,0)
    f00 = f00 + fy_i

!        fy_1 = A(1,0) + A(1,1)*dy + A(1,2)*(dy*dy) + A(1,3)*((dy*dy)*dy)
    fy_i = ((A(1,3)*dy + A(1,2))*dy + A(1,1))*dy + A(1,0)
    f00 = f00 + dx*fy_i
          

!        fy_2 = A(2,0) + A(2,1)*dy + A(2,2)*(dy*dy) + A(2,3)*((dy*dy)*dy)
    fy_i = ((A(2,3)*dy + A(2,2))*dy + A(2,1))*dy + A(2,0)
    f00 = f00 + (dx*dx)*fy_i

!        fy_i = A(3,0) + A(3,1)*dy + A(3,2)*(dy*dy) + A(3,3)*((dy*dy)*dy)
    fy_i = ((A(3,3)*dy + A(3,2))*dy + A(3,1))*dy + A(3,0)
    f00 = f00 + (dx*dx)*dx*fy_i

#endif
    psi_interpol_00 = f00
    return
  end function psi_interpol_00


end module bicub_mod_gpu


#undef RC_CUB_GPU 
#undef ZC_CUB_GPU 
#undef ACOEF_CUB_GPU 
