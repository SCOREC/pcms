module ptl_module_gpu
use dimensions_mod_gpu
use precision_mod_gpu
use ptl_module, only :              &
   ptl_nphase_host => ptl_nphase,   &
   ptl_nphase2_host => ptl_nphase2, &
   ptl_nconst_host => ptl_nconst,   &
   ptl_nsp_max_host => ptl_nsp_max, &
   pir_host => pir,                 &
   piz_host => piz,                 &
   pip_host => pip,                 &
   pirho_host => pirho,             & 
   pim_host => pim,                 &
   piw1_host => piw1,               &
   piw2_host => piw2,               &
   piw0_host => piw0

use util_mod_gpu, only : get_gpu_streamid
use cudafor
implicit none

!  logical, parameter :: use_ptl_block_copy = .false.

  integer, parameter :: ptl_nphase = ptl_nphase_host
  integer, parameter :: ptl_nphase2 = ptl_nphase2_host
  integer, parameter :: ptl_nconst = ptl_nconst_host
  integer, parameter :: ptl_nsp_max = ptl_nsp_max_host

  integer, parameter :: pir = pir_host
  integer, parameter :: piz = piz_host
  integer, parameter :: pip = pip_host
  integer, parameter :: pirho = pirho_host
  integer, parameter :: pim = pim_host
  integer, parameter :: piw1 = piw1_host
  integer, parameter :: piw2 = piw2_host
  integer, parameter :: piw0 = piw0_host

!  integer, parameter :: spe_num = spe_num_host
!  integer, parameter :: spe_maxnum = spe_maxnum_host
!  integer, parameter :: spe_type = spe_type_host

  integer (kind=work_p),allocatable :: sp_ptl_gid(:) !!! new CPU arrays
  real (kind=work_p),allocatable :: sp_ptl_ph(:,:) !!! new CPU arrays
  real (kind=work_p),allocatable :: sp_ptl_ct(:,:) !!! new CPU arrays
  real (kind=work_p),allocatable :: sp_phase0(:,:) !!! new CPU arrays

#ifdef USE_GPU_ASYNC
! --------------------------------------------------
! asynchronous transfer by cudaMemcpyAsync requires
! non-pageable pinned memory
! --------------------------------------------------
  attributes(pinned) :: sp_ptl_gid, sp_ptl_ph, sp_ptl_ct, sp_phase0
#endif

  logical, parameter :: use_update_phase0 = .true.

  integer :: ptl_isp, ptl_nsp
  integer ::  num_gpu,maxnum_gpu,type_gpu
  integer(8),allocatable :: ptl_gid_gpu(:)
  real (kind=work_p),allocatable :: ptl_ph_gpu(:,:)
  real (kind=work_p),allocatable :: ptl_ct_gpu(:,:)
  real (kind=work_p),allocatable :: phase0_gpu(:,:)
  real (kind=work_p),allocatable :: rhoi_gpu(:)

! space to store the last known triangle number for a particle
! Not allocated unless USE_TR_CHECK is defined!
!! NB: tr_save_gpu is being used differently than tr_save
!! on the host. Here it stores the current triangle of the
!! particle, not where it was last timestep.

  integer, allocatable :: tr_save_gpu(:)
  attributes(device) :: tr_save_gpu


  real(kind=work_p), dimension(0:ptl_nsp_max) :: ptl_mass, ptl_c_m, ptl_c2_2m, ptl_charge
  logical,dimension(0:ptl_nsp_max) :: ptl_deltaf_sp
  integer,parameter :: idebug=0


  attributes(device) :: ptl_gid_gpu, ptl_ph_gpu, ptl_ct_gpu, phase0_gpu, rhoi_gpu
  attributes(constant) :: num_gpu,maxnum_gpu,type_gpu
  attributes(constant) :: ptl_mass, ptl_c_m, ptl_c2_2m, ptl_charge, ptl_deltaf_sp
  attributes(constant) :: ptl_isp, ptl_nsp

  contains

  attributes(host) &
  subroutine update_device_species_type(sp, gpu_ibegin,gpu_iend)
  use sml_module, only : sml_mype
  use ptl_module, only :  &
      species_type              

implicit none

  type(species_type),intent(in) :: sp
  integer, intent(in) :: gpu_ibegin, gpu_iend !!! isp=0 electron, isp=1 ion

!  integer :: lb1,ub1,lb2,ub2, llb,uub
  integer :: ierr
  integer :: i, nn, nc, np , iphase
#ifdef USE_GPU_EMU
  logical, parameter :: use_cudaMemcpy = .false.
#else
  logical, parameter :: use_cudaMemcpy = .true.
#endif
  logical :: isvalid

  integer :: icount,streamid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!! convert the ptl_type data structure !!!


  streamid = get_gpu_streamid()


  allocate(sp_ptl_ct(gpu_ibegin:gpu_iend,ptl_nconst), &
           sp_ptl_ph(gpu_ibegin:gpu_iend,ptl_nphase), &
           sp_ptl_gid(gpu_ibegin:gpu_iend), stat=ierr )
 

  do nn = gpu_ibegin, gpu_iend
     sp_ptl_gid(nn)=sp%ptl(nn)%gid
     do np = 1, ptl_nphase
        sp_ptl_ph(nn,np) = sp%ptl(nn)%ph(np)
     enddo
     do nc = 1, ptl_nconst
        sp_ptl_ct(nn,nc) = sp%ptl(nn)%ct(nc)
     enddo
  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!! copy data from CPU to GPU !!!

  if (.not.allocated(ptl_gid_gpu)) then
     allocate( ptl_gid_gpu( gpu_ibegin:gpu_iend), stat=ierr)
     call assert(ierr.eq.0,'alloc(ptl_gid_gpu)',ierr)
  endif

  if (.not.allocated( ptl_ph_gpu)) then
    allocate( ptl_ph_gpu(  gpu_ibegin:gpu_iend,ptl_nphase),stat=ierr)
    call assert(ierr.eq.0,'alloc(ptl_ph_gpu)',ierr)
  endif

  if (.not.allocated(ptl_ct_gpu)) then
    allocate(  ptl_ct_gpu(  gpu_ibegin:gpu_iend,ptl_nconst),stat=ierr)
    call assert(ierr.eq.0,'alloc(ptl_ct_gpu)',ierr)
  endif

  if (.not.allocated(phase0_gpu)) then
    allocate( phase0_gpu(  gpu_ibegin:gpu_iend,ptl_nphase),stat=ierr)
    call assert(ierr.eq.0,'alloc(phase0_gpu)',ierr)
  endif

#ifdef USE_TR_CHECK
! Space is needed for tr_save when using the tr_check optimization,
! but it's probably not worth copying it back and forth.
  if (.not.allocated(tr_save_gpu)) then
    allocate( tr_save_gpu(  gpu_ibegin:gpu_iend),stat=ierr)
    call assert(ierr.eq.0,'alloc(tr_save_gpu)',ierr)
  endif
#endif


!!  allocate( ptl_gid_gpu( gpu_ibegin:gpu_iend),                      &
!!            ptl_ph_gpu(  gpu_ibegin:gpu_iend,ptl_nphase),           &
!!            ptl_ct_gpu(  gpu_ibegin:gpu_iend,ptl_nconst),           &
!!            phase0_gpu(  gpu_ibegin:gpu_iend,ptl_nphase)  )

  num_gpu = sp%num
  maxnum_gpu = sp%maxnum
  type_gpu = sp%type

     isvalid = size(ptl_gid_gpu).eq.size(sp_ptl_gid)
     call assert(isvalid,'invalid size(ptl_gid_gpu)',size(ptl_gid_gpu))

     isvalid = size(ptl_ph_gpu).eq.size(sp_ptl_ph)
     call assert(isvalid,'invalid size(ptl_ph_gpu)',size(ptl_ph_gpu))

     isvalid = size(ptl_ct_gpu).eq.size(sp_ptl_ct)
     call assert( isvalid,'invalid size(ptl_ct_gpu)',size(ptl_ct_gpu))

 if(use_cudaMemcpy) then


     icount = size(ptl_gid_gpu)
#ifdef USE_GPU_ASYNC
     ierr = cudaMemcpyAsync( ptl_gid_gpu, sp_ptl_gid,  &
                icount, cudaMemcpyHostToDevice,streamid)
     call assert(ierr.eq.0,'cudaMemcpyAsync(ptl_gid_gpu)',ierr)
#else
     ierr = cudaMemcpy( ptl_gid_gpu, sp_ptl_gid,  &
                icount, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(ptl_gid_gpu)',ierr)
#endif




     icount = size(ptl_ph_gpu)

#ifdef USE_GPU_ASYNC
     ierr = cudaMemcpyAsync( ptl_ph_gpu, sp_ptl_ph,    &
               icount, cudaMemcpyHostToDevice,streamid)
     call assert(ierr.eq.0,'cudaMemcpyAsync(ptl_ph_gpu)',ierr)
#else
     ierr = cudaMemcpy( ptl_ph_gpu, sp_ptl_ph,    &
               icount, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(ptl_ph_gpu)',ierr)
#endif




     icount = size(ptl_ct_gpu)
#ifdef USE_GPU_ASYNC
     ierr = cudaMemcpyAsync( ptl_ct_gpu, sp_ptl_ct,    &
               icount, cudaMemcpyHostToDevice,streamid)
     call assert(ierr.eq.0,'cudaMemcpyAsync(ptl_ct_gpu)',ierr)
#else
     ierr = cudaMemcpy( ptl_ct_gpu, sp_ptl_ct,    &
               icount, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(ptl_ct_gpu)',ierr)
#endif

!     if(sp%type/=0) then
!        allocate(rhoi_gpu(gpu_ibegin:gpu_iend))
!        ierr = cudaMemcpyAsync( rhoi_gpu, sp%rhoi,    &
!                  icount_gid, cudaMemcpyHostToDevice,streamid)
!     endif
 else
     ptl_gid_gpu(gpu_ibegin:gpu_iend) = sp_ptl_gid(gpu_ibegin:gpu_iend)
     ptl_ph_gpu(gpu_ibegin:gpu_iend,1:ptl_nphase) = sp_ptl_ph(gpu_ibegin:gpu_iend,1:ptl_nphase)
     ptl_ct_gpu(gpu_ibegin:gpu_iend,1:ptl_nconst) = sp_ptl_ct(gpu_ibegin:gpu_iend,1:ptl_nconst)

!     if(type_gpu/=0) then
!        allocate(rhoi_gpu(gpu_ibegin:gpu_iend))
!        rhoi_gpu(gpu_ibegin:gpu_iend) = sp%rhoi(gpu_ibegin:gpu_iend)
!     endif

 endif

#ifdef USE_GPU_ASYNC
   ierr = cudaStreamSynchronize( streamid )
   call assert(ierr.eq.0,'cudStreamSync',ierr)
#endif
  
  deallocate(sp_ptl_gid, stat=ierr)
  call assert(ierr.eq.0,'dealloc(sp_ptl_gid)',ierr)

  deallocate(sp_ptl_ph, stat=ierr)
  call assert(ierr.eq.0,'dealloc(sp_ptl_ph)',ierr)

  deallocate(sp_ptl_ct,stat=ierr)
  call assert(ierr.eq.0,'dealloc(sp_ptl_ct)',ierr)



!  ---------------------
!  update phase0 array?
!  ---------------------

  if (use_update_phase0) then
    allocate( sp_phase0(gpu_ibegin:gpu_iend,1:ptl_nphase),stat=ierr)
    call assert(ierr.eq.0,'alloc(sp_phase0)',ierr)

    do nn=gpu_ibegin,gpu_iend
     sp_phase0(nn,1:ptl_nphase) = sp%phase0(1:ptl_nphase,nn)
    enddo

    if (use_cudaMemcpy) then
       icount = size(phase0_gpu)
#ifdef USE_GPU_ASYNC
       ierr = cudaMemcpyAsync( phase0_gpu, sp_phase0, &
                 icount,  cudaMemcpyHostToDevice,streamid)
       call assert(ierr.eq.0,'cudaMemcpyAsync(phase0_gpu)',ierr)
#else
       ierr = cudaMemcpy( phase0_gpu, sp_phase0, &
                 icount,  cudaMemcpyHostToDevice)
       call assert(ierr.eq.0,'cudaMemcpy(phase0_gpu)',ierr)
#endif

    else
      phase0_gpu = sp_phase0
    endif


#ifdef USE_GPU_ASYNC
   ierr = cudaStreamSynchronize( streamid )
   call assert(ierr.eq.0,'cudStreamSync',ierr)
#endif



    deallocate( sp_phase0,stat=ierr)
    call assert(ierr.eq.0,'dealloc(sp_phase0)',ierr)
  endif


  return
  end subroutine update_device_species_type


  attributes(host) &
  subroutine update_host_species_type(sp,gpu_ibegin,gpu_iend)
  use sml_module, only : sml_mype

   use ptl_module, only :   species_type,  ptl_type

  integer, intent(in) ::  gpu_ibegin, gpu_iend
  type(species_type),intent(inout) :: sp
  integer :: ierr, stat
!  logical :: isvalid
  integer :: i, nn, np, nc
  integer :: icount_gid,icount_ph,icount_ct
  integer, parameter :: idebug = 0

#ifdef USE_GPU_EMU
  logical, parameter :: use_cudaMemcpy = .false.
#else
  logical, parameter :: use_cudaMemcpy = .true.
#endif
  integer :: icount, streamid


  streamid = get_gpu_streamid()

  allocate(sp_ptl_ct(gpu_ibegin:gpu_iend,ptl_nconst), &
           sp_ptl_ph(gpu_ibegin:gpu_iend,ptl_nphase), &
           sp_ptl_gid(gpu_ibegin:gpu_iend), stat=ierr )
  call assert(ierr.eq.0,'alloc(sp_ptl_ct)',ierr)

  icount_ph = ptl_nphase*(gpu_iend-gpu_ibegin+1)
  icount_gid = gpu_iend-gpu_ibegin+1
  icount_ct = ptl_nconst*(gpu_iend-gpu_ibegin+1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!! copy data from GPU to CPU !!!
  
  sp%num = num_gpu
  sp%maxnum = maxnum_gpu
  sp%type = type_gpu


!  sp%tr_save(gpu_ibegin:gpu_iend) = tr_save_gpu(gpu_ibegin:gpu_iend)
!  sp%p_save(1:3,gpu_ibegin:gpu_iend) = p_save_gpu(1:3,gpu_ibegin:gpu_iend)

 if(use_cudaMemcpy) then

     icount = size(sp_ptl_gid)
#ifdef USE_GPU_ASYNC
     ierr = cudaMemcpyAsync( sp_ptl_gid, ptl_gid_gpu,    &
               icount, cudaMemcpyDeviceToHost,streamid)
     call assert( ierr.eq.0,'cudaMemcpyAsync(sp_ptl_gid)', ierr )
#else
     ierr = cudaMemcpy( sp_ptl_gid, ptl_gid_gpu,    &
               icount, cudaMemcpyDeviceToHost)
     call assert( ierr.eq.0,'cudaMemcpy(sp_ptl_gid)', ierr )
#endif



     icount = size(sp_ptl_ph)
#ifdef USE_GPU_ASYNC
     ierr = cudaMemcpyAsync( sp_ptl_ph, ptl_ph_gpu,    &
               icount, cudaMemcpyDeviceToHost,streamid)
     call assert( ierr.eq.0,'cudaMemcpyAsync(sp_ptl_ph)', ierr )
#else
     ierr = cudaMemcpy( sp_ptl_ph, ptl_ph_gpu,    &
               icount, cudaMemcpyDeviceToHost)
     call assert( ierr.eq.0,'cudaMemcpy(sp_ptl_ph)', ierr )
#endif




     icount = size(sp_ptl_ct)
#ifdef USE_GPU_ASYNC
     ierr = cudaMemcpyAsync( sp_ptl_ct, ptl_ct_gpu,    &
               icount, cudaMemcpyDeviceToHost, streamid)
     call assert( ierr.eq.0,'cudaMemcpyAsync(sp_ptl_ct)', ierr )
#else
     ierr = cudaMemcpy( sp_ptl_ct, ptl_ct_gpu,    &
               icount, cudaMemcpyDeviceToHost )
     call assert( ierr.eq.0,'cudaMemcpy(sp_ptl_ct)', ierr )
#endif

!     if(sp%type/=0) then
!        ierr = cudaMemcpyAsync( sp%rhoi, rhoi_gpu,    &
!                  icount_gid, cudaMemcpyDeviceToHost,streamid)
!     endif
 else
     sp_ptl_gid = ptl_gid_gpu
     sp_ptl_ph = ptl_ph_gpu
     sp_ptl_ct = ptl_ct_gpu

!     if(sp%type/=0) then
!        sp%rhoi(gpu_ibegin:gpu_iend) = rhoi_gpu(gpu_ibegin:gpu_iend)
!     endif
 endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!! convert the simple array to ptl_type data structure !!!
#ifdef USE_GPU_ASYNC
   ierr = cudaStreamSynchronize( streamid )
   call assert(ierr.eq.0,'cudStreamSync',ierr)
#endif

  do nn = gpu_ibegin, gpu_iend
     sp%ptl(nn)%gid = sp_ptl_gid(nn)
     do np = 1,ptl_nphase
        sp%ptl(nn)%ph(np) = sp_ptl_ph(nn,np)
     enddo
     do nc = 1,ptl_nconst
        sp%ptl(nn)%ct(nc) = sp_ptl_ct(nn,nc)
     enddo
  enddo

  deallocate(ptl_gid_gpu, ptl_ph_gpu, ptl_ct_gpu,stat=ierr) 
  call assert(ierr.eq.0,'dealloc(ptl_gid_gpu)',ierr)

  deallocate(sp_ptl_gid, sp_ptl_ph, sp_ptl_ct,stat=ierr )
  call assert(ierr.eq.0,'dealloc(sp_ptl_gid)',ierr)

  if (use_update_phase0) then

    allocate( sp_phase0(gpu_ibegin:gpu_iend,1:ptl_nphase), stat=ierr)
    call assert(ierr.eq.0,'alloca(sp_phase0)',ierr)


    if (use_cudaMemcpy) then
       icount = size(sp_phase0)
#ifdef USE_GPU_ASYNC
       ierr = cudaMemcpyAsync( sp_phase0, phase0_gpu,  &
                 icount, cudaMemcpyDeviceToHost, streamid)
       call assert(ierr.eq.0,'cudaMemcpyAsync(sp_phase0)',ierr)
#else
       ierr = cudaMemcpy( sp_phase0, phase0_gpu,  &
                 icount, cudaMemcpyDeviceToHost )
       call assert(ierr.eq.0,'cudaMemcpy(sp_phase0)',ierr)
#endif
    else 
       sp_phase0 = phase0_gpu
    endif

#ifdef USE_GPU_ASYNC
   ierr = cudaStreamSynchronize( streamid )
   call assert(ierr.eq.0,'cudStreamSync',ierr)
#endif

    do nn=gpu_ibegin,gpu_iend
      sp%phase0(1:ptl_nphase, nn) = sp_phase0(nn,1:ptl_nphase)
    enddo

    deallocate( sp_phase0, stat=ierr)
    call assert(ierr.eq.0,'dealloc(sp_phase0)',ierr)

  endif

  deallocate(phase0_gpu,stat=ierr)
  call assert(ierr.eq.0,'dealloc(phase0_gpu)',ierr)

#ifdef USE_TR_CHECK
  deallocate(tr_save_gpu,stat=ierr)
  call assert(ierr.eq.0,'dealloc(tr_save_gpu)',ierr)
#endif

  return
  end subroutine update_host_species_type


  attributes(host) &
  subroutine update_device_ptl( )
  use ptl_module, only : &
    ptl_mass_host => ptl_mass, &
    ptl_c_m_host => ptl_c_m, &
    ptl_charge_host => ptl_charge, &
    ptl_c2_2m_host => ptl_c2_2m, &
    ptl_deltaf_sp_host => ptl_deltaf_sp, &
    ptl_isp_host => ptl_isp, &
    ptl_nsp_host => ptl_nsp


    integer :: lb,ub
    logical :: isvalid

    lb = lbound(ptl_mass_host,1)
    ub = ubound(ptl_mass_host,1)
    isvalid = (lbound(ptl_mass,1) <= lb) 
    call assert(isvalid, &
       'ptl_module_gpu: invalid lb ', lb)
    isvalid = (ub <= ubound(ptl_mass,1) )
    call assert(isvalid, &
       'ptl_module_gpu: invalid ub ', ub)

    ptl_mass(lb:ub)   = ptl_mass_host(lb:ub)
    ptl_c_m(lb:ub)    = ptl_c_m_host(lb:ub)
    ptl_c2_2m(lb:ub)  = ptl_c2_2m_host(lb:ub)
    ptl_charge(lb:ub) = ptl_charge_host(lb:ub)
    ptl_deltaf_sp(lb:ub)   = ptl_deltaf_sp_host(lb:ub)
    ptl_isp = ptl_isp_host
    ptl_nsp = ptl_nsp_host
  return
  end subroutine update_device_ptl

end module ptl_module_gpu
