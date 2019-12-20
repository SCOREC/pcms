module one_d_cub_mod_gpu
use precision_mod_gpu
use cudafor

implicit none
   integer, parameter :: ndeg = 3
   integer, parameter :: max_npsi = eq_mpsi-1

   real (kind=work_p),allocatable,target  :: one_d_cub_psi_acoef_gpu(:,:)
   real (kind=work_p),pointer,texture :: tex_acoef(:,:)
   real (kind=work_p),constant :: one_d_cub_dpsi_inv_gpu, one_d_cub_psimin_gpu

   attributes(device)  :: one_d_cub_psi_acoef_gpu
   private :: tex_acoef
contains

   attributes(host) &
   subroutine update_device_one_d()
   use assert_mod
   use precision_mod_gpu
   use one_d_cub_mod, only : &
       one_d_cub_psi_acoef,  one_d_cub_dpsi_inv, one_d_cub_psimin
   use cudafor

implicit none

   integer :: npsi
   integer :: icount_one_d_acoef
   integer :: ierr, i, j
   integer :: lb1,ub1, lb2,ub2
   logical :: isvalid
   integer, parameter :: idebug = 0
#ifdef USE_GPU_EMU
   logical, parameter :: use_cudaMemcpy = .false.
#else
   logical, parameter :: use_cudaMemcpy = .true.
#endif



   if (.not.allocated(one_d_cub_psi_acoef_gpu)) then
     lb1 = lbound(one_d_cub_psi_acoef,1)
     ub1 = ubound(one_d_cub_psi_acoef,1)
     lb2 = lbound(one_d_cub_psi_acoef,2)
     ub2 = ubound(one_d_cub_psi_acoef,2)
     allocate( one_d_cub_psi_acoef_gpu( lb1:ub1, lb2:ub2 ),stat=ierr )
     call assert(ierr.eq.0,'alloc(one_d_cub_psi_acoef_gpu)',ierr)

     one_d_cub_psi_acoef_gpu = 0
   endif


   isvalid = size(one_d_cub_psi_acoef).eq.size(one_d_cub_psi_acoef_gpu)
   call assert(isvalid,                                                  &
     'invalid size(one_d_cub_psi_acoef_gpu',                             &
     size(one_d_cub_psi_acoef_gpu))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!! copy data from CPU to GPU !!!

   
   one_d_cub_dpsi_inv_gpu = one_d_cub_dpsi_inv
   one_d_cub_psimin_gpu = one_d_cub_psimin

   if(use_cudaMemcpy) then
        icount_one_d_acoef = size( one_d_cub_psi_acoef )
        ierr = cudaMemcpy( one_d_cub_psi_acoef_gpu,                 &
                           one_d_cub_psi_acoef,                          &
                           icount_one_d_acoef,                           &
                           cudaMemcpyHostToDevice )

        call assert(ierr.eq.0,                                           &
               'cudaMemcpy(one_d_cub_psi_acoef_gpu',ierr)
   else

     one_d_cub_psi_acoef_gpu = one_d_cub_psi_acoef

   endif
   
  tex_acoef => one_d_cub_psi_acoef_gpu   
   return
   end subroutine update_device_one_d


   attributes(device) &
   subroutine I_interpol_wo_pspline(psi, ideriv, ivalue)
   use precision_mod_gpu
      implicit none

      real (kind=work_p), intent(in) :: psi
      integer, intent(in) :: ideriv
      real (kind=work_p), intent(inout) :: ivalue

      real (kind=work_p) :: pn, wp
      integer :: ip

      pn=psi*one_d_cub_dpsi_inv_gpu
      ip=floor(pn)+1
      ip=min(max(ip,1),ubound(one_d_cub_psi_acoef_gpu,2))
      wp=pn-real(ip-1,kind=work_p)

      if (ideriv==0) then
         ivalue=tex_acoef(0,ip)+(tex_acoef(1,ip)+(tex_acoef(2,ip)+tex_acoef(3,ip)*wp)*wp)*wp
      elseif (ideriv==1) then
         ivalue=(tex_acoef(1,ip)+(2.0_work_p*tex_acoef(2,ip)+3.0_work_p*tex_acoef(3,ip)*wp)*wp)*one_d_cub_dpsi_inv_gpu
      endif
!         print *, 'ideriv in I_interpol_wo_pspline should be 0 or 1'
!         stop
   end subroutine I_interpol_wo_pspline

end module one_d_cub_mod_gpu
