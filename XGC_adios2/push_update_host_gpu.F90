      subroutine push_update_host_gpu(  sp, &
                                       psn, diag_on, &
                                       gpu_ibegin, gpu_iend )
      use sml_module, only : sml_neutral_use_ion_loss
      use psn_class, only :  psn_type
      use ptl_module, only :  species_type
      use ptl_module_gpu, only : &
         update_host_species_type

      use diag_module_gpu, only : &
          update_host_diag 

      use psn_class_gpu, only : &
          update_host_psn_type

      use neu_module_gpu, only : &
          update_host_neu

      implicit none
      logical, intent(in) :: diag_on
      integer, intent(in) :: gpu_ibegin, gpu_iend

      type(psn_type) :: psn
      type(species_type) :: sp


      integer, parameter :: idebug = 0
      ! need to copy 
      ! ptl(1:sp%num)%ph 
      ! ptl(1:sp%num)%ct
      ! particle-triangle data %itr, %p(3) for E-field calculation
      ! psn%E_rho_ff
      !
      if (idebug >= 1) then
!        if (sml_mype == 0) print*,'before update_host_sml_type '
!        if (sml_mype == 0) print *,'diag_on = ', diag_on
      endif
      ! call update_host_sml()

      if (idebug >= 1)  then
!        if (sml_mype == 0) print*,'before update_host_psn_type '
      endif
      call update_host_psn_type(psn)

      if (idebug >= 1)  then
!        if (sml_mype == 0) print*,'before update_host_neu '
      endif
      if(.not. sml_neutral_use_ion_loss) then
        call update_host_neu()
      endif


      if (idebug >= 1) then
!        if (sml_mype == 0) print*,'before update_host_species_type '
      endif
      call update_host_species_type(sp, &
                                    gpu_ibegin,gpu_iend)


!      call update_host_grid_type()

      !if (diag_on) then -- diag_on is true only when ipc==1
        if (idebug >= 1) then
!          if (sml_mype == 0) print*,'before update_host_diag '
        endif
        call update_host_diag()
      !endif


      return
      end subroutine push_update_host_gpu
