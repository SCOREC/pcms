!    ----------------------------------------
!    calculate E_phi_ff(1:3, 0:1, inode, iphi)
!    ----------------------------------------
      attributes(device)                                                 &
     &subroutine calc_E_phi_ff_gpu(                                      &
     &                          inode, iphi, E_phi_ff)
      use psn_class_gpu, only : psn_pot_phi_real,                        &
     &          psn_ff_hdp_tr, psn_ff_hdp_p,                             &
     &          psn_ff_1dp_tr, psn_ff_1dp_p,                             &
     &          psn_bfollow_1_dx

      use grid_class_gpu, only :  grid_max_v, grid_nnode,                &
     &       grid_ntriangle, grid_gradx,grid_grady, grid_v_node,         &
     &       grid_num_v_node

      use sml_module_gpu, only : sml_bt_sign 
      implicit none


      integer, intent(in) :: inode
      integer, intent(in) :: iphi
      real(8), intent(inout) :: E_phi_ff(1:3,0:1)


!     ---------------
!     local variables
!     ---------------
      integer :: max_v
      integer :: nphi, iphip1,iphip2, iphim1 
      real(8) :: E_r_ff(0:1), E_z_ff(0:1), E_para(0:1)


      integer, parameter :: idebug = 0

      integer :: INC, DEC
      INC(iphi,nphi) = mod( ((iphi-1)+nphi) +1,nphi)+1  
      DEC(iphi,nphi) = mod( ((iphi-1)+nphi) -1,nphi)+1

      E_phi_ff(:,:) = 0

      ! ---------------------------------------------------------
      ! assume dimension of psn_pot_phi_real(1:grid%nnode, 0:nphim1 )
      ! ---------------------------------------------------------
      nphi = size( psn_pot_phi_real, dim=2)

      max_v = grid_max_v


      if (lbound(psn_pot_phi_real,2) .eq. 0) then
!       ---------------------------
!       note 0 <= iphi  <= (nphi-1)
!       ---------------------------

        iphip1 = mod(iphi+1,nphi)
        iphip2 = mod(iphip1 + 1, nphi)
        iphim1 = mod((iphi-1) + nphi, nphi )
      else
!       ----------------------
!       note 1 <= iphi <= nphi
!       ----------------------
        iphip1 = INC(iphi,nphi)
        iphip2 = INC(iphip1,nphi)
        iphim1 = DEC(iphi,nphi)
      endif



!     ----------------------------------
!     compute (E_r_ff, E_z_ff) from (E_r_real,E_z_real)
!     ----------------------------------

      call calc_Er_Ez_gpu(grid_nnode,grid_ntriangle,                    &
     &        max_v,                                                    &
     &        psn_ff_hdp_tr,psn_ff_hdp_p,                               &
     &        psn_pot_phi_real(:,iphi), psn_pot_phi_real(:,iphip1),     &
     &        inode, E_r_ff, E_z_ff )

!     --------------
!     compute E_para
!     --------------


          call calc_E_para_gpu( grid_nnode, grid_ntriangle,             &
     &       psn_ff_hdp_tr,psn_ff_hdp_p,                                &
     &      psn_ff_1dp_tr, psn_ff_1dp_p,                                &
     &      nphi, iphi, psn_pot_phi_real,                               &
     &      sml_bt_sign, psn_bfollow_1_dx,                              &
     &      inode, E_para )

!         ----------------
!         final assignment
!         ----------------
          E_phi_ff(1,0:1) = E_r_ff(0:1)
          E_phi_ff(2,0:1) = E_z_ff(0:1)
          E_phi_ff(3,0:1) = E_para(0:1)

          return
          end subroutine calc_E_phi_ff_gpu
