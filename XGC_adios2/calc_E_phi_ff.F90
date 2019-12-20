!    ----------------------------------------
!    calculate E_phi_ff(1:3, 0:1, inode, iphi)
!    ----------------------------------------
      subroutine calc_E_phi_ff( grid, psn,                              &
     &                          inode, iphi, E_phi_ff)
      use psn_class, only : psn_type
      use grid_class, only : grid_type
      use sml_module, only : sml_bt_sign 
      implicit none

      type(grid_type), intent(in) :: grid
      type(psn_type), intent(in) :: psn

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
      ! assume dimension of psn%pot_phi_real(1:grid%nnode, 0:nphim1 )
      ! ---------------------------------------------------------
      nphi = size( psn%pot_phi_real, dim=2)

      max_v = size(grid%v_node,dim=1)


      if (lbound(psn%pot_phi_real,2) .eq. 0) then
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

      call calc_Er_Ez(grid%nnode,grid%ntriangle, grid%nd,               &
     &        max_v, grid%num_v_node, grid%v_node,                      &
     &        grid%gradx, grid%grady,                                   &
     &        psn%ff_hdp_tr,psn%ff_hdp_p,                               &
     &        psn%pot_phi_real(:,iphi), psn%pot_phi_real(:,iphip1),     &
     &        inode, E_r_ff, E_z_ff )

!     --------------
!     compute E_para
!     --------------


          call calc_E_para( grid%nnode, grid%ntriangle, grid%nd,        &
     &       psn%ff_hdp_tr,psn%ff_hdp_p,                                &
     &      psn%ff_1dp_tr, psn%ff_1dp_p,                                &
     &      nphi, iphi, psn%pot_phi_real,                               &
     &      sml_bt_sign, psn%bfollow_1_dx,                              &
     &      inode, E_para )

!         ----------------
!         final assignment
!         ----------------
          E_phi_ff(1,0:1) = E_r_ff(0:1)
          E_phi_ff(2,0:1) = E_z_ff(0:1)
          E_phi_ff(3,0:1) = E_para(0:1)

          return
          end subroutine calc_E_phi_ff
