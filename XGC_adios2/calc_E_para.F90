      subroutine calc_E_para( grid_nnode, grid_ntriangle, grid_nd,      &
     &              ff_hdp_tr, ff_hdp_p,                                &
     &              ff_1dp_tr, ff_1dp_p,                                &
     &              nphi, iphi, dpot,                                   &
     &              sml_bt_sign, psn_bfollow_1_dx,                      &
     &                        inode, E_para )
      implicit none
      integer, intent(in) :: grid_nnode
      integer, intent(in) :: grid_ntriangle
      integer, intent(in) :: grid_nd(3,grid_ntriangle)

      integer, intent(in) :: ff_hdp_tr(grid_nnode,0:1)
      real(8), intent(in) :: ff_hdp_p(3,grid_nnode,0:1)

      integer, intent(in) :: ff_1dp_tr(grid_nnode,0:1)
      real(8), intent(in) :: ff_1dp_p(3,grid_nnode,0:1)

      integer, intent(in) :: nphi
      integer, intent(in) :: iphi
      real(8), intent(in) :: dpot(grid_nnode,0:(nphi-1))

      real(8), intent(in) :: sml_bt_sign
      real(8), intent(in) :: psn_bfollow_1_dx(grid_nnode)

      integer, intent(in) :: inode
      real(8), intent(inout) :: E_para(0:1)

      integer :: dir, i, j, vj, trig
      real(8) :: wt, E_mid, E_left, E_right
      real(8) :: potm1, pot0, pot1, pot2
      real(8) :: E_para_right_vj, E_para_left_vj


      integer :: vi, jphi,jdir 
      real(8) :: E_mid_0, E_mid_p1, E_mid_m1
      real(8) :: dpot_ff
      integer :: DEC, INC

!     ----------------
!     inline functions
!     ----------------
      INC(i,nphi) = mod(i+1+nphi,nphi)
      DEC(i,nphi) = mod(i-1+nphi,nphi)
      
      dpot_ff(vi,jphi,jdir) = sum( ff_hdp_p(1:3,vi,jdir) *                         &
     &                  dpot( grid_nd(1:3,ff_hdp_tr(vi,jdir)), jphi ) )


      E_mid_0(vi,jphi) = -sml_bt_sign*(dpot_ff(vi,INC(jphi,nphi),1)-        &
     &        dpot_ff(vi,jphi,0) )*psn_bfollow_1_dx(vi)*2.0d0

      E_mid_p1(vi,jphi) = E_mid_0(vi,INC(jphi,nphi))
      E_mid_m1(vi,jphi) = E_mid_0(vi,DEC(jphi,nphi))
!     ----------------------------
!     compute the E_para component
!     ----------------------------
!
!   ================
!   original code
!
! E-parallel
! pot(:,0,1) is tmp variable
! obtain E_parallel at half position  --  1_dx is ~ 2 R dphi
!            pot(:,0,1)= -sml_bt_sign*(pot(:,1,0)-pot(:,0,0))*
!                 psn%bfollow_1_dx(:)*2D0  ! ### replace bfollow_1_dx
!
!
!  ! send left recv right
!  idest=mod(sml_intpl_mype-1+sml_intpl_totalpe,sml_intpl_totalpe)
!  isendtag=sml_intpl_mype
!  isource=mod(sml_intpl_mype+1,sml_intpl_totalpe)
!  irecvtag=isource
!  call mpi_sendrecv(pot(:,0,1),grid%nnode,MPI_REAL8,  idest,isendtag, &
!       E_para(:,1),grid%nnode,MPI_REAL8,isource,irecvtag, &
!       sml_intpl_comm,istatus,ierr)
!  
!  ! send right recv left
!  idest=mod(sml_intpl_mype+1,sml_intpl_totalpe)
!  isendtag=sml_intpl_mype
!  isource=mod(sml_intpl_mype-1+sml_intpl_totalpe,sml_intpl_totalpe)
!  irecvtag=isource
!  call mpi_sendrecv(pot(:,0,1),grid%nnode,MPI_REAL8,  idest,isendtag, &
!       E_para(:,0),grid%nnode,MPI_REAL8,isource,irecvtag, &
!       sml_intpl_comm,istatus,ierr)
!  call t_stopf("GET_POT_SR")
!  
!  ! require converting of E_para(:,0) and E_perp(:,1) ????
!  call t_startf("GET_POT_CNVRT")
!  call cnvt_grid_real2ff(grid,psn%ff_1dp_tr,psn%ff_1dp_p,E_para(:,:),pot(:,:,2)) ! ???
!  E_para(:,0)=0.5D0*(pot(:,0,2)+pot(:,0,1))
!  E_para(:,1)=0.5D0*(pot(:,1,2)+pot(:,0,1))
!   ================

      E_para(:) = 0

      i = inode
!     -----------------------------
!     compute the left contribution,
!     virtual plane at dphi=-0.5
!
!     perform action of
!     cnvt_grid_real2ff(grid,  tr=>psn%ff_1dp_tr,   p=>psn%ff_1dp_p....)
!     -----------------------------
      E_left = 0
      dir = 0
      trig = ff_1dp_tr(i,dir)
      do j=1,3
        vj = grid_nd(j, trig)
        wt = ff_1dp_p(j,i,dir)


        E_para_left_vj = E_mid_m1(vj,iphi)
        E_left = E_left + wt*E_para_left_vj
      enddo

!     -----------------------------
!     compute the right contribution,
!     virtual plane at dphi=1.5
!
!     perform action of
!     cnvt_grid_real2ff(grid,  tr=>psn%ff_1dp_tr,   p=>psn%ff_1dp_p....)
!     -----------------------------
      E_right = 0
      dir = 1
      trig = ff_1dp_tr(i,dir)
      do j=1,3
        vj = grid_nd(j, trig)
        wt = ff_1dp_p(j,i,dir)


        E_para_right_vj = E_mid_p1(vj,iphi)
        E_right = E_right + wt*E_para_right_vj
      enddo
!     ----------------------------
!     compute the mid contribution,
!     virtual plane at dphi=0.5
!     ----------------------------
      E_mid = E_mid_0(inode,iphi)

      E_para(0) = 0.5d0*(E_left  + E_mid)
      E_para(1) = 0.5d0*(E_right + E_mid)

      return
      end subroutine calc_E_para


