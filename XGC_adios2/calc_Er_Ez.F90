      subroutine calc_Er_Ez(nnode, ntriangle, nd,                       &
     &    max_v, num_v_node, v_node, gx, gy,                            &
     &    tr, p,                                                        &
     &    pot0, pot1,                                                   &
     &    inode, E_r_ff, E_z_ff )
      implicit none

      integer, intent(in) :: nnode
      integer, intent(in) :: ntriangle
      integer, intent(in) :: nd(3,ntriangle)
      integer, intent(in) :: max_v
      integer, intent(in) :: num_v_node(1:nnode)
      integer, intent(in) :: v_node(max_v,nnode)
      real(8), intent(in) :: gx(max_v,nnode)
      real(8), intent(in) :: gy(max_v,nnode)

      integer, intent(in) :: tr(nnode,0:1)
      real(8), intent(in) :: p(3,nnode,0:1)

      real(8), intent(in) :: pot0(nnode)
      real(8), intent(in) :: pot1(nnode)


      integer, intent(in) :: inode
      real(8), intent(inout) :: E_r_ff(0:1)
      real(8), intent(inout) :: E_z_ff(0:1)

!     ----------------------------------
!     compute (E_r_ff, E_z_ff) from (E_r_real,E_z_real)
!     ----------------------------------
      integer :: i,j, dir,itr, vj 
      real(8) :: wt, gradx, grady
      real(8) :: E_r_real_dir, E_z_real_dir

      E_r_ff(:) = 0
      E_z_ff(:) = 0

      i = inode
      do dir=0,1
        itr = tr(i,dir)

        do j=1,3
          vj = nd(j,itr)
          wt = p(j,i,dir)

!         -------------------------
!         compute E_r_real,E_z_real
!         -------------------------
          

          if (dir == 1) then
           call calc_gradxy( nnode, max_v,                              &
     &          num_v_node, v_node,                                     &
     &          gx, gy,                                                 &
     &          pot1,                                                   &
     &          vj,  gradx, grady )
          else
           call calc_gradxy( nnode, max_v,                              &
     &          num_v_node, v_node,                                     &
     &          gx, gy,                                                 &
     &          pot0,                                                   &
     &          vj,  gradx, grady )
          endif


          E_r_real_dir = gradx
          E_z_real_dir = grady

!         ----------------------
!         compute E_r_ff, E_z_ff
!         ----------------------


          E_r_ff(dir) = E_r_ff(dir) + wt*E_r_real_dir
          E_z_ff(dir) = E_z_ff(dir) + wt*E_z_real_dir
        enddo
      enddo

      return
      end subroutine calc_Er_Ez
