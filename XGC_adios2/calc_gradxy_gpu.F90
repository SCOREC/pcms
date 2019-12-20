
      attributes(device)                                                &
     &subroutine calc_gradxy_gpu( nnode, max_v,                         &
     &             pot,                                                 &
     &             vi,  gradx, grady )

      use grid_class_gpu, only :                                        &
     &   num_v_node => grid_num_v_node,                                 &
     &   v_node => grid_v_node,                                         &
     &   gx => grid_gradx,                                              &
     &   gy => grid_grady
      implicit none
!     -----------------------------------------------
!     compute the gradient at a vertex  of a triangle
!     -----------------------------------------------
      integer,intent(in) :: nnode, max_v, vi
      real(8), dimension(1:nnode), intent(in) :: pot
      real(8),intent(inout) :: gradx, grady

!     ---------------
!     local variables
!     ---------------

      integer :: j, vj

      gradx = 0
      grady = 0
      do j=1,num_v_node(vi)
         vj = v_node(j,vi)
         gradx = gradx + pot(vj) * gx(j, vi)
         grady = grady + pot(vj) * gy(j, vi)
      enddo
      return
      end subroutine calc_gradxy_gpu
