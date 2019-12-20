      subroutine calc_gradxy( nnode, max_v,                             &
     &             num_v_node, v_node, gx, gy, pot,                     &
     &             vi,  gradx, grady )
      implicit none
!     -----------------------------------------------
!     compute the gradient at a vertex  of a triangle
!     -----------------------------------------------
      integer,intent(in) :: nnode, max_v, vi
      integer, dimension(1:nnode), intent(in) :: num_v_node
      integer, dimension(1:max_v,1:nnode), intent(in) :: v_node
      real(8), dimension(1:max_v,1:nnode), intent(in) :: gx
      real(8), dimension(1:max_v,1:nnode), intent(in) :: gy
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
      end subroutine calc_gradxy
