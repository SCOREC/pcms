      subroutine gen_gradxy(nnode,ntr,max_t,max_v,                       &
     &      nd, mapping, tr_area, num_t_node, tr_node,                   &
     &      num_v_node, v_node,  gx, gy, do_write )
      implicit none
!     ---------------------------------------------------------
!     generate weights to estimate the gradient at vertices
!     compute gradient at center of triangle
!     then perform weighted average based on areas of triangles
!     ---------------------------------------------------------
      integer, parameter :: idebug = 0
      integer, parameter :: dp = selected_real_kind(15,100)

      integer :: nnode, ntr, max_t, max_v
      integer, dimension(1:nnode) :: num_t_node
      integer, dimension(1:max_t, 1:nnode) :: tr_node
      integer, dimension(1:3, 1:ntr) :: nd
      real(kind=dp), dimension(1:2,1:3,1:ntr) :: mapping
      real(kind=dp), dimension(1:ntr) :: tr_area

      real(kind=dp), dimension(1:max_v,1:nnode) :: gx, gy

      integer, dimension(1:nnode) :: num_v_node
      integer, dimension(1:max_v,1:nnode) :: v_node
      logical :: do_write


      intent(in) :: do_write
      intent(in) :: nnode, ntr, max_t, max_v
      intent(in) :: mapping, nd, tr_area
      intent(in) :: num_t_node, tr_node
      intent(in) :: num_v_node, v_node

      intent(out) :: gx, gy

!      ---------------
!      local variables
!      ---------------
       real(kind=dp), allocatable, dimension(:) :: pot
       real(kind=dp), allocatable, dimension(:,:) :: E_tr
       real(kind=dp), allocatable, dimension(:,:) :: E_v
       real(kind=dp) :: abserrx, abserry, max_abserrx,max_abserry
       real(kind=dp) :: max_Ex, max_Ey, Ex, Ey
       real(kind=dp) :: area_sum
       real(kind=dp), dimension(2) :: E
       real(kind=dp) :: weight, dp1, dp2, m11,m12,m21,m22
       integer, allocatable, dimension(:) :: ip
       integer :: i,j,k,vj, itr, ipv1,ipv2,ipv3
       integer, dimension(3) :: nodes


      gx(:,:) = 0
      gy(:,:) = 0

      allocate( ip(1:nnode) )
      ip(1:nnode) = 0

      do i=1,nnode
        area_sum = 0
        do j=1,num_t_node(i)
          itr = tr_node(j,i)
          area_sum = area_sum + tr_area(itr)
        enddo

        do j=1,num_v_node(i)
          vj = v_node(j,i)
          ip(vj) = j
        enddo

        do j=1,num_t_node(i)
          itr = tr_node(j, i)
          nodes(1:3) = nd(1:3, itr)
!         ---------------------------------------------------------
!         dp1 = pot( nodes(1) ) - pot( nodes(3) )
!         dp2 = pot( nodes(2) ) - pot( nodes(3) )
!         E_tr(1:2, itr) = -(dp1*mapping(1,1:2,itr) + dp2*mapping(2,1:2,itr))
!         or
!         E_tr(1,itr) = (-dp1)*mapping(1,1) + (-dp2)*mapping(2,1)
!         E_tr(2,itr) = (-dp1)*mapping(1,2) + (-dp2)*mapping(2,2)
!         ---------------------------------------------------------
          m11 = mapping(1,1,itr)
          m12 = mapping(1,2,itr)
          m21 = mapping(2,1,itr)
          m22 = mapping(2,2,itr)

          ipv1 = ip(nodes(1))
          ipv2 = ip(nodes(2))
          ipv3 = ip(nodes(3))

          weight = tr_area(itr)/area_sum
          gx(ipv1,i) = gx(ipv1,i) + weight*(-m11)
          gx(ipv2,i) = gx(ipv2,i) + weight*(-m21)
          gx(ipv3,i) = gx(ipv3,i) + weight*( m11 + m21 )

          gy(ipv1,i) = gy(ipv1,i) + weight*(-m12)
          gy(ipv2,i) = gy(ipv2,i) + weight*(-m22)
          gy(ipv3,i) = gy(ipv3,i) + weight*(m12 + m22)
         enddo
       enddo

       deallocate(ip)

!      ------------
!      double check
!      ------------
       if (idebug >= 1) then
       allocate( pot(1:nnode) )
       allocate( E_tr(1:2,1:ntr) )
       allocate( E_v(1:2,1:nnode) )

       call random_number(pot)
       pot = 2.0*pot - 1.0

       E_tr(:,:) = 0
       E_v(:,:) = 0

!      --------------------------
!      gradient in each triangle
!      --------------------------
       do i=1,nnode
         do j=1,num_t_node(i)
           itr = tr_node(j,i)
           nodes(1:3) = nd(1:3, itr)

           dp1 = pot( nodes(1) ) - pot( nodes(3) )
           dp2 = pot( nodes(2) ) - pot( nodes(3) )
           E_tr(1:2,itr) = -(dp1 * mapping(1,1:2,itr) +                 &
     &                        dp2 * mapping(2,1:2,itr) )

          enddo
       enddo

!      -------------------------------
!      area weighted average at vertex
!      -------------------------------
       do i=1,nnode
         E(1:2) = 0
         area_sum = 0
         do j=1,num_t_node(i)
           itr = tr_node(j,i)
           area_sum = area_sum + tr_area(itr)

           E(1:2) = E(1:2) + tr_area(itr)*E_tr(1:2,itr)
         enddo
         E_v(1:2,i) = E(1:2) / area_sum
       enddo

       do i=1,nnode

!          k = num_v_node(i)
!          Ex = sum( gx(1:k,i) * pot( v_node(1:k,i) ) )
!          Ey = sum( gy(1:k,i) * pot( v_node(1:k,i) ) )

         call calc_gradxy(  nnode, max_v,                               &
     &             num_v_node,v_node,gx,gy,pot,                         &
     &             i, Ex, Ey )


         abserrx = abs(Ex - E_v(1,i))
         max_abserrx = max( max_abserrx, abserrx)
         abserry = abs(Ey - E_v(2,i))
         max_abserry = max( max_abserry, abserry)
       
         max_Ex = max( max_Ex, abs(E_v(1,i)) )
         max_Ey = max( max_Ey, abs(E_v(2,i)) )
       enddo

       if (do_write) then
       write(*,*) 'max_abserrx, max_Ex ', max_abserrx, max_Ex
       write(*,*) 'max_abserry, max_Ey ', max_abserry, max_Ey
       endif


!      --------
!      clean up
!      --------
       deallocate( pot )
       deallocate( E_tr )
       deallocate( E_v)

       endif

       return
       end subroutine gen_gradxy
