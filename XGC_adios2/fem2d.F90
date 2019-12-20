
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     helm2dElem 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      subroutine helm2dElem(alpha,beta,ul,xl,ss,pp,isw)
!
!     Two dimensional (plane) Linear Thermal Element: -alpha del^2 u + beta u
!
!-----[--.----+----.----+----.-----------------------------------------]        
!
!     Inputs:
!     alpha      - -alpha del^2 u + beta u
!     beta       - 
!     ul(3)      - potential at nodes
!     xl(2,3)    - 2D coordinate of nodes
!     ss(3,2)    - isw==3: stores B field
!     isw        - 1, matrix (only); 2, helm op (only); 3, b dot grad op.
!
!     Outputs:
!     ss     - isw 1: element matrix. isw 3: area of nodes
!     pp(3)  - f(u) (isw== 2 or 3)
!-----[--.----+----.----+----.-----------------------------------------]
      implicit  none
      integer,intent(in)::isw
      real*8,intent(in)::ul(*),alpha,beta,xl(2,*)
      real*8::ss(3,3) ! isw: 1,tang mat (out); 2 not used; 3 b(3,2) (in), area(3,1) (out)
      real*8,intent(out):: pp(3)
      !
      integer :: i,j,l,lint,i1,j1
      real*8 :: xsj,a1,a2,a4,shj,ttot

      real*8 :: dd(2,2),shp(3,3),el(4,3),gradpot(2),pot,bb(2),b(3,2)

      if( isw == 1 ) ss = 0d0 ! 2: not used; 3: in/out
      if( isw == 2 .or. isw == 3 ) pp = 0.d0 ! output vector, init

      if (alpha/=0.d0 .or. isw == 3) then
         ! copy B into temparary B from ss, will clobber ss with area
         if( isw == 3 ) then
            do i=1,3
               b(i,1) = ss(i,1)   ! grap B
               b(i,2) = ss(i,2)
               ss(i,1) = 0d0 ! return area
            enddo
         end if
         ! one point integration for del^2
         l    =  1
         call new_tint2d(l,lint,el) ! area coordiante for point 
         if (lint/=1) stop 'lint != 1'
         do l = 1,lint
            !     get FE shape functions
            call new_trishp(el(1,l),xl,xsj,shp)
            xsj = xsj*el(4,l) ! area of point in this elem
            ! sets material tensor 'dd' with alpha (isw==1). computes grad (isw==2 & 3)
            call new_thfx2d(alpha,ul,shp,gradpot,pot,dd)
            if( isw==1 .or. isw==2 ) then
               !     Normal element evaluation
               do j = 1,3 ! nel
                  a1 = (dd(1,1)*shp(1,j) + dd(1,2)*shp(2,j))*xsj
                  a2 = (dd(2,1)*shp(1,j) + dd(2,2)*shp(2,j))*xsj
                  !a4 = beta*shp(3,j)*xsj 
                  !     Consistent rate and conductivity terms
                  if( isw == 1 ) then
                     do i = 1,3
                        ss(i,j) = ss(i,j) + a1*shp(1,i) + a2*shp(2,i) !+ a4*shp(3,i)
                     end do
                  else ! if( isw == 2 ) then
                     pp(j) = pp(j) + a1*gradpot(1) + a2*gradpot(2) !+ a4*pot
                  endif
               end do
            else if( isw == 3 ) then 
               stop 'isw == 3 not used???'
               ! b dot grad, return B dot grad(U) in pp, should this use 3 point integration???
               do j = 1,3
                  bb(1) = b(j,1)*shp(3,j)*xsj ! B is stored in first two columns of ss (yuck)
                  bb(2) = b(j,2)*shp(3,j)*xsj
                  pp(j)   = pp(j) + bb(1)*gradpot(1) + bb(2)*gradpot(2)  
                  ss(j,1) = ss(j,1) + xsj*shp(3,j) ! return area in ss (yuck)
               end do
            else
               stop 'bad isw'
            endif
         end do
      end if
      ! add mass 
      if (beta/=0.d0 .and. isw /= 3) then
         l     = -3                         ! 3-pt formula
         call new_tint2d(l,lint,el) ! get area coordinates in 'el'
         if (lint/=3) stop 'lint != 3'
         do l = 1,lint
            call new_trishp(el(1,l),xl,xsj,shp) ! uses el(1:3,l)
            xsj = el(4,l)*xsj ! area of point in this elem/integration point
            if (isw == 2) then 
               pot = 0.0d0
               do i = 1,3 ! nel
                  pot = pot + shp(3,i)*ul(i)
               end do
            end if
            do j = 1,3 ! nel
               a4 = beta*shp(3,j)*xsj
               if( isw == 1 ) then
                  do i = 1,3 ! nel
                     ss(i,j) = ss(i,j) + a4*shp(3,i)
                  end do
               else ! if( isw == 2 ) then
                  pp(j) = pp(j) + a4*pot
               end if
            end do
         end do ! l
      end if
      end subroutine helm2dElem
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     bdotgradelem
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      subroutine bdotgradelem(ul,xl,bb,pp,area)
!
!     Two dimensional (plane) Linear Thermal Element B dot grad operator
!
!-----[--.----+----.----+----.-----------------------------------------]
!     Inputs:
!     ul(3)      - potential at nodes
!     xl(2,3)    - 2D coordinate of nodes
!     bb(3,2)     - isw==3: stores B field
!
!     Outputs:
!     pp(3)   - b dot grad (u) 
!     area(3) - area of nodes
!-----[--.----+----.----+----.-----------------------------------------]
      implicit  none
      real*8,intent(in)::ul(3),xl(2,3),bb(3,2) 
      real*8,intent(out):: pp(3),area(3)
      !
      integer :: l,lint,j
      real*8 :: xsj

      real*8 :: dum(2,2),shp(3,3),el(4,3),gradpot(2),pot,bt(2),alpha

      pp = 0.d0 ! output vector, init
      area = 0.d0
      alpha = 1.d300 ! dummy
      ! 3 point integration
      l    =  -3
      call new_tint2d(l,lint,el) ! area coordiante for point 
      do l = 1,lint
         !     get FE shape functions
         call new_trishp(el(1,l),xl,xsj,shp)
         xsj = xsj*el(4,l) ! area of point in this elem
         ! sets material tensor 'dum' with alpha
         call new_thfx2d(alpha,ul,shp,gradpot,pot,dum)
         ! return b dot grad in pp, and area
         do j = 1,3
            bt(1) = bb(j,1)*shp(3,j)*xsj
            bt(2) = bb(j,2)*shp(3,j)*xsj
            pp(j)   = pp(j) + bt(1)*gradpot(1) + bt(2)*gradpot(2)
            area(j) = area(j) + xsj*shp(3,j)
         end do
      end do
    end subroutine bdotgradelem
!-----[--.----+----.----+----.-----------------------------------------]
!  new_thfx2d
!-----[--.----+----.----+----.-----------------------------------------]
      subroutine new_thfx2d(alpha,ul,shp,gradpot,pot,dd)
!     Compute thermal gradient and flux
      implicit  none
      real*8,intent(in)::ul(3),alpha,shp(3,*)
      real*8,intent(out)::dd(2,2),gradpot(2),pot
      !
      integer   i
      !
      gradpot(1) = 0.0d0
      gradpot(2) = 0.0d0
      pot     = 0.0d0
      do i = 1,3
         gradpot(1) = gradpot(1) + shp(1,i)*ul(i)
         gradpot(2) = gradpot(2) + shp(2,i)*ul(i)
         pot     = pot     + shp(3,i)*ul(i)
      end do
      !
      dd(1,1) = alpha
      dd(2,2) = alpha
      dd(1,2) = 0.d0    
      dd(2,1) = 0.d0    

      end subroutine new_thfx2d

!-----[--.----+----.----+----.-----------------------------------------]
      subroutine new_tint2d(l,lint,el)
!-----[--.----+----.----+----.-----------------------------------------]
!     Purpose: Set gauss points and weights for triangular elements
      
!     Inputs:
!     l       - Number of gauss points indicator
      
!     Outputs:
!     lint    - Total number of points
!     el(4,*) - Area coordinate points and weights for quadrature
!-----[--.----+----.----+----.-----------------------------------------]
      implicit  none
            
      integer   l,lint
      real*8    el(4,*),r0,r1,r2,ww,eta
      
      data ww, eta / 0.3333333333333333d0 , 0.1666666666666667d0 /
      
!     1-point gauss integration
      
      if(l==1) then
         el(1,1) = ww
         el(2,1) = ww
         el(3,1) = ww
         el(4,1) = 1.d0
         lint    = 1
      elseif(l==-3) then
        el(1,1) = 1.0d0 - ww
        el(2,1) = eta
        el(3,1) = eta
        el(4,1) = ww

        el(1,2) = eta
        el(2,2) = 1.0d0 - ww
        el(3,2) = eta
        el(4,2) = ww

        el(1,3) = eta
        el(2,3) = eta
        el(3,3) = 1.0d0 - ww
        el(4,3) = ww

        lint    = 3
      else
        write(  *,2000) l
        lint    = -1
      endif
      
!     Format

 2000 format(' *ERROR* NEW_TINT2D: Wrong quadrature, l =',i3)
      
      end subroutine new_tint2d

!-----[--.----+----.----+----.-----------------------------------------]
      subroutine new_trishp(el,xl,xsj,shp)
!-----[--.----+----.----+----.-----------------------------------------]
!      Purpose: Triangular shape function routine

!     Type:  |iord| = 1:  Linear  three-node
!     |iord| = 2:  Quadratic six-node
!     |iord| = 3:  Quadratic seven-node
!     |iord| = 4:  Quadratic + 3 bubbles (Zienkiewicz/Lefebre)
      
!     iord  > 0:  Mid-side and center node are  global  coords.
!     iord  < 0:  Mid-side and center node heirarchical coords.
      
!     Inputs:
!     el(3)     - Area coordinates for point
!     xl(ndm,*) - Nodal coordinates for element
!     ndm       - Spatial dimension of mesh
!     iord      - Order of shape functions (see above)
      
!     Outputs:
!     xsj       - Jacobian determinant at point
!     shp(3,*)  - Shape functions and derivatives at point
!-----[--.----+----.----+----.-----------------------------------------]
      implicit  none
      
      integer   i
      real*8    xsj,xsjr
      real*8    x1,x2,x3,y1,y2,y3
      real*8    el(3),xl(2,*),shp(3,*)
      
!     Form Jacobian terms
      
      x1 = xl(1,1)
      x2 = xl(1,2)
      x3 = xl(1,3)
      
      y1 = xl(2,1)
      y2 = xl(2,2)
      y3 = xl(2,3)
 
      xsj  = x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2) ! det(A)
      xsjr = 1.0d0
      if(xsj.ne.0.0d0) then
         xsjr = 1.0d0/xsj
      else 
         stop 'xsj==0'
      endif
      xsj  = 0.5d0*xsj ! area (out)

!     Specify shape functions and their derivatives
      
      shp(1,1) = (y2-y3)*xsjr
      shp(2,1) = (x3-x2)*xsjr
      shp(3,1) = el(1)
      
      shp(1,2) = (y3-y1)*xsjr
      shp(2,2) = (x1-x3)*xsjr
      shp(3,2) = el(2)
      
      shp(1,3) = (y1-y2)*xsjr
      shp(2,3) = (x2-x1)*xsjr
      shp(3,3) = el(3)
      
    end subroutine new_trishp
      
