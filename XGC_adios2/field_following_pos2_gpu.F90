attributes(device) &
#ifndef POS2_ORG
subroutine field_following_pos2_gpu(x_org,phi_org,phi_dest,x_dest)
  use precision_mod_gpu
  implicit none
  real (kind=work_p), intent(in) :: x_org(2),phi_org, phi_dest
  real (kind=work_p), intent(out) :: x_dest(2)
  real (kind=work_p) :: phi, x(2),dphi
  real (kind=work_p) :: hh,h6,dx1(2),dx2(2),dx3(2),x_tmp(2)
  ! FIXME : This *should* be a parameter, but it generates considerably slower
  !         GPU code, for some reason.
  integer :: sml_ff_step=1
  integer :: i

  dphi=(phi_dest-phi_org)/real(sml_ff_step)

  ! 0 step
  phi=phi_org
  x=x_org

  do i=1, sml_ff_step
     ! get first derivative
     call derivs_gpu(x,phi,dx1)

     if( .false. ) then ! first order calculation
        x = x + dx1*dphi
     else if( .true. ) then  ! second order calculation - rk2        
        ! obtain mid point
        hh=dphi*0.5D0

        x_tmp = x + dx1*hh

        ! get new derivative
        call derivs_gpu(x_tmp,phi+hh,dx2)

        ! advance one step using mid-point derivative
        x = x + dx2*dphi
     else
        ! 4th Order Calculation - rk4
        !
        hh=dphi*0.5D0
        h6=dphi/6D0
        ! derivative 1 (from x)- obtained already
        x_tmp=x + hh*dx1             ! yt=y+hh*dydx      : yt -> x_tmp
        !
        ! derivative 2 (from x_tmp) 
        call derivs_gpu(x_tmp,phi+hh,dx2)! dyt from yt       : dyt -> dx2
        x_tmp=x + hh*dx2             ! yt=y+hh*dyt       : yt -> x_tmp
        !
        ! derivative 3 (from x_tmp)
        call derivs_gpu(x_tmp,phi+hh,dx3)! dym from yt       : dym -> dx3 
        x_tmp=x + dphi*dx3           ! yt=y + h*dym      : yt -> x_tmp
        dx3 = dx2 + dx3              ! dym = dyt + dym   : dym -> dx3 , dyt -> dx2
        !
        ! derivative 4 (from x_tmp)
        call derivs_gpu(x_tmp,phi+dphi,dx2)    ! dyt from yt       : dyt -> dx2, yt -> x_tmp        
        x = x + h6 * (dx1 + dx2 + 2D0*dx3) ! yout = y + h6* (dydx+dyt+2D0*dym)  
     endif
     phi=phi+dphi
  enddo
  x_dest=x

end subroutine field_following_pos2_gpu

#else

subroutine field_following_pos2_gpu(x_org,phi0,phi,x_dest)
  use precision_mod_gpu
  implicit none
  real (kind=work_p), intent(in) :: x_org(2), phi,phi0
  real (kind=work_p), intent(out) :: x_dest(2)
  real (kind=work_p) :: B(2),Bphi, x_mid(2), dphi
  real (kind=work_p) :: hh,h6,dx1(2),dx2(2),dx3(2),x_tmp(2)


  dphi=phi-phi0

  call bvec_interpol_gpu(x_org(1),x_org(2),phi0,b(1),b(2),bphi)

  if( .false. ) then ! first order calculation
     x_dest(:)=x_org(:) + b(:)/bphi*(x_org(1)*dphi)
  else if( .true. ) then
     ! second order calculation - rk2
     ! obtain mid point
     hh=dphi*0.5D0
     x_mid(:)=x_org(:) + b(:)/bphi*(x_org(1)*hh)

     ! get new derivative
     call bvec_interpol_gpu(x_mid(1),x_mid(2),phi0+hh,b(1),b(2),bphi)

     ! advance one step using mid-point derivative
     x_dest(:)=x_org(:) + b(:)/bphi*(x_mid(1)*dphi)
  else
     ! 4th Order Calculation - rk4
     !
     hh=dphi*0.5D0
     h6=dphi/6D0
     ! derivative 1 (from x_org)- obtained already
     dx1=b/bphi*x_org(1)     ! dydx from y       : y-> x_org , dx1 -> dydx
     x_tmp=x_org + hh*dx1    ! yt=y+hh*dydx      : yt -> x_tmp
     !
     ! derivative 2 (from x_tmp) 
     call bvec_interpol_gpu(x_tmp(1),x_tmp(2),phi0+hh,b(1),b(2),bphi)
     dx2=b/bphi*x_tmp(1)     ! dyt from yt       : dyt -> dx2
     x_tmp=x_org + hh*dx2    ! yt=y+hh*dyt       : yt -> x_tmp
     !
     ! derivative 3 (from x_tmp)
     call bvec_interpol_gpu(x_tmp(1),x_tmp(2),phi0+hh,b(1),b(2),bphi)
     dx3=b/bphi*x_tmp(1)     ! dym from yt       : dym -> dx3 
     x_tmp=x_org + dphi*dx3     ! yt=y + h*dym      : yt -> x_tmp
     dx3 = dx2 + dx3         ! dym = dyt + dym   : dym -> dx3 , dyt -> dx2
     !
     ! derivative 4 (from x_tmp)
     call bvec_interpol_gpu(x_tmp(1),x_tmp(2),phi0+dphi,b(1),b(2),bphi)
     dx2=b/bphi*x_tmp(1)    ! dyt from yt       : dyt -> dx2, yt -> x_tmp
     x_dest=x_org+h6* ( dx1 + dx2 + 2D0*dx3)     ! yout = y + h6* (dydx+dyt+2D0*dym)  : yout -> x_dest, dydx -> dx1, dyt-> dx2 , dym -> dx3
  endif

  !x_dest(:)=x_org(:)
end subroutine field_following_pos2_gpu
#endif


