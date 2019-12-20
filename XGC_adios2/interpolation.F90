!****************************************************************************
! retrun I
!
! first created : 2000/08/26
! last modified : 2002/02/08
!adopted from xorbit
!
!****************************************************************************
real (kind=8)  function I_interpol(in_psi, ideriv,region)
    use EZspline_obj
    use EZspline
    use itp_module, only : spl_psi,itp_min_psi,itp_max_psi
    use eq_module, only : eq_x_psi
    use sml_module, only : sml_bt_sign
#ifdef USE_ONE_D_I_CUB_MOD
    use one_d_cub_mod
#endif

    implicit none
    integer , intent(in) :: ideriv,region
    real (kind=8) , intent(in) :: in_psi
    real (kind=8)  :: psi
    integer :: ier
    real (kind=8) :: r8value


!    sign=-1D0 !-1D0 : cocurrent, 1D0 :counter current 
!    sml_bt_sign can be changed in setup.f90 2002/02/08
    
    if(region == 3 ) then
       
       psi=min(eq_x_psi,itp_max_psi) ! for itp_max_psi < eq_x_psi case 2002/01/22
       if(ideriv == 0) then
#ifdef USE_ONE_D_I_CUB_MOD
          call I_interpol_wo_pspline(psi, ideriv, r8value)
#else
          call EZspline_interp(spl_psi,psi,r8value,ier)
          call EZspline_error(ier)
#endif
          I_interpol=sml_bt_sign*r8value
       else
          I_interpol=0
       endif
       
    else
       
       psi = in_psi
       if(psi < itp_min_psi) then
          if(psi < itp_min_psi - 1D-4) then
             print *, 'psi range exceeded',psi
             !call err_count
          endif
	  psi=itp_min_psi
       elseif(psi > itp_max_psi) then
          psi=itp_max_psi ! I is constant outside of edge
       endif
#ifdef USE_ONE_D_I_CUB_MOD
       call I_interpol_wo_pspline(psi, ideriv, r8value)
#else
       call EZspline_derivative(spl_psi,ideriv,psi,r8value,ier)
       call EZspline_error(ier)
#endif
       I_interpol=sml_bt_sign*r8value
       
    endif
end function
!**************************************************************
! B-field interpolation
!**************************************************************

real (kind=8) function b_interpol(r,z,phi)
  use sml_module, only: sml_cylindrical
  use eq_module
#ifdef USE_BICUB_MOD
  use bicub_mod
#endif
  implicit none
  real (kind=8) , intent(in) :: r,z,phi
  real (kind=8)              :: psi,dpsi_dr,dpsi_dz,br,bz,bphi,fi
  real (kind=8)              :: ripp,dripp_dr,dripp_dz
  real (kind=8) , external   :: I_interpol,psi_interpol
  real (kind=8)              :: r0

  r0=eq_axis_r

#ifdef USE_BICUB_MOD
  call bicub_interpol(psi_bicub,r,z,psi,dpsi_dr,dpsi_dz)
#else
  psi     = psi_interpol(r,z,0,0)
  dpsi_dr = psi_interpol(r,z,1,0)
  dpsi_dz = psi_interpol(r,z,0,1)
#endif
  
  if(psi<eq_x_psi .AND. z<eq_x_z) then
     fi=I_interpol(psi,0,3)
  else
     fi=I_interpol(psi,0,1)
  endif

  if (sml_cylindrical) then
    br=- dpsi_dz / r0   ! sign ignored   sml_bp_sign
    bz= dpsi_dr / r0
    bphi=fi / r0
  else
    br=- dpsi_dz / r   ! sign ignored   sml_bp_sign
    bz= dpsi_dr / r
    bphi=fi / r
  endif
  
  b_interpol= sqrt(br**2+bz**2+bphi**2)

end function b_interpol

real (kind=8) function b_interpol_sym(r,z)  ! axi-symmetric b interpolation
  use sml_module, only: sml_cylindrical
  use eq_module
#ifdef USE_BICUB_MOD
  use bicub_mod
#endif
  implicit none
  real (kind=8) , intent(in) :: r,z
  real (kind=8)              :: psi,dpsi_dr,dpsi_dz,br,bz,bphi,fi
  real (kind=8)              :: ripp,dripp_dr,dripp_dz
  real (kind=8) , external   :: I_interpol,psi_interpol
  real (kind=8)              :: r0

  r0=eq_axis_r

#ifdef USE_BICUB_MOD
 call bicub_interpol(psi_bicub,r,z,psi,dpsi_dr,dpsi_dz)
#else
  psi     = psi_interpol(r,z,0,0)
  dpsi_dr = psi_interpol(r,z,1,0)
  dpsi_dz = psi_interpol(r,z,0,1)
#endif
  
  if(psi<eq_x_psi .AND. z<eq_x_z) then
     fi=I_interpol(psi,0,3)
  else
     fi=I_interpol(psi,0,1)
  endif
  
  if (sml_cylindrical) then
    br=- dpsi_dz / r0   ! sign ignored   sml_bp_sign
    bz= dpsi_dr / r0
    bphi=fi / r0
  else
    br=- dpsi_dz / r   ! sign ignored   sml_bp_sign
    bz= dpsi_dr / r
    bphi=fi / r
  endif

  b_interpol_sym= sqrt(br**2+bz**2+bphi**2)
  
end function b_interpol_sym

! return B_phi
subroutine bphi_interpol_rzpsi(r,z,phi,psi,bphi)
  use sml_module, only: sml_cylindrical
  use eq_module 
  implicit none
  real (kind=8) , intent(in) :: r,z,phi,psi
  real (kind=8)              :: bphi,fi
  real (kind=8)              :: ripp, dripp_dr, dripp_dz
  real (kind=8) , external   :: I_interpol
  real (kind=8)              :: r0

  r0=eq_axis_r

  if(psi<eq_x_psi .AND. z<eq_x_z) then
     fi=I_interpol(psi,0,3)
  else
     fi=I_interpol(psi,0,1)
  endif

  if (sml_cylindrical) then
    bphi=fi / r0
  else
    bphi=fi / r
  endif
  
end subroutine bphi_interpol_rzpsi


subroutine bvec_interpol(r,z,phi,br,bz,bphi)
  use eq_module 
  use sml_module, only : sml_bp_sign, sml_cylindrical
#ifdef USE_BICUB_MOD
  use bicub_mod
#endif
  implicit none
  real (kind=8) , intent(in) :: r,z,phi
  real (kind=8)              :: psi,dpsi_dr,dpsi_dz,br,bz,bphi,fi
  real (kind=8)              :: ripp, dripp_dr, dripp_dz
  real (kind=8) , external   :: I_interpol,psi_interpol
  real (kind=8)              :: r0

  r0=eq_axis_r

#ifdef USE_BICUB_MOD
 call bicub_interpol(psi_bicub,r,z,psi,dpsi_dr,dpsi_dz)
#else
 psi     = psi_interpol(r,z,0,0)
 dpsi_dr = psi_interpol(r,z,1,0)
 dpsi_dz = psi_interpol(r,z,0,1)
#endif
 
 if(psi<eq_x_psi .AND. z<eq_x_z) then
    fi=I_interpol(psi,0,3)
 else
    fi=I_interpol(psi,0,1)
 endif

 if (sml_cylindrical) then
   br=- dpsi_dz / r0 * sml_bp_sign
   bz= dpsi_dr / r0  * sml_bp_sign
   bphi=fi / r0
 else
   br=- dpsi_dz / r * sml_bp_sign
   bz= dpsi_dr / r  * sml_bp_sign
   bphi=fi / r
 endif
 
end subroutine bvec_interpol


real (kind=8) function psi_interpol_pspline(r,z,r_der,z_der)
  use ITP_module
  use EZspline_obj
  use EZspline
  
  implicit none
  real (kind=8), intent(in) :: r, z
  integer , intent(in) :: r_der, z_der
  
  integer :: ier
  real (kind=8) :: r8value
  
  
#if defined(USE_EXPLICIT)
  call EZspline_interp(spl(r_der,z_der),r,z,r8value,ier)
#else
  call EZspline_derivative(spl(0,0),r_der,z_der,r,z,r8value,ier)
#endif
  call EZspline_error(ier)
  
  if (r_der==0 .and. z_der==0) then
     psi_interpol_pspline = max(1D-99,r8value)
  else
     psi_interpol_pspline = r8value
  endif
end function psi_interpol_pspline


!****************************************************************************
! retrun psi and dpsi as a function of r, z
!
! first created : 2000/10/19
! last modified : 2000/10/19
!****************************************************************************

real (kind=8) function psi_interpol(r_in,z_in,r_der,z_der)
  use itp_module
#ifdef USE_BICUB_MOD
  use bicub_mod
#endif

#ifdef PROFILE_DEBUG
  use perf_mod, only : t_startf, t_stopf
#endif

  implicit none
  real (kind=8), intent(in) :: r_in, z_in
  integer , intent(in) :: r_der, z_der
!  real (kind=8) :: r,z
  real (kind=8) ,external:: psi_interpol_pspline
#ifdef USE_BICUB_MOD
  real (8) :: dpsi_dr, dpsi_dz, psi
#endif

#ifdef PROFILE_DEBUG
  character*17, save ::  str(0:3,0:3)
  logical, save :: isfirst = .true.

  if (isfirst) then
  isfirst = .false.
  str(0,0) = 'psi_interpol(0,0)'
  str(1,0) = 'psi_interpol(1,0)'
  str(2,0) = 'psi_interpol(2,0)'
  str(3,0) = 'psi_interpol(3,0)'

  str(0,1) = 'psi_interpol(0,1)'
  str(1,1) = 'psi_interpol(1,1)'
  str(2,1) = 'psi_interpol(2,1)'
  str(3,1) = 'psi_interpol(3,1)'

  str(0,2) = 'psi_interpol(0,2)'
  str(1,2) = 'psi_interpol(1,2)'
  str(2,2) = 'psi_interpol(2,2)'
  str(3,2) = 'psi_interpol(3,2)'
 
  str(0,3) = 'psi_interpol(0,3)'
  str(1,3) = 'psi_interpol(1,3)'
  str(2,3) = 'psi_interpol(2,3)'
  str(3,3) = 'psi_interpol(3,3)'
  endif
 

  call t_startf( str(r_der,z_der) )

#endif

#ifdef USE_BICUB_MOD
  call bicub_interpol(psi_bicub,r_in,z_in,psi,dpsi_dr,dpsi_dz)
  if(r_der==0 .and. z_der==0) then
     psi_interpol=psi
  elseif(r_der==1 .and. z_der==0) then
     psi_interpol=dpsi_dr
  elseif(r_der==0 .and. z_der==1) then
     psi_interpol=dpsi_dz
  else
     print *, 'Do not use second derivative of psi_interpol with USE_BICUB_MOD.'
     print *, 'Replace it with psi_interpol_bicub(r,z,psi,dpsi_dr,dpsi_dz,d2psi_drdz,d2psi_dr2,d2psi_dz2).'
     stop 
  endif

  if (r_der==0 .and. z_der==0) then
     psi_interpol = max(1D-99, psi_interpol)
  endif
#else
  psi_interpol=psi_interpol_pspline(r_in,z_in,r_der,z_der)
#endif
  
#ifdef PROFILE_DEBUG
 call t_stopf( str(r_der,z_der) )
#endif
  
end function psi_interpol
  



!****************************************************************************
! initialization of interpolation routines
!
! first created : 2000/08/26
! last modified : 2003/01/28
! adopted from xorbit
!****************************************************************************
subroutine init_interpolation
  use eq_module
  use itp_module
  use EZspline_obj
  use EZspline
#if defined(USE_BICUB_MOD)
    use bicub_mod
#endif
#ifdef USE_ONE_D_I_CUB_MOD
    use one_d_cub_mod
#endif

    !use numerical_libraries
    implicit none
    integer :: i,j, max_mpsi, max_nearx_mr , max_mr ,region
    real (kind=8)  :: psi,theta,r,z

    real (kind=8) ,external:: psi_interpol_pspline
    integer, parameter :: idebug = 0
    integer, dimension(2) :: BCS1, BCS2
    integer :: ier,  ir,iz,  r_der, z_der
    real (kind=8), dimension(eq_mr,eq_mz) :: rzgrid
    real (kind=8) :: rr,zz, abserr,relerr 
    real (kind=8) :: r8value, maxrelerr, maxabserr
    real (kind=8), parameter  ::  one = 1.0d0, zero = 0.0d0

    max_mpsi=eq_mpsi
    max_mr=eq_mr

    !set valid psi, theta range for interpolation
    itp_min_psi=eq_psi_grid(1)
    itp_max_psi=eq_psi_grid(itp_mpsi)
        
    !------------------------------------------------------------

    !initializing interpolation routines
    
    ! for psi grid 
    !common - psi
!    print *,eq_psi_grid(1:2)
!    print *, 'psi_grid'
!    do i=1, eq_mpsi
!       write(111,*) i, eq_psi_grid(i)
!    enddo
!    close(111)
    

    BCS1 = 0 ! not a knot
    BCS2 = 0 ! not a knot
    call EZspline_init(spl_psi,eq_mpsi,BCS1,ier)
    call EZspline_error(ier)

    spl_psi%x1 = eq_psi_grid

    call EZspline_setup(spl_psi,eq_I,ier)
    call EZspline_error(ier)

    !all rz mode interpolation

    ! for r-z grid  -- additional term added for numerical safety.
!    print *, 'rz grid'
    itp_min_r=eq_rgrid(1)! *(1D0-1D-5) + 1D-5*eq_rgrid(2)
    itp_max_r=eq_rgrid(eq_mr)! *(1D0-1D-5)+ 1D-5 *eq_rgrid(eq_mr-1)
    itp_min_z=eq_zgrid(1) !*(1D0-1D-5) + 1D-5*eq_zgrid(2)
    itp_max_z=eq_zgrid(eq_mz) !*(1D0-1D-5) + 1D-5*eq_zgrid(eq_mz-1)


    BCS1 = 0 ! not a knot
    BCS2 = 0 ! not a knot

     call EZspline_init(spl(0,0),eq_mr,eq_mz,BCS1,BCS2,ier)
     call EZspline_error(ier)
     spl(0,0)%x1 = eq_rgrid
     spl(0,0)%x2 = eq_zgrid


    !debug
    !do i=1, eq_mr
    !   do j=1, eq_mz
    !      write(800,800) eq_rgrid(i),eq_zgrid(j),eq_psi_rz(i,j)
    !   enddo
    !   write(800,*) ' '
    !enddo
    !do i=1, eq_mr+itp_korder_rz
    !   write(800,900) itp_r_KNOT(i)
    !enddo
    !do i=1, eq_mz+itp_korder_rz
    !   write(800,900) itp_z_KNOT(i)
    !enddo
900 format (e19.13)    
800 format (e19.13,' ',e19.13,' ',e19.13)

    !for psi

    if (idebug.ge.1) then
       write(*,*) 'mindr ',minval(eq_rgrid(2:eq_mr)-eq_rgrid(1:(eq_mr-1)))
       write(*,*) 'maxdr ',maxval(eq_rgrid(2:eq_mr)-eq_rgrid(1:(eq_mr-1)))
       write(*,*) 'mindz ',minval(eq_zgrid(2:eq_mz)-eq_zgrid(1:(eq_mz-1)))
       write(*,*) 'maxdz ',maxval(eq_zgrid(2:eq_mz)-eq_zgrid(1:(eq_mz-1)))
    endif
    
    call EZspline_setup(spl(0,0),eq_psi_rz,ier)
    call EZspline_error(ier)

#ifdef USE_BICUB_MOD
      call setup_bicub(psi_bicub,eq_mr,eq_rgrid,eq_mz,eq_zgrid,psi_interpol_pspline)
#endif
#ifdef USE_ONE_D_I_CUB_MOD
      call setup_psi_one_d_cub(eq_mpsi, eq_psi_grid)
#endif
end subroutine

