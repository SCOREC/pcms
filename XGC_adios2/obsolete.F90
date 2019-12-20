subroutine f_collision_single_sp_org(grid,st, &
      col_f_mat,col_f_ksp,col_f_vecb,col_f_vecx)
  use sml_module
  use grid_class
  use col_module,only : col_pin, col_pout
  use col_f_module
  use f0_module
  use ptl_module
  use eq_module, only : eq_axis_b
  use perf_monitor
  implicit none

#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscmat.h>
#include <finclude/petscksp.h>
#include <finclude/petscpc.h>

  Mat :: col_f_mat
  KSP :: col_f_ksp
  Vec :: col_f_vecb
  Vec :: col_f_vecx

  type(grid_type) :: grid
  integer, intent(in) :: st ! sp%type of ion - specise for collision
  integer :: node, imu, ivp
  real (8) :: den, ti_ev, te_ev, massi, masse, gammac(4)
  real (8) :: f(0:f0_nmu,-f0_nvp:f0_nvp), vol(0:f0_nmu)
  real (8), dimension(0:f0_nmu,-f0_nvp:f0_nvp) :: df, df2, f_corr
  real (8) :: vth, lz, deltay, deltax
  integer :: vpic_ierr
  real (8) :: conv_factor, dum(-f0_nvp:f0_nvp)
  real (8) :: smu_n, pi4bb0vth3_dd
  ! For advanced correction -->
  integer :: is, is2, ie, ie2, i, j, wi
  real (kind=8) :: sigma, val, norm, dist, mu_vol

  df=0.D0
  df2=0.D0
  f_corr=0.D0
  dum=0.D0

  ! for all node, no_electron
  do node=f0_inode1, f0_inode2
     ! check simulation range
     !
     if(col_pin > grid%psi(node) .or. col_pout < grid%psi(node) ) cycle

     ti_ev=f0_t_ev(node,st)
     !te_ev=psn%tempe_ev(node)
     te_ev=ti_ev ! dummay variable for single species
     den=f0_den(node)
     massi=ptl_mass(1)
     !masse=ptl_mass(0)
     masse=massi  ! dummy variable for single species

     ! 1. get lambda gamma
!     call t_startf("COL_F_LAMBDA_GAMMA")
     call col_f_lambda_gamma(den,ti_ev,te_ev,massi,masse,gammac)
!     call t_stopf( "COL_F_LAMBDA_GAMMA")

     ! 2. perform collision
     !*** this routine work only when f0_imu1==0 and f0_imu2==f0_nmu.
     !***** need mu decomposition

     ! unit conversion
     ! f_col is defined in d^3 v space
     ! f_xgc is defined in sqrt(m/(2pi B0 e^2)) B dsmu dv|| space
     ! conversion factor = {m/(4pi e)} sqrt(m/(2pi T)) /sqrt(mu_n)
     ! conversion factor = {m/(4pi e)} sqrt(m/(2pi T)) /sqrt(mu_n)     --> wrong
     ! conversion factor = 1/sqrt(2 * (2pi)^3 * B0 * mu * e^2 / m^3)  or
     !                   = 1/sqrt((2pi)^3 T_ev e^3 / m^3)   --> 8/5/2013 note
     ! Volume elements are
     ! 4 Pi (B/B0) vth^3 sqrt(mu_n) dsqrt(mu_n) dv||_n
     ! vth=sqrt(T/m)

     !conv_factor=(massi/(2D0*sml_2pi*sml_e_charge))*sqrt(massi/(sml_2pi*ti_ev)) !!--> This seems to be wrong
     !conv_factor=1D0/sqrt(2D0*ti_ev*(sml_2pi * sml_e_charge / massi)**3)  !! --> wrong, too!
     conv_factor=0.5D0/sqrt(ti_ev*(sml_2pi * sml_e_charge / massi)**3)

#ifdef COL_F_NAN_CHECK
     if(.not. (f0_f(2,node,1,st) > 1D0 .or. f0_f(2, node, 1, st) < 2D0) ) then
        print *, sml_mype, '] NAN FOUND f0_f ', f0_f(2, node,1, st), node, conv_factor
        stop
     endif
#endif

     ! prepare local f
     do imu=0, f0_nmu
        if(imu==0) then
           smu_n=f0_dsmu/f0_mu0_factor
        else
           smu_n=f0_dsmu*real(imu,8)
        endif

        ! Local f with basic correction:
        ! Simply cut off all negative values
        f(imu,:)=max(f0_f(:,node,imu,st),0.D0)*conv_factor/smu_n

        if (f0_col_advanced_correction) then
          ! Instead of just truncating negative values of f, gaussian smoothing
          do ivp=-f0_nvp,f0_nvp
            if (f0_f(ivp,node,imu,st) .ge. 0.D0) cycle
            ! Determine correction for f0_f
            ! Determine well-behaved value by averaging over neighbours
            ! Stencil (e.g.):
            ! 1 1 1
            ! 1 0 1
            ! 1 1 1
            !rh if (imu==0) then
            !rh   ! More aggressive smoothing at mu=0
            !rh   wi=f0_smooth_df_width+1
            !rh   sigma=f0_smooth_df_sigma*1.5D0
            !rh else
              wi=f0_smooth_df_width
              sigma=f0_smooth_df_sigma
            !rh endif
            is2=max(-f0_nvp,ivp-wi)
            ie2=min(f0_nvp,ivp+wi)
            is=max(f0_imu1,imu-wi)
            ie=min(f0_imu2,imu+wi)
            val=0.D0
            norm=0.D0
            do i=is,ie
              do j=is2,ie2
                if (i==imu .and. j==ivp) cycle  !! Comment this line to take the point with f<0 into account
                if(i==0 .or. i==f0_nmu) then
                  mu_vol=0.5D0
                else
                  mu_vol=1D0
                endif
                dist=sqrt((real(i-imu,8))**2+(real(j-ivp,8))**2)
                val=val+max(f0_f(j,node,i,st),0.D0)*exp(-(dist/sigma)**2)*mu_vol
                norm=norm+exp(-(dist/sigma)**2)*mu_vol
              enddo
            enddo
            f_corr(imu,ivp)=val/norm-f0_f(ivp,node,imu,st)
            f(imu,ivp)=val/norm*conv_factor/smu_n
          enddo
        else
          f_corr(imu,:)=-min(f0_f(:,node,imu,st),0.D0)
        endif

     enddo

     ! Smooth f at high energies before passing it to the collision operator
     ! The hope here is that this is not bad for the result because we average the
     ! particles over several time steps anyway to obtain f0_f
!rh Shifted to f_source in f0module.F90
!rh     if (f0_smooth_df_on) then
!rh        call smooth_f2(f,df2)
!rh        df2=0.D0
!rh     endif

     vth=sqrt(ti_ev*sml_ev2j/massi)       !  SI, v_parallel
     pi4bb0vth3_dd=2D0*sml_2pi*f0_B_B0(node)*vth*vth*vth*f0_dsmu*f0_dvp

     !prepare phase space volume
     smu_n=f0_dsmu/f0_mu0_factor
     vol(0)=0.5D0*pi4bb0vth3_dd*smu_n

     smu_n=f0_dsmu*(real(imu,8)-0.5D0)
     vol(f0_nmu)=0.5D0*pi4bb0vth3_dd*smu_n

     do imu=0+1, f0_nmu-1
        smu_n=f0_dsmu*real(imu,8)
        vol(imu)=pi4bb0vth3_dd*smu_n
     enddo


     !v-grid information

     lz=-f0_vp_max*vth ! minimum
     deltay = f0_dsmu*vth*sqrt(2D0*f0_b_b0(node)) ! perp spacing
     deltax = f0_dvp *vth                  ! parallel spacing - "ES : please confirm this value"

     call t_startf("COL_F_CORE_S")
#ifdef COL_F_NAN_CHECK
     if(.not. (f(1,2) > 1D0 .or. f(1,2) < 2D0) ) then
        print *, sml_mype, '] NAN FOUND f ', f(1,2), node, smu_n
        stop
     endif
#endif
     call col_f_core_s(st, lz, deltay, deltax,&
          f, vol, gammac(1),df, vpic_ierr,    &
          col_f_mat,col_f_ksp,col_f_vecb,col_f_vecx)  ! df is variation of f - return value
     call t_stopf( "COL_F_CORE_S")

!rh     ! Smooth df before multiplying conversion factor
!rh     if (f0_smooth_df_on) then
!rh        call smooth_f2(df,df2)
!rh     endif

     ! apply modification on f0g
     do imu=0, f0_nmu
        if(imu==0) then
           smu_n=f0_dsmu/f0_mu0_factor
        else
           smu_n=f0_dsmu*real(imu,8)
        endif
        ! Update correction
        dum(:)=min(f0_f(:,node,imu,st)+f_corr(imu,:)+(df(imu,:)+df2(imu,:))/conv_factor*smu_n,0.D0)
        f_corr(imu,:)=f_corr(imu,:)-dum(:)+1.D-10
        ! Add result of collisions to separate array for diagnosis
        f0_df0g(:,node,imu,st) = f0_df0g(:,node,imu,st) + df(imu,:)/conv_factor*smu_n
!rh        if (f0_smooth_df_on .and. f0_smooth_df_conserve) then
!rh          f0_df0g3(:,node,imu,st) = df2(imu,:)/conv_factor*smu_n
!rh        endif
        if (f0_f_correction) then
          ! Add correction to f0_f0g -->
          f0_df0g(:,node,imu,st)=f0_df0g(:,node,imu,st) + f_corr(imu,:)
          ! Add correction to separate array for diagnosis -->
          f0_df0g2(:,node,imu,st)=f0_df0g2(:,node,imu,st) + f_corr(imu,:)
        endif
     enddo
  enddo

!rh  if (f0_smooth_df_on) then
!rh     ! Smooth f0_df0g
!rh     call smooth_f(st)
!rh  endif

end subroutine f_collision_single_sp_org

! electron and ion collision
subroutine f_collision_multi_sp_org(grid,st0,st1, &
             col_f_mat,col_f_ksp,col_f_vecb, col_f_vecx)
  use sml_module
  use grid_class
  use col_f_module
  use f0_module
  use ptl_module
  use perf_monitor
  implicit none
#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscmat.h>
#include <finclude/petscksp.h>
#include <finclude/petscpc.h>

  Mat :: col_f_mat
  KSP :: col_f_ksp
  Vec :: col_f_vecb
  Vec :: col_f_vecx

  type(grid_type) :: grid
  integer, intent(in) :: st0, st1  ! st0 : electron sp%type, st1: ion sp%type
  integer :: node, imu, ivp, isp, nsp, st(2)
  real (8) :: den, ti_ev, te_ev, massi, masse, gammac(4)
  real (8), dimension(0:f0_nmu,-f0_nvp:f0_nvp,st0:st1) :: f, df, df2
  real (8), dimension(0:f0_nmu) :: voli, vole
  real (8) :: vth
  real (8) :: lzi, deltayi, deltaxi
  real (8) :: lze, deltaye, deltaxe
  integer :: vpic_ierr
  real (8) :: conv_factor(0:1), dum(-f0_nvp,f0_nvp)
  real (8) :: smu_n, pi4bb0vth3_dd, vthi, vthe
  ! For advanced correction -->
  integer :: is, is2, ie, ie2, i, j, wi
  real (kind=8) :: sigma, val, norm, dist, mu_vol

  if (st0 .ne. 0 .or. st1 .ne. 1) then
    print *,'General species not implemented properly yet in f_collision_multi_sp'
    stop
  endif

  nsp=2
  st(1)=st0
  st(2)=st1

  df=0.D0
  df2=0.D0
  dum=0.D0

  ! for all node, no_electron
  do node=f0_inode1, f0_inode2
     ! check simulation range
#ifndef NO_COL_AXIS
     if(sml_inpsi > grid%psi(node) .or. sml_outpsi < grid%psi(node)) cycle
#else
     if(sml_inpsi+ 0.02 > grid%psi(node) .or. sml_outpsi -0.001 < grid%psi(node)) cycle
#endif
     ti_ev=f0_t_ev(node,st1)
     te_ev=f0_t_ev(node,st0)
     den=f0_den(node)
     massi=ptl_mass(st1)
     masse=ptl_mass(st0)


     ! 1. get lambda gamma
     call t_startf("COL_F_LAMBDA_GAMMA")
     call col_f_lambda_gamma(den,ti_ev,te_ev,massi,masse,gammac)
     call t_stopf( "COL_F_LAMBDA_GAMMA")

     ! 2. perform collision
     !*** this routine work only when f0_imu1==0 and f0_imu2==f0_nmu.
     !***** need mu decomposition

     ! unit conversion
     ! f_col is defined in d^3 v space
     ! f_xgc is defined in sqrt(m/(2pi B0 e^2)) B dsmu dv|| space
     ! conversion factor = {m/(4pi e)} sqrt(m/(2pi T)) /sqrt(mu_n)
     ! Volume elements are
     ! 4 Pi (B/B0) vth^3 sqrt(mu_n) dsqrt(mu_n) dv||_n
     ! vth=sqrt(T/m)

     !rh Are these correct (compare to f_collision_single_sp)?
     !rh wrong ---> conv_factori=(massi/(2D0*sml_2pi*sml_e_charge))*sqrt(massi/(sml_2pi*ti_ev))
     !rh wrong ---> conv_factore=(masse/(2D0*sml_2pi*sml_e_charge))*sqrt(masse/(sml_2pi*te_ev))

     !!conv_factori=(massi/(sml_2pi*sml_e_charge))**1.5D0
     !!conv_factore=(masse/(sml_2pi*sml_e_charge))**1.5D0

     !rh These are correct
     conv_factor(1)=0.5D0/sqrt(ti_ev*(sml_2pi * sml_e_charge / massi)**3)
     conv_factor(0)=0.5D0/sqrt(te_ev*(sml_2pi * sml_e_charge / masse)**3)


     ! prepare local f
     do isp=1,nsp
       do imu=0, f0_nmu
          if(imu==0) then
             smu_n=f0_dsmu/f0_mu0_factor
          else
             smu_n=f0_dsmu*real(imu,8)
          endif

          ! Local f with basic correction:
          ! Simply cut off all negative values
          f(imu,:,st(isp))=max(f0_f(:,node,imu,st(isp)),0.D0)*conv_factor(st(isp))/smu_n

          if (f0_col_advanced_correction) then
            ! Instead of just truncating negative values of f, gaussian smoothing
            do ivp=-f0_nvp,f0_nvp
              if (f0_f(ivp,node,imu,st(isp)) .ge. 0.D0) cycle
              ! Determine correction for f0_f
              ! Determine well-behaved value by averaging over neighbours
              ! Stencil (e.g.):
              ! 1 1 1
              ! 1 0 1
              ! 1 1 1
              !rh if (imu==0) then
              !rh   ! More aggressive smoothing at mu=0
              !rh   wi=f0_smooth_df_width+1
              !rh   sigma=f0_smooth_df_sigma*1.5D0
              !rh else
                wi=f0_smooth_df_width
                sigma=f0_smooth_df_sigma
              !rh endif
              is2=max(-f0_nvp,ivp-wi)
              ie2=min(f0_nvp,ivp+wi)
              is=max(f0_imu1,imu-wi)
              ie=min(f0_imu2,imu+wi)
              val=0.D0
              norm=0.D0
              do i=is,ie
                do j=is2,ie2
                  if (i==imu .and. j==ivp) cycle  !! Comment this line to take the point with f<0 into account
                  if(i==0 .or. i==f0_nmu) then
                    mu_vol=0.5D0
                  else
                    mu_vol=1D0
                  endif
                  dist=sqrt((real(i-imu,8))**2+(real(j-ivp,8))**2)
                  val=val+max(f0_f(j,node,i,st(isp)),0.D0)*exp(-(dist/sigma)**2)*mu_vol
                  norm=norm+exp(-(dist/sigma)**2)*mu_vol
                enddo
              enddo
              f(imu,ivp,st(isp))=val/norm*conv_factor(st(isp))/smu_n
            enddo
          endif

       enddo
!rh Shifted to f_source in f0module.F90
!rh       if (f0_smooth_df_on) then
!rh         ! Smooth f at high energies before passing it to the collision operator
!rh         ! The hope here is that this is not bad for the result because we average the
!rh         ! particles over several time steps anyway to obtain f0_f
!rh         call smooth_f2(f(:,:,st(isp)),df2)
!rh         df2=0.D0
!rh       endif
     enddo


     vthe=sqrt(te_ev*sml_ev2j/masse)       !  SI, v_parallel
     vthi=sqrt(ti_ev*sml_ev2j/massi)

     !prepare phase space volume
     ! - electron
     vth=vthe
     pi4bb0vth3_dd=2D0*sml_2pi*f0_B_B0(node)*vth*vth*vth*f0_dsmu*f0_dvp

     smu_n=f0_dsmu/f0_mu0_factor
     vole(0)=0.5D0*pi4bb0vth3_dd*smu_n

     smu_n=f0_dsmu*(real(imu,8)-0.5D0)
     vole(f0_nmu)=0.5D0*pi4bb0vth3_dd*smu_n

     do imu=0+1, f0_nmu-1
        smu_n=f0_dsmu*real(imu,8)
        vole(imu)=pi4bb0vth3_dd*smu_n
     enddo

     ! - ion
     vth=vthi
     pi4bb0vth3_dd=2D0*sml_2pi*f0_B_B0(node)*vth*vth*vth*f0_dsmu*f0_dvp

     smu_n=f0_dsmu/f0_mu0_factor
     voli(0)=0.5D0*pi4bb0vth3_dd*smu_n

     smu_n=f0_dsmu*(real(imu,8)-0.5D0)
     voli(f0_nmu)=0.5D0*pi4bb0vth3_dd*smu_n

     do imu=0+1, f0_nmu-1
        smu_n=f0_dsmu*real(imu,8)
        voli(imu)=pi4bb0vth3_dd*smu_n
     enddo


     !v-grid information
     !ion
     vth=vthi       !  SI, v_parallel


     lzi=-f0_vp_max*vth ! minimum
     deltayi = f0_dsmu*vth*sqrt(2D0*f0_b_b0(node)) ! perp spacing
     deltaxi = f0_dvp*vth                  ! parallel spacing - "ES : please confirm this value"

     !electron
     vth=vthe       !  SI, v_parallel


     lze=-f0_vp_max*vth ! minimum
     deltaye = f0_dsmu*vth*sqrt(2D0*f0_b_b0(node)) ! perp spacing
     deltaxe = f0_dvp*vth                  ! parallel spacing - "ES : please confirm this value"



     ! collision in action
     ! ES : please exchange parameters for col_f_core_m below
     !dfi and dfe have been added
     call t_startf("COL_F_CORE_M")
     call col_f_core_m(st1, lzi, deltayi, deltaxi, f(:,:,st1), voli, &
                       st0, lze, deltaye, deltaxe, f(:,:,st0), vole, &
                       gammac, df(:,:,st1), df(:,:,st0), node, vpic_ierr, &
                       col_f_mat, col_f_ksp, col_f_vecb, col_f_vecx)
     call t_stopf( "COL_F_CORE_M")

     do isp=1,nsp
       do imu=0, f0_nmu
          if(imu==0) then
             smu_n=f0_dsmu/f0_mu0_factor
          else
             smu_n=f0_dsmu*real(imu,8)
          endif
          f0_df0g(:,node,imu,st(isp)) = f0_df0g(:,node,imu,st(isp)) + df(imu,:,st(isp))/conv_factor(st(isp))*smu_n
       enddo
     enddo

  enddo

end subroutine f_collision_multi_sp_org


