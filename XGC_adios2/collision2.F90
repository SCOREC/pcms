! momentum conserving collision
! based on Plasma Phys. Control. Fusion, W X Wang, N nakjima, M Okamoto 41, 1091 (1999)
subroutine conserving_collision(sp,istep,flag)
  use sml_module
  use col_module
  use ptl_module
  use omp_module, only : split_indices
  use eq_module,only : eq_axis_r, eq_axis_z, is_rgn12, eq_x_psi, eq_x_z, eq_x2_z
  use perf_monitor
  implicit none
  type(species_type):: sp
  integer, intent(in) :: flag
  integer :: ith, i_beg(sml_nthreads), i_end(sml_nthreads)

  integer :: isp, istep
  integer :: i, j, l, iflag2

  real (8), allocatable :: col_c_org_vsum(:,:,:,:,:) 
  real (8), allocatable :: col_c_new_vsum(:,:,:,:,:) 
  real (8), allocatable :: col_c_delta(:,:,:,:) 
  real (8), allocatable :: col_c_matrix(:,:,:,:,:), dumv(:,:,:,:)
  real (8), allocatable :: org_avg(:,:,:,:,:), dum_avg(:,:,:,:)

  real (8) :: dt, c_m_sp, c2_2m_sp
  real (8) :: r,z, phi, theta, psi, b
  real (8) :: dnb, ti_ev
  real (8) :: vp, up, rho, rho_m, mu, w, deltaw, w1
  real (8) :: accel_factor, ekmin, ekin, pitch, s_mass, s_charge
  real (8) :: v(3), b_vec(3), tmp(2)
  real (8) :: two_therm_en, v_square, y, sqrty, dphidy, phiy
  real (kind=8),external :: psi_interpol, b_interpol, derf
#ifdef COL_NO_EDGE_CONSERV
  integer :: edge_no_conserv_j_start
#endif


!pw  call t_startf("CONSERV_COL_INIT")
#ifdef COL_NO_EDGE_CONSERV
  edge_no_conserv_j_start=col_2_m-col_no_conserv_m_num+1
#endif
  allocate(org_avg(3,col_2_mtheta,0:col_2_m,ptl_isp:ptl_nsp,sml_nthreads),&
           dum_avg(3,col_2_mtheta,0:col_2_m,ptl_isp:ptl_nsp))
  allocate(col_c_org_vsum(3,col_2_mtheta,0:col_2_m,ptl_isp:ptl_nsp,sml_nthreads),&
           col_c_new_vsum(3,col_2_mtheta,0:col_2_m,ptl_isp:ptl_nsp,sml_nthreads),&
           col_c_delta(3,col_2_mtheta,0:col_2_m,ptl_isp:ptl_nsp))
  
  if(col_2_iteration_method==2) then
     allocate(col_c_matrix(9,col_2_mtheta,0:col_2_m,ptl_isp:ptl_nsp,sml_nthreads),&
              dumv(9,col_2_mtheta,0:col_2_m,ptl_isp:ptl_nsp))
     col_c_matrix=0D0
  endif

  dt = sml_dt * col_period
  iflag2= 1 + col_en_col_on*2
  org_avg=0D0
  col_c_org_vsum=0D0
  col_c_new_vsum=0D0

  tmp=0D0
!pw  call t_stopf("CONSERV_COL_INIT")

  do isp=ptl_isp , ptl_nsp
      c_m_sp=ptl_c_m(isp)
      c2_2m_sp=ptl_c2_2m(isp)
      s_mass=sp%mass
      s_charge=sp%charge
      ekmin = 1.d-3 * sml_ev2j  ! 1 eV minimum temperature for collision
      call split_indices(sp%num, sml_nthreads, i_beg, i_end)

!$OMP PARALLEL DO &
!$OMP PRIVATE( ITH, I, R, Z, PHI, PSI, THETA, J, L, B, RHO, MU, &
!$OMP&         DNB, TI_EV, UP, RHO_M, W, V, EKIN, PITCH, ACCEL_FACTOR, &
!$OMP&         TWO_THERM_EN, V_SQUARE, Y, SQRTY, DPHIDY, PHIY, TMP)
      do ith=1,sml_nthreads
         call t_startf("CONSERV_COL_LOOP1")
         do i=i_beg(ith),i_end(ith)
             r=sp%ptl(i)%ph(1)
             z=sp%ptl(i)%ph(2)
             phi=sp%ptl(i)%ph(3)
             psi=psi_interpol(r,z,0,0)

             if(col_2_pin<psi .AND. psi<col_2_pout .AND. is_rgn12(r,z,psi) .AND. z>eq_x_z .and. z<eq_x2_z) then
                theta=acos((r-eq_axis_r)/sqrt((r-eq_axis_r)**2 + (z-eq_axis_z)**2))
                if(z < 0D0) then
                   theta= sml_2pi - theta
                endif

                j=nint((psi-col_2_pin)*col_2_inv_dp)
                l=mod(nint(theta*col_2_inv_dtheta),col_2_mtheta)+1
                b=b_interpol(r,z,phi)
                rho=sp%ptl(i)%ph(pirho)
                mu=sp%ptl(i)%ct(pim)

                call background_ion_profile(theta,r,z,psi,dnb,ti_ev,up)

                if (col_moving_frame) then
                   rho_m=rho - up/(c_m_sp*b)
                else
                   rho_m=rho
                endif

                if(ptl_deltaf_sp(sp%type) .and. col_varying_bg==1) then 
                    print *,'deltaf collision is not implemented yet.'
                    stop
!                    w=ptl%ion%phase(6,i)*ptl%ion%phase(8,i)
                else
                    w=sp%ptl(i)%ct(piw0)
                endif

                !quantity before test ptl collision
                v(1)=1D0    !weight
                v(3)=c_m_sp*rho_m*b     !velocity in moving frame
                v(2)=0.5D0*s_mass*v(3)**2 + mu*b    !energy in moving frame
                col_c_org_vsum(:,l,j,isp,ith)=col_c_org_vsum(:,l,j,isp,ith)+v(:)*w

                !test ptl collision
                call rho_mu_to_ev_pitch2(rho_m,mu,b,ekin,pitch,sp%type) ! ekin is in eV unit
                ekin=ekin*sml_ev2j         ! conversion to SI unit
                ekin = max(ekmin, ekin)        ! 2002/09/17 added for preventing small energy
                
                if(col_accel) then
                   call col_accel_factor(psi,accel_factor)
                else
                   accel_factor=1D0
                endif
                
                call scatr_one(ekin, pitch, s_mass, s_charge, dnb, ti_ev, ptl_mass_au, ptl_charge_eu,&
                 accel_factor, dt, ekmin, iflag2)
                 
                call ev_pitch_to_rho_mu2(ekin*sml_j2ev, pitch, b, rho_m, mu,sp%type)

                if (col_moving_frame) then
                   sp%ptl(i)%ph(pirho)=rho_m + up/(c_m_sp*b)
                else
                   sp%ptl(i)%ph(pirho)=rho_m
                endif
                sp%ptl(i)%ct(pim)=mu

                !quantity after test ptl collision
                v(3)=c_m_sp*rho_m*b     !velocity in moving frame
                v(2)=0.5D0*s_mass*v(3)**2 + mu*b    !energy in moving frame
                col_c_new_vsum(:,l,j,isp,ith)=col_c_new_vsum(:,l,j,isp,ith)+v(:)*w

                !coefficient calculation
                two_therm_en = 2D0*ti_ev*sml_ev2j/s_mass
                v_square=2D0/s_mass*v(2)
                y=v_square/two_therm_en
                sqrty = sqrt(y)

                dphidy=2D0/sml_sqrtpi*exp(-y)*sqrty
                phiy=derf(sqrty)-dphidy

                col_org_nuE(i)=(phiy-dphidy)/sqrty**3
                col_org_nuD(i)=((1D0-0.5D0/y)*phiy+dphidy)/sqrty**3

                !org_nuE_v4_avg
                org_avg(1,l,j,isp,ith)=org_avg(1,l,j,isp,ith)+w*col_org_nuE(i)*(v_square**2)
                !org_nuD_vp2_avg
                org_avg(2,l,j,isp,ith)=org_avg(2,l,j,isp,ith)+w*col_org_nuD(i)*(v(3)**2)
                !org_v2_avg
                org_avg(3,l,j,isp,ith)=org_avg(3,l,j,isp,ith)+w*v_square

                !save ptl information
                col_org_j(i)=j
                col_org_l(i)=l
                col_org_vp_m(i)=v(3)
                col_org_v2(i)=v_square
                col_gid(i)=.true.
             else
                col_gid(i)=.false.
             endif
          enddo
          call t_stopf("CONSERV_COL_LOOP1")
      enddo
  enddo

  do ith=2,sml_nthreads
     col_c_org_vsum(:,:,:,:,1)=col_c_org_vsum(:,:,:,:,1)+col_c_org_vsum(:,:,:,:,ith)
     col_c_new_vsum(:,:,:,:,1)=col_c_new_vsum(:,:,:,:,1)+col_c_new_vsum(:,:,:,:,ith)
     org_avg(:,:,:,:,1)=org_avg(:,:,:,:,1)+org_avg(:,:,:,:,ith)
  enddo

  call t_startf("CONSERV_COL_RED1")
  call my_mpi_allreduce(col_c_org_vsum(:,:,:,:,1),dum_avg,3*col_2_mtheta*(col_2_m+1)*(ptl_nsp-ptl_isp+1))
  col_c_org_vsum(:,:,:,:,1)=dum_avg
  call my_mpi_allreduce(col_c_new_vsum(:,:,:,:,1),dum_avg,3*col_2_mtheta*(col_2_m+1)*(ptl_nsp-ptl_isp+1))
  col_c_new_vsum(:,:,:,:,1)=dum_avg
  col_c_delta=col_c_new_vsum(:,:,:,:,1)-col_c_org_vsum(:,:,:,:,1)

  call my_mpi_allreduce(org_avg(:,:,:,:,1),dum_avg,3*col_2_mtheta*(col_2_m+1)*(ptl_nsp-ptl_isp+1))
  do i=1,3
     org_avg(i,:,:,:,1)=dum_avg(i,:,:,:)/col_c_org_vsum(1,:,:,:,1)
  enddo
  call t_stopf("CONSERV_COL_RED1")

  do isp=ptl_isp , ptl_nsp
      s_mass=sp%mass
      call split_indices(sp%num, sml_nthreads, i_beg, i_end)

!$OMP PARALLEL DO &
!$OMP PRIVATE( ITH, I, J, L, VP, V_SQUARE, V )
      do ith=1,sml_nthreads
!pw      call t_startf("CONSERV_COL_LOOP2")
         do i=i_beg(ith),i_end(ith)
            if(col_gid(i)) then
               j=col_org_j(i)
               l=col_org_l(i)

               vp=col_org_vp_m(i)
               v_square=col_org_v2(i)

               col_coeff_a(i) = -col_org_nuD(i)*vp/org_avg(2,l,j,isp,1)
               col_coeff_b(i) = -col_org_nuE(i)*v_square/org_avg(1,l,j,isp,1)
               col_coeff_d(i) = col_org_nuE(i)*v_square*org_avg(3,l,j,isp,1)&
                                /org_avg(1,l,j,isp,1)- 1d0

               if(col_2_iteration_method==2) then
                  ! ES: Matrix coefficient : memory access order in Fortran
                  v(1) = 1D0
                  v(2) = 0.5D0 * s_mass * v_square
                  v(3) = vp
         
                  col_c_matrix(1:3,l,j,isp,ith) = col_c_matrix(1:3,l,j,isp,ith) +&
                                              col_coeff_a(i) * v(:)
                  col_c_matrix(4:6,l,j,isp,ith) = col_c_matrix(4:6,l,j,isp,ith) +&
                                              col_coeff_b(i) * v(:)
                  col_c_matrix(7:9,l,j,isp,ith) = col_c_matrix(7:9,l,j,isp,ith) +&
                                              col_coeff_d(i) * v(:)
               endif
            endif
         enddo
!pw      call t_stopf("CONSERV_COL_LOOP2")
      enddo
  enddo

  if(col_2_iteration_method==1) then
     !iterative way
     if (sml_mype==0) then
        print *,'iterative collision is not available now.'
     endif
     stop
  elseif(col_2_iteration_method==2) then   !matrix solving from xgc1p

     do ith=2,sml_nthreads
        col_c_matrix(:,:,:,:,1)=col_c_matrix(:,:,:,:,1)+col_c_matrix(:,:,:,:,ith)
     enddo
     !************* ES: Reduce matrix coefficients ********************
     call t_startf("CONSERV_COL_RED2")
     call my_mpi_allreduce(col_c_matrix(:,:,:,:,1),dumv,9*col_2_mtheta*(col_2_m+1)*(ptl_nsp-ptl_isp+1))
     col_c_matrix(:,:,:,:,1) = dumv 
     call t_stopf("CONSERV_COL_RED2")

     !************* Solve linear equation for each grid ***************  
     ! ES: coorperation among processors(?)-reasonable?
     ! load allocation
     ! tmpi1 = mod(col_c_ntheta*(col_c_nr+1), sml_totalpe)
     ! tmpi2 = col_c_ntheta*(col_c_nr+1)/sml_totalpe

!pw  call t_startf("CONSERV_COL_LINSOLV")
     ! Ax=B
     ! OpenMP can be utilized
     do isp=ptl_isp, ptl_nsp
        do i=0, col_2_m
           do j=1, col_2_mtheta
              v = 0d0
              b_vec(:) = col_c_delta(:,j,i,isp)
      
              call linsol3(col_c_matrix(:,j,i,isp,1),b_vec,v)
              ! correction factor is saved in 'dumv' for memory saving
              dumv(1:3,j,i,isp) = v(:) 
           enddo
        enddo    
     enddo
!pw  call t_stopf("CONSERV_COL_LINSOLV")
   
     !************* Change weight *************************************
     do isp=ptl_isp , ptl_nsp       
        ! for reductions, divide particles among OpenMP threads
        ! and calculate contributions for each subset
        ! of particles independently. This will introduce a
        ! roundoff difference in the results for different numbers
        ! of threads.
        s_mass=sp%mass
     
!$OMP PARALLEL DO &
!$OMP PRIVATE( ITH, I, J, L, &
!$OMP& W, DELTAW )
        do ith=1,sml_nthreads
!pw        call t_startf("CONSERV_COL_LOOP3")
           do i=i_beg(ith),i_end(ith)
              if(col_gid(i)) then
                 j=col_org_j(i)
                 l=col_org_l(i)

                if(ptl_deltaf_sp(sp%type) .and. col_varying_bg==1) then
                    print *,'deltaf collision with varying bg is not implemented yet.'
                    stop
!                    w=ptl%ion%phase(6,i)*ptl%ion%phase(8,i)
                elseif(ptl_deltaf_sp(sp%type) .and. col_varying_bg==0) then
                    !Julien, copy paste to have only one "if then else" concerning the df / full-f switch
                    w=sp%ptl(i)%ct(piw0)
                    w1=sp%ptl(i)%ph(piw1)*w


                    !update the weight!!!
                    deltaw =  dumv(1,l,j,isp) * col_coeff_a(i) &
                            + dumv(2,l,j,isp) * col_coeff_b(i) &
                            + dumv(3,l,j,isp) * col_coeff_d(i)
#ifdef COL_NO_EDGE_CONSERV
                    if(j<edge_no_conserv_j_start) then
#endif
                      !Julien, f = f0 + df, one needs to alter both weights as f contains the df.
                      w  =  w - deltaw
                      w1 = w1 - deltaw
#ifdef COL_NO_EDGE_CONSERV
                    endif
#endif

                   if(w>0D0) then
                       sp%ptl(i)%ct(piw0)=w
                       sp%ptl(i)%ph(piw1)=w1/w
                    else
                       w  =   w + deltaw
                       w1 = (w1 + deltaw)/w
                       sp%ptl(i)%ct(piw0)=0.1*w
                       sp%ptl(i)%ph(piw1)=w1 - 0.9
                    endif
                else ! Usual full-f scheme

                    w=sp%ptl(i)%ct(piw0)

                    !update the weight!!!
                    deltaw =  dumv(1,l,j,isp) * col_coeff_a(i) &
                            + dumv(2,l,j,isp) * col_coeff_b(i) &
                            + dumv(3,l,j,isp) * col_coeff_d(i)
#ifdef COL_NO_EDGE_CONSERV
                    if(j<edge_no_conserv_j_start) then
#endif
                       w = w - deltaw
#ifdef COL_NO_EDGE_CONSERV
                    endif
#endif
                    if(w>0D0) then
                       sp%ptl(i)%ct(piw0)=w
                    else
                       sp%ptl(i)%ct(piw0)=0.1*(w+deltaw)
                    endif
                endif

              endif
           enddo
        enddo
!pw     call t_stopf("CONSERV_COL_LOOP3")
     enddo
  endif

  deallocate(org_avg,dum_avg, col_c_org_vsum, col_c_new_vsum, col_c_delta)
  if(col_2_iteration_method==2) then
     deallocate(col_c_matrix,dumv)
  endif

end subroutine conserving_collision

! 3 by 3 linear solver a*x=b
subroutine linsol3(mata,bs,x)
 real(8) :: det
 real(8) :: a(9),b(3),x(3)
 real(8) :: mata(9), bs(3)
 !real(8) :: normr1, normr2, normr3
 !real(8) :: normc1, normc2, normc3 
 a = mata
 b = bs

 !normr1 = a(1)
 !normr2 = a(2)
 !normr3 = a(3)
 !a(1) = 1d0
 !a(2) = 1d0
 !a(3) = 1d0
 !b(1) = b(1)/normr1
 !b(2) = b(2)/normr2
 !b(3) = b(3)/normr3
 !a(4) = a(4)/normr1
 !a(5) = a(5)/normr2
 !a(6) = a(6)/normr3
 !a(7) = a(7)/normr1
 !a(8) = a(8)/normr2
 !a(9) = a(9)/normr3

 !normc2 = a(4)
 !normc3 = a(7)
 !a(4) = 1d0
 !a(7) = 1d0
 !a(5) = a(5)/normc2
 !a(6) = a(6)/normc2
 !a(8) = a(8)/normc3
 !a(9) = a(9)/normc3
 
! Can LU decomposition make difference??--a??
! Simple normalization doens't change accuracy of solution
 det = a(1)*(a(5)*a(9)-a(6)*a(8)) + a(2)*(a(6)*a(7)-a(4)*a(9)) + a(3)*(a(4)*a(8)-a(5)*a(7))
 x(1) = ( b(1)*(a(5)*a(9)-a(6)*a(8)) + b(2)*(a(6)*a(7)-a(4)*a(9)) + b(3)*(a(4)*a(8)-a(5)*a(7)) )/det
 x(2) = -( b(1)*(a(2)*a(9)-a(3)*a(8)) + b(2)*(a(3)*a(7)-a(1)*a(9)) + b(3)*(a(1)*a(8)-a(2)*a(7)) )/det
 x(3) = ( b(1)*(a(2)*a(6)-a(3)*a(5)) + b(2)*(a(3)*a(4)-a(1)*a(6)) + b(3)*(a(1)*a(5)-a(2)*a(4)) )/det

 !x(2) = x(2)/normc2
 !x(3) = x(3)/normc3
end subroutine linsol3


