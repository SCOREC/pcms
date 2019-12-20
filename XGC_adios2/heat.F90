subroutine source(istep,sp)
  use sml_module
  use src_module
  use ptl_module
  use eq_module
  use perf_monitor
  use omp_module, only : split_indices
  implicit none
  type(species_type) :: sp
  integer, intent(in) :: istep
  integer, parameter :: N=5
  real (8) :: c_m, mass
  integer ::  i, j
  real (8) :: r, z, phi, rho, mu, w0, psi, factor, br,bz,bphi,b, vp, v(N)
  real (8), dimension(src_narea) :: pimw, popw, pi, po, width
  real (8), dimension(src_narea) :: del_en, del_mo, alpha, beta, tmp, dToverT
  real (8), dimension(sp%num) :: B_info

  real (8) :: sumall(N,src_narea)
  logical, parameter :: iterative=.true.
  !omp
  real (8) :: l_sum(N,src_narea,sml_nthreads)
  integer :: ith, i_beg(sml_nthreads), i_end(sml_nthreads) ! omp
  ! external function
  real (8), external :: psi_interpol, b_interpol

  call split_indices(sp%num, sml_nthreads, i_beg, i_end)
  

  c_m=sp%c_m
  mass=sp%mass

  pi = src_pin
  po = src_pout
  width = src_decay_width
  pimw = pi - width
  popw = po + width

!$OMP PARALLEL DO &
!$OMP PRIVATE( ITH, I, R, Z, PHI, RHO, MU, W0, PSI, J, FACTOR, Br,Bz,BPhi,B, VP, V)
  do ith=1,sml_nthreads
     call t_startf("SOURCE_LOOP1") 
     l_sum(:,:,ith)=0D0
     
     do i=i_beg(ith),i_end(ith)
        if(sp%ptl(i)%gid > 0) then
           r=sp%ptl(i)%ph(1)
           z=sp%ptl(i)%ph(2)
           phi=sp%ptl(i)%ph(3)
           rho=sp%ptl(i)%ph(4)
           mu=sp%ptl(i)%ct(1)
           w0=sp%ptl(i)%ct(2)
           psi=psi_interpol(r,z,0,0)
!           psi= sp%psi(i)

           do j=1, src_narea
              if(pimw(j) < psi .AND. psi < popw(j) .and. is_rgn12(r,z,psi) ) then
                 if( pi(j) < psi .and. psi < po(j) ) then
                    factor=1D0
                 elseif( pimw(j) < psi .and. psi < pi(j) ) then
                    factor= (psi-pimw(j))/width(j)
                 else
                    factor= (popw(j)-psi)/width(j)
                 endif
                 
                 call bvec_interpol(r,z,phi,br,bz,bphi)
                 B=sqrt(br*br+bz*bz+bphi*bphi)
                 B_info(i)=B
                 vp=c_m*rho*B
                 

                 v(1)=0.5D0*mass*vp*vp + mu*B  ! energy
                 v(2)=vp                              ! parallel v
                 v(3)=factor                          ! weight sum * factor^2
                 v(4)=vp*r*Bphi/B
                 v(5)=r*Bphi/B


                 l_sum(:,j,ith)=l_sum(:,j,ith) + v(:)*factor*w0   ! full weight
              endif
              
           enddo
        endif
     enddo
     call t_stopf("SOURCE_LOOP1") 
  enddo
     
  ! omp summation
  do ith=2, sml_nthreads
     l_sum(:,:,1)=l_sum(:,:,1)+l_sum(:,:,ith)
  enddo
  
  ! mpi summation  
  call t_startf("SOURCE_RED")
  call my_mpi_allreduce(l_sum(:,:,1), sumall, N*src_narea)
  call t_stopf("SOURCE_RED")
  
  
  ! find alpha and beta -- iteratively or directly(=default)
  del_en(:)=src_heat_power(:) * sml_dt * real(src_period)
  del_mo(:)=src_torque(:)     * sml_dt * real(src_period)
  
  if(iterative) then
     alpha=0D0
     beta=0D0
     do i=1, src_niter
        tmp=(del_en - mass*sqrt(1D0+alpha)*beta*sumall(2,:) - 0.5D0*mass*beta*beta*sumall(3,:) )/sumall(1,:)
        beta =(del_mo - mass*(sqrt(1D0+alpha)-1D0)*sumall(4,:) ) / (sumall(5,:)*mass)
        alpha=tmp
     enddo
     
     alpha=max(alpha,-0.99D0)
  else
     alpha=del_mo
     dToverT=del_en/sumall(1,:)
     beta= sqrt(1D0+dToverT)-1D0
  endif

  if(iterative) then     
!$OMP PARALLEL DO &
!$OMP PRIVATE( ITH, I, R, Z, PSI, J, FACTOR, B)
     do ith=1,sml_nthreads
        call t_startf("SOURCE_LOOP2") 
        do i=i_beg(ith),i_end(ith)
           if(sp%ptl(i)%gid > 0)  then

              r=sp%ptl(i)%ph(1)
              z=sp%ptl(i)%ph(2)
              psi=psi_interpol(r,z,0,0)
              !psi= sp%psi(i)
              
              do j=1, src_narea
 
                 
                 if(pimw(j) < psi .AND. psi < popw(j) .and. is_rgn12(r,z,psi) ) then
                    if( pi(j) < psi .and. psi < po(j) ) then
                       factor=1D0
                    elseif( pimw(j) < psi .and. psi < pi(j) ) then
                       factor= (psi-pimw(j))/width(j)
                    else
                       factor= (popw(j)-psi)/width(j)
                    endif
                    B=B_info(i)
                    
                    ! adjust particle velocity
                    sp%ptl(i)%ct(1)=sp%ptl(i)%ct(1)*(1D0+alpha(j)*factor)    ! mu change
                    sp%ptl(i)%ph(4)=sp%ptl(i)%ph(4)*sqrt(1D0+alpha(j)*factor) + beta(j)*factor/(B*c_m) ! vp or rho change
                 endif

              enddo
           endif
        enddo
        call t_stopf("SOURCE_LOOP2") 
     enddo
  else !(.not. iterative)
!$OMP PARALLEL DO &
!$OMP PRIVATE( ITH, I, R, Z, PSI, J, PIMW, POPW, FACTOR)
     do ith=1,sml_nthreads
        call t_startf("SOURCE_LOOP2") 
        do i=i_beg(ith),i_end(ith)           
           if(sp%ptl(i)%gid > 0) then

              r=sp%ptl(i)%ph(1)
              z=sp%ptl(i)%ph(2)
              psi=psi_interpol(r,z,0,0)
              !psi= sp%psi(i)

              do j=1, src_narea
                 if(pimw(j) < psi .AND. psi < popw(j) .and. is_rgn12(r,z,psi) ) then
                    if( pi(j) < psi .and. psi < po(j) ) then
                       factor=1D0
                    elseif( pimw(j) < psi .and. psi < pi(j) ) then
                       factor= (psi-pimw(j))/width(j)
                    else
                       factor= (popw(j)-psi)/width(j)
                    endif
              
                    ! adjust particle velocity
                    sp%ptl(i)%ct(1)=sp%ptl(i)%ct(1)*(1D0+factor*dToverT(j))
                    sp%ptl(i)%ph(4)=sp%ptl(i)%ph(4)+factor*(alpha(j)+beta(j)*(sp%ptl(i)%ph(4)-sumall(4,j)/sumall(5,j)))
                 endif
              enddo
           endif           
        enddo
        call t_stopf("SOURCE_LOOP2") 
     enddo
  endif !(.not. iterative)
  
end subroutine source


