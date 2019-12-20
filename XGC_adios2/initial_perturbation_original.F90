module initial_perturbation

IMPLICIT NONE

real(8) :: dR, dZ, sigma, rho0, dN0, Rmin, Zmin,pertw,cos_factor,sigma_n,sigma_chi,bt_sign,c1,c2,R0,Z0,B0;
logical :: radial_gaussian,circular_ad_hoc;
real(8), dimension(:)  , allocatable :: nys
real(8), dimension(:,:), allocatable :: filter_ij,qs_ij,chi_ij,rho_ij
integer :: nr,nz,nky,delta_nky,poloidal_envelope,nkymax
logical :: FIRST
contains

subroutine ip_initialize()
  use sml_module, only:sml_bt_sign,sml_nphi_total,sml_wedge_n,sml_mype

  logical :: nzero

  namelist /perturbation/ nr,nz,dR,dZ,sigma,rho0,nky,delta_nky,nzero,Rmin,Zmin,pertw,&
                          cos_factor,radial_gaussian,poloidal_envelope,nkymax,       &
                          sigma_n,sigma_chi,bt_sign,circular_ad_hoc,c1,c2,Z0,R0,B0
  
  namelist /perturbation_data/ filter_ij,qs_ij,chi_ij,rho_ij
  FIRST=.true.
  !Default values
  cos_factor= 1D0
  pertw     = 1D-3
  nky       = sml_nphi_total/2
  delta_nky = sml_wedge_n
  nkymax=delta_nky*nky
  nzero     = .false.
  poloidal_envelope=0
  sigma_n=24
  sigma_chi=3.14
  circular_ad_hoc=.false.
  R0=1.7D0;
  Z0=0D0;
  c1=0.854D0;
  c2=13.6893D0; 
    
  if(sml_bt_sign>0)then
    bt_sign=-1D0;
  else
    bt_sign=+1D0;
  endif

  open(unit=20,file='perturbation.in', status='old',action='read')
  READ(20, NML=perturbation)
  close(unit=20)

if(sml_mype.eq.0)then
  print *, 'nr=',nr
  print *, 'nz=',nz
  print *, 'dR=',dR
  print *, 'dZ=',dZ
  print *, 'sigma=',sigma
  print *, 'rho0=',rho0
  print *, 'nky=',nky
  print *, 'delta_nky=',delta_nky
  print *, 'nzero=',nzero
  print *, 'Rmin=',Rmin
  print *, 'Zmin=',Zmin
  print *, 'pertw=',pertw
  print *, 'cos_factor=',cos_factor
  print *, 'radial_gaussian=',radial_gaussian
  print *, 'poloidal_envelope=',poloidal_envelope
  print *, 'nkymax=',nkymax
  print *, 'sigma_n=',sigma_n
  print *, 'sigma_chi=',sigma_chi 
endif

  if(nzero)then
    dN0=1D0
  else
    dN0=0D0
  endif
  !print *,nr,nz,dR,dZ,sigma,rho0
  if(.not.circular_ad_hoc)then 
    allocate(filter_ij(nr,nz))
    allocate(qs_ij(nr,nz))
    allocate(chi_ij(nr,nz))
    allocate(rho_ij(nr,nz))
  
    open(unit=20,file='perturbation.in', status='old',action='read')
    READ(20, NML=perturbation_data)
    close(unit=20)
  endif
end subroutine

subroutine ip_eval(dN,R,Z,varphi)
  use sml_module, only:sml_mype
  real(8), intent(inout) :: dN
  real(8), intent(in)  :: R,Z,varphi

  integer :: i,j,k
  real(8) :: wi,wj,qs,chi,rho,di,dj,tta,eps

  if(circular_ad_hoc)then 

    rho=sqrt((R-R0)**2+(Z-Z0)**2);!r
    tta=atan2((Z-Z0),(R-R0));
    eps=rho/R0;
    chi=2*atan(sqrt((1-eps)/(1+eps))*tan(0.5*tta));
    dN=dN0; !Potential contribution from n=0 mode
    qs=bt_sign*(c1+c2*rho**2);
    if(FIRST.AND.(0.18.lt.rho.and.rho.lt.0.22))then
      print *,'qs(r)',R,R0,Z,Z0,rho,rho0,tta,chi,qs,dble(k*delta_nky)*qs
      FIRST=.false.;
    endif
    do k=1,nky
      !dN=dN+2D0*cos(dble(k*delta_nky)*(qs*chi-varphi));
      if(poloidal_envelope==0)then
        dN=dN+ cos_factor*cos(dble(k*delta_nky)*(qs*chi-varphi));
      elseif(poloidal_envelope==1)then
        dN=dN+ cos_factor*cos(dble(k*delta_nky)*(qs*chi-varphi))/sqrt(dble(k));
      elseif(poloidal_envelope==2)then
        dN=dN+ cos_factor*cos(dble(k*delta_nky)*(qs*chi-varphi))/dble(k);
      elseif(poloidal_envelope==3)then
        dN=dN+ cos_factor*cos(dble(k*delta_nky)*(qs*chi-varphi))*exp(-0.5*(dble(k)/sigma_n)**2);
      elseif(poloidal_envelope==4)then
        dN=dN+ cos_factor*cos(dble(k*delta_nky)*(qs*chi-varphi)) &
               *exp(-0.5*(dble(k)/sigma_n)**2)*exp(-0.5*(chi/sigma_chi)**2);
      endif
    enddo
    if(radial_gaussian)then
      dN=exp(-0.5*((rho-rho0)/sigma)**2) * dN !* pertw
    endif
!  if(sml_mype.eq.0) print *,'[initial_perturbation.F90]',rho,tta,eps,chi,qs,dN
  else

  di=(Z-Zmin)/dZ;
  i=1+floor(di);
  wi=1D0-(di-dble(i-1));
  
  dj=(R-Rmin)/dR;
  j=1+floor(dj);
  wj=1D0-(dj-dble(j-1));
  
  !print *,'i,j,filter',i,j,wi,wj,(filter_ij(i,j)+filter_ij(i+1,j)+filter_ij(i,j+1)+filter_ij(i+1,j+1))>0.5

  if( 0.5>(filter_ij(i,j)+filter_ij(i+1,j)+filter_ij(i+1,j+1)+filter_ij(i+1,j+1)))then
    qs  = wj*(wi* qs_ij(i,j)+(1D0-wi)* qs_ij(i+1,j))+(1D0-wj)*(wi* qs_ij(i+1,j)+(1D0-wi)* qs_ij(i+1,j+1));
    chi = wj*(wi*chi_ij(i,j)+(1D0-wi)*chi_ij(i+1,j))+(1D0-wj)*(wi*chi_ij(i+1,j)+(1D0-wi)*chi_ij(i+1,j+1));
    rho = wj*(wi*rho_ij(i,j)+(1D0-wi)*rho_ij(i+1,j))+(1D0-wj)*(wi*rho_ij(i+1,j)+(1D0-wi)*rho_ij(i+1,j+1));
    !print *,'(qs,chi,rho)',qs,chi,rho
    dN=dN0; !Potential contribution from n=0 mode
    qs=bt_sign*qs
    do k=1,nky
      !dN=dN+2D0*cos(dble(k*delta_nky)*(qs*chi-varphi));
      if(poloidal_envelope==0)then
        dN=dN+ cos_factor*cos(dble(k*delta_nky)*(qs*chi-varphi));
      elseif(poloidal_envelope==1)then
        dN=dN+ cos_factor*cos(dble(k*delta_nky)*(qs*chi-varphi))/sqrt(dble(k));
      elseif(poloidal_envelope==2)then
        dN=dN+ cos_factor*cos(dble(k*delta_nky)*(qs*chi-varphi))/dble(k);
      elseif(poloidal_envelope==3)then
        dN=dN+ cos_factor*cos(dble(k*delta_nky)*(qs*chi-varphi))*exp(-0.5*(dble(k)/sigma_n)**2);
      elseif(poloidal_envelope==4)then
        dN=dN+ cos_factor*cos(dble(k*delta_nky)*(qs*chi-varphi)) &
               *exp(-0.5*(dble(k)/sigma_n)**2)*exp(-0.5*(chi/sigma_chi)**2);
      endif
    enddo
    if(radial_gaussian)then
      dN=exp(-0.5*((rho-rho0)/sigma)**2) * dN !* pertw
    endif
  else
    dN=0D0
  endif
  !dN=(filter_ij(i,j)+filter_ij(i+1,j)+filter_ij(i,j+1)+filter_ij(i+1,j+1))
  endif
end subroutine

subroutine ip_test()

  integer :: i,j
  
  real(8) :: R,Z,dN
  
  open(unit=212,file='perturbation.out', status='replace',action='write')
  do i=1,nr-1
    R=Rmin+dble(i-1)*dR
    do j=1,nz-1
      Z=Zmin+dble(j-1)*dZ
      call ip_eval(dN,R,Z,0D0)
      !write(unit=212,fmt='(I5,I5,ES14.7)') i,j,dN
      write(unit=212,fmt='(ES14.7,1x,ES14.7,1x,ES14.7)') R,Z,dN
    enddo
  enddo
  close(unit=212)

end

subroutine ip_destroy()
  if(allocated(filter_ij)) deallocate(filter_ij)
  if(allocated(qs_ij)) deallocate(qs_ij)
  if(allocated(chi_ij)) deallocate(chi_ij)
  if(allocated(rho_ij)) deallocate(rho_ij)
end subroutine


end module
