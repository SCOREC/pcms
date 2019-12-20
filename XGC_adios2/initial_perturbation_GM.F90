#define READ_IN
#undef DEBUG_MAP
#undef TRANSPOSE
#undef DOUBLE_READ
module initial_perturbation

IMPLICIT NONE

real(8) :: dR, dZ, sigma, rho0, dN0, Rmin, Zmin,pertw,cos_factor,sigma_n,sigma_chi,bt_sign,phi_sign,c1,c2,R0,Z0,B0;
logical :: radial_gaussian,circular_ad_hoc
real(8), dimension(:)  , allocatable :: nys
real(8), dimension(:,:), allocatable :: filter_ij,qs_ij,chi_ij,rho_ij
integer :: nr,nz,nky,delta_nky,nkymax, poloidal_envelope
logical :: FIRST
logical :: GM_perturbation
integer :: NPSI,nz_GENE
real(8), dimension(:), allocatable :: rho_tor,psi_grid,q,R_grid,Z_grid
real(8), dimension(:,:), allocatable :: chi_map,psi_map,filter
real(8) :: RminJD, ZminJD, dZJD, drJD
contains

subroutine ip_initialize()
  use sml_module, only:sml_bt_sign,sml_nphi_total,sml_wedge_n,sml_mype, sml_comm_null, sml_comm
  use adios_read_mod
  include 'mpif.h'

  logical :: nzero

  integer :: ir, iz, ierr
  character(100) :: fname
  integer :: adios_err
  integer(8) :: adios_handle, bb_sel
  integer(8), dimension(2) :: bounds, counts

     ! change second index afterswitching node order in XGC
  namelist /perturbation/ nr,nz,dR,dZ,sigma,rho0,nky,delta_nky,nzero,Rmin,Zmin,pertw,&
                          cos_factor,radial_gaussian,poloidal_envelope,nkymax,       &
                          sigma_n,sigma_chi,bt_sign,phi_sign,circular_ad_hoc,c1,c2,Z0,R0,B0,  &
                          GM_perturbation,npsi, nz_GENE
  
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
  GM_perturbation=.true.
  R0=1.7D0;
  Z0=0D0;
  c1=0.854D0;
  c2=13.6893D0; 
    
  if(sml_bt_sign>0)then
    bt_sign=-1D0;
  else
    bt_sign=+1D0;
  endif

  phi_sign=1d0

  open(unit=20,file='perturbation.in', status='old',action='read')
  READ(20, NML=perturbation)
  close(unit=20)

if(sml_mype.eq.0)then
  print *, 'nr=',nr
  print *, 'nz=',nz
  print *, 'npsi=',npsi
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
  print *, 'GM_perturbation=',GM_perturbation
endif

  if(nzero)then
    dN0=1D0
  else
    dN0=0D0
  endif
  !print *,nr,nz,dR,dZ,sigma,rho0
  if ((.not.circular_ad_hoc).and.(.not.GM_perturbation))then 
    allocate(filter_ij(nr,nz))
    allocate(qs_ij(nr,nz))
    allocate(chi_ij(nr,nz))
    allocate(rho_ij(nr,nz))
  
    open(unit=20,file='perturbation.in', status='old',action='read')
    READ(20, NML=perturbation_data)
    close(unit=20)
  endif

#ifdef READ_IN
  if (GM_perturbation) then
#ifdef DOUBLE_READ
     allocate(filter_ij(2000, 2000))
     allocate(qs_ij(2000, 2000))
     allocate(chi_ij(2000, 2000))
     allocate(rho_ij(2000, 2000))

     open(unit=20,file='perturbation.in', status='old',action='read')
     READ(20, NML=perturbation_data)
     close(unit=20)
     drJD=dr
     dzJD=dz
     RminJD=Rmin
     ZminJD=Zmin
#endif
     allocate(chi_map(NR,NZ))
     allocate(psi_map(NR,NZ))
     allocate(filter(NR,NZ))
     allocate(rho_tor(NPSI),psi_grid(NPSI), q(NPSI))
     allocate(R_grid(NR))
     allocate(Z_grid(NZ))

     if (sml_mype.eq.0)then
        fname="ogyropsi_init_cond.bp"

        call adios_read_open_file(adios_handle, fname, 0,MPI_COMM_SELF, adios_err)

        bounds(1) = 0
        bounds(2) = 0
        counts(1) = int(NPSI, kind=8)
        counts(2) = 1
        call adios_selection_boundingbox(bb_sel, 1, bounds, counts)
        call adios_schedule_read(adios_handle, bb_sel, "/data/var1d/rho_tor", 0, 1, rho_tor, adios_err)
        call adios_schedule_read(adios_handle, bb_sel, "/data/var1d/q", 0, 1, q, adios_err)
        call adios_schedule_read(adios_handle, bb_sel, "/data/grid/PSI", 0, 1, psi_grid, adios_err)
        call adios_perform_reads(adios_handle, adios_err)
        call adios_selection_delete(bb_sel)

        counts(1) = int(NR, kind=8)
        call adios_selection_boundingbox(bb_sel, 1, bounds, counts)
        call adios_schedule_read(adios_handle, bb_sel, "/data/var1d/rmesh", 0, 1, R_grid, adios_err)
        call adios_perform_reads(adios_handle, adios_err)
        call adios_selection_delete(bb_sel)

        counts(1) = int(NZ, kind=8)
        call adios_selection_boundingbox(bb_sel, 1, bounds, counts)
        call adios_schedule_read(adios_handle, bb_sel, "/data/var1d/zmesh", 0, 1, Z_grid, adios_err)
        call adios_perform_reads(adios_handle, adios_err)
        call adios_selection_delete(bb_sel)

        counts(1) = int(NR, kind=8)
        counts(2) = int(NZ, kind=8)
        call adios_selection_boundingbox(bb_sel, 2,  bounds, counts)
        call adios_schedule_read(adios_handle, bb_sel, "/data/var2d/psiRZ", 0, 1, psi_map, adios_err)
        call adios_schedule_read(adios_handle, bb_sel, "/data/var2d/chiRZ", 0, 1, chi_map, adios_err)
        call adios_perform_reads(adios_handle, adios_err)
        call adios_selection_delete(bb_sel)

        call adios_read_close(adios_handle, adios_err)

!        psi_map=transpose(psi_map)
!        chi_map=transpose(chi_map)

     endif
     
     call mpi_bcast(psi_grid,npsi,mpi_real8,0,sml_comm,ierr) 
     call mpi_bcast(q,npsi,mpi_real8,0,sml_comm,ierr) 
     call mpi_bcast(rho_tor,npsi,mpi_real8,0,sml_comm,ierr) 
     call mpi_bcast(R_grid,nr,mpi_real8,0,sml_comm,ierr) 
     call mpi_bcast(Z_grid,nz,mpi_real8,0,sml_comm,ierr) 
     call mpi_bcast(psi_map,nr*nz,mpi_real8,0,sml_comm,ierr)
     call mpi_bcast(chi_map,nr*nz,mpi_real8,0,sml_comm,ierr)

     do ir=1,NR
        do iz=1,NZ
           if (psi_map(ir,iz).le.1.0) then
              filter(ir,iz)=0.0D0
           else
              filter(ir,iz)=1.0D0
           endif
        end do
     end do

     dZ=Z_grid(2)-Z_grid(1)     
     dR=R_grid(2)-R_grid(1)
     Zmin=Z_grid(1)
     Rmin=R_grid(1)

     if (sml_mype.eq.0) then
        print*, "importing from "// fname
        print*, "Zmin", Z_grid(1)
        print*, "Rmin", R_grid(1)
        print*, "Zmax", Z_grid(NZ)
        print*, "Rmax", R_grid(NR)
        print*, "dZ", dZ
        print*, "dR", dR
     endif

  endif
#endif

end subroutine

subroutine ip_eval(dN,R,Z,varphi)
  use sml_module, only:sml_mype,sml_nphi_total,sml_wedge_n,sml_mype,sml_2pi
  real(8), intent(inout) :: dN
  real(8), intent(in)  :: R,Z,varphi

  integer :: i,j,k
  real(8) :: wi,wj,qs,chi,rho,di,dj,tta,eps

  integer :: ind_l_tmp, ind_h_tmp, tmp_ind, y_res, ind_l, ind_u, ind_tmp
  real(8) :: y_cut, dy, dl, du, w_chi
  real(8) :: mychi, mypsi
  logical, external :: is_nan

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
  elseif (.not.GM_perturbation) then

  di=(Z-Zmin)/dZ;
  i=1+floor(di);
  wi=1D0-(di-dble(i-1));
  
  dj=(R-Rmin)/dR;
  j=1+floor(dj);
  wj=1D0-(dj-dble(j-1));
  
  !print *,'i,j,filter',i,j,wi,wj,(filter_ij(i,j)+filter_ij(i+1,j)+filter_ij(i,j+1)+filter_ij(i+1,j+1))>0.5

  if( 0.5>(filter_ij(i,j)+filter_ij(i+1,j)+filter_ij(i,j+1)+filter_ij(i+1,j+1)))then
    qs  = wj*(wi* qs_ij(i,j)+(1D0-wi)* qs_ij(i+1,j))+(1D0-wj)*(wi* qs_ij(i,j+1)+(1D0-wi)* qs_ij(i+1,j+1));
    chi = wj*(wi*chi_ij(i,j)+(1D0-wi)*chi_ij(i+1,j))+(1D0-wj)*(wi*chi_ij(i,j+1)+(1D0-wi)*chi_ij(i+1,j+1));
    rho = wj*(wi*rho_ij(i,j)+(1D0-wi)*rho_ij(i+1,j))+(1D0-wj)*(wi*rho_ij(i,j+1)+(1D0-wi)*rho_ij(i+1,j+1));
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
    endif

  else if (GM_perturbation) then
#ifdef READ_IN

#ifdef DOUBLE_READ
     di=(Z-ZminJD)/dZJD;
     i=1+floor(di);
     wi=1D0-(di-dble(i-1));

     dj=(R-RminJD)/dRJD;
     j=1+floor(dj);
     wj=1D0-(dj-dble(j-1));


     if( 0.5>(filter_ij(i,j)+filter_ij(i+1,j)+filter_ij(i,j+1)+filter_ij(i+1,j+1)))then
        qs  = wj*(wi* qs_ij(i,j)+(1D0-wi)* qs_ij(i+1,j))+(1D0-wj)*(wi* qs_ij(i,j+1)+(1D0-wi)* qs_ij(i+1,j+1));
        chi = wj*(wi*chi_ij(i,j)+(1D0-wi)*chi_ij(i+1,j))+(1D0-wj)*(wi*chi_ij(i,j+1)+(1D0-wi)*chi_ij(i+1,j+1));
        rho = wj*(wi*rho_ij(i,j)+(1D0-wi)*rho_ij(i+1,j))+(1D0-wj)*(wi*rho_ij(i,j+1)+(1D0-wi)*rho_ij(i+1,j+1));


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
    	if (sml_mype.eq.0) print *, 'txt-map (R,Z,qs,chi,rho)',R,Z,qs,chi,rho, dN
     endif
#endif

#ifndef TRANSPOSE
     di=(R-Rmin)/dR;
     i=1+floor(di);
     wi=1D0-(di-dble(i-1));
  
     dj=(Z-Zmin)/dZ;
     j=1+floor(dj);
     wj=1D0-(dj-dble(j-1));
#else
     di=(Z-Zmin)/dZ;
     i=1+floor(di);
     wi=1D0-(di-dble(i-1));

     dj=(R-Rmin)/dR;
     j=1+floor(dj);
     wj=1D0-(dj-dble(j-1));
#endif
     if (i .gt. NR) then
        if (sml_mype.eq.0) print *, "ERROR", R,i,wi, Rmin, dR
     endif
     if (j .gt. NZ) then
        if (sml_mype.eq.0) print *, "ERROR", Z,j,wj, Zmin, dZ
     endif

     if ((i.le.NR).and.(j.le.NZ)) then
        if( 0.5>(filter(i,j)+filter(i+1,j)+filter(i,j+1)+filter(i+1,j+1)))then
           mychi =      wj*(wi*chi_map(i,j)+(1D0-wi)*chi_map(i+1,j))+ &
                (1D0-wj)*(wi*chi_map(i,j+1)+(1D0-wi)*chi_map(i+1,j+1))

           mypsi =       wj*(wi*psi_map(i,j)+(1D0-wi)*psi_map(i+1,j))+ &
                (1D0-wj)*(wi*psi_map(i,j+1)+(1D0-wi)*psi_map(i+1,j+1))

           if (is_nan(mychi)) then
              print *, 'Found a NaN for (R,Z)',R,Z
           endif

           if (mypsi.gt.1) then
              print *, 'wrong psi for (R,Z)',R,Z,mypsi,wi,wj,i,j
           endif

           ind_tmp=minloc(abs(psi_grid-mypsi),1)

           !not wokirng at LCFS
           if (psi_grid(ind_tmp).le.mypsi) then
              ind_l=ind_tmp
              ind_u=ind_tmp+1
              dl=psi_grid(ind_u)-mypsi
              du=mypsi-psi_grid(ind_l)
           else
              ind_l=ind_tmp-1
              ind_u=ind_tmp
              dl=psi_grid(ind_u)-mypsi
              du=mypsi-psi_grid(ind_l)
           endif

           qs=(q(ind_l)*dl+q(ind_u)*du)/(psi_grid(ind_u)-psi_grid(ind_l))
           rho=(rho_tor(ind_l)*dl+rho_tor(ind_u)*du)/(psi_grid(ind_u)-psi_grid(ind_l))

#ifdef DEBUG_MAP
           if (sml_mype.eq.0.or.sml_mype.eq.2) print *, '(R,Z,psi,chi,q,rho,mype)',R,Z,mypsi,mychi,qs,rho,i,j,psi_map(i,j),chi_map(i,j),sml_mype
           !        if (sml_mype.eq.0) print *, '(wi,wj,ind_l,ind_u,dl,du)',wi,wj,ind_l,ind_u,dl,du
#endif


#else
           di=(Z-Zmin)/dZ;
           i=1+floor(di);
           wi=1D0-(di-dble(i-1));

           dj=(R-Rmin)/dR;
           j=1+floor(dj);
           wj=1D0-(dj-dble(j-1));

           if( 0.5>(filter_ij(i,j)+filter_ij(i+1,j)+filter_ij(i,j+1)+filter_ij(i+1,j+1)))then
              qs  = wj*(wi* qs_ij(i,j)+(1D0-wi)* qs_ij(i+1,j))+(1D0-wj)*(wi* qs_ij(i,j+1)+(1D0-wi)* qs_ij(i+1,j+1));
              mychi = wj*(wi*chi_ij(i,j)+(1D0-wi)*chi_ij(i+1,j))+(1D0-wj)*(wi*chi_ij(i,j+1)+(1D0-wi)*chi_ij(i+1,j+1));
              rho = wj*(wi*rho_ij(i,j)+(1D0-wi)*rho_ij(i+1,j))+(1D0-wj)*(wi*rho_ij(i,j+1)+(1D0-wi)*rho_ij(i+1,j+1));
#endif
              if (mychi .gt. sml_2pi/2*(1-sml_2pi/2/dble(nz_GENE)) ) then
                 w_chi=1.0!sml_2pi/2-chi
              else
                 w_chi=1.0
              end if

              !print *,'(qs,chi,rho)',qs,chi,rho
              dN=dN0; !Potential contribution from n=0 mode
              qs=bt_sign*qs

              dy=sml_2pi/sml_nphi_total/sml_wedge_n
              y_res=sml_nphi_total

              y_cut=(qs*mychi-phi_sign*varphi)/dy    

              ! #ifdef DEBUG_MAP
              !        if (sml_mype.eq.0) print *, '(qs,mychi,varphi,dy,y_cut)',qs,mychi,varphi,dy,y_cut
              ! #endif
              y_cut=mod(mod(y_cut,real(y_res,kind=8))+real(y_res,kind=8), real(y_res,kind=8))

              tmp_ind=int(y_cut)
              ind_l_tmp=mod(mod(tmp_ind,y_res)+y_res,y_res)
              ind_h_tmp=mod(mod(tmp_ind+1,y_res)+y_res,y_res)

              do k=1,nky
                 dN=dN + w_chi*( &
                      &     cos(dble(k*delta_nky)*dy*dble(ind_h_tmp))*(y_cut-dble(tmp_ind)) &
                      &   + cos(dble(k*delta_nky)*dy*dble(ind_l_tmp))*(1.0D0-(y_cut-dble(tmp_ind)))) !&
                 !             &   + (1.0-w_chi)*( &
                 !             &     cos(dble(k*delta_nky*(1+sml_2pi*qs))*dy*ind_h_tmp)*(y_cut-tmp_ind) &
                 !             &   + cos(dble(k*delta_nky*(1+sml_2pi*qs))*dy*ind_l_tmp)*(1.0-(y_cut-tmp_ind)))
                 ! #ifdef DEBUG_MAP
                 !              if (sml_mype.eq.0) print *,'(y_cut,ind_h_tmp,ind_l_tmp,tmp_ind,dN)',y_cut,ind_h_tmp,ind_l_tmp,tmp_ind,dN
                 ! #endif
              end do

              if (radial_gaussian) then
                 dN=exp(-0.5*((rho-rho0)/sigma)**2) * dN !* pertw
              endif
              if (poloidal_envelope==1) then
                 dN=exp(-0.5*(mychi/sigma_chi)**2) * dN
              endif
           else
              dN=0d0
           endif
#ifdef DOUBLE_READ
           if (sml_mype.eq.0) print *, 'che-map (R,Z,qs,chi,rho)',R,Z,qs,mychi,rho,dN
#endif
        else
           dN=0D0
        endif
     else
        dN=0D0
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
#ifdef READ_IN
  if (GM_perturbation) then
     deallocate(chi_map,psi_map,filter)
     deallocate(rho_tor,psi_grid, q)
     deallocate(R_grid,Z_grid)
  endif
#endif

end subroutine


end module
