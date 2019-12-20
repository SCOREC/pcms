attributes(device) &
subroutine bounce_gpu(new_phase,old_phase,rtn)
  use eq_module_gpu , only : eq_axis_r,eq_x_psi
  use ptl_module_gpu, only : ptl_nphase, piw1, piw2
  use sml_module_gpu
  use bnc_module_gpu, only :  bnc_max_r, bnc_min_r
  use precision_mod_gpu
  implicit none
  real (kind=work_p),intent(inout) :: new_phase(ptl_nphase),old_phase(ptl_nphase)
  integer, intent(out) :: rtn
  integer :: inner
  real (kind=work_p) :: r, z, psi,sign
  real (kind=work_p) :: b , br, bz, bphi,mid,z1,z2,z_pmin
  integer :: i,j,count
  integer, parameter :: JMAX=5
  !real (kind=8), external :: psi_interpol,z_psi_min,I_interpol,B_interpol_sym
  !  real (kind=8), parameter :: ZTOL=1D-6, NPSI_TOL=1D-5
  !  real (kind=8), parameter :: ZTOL=1D-10, NPSI_TOL=1D-5
  real (kind=work_p), parameter :: ZTOL=1E-6_work_p, NPSI_TOL=1E-4_work_p, BTOL=1E-8_work_p, PTOL=1E-8_work_p
  real (kind=work_p) :: psitol
  real (kind=work_p) :: delz,deltab2,deltap2
  real (kind=work_p) :: RBt,oldB,oldP,newB,newP,deltab,deltap,sign2,deltar,deltaz,zh,rh,dz1,dz2,dz3
  real (kind=work_p) :: dpsi_dr,dpsi_dz,d2psi_d2r,d2psi_drdz,d2psi_d2z,db_dr,db_dz,denom
  logical :: use_current, diverge
  psitol=NPSI_TOL*eq_x_psi
  rtn=0

  ! reverse angle --   z -> -z, phase variable -> old phase variable except z
  !if(sml_concentric) then
  if(.true.) then
     psi=psi_interpol_gpu(new_phase(1),new_phase(2),0,0)
!     if(psi<=sml_inpsi .or. sml_bounce==2 .and. psi>= sml_outpsi) then
     if( (sml_bounce==1 .or. sml_bounce==2) .and. psi < sml_inpsi .or. (sml_bounce==2 .or. sml_bounce==3) .and. psi > sml_outpsi ) then
        new_phase(1:ptl_nphase)=old_phase(1:ptl_nphase)
        new_phase(2)=-old_phase(2)
        if(sml_bounce_zero_weight==1) new_phase(piw1:piw2)=0.0_work_p
     endif
  else 
     ! bounce at inner (& outer )boundary
     ! Finding z value which has same psi and r of old particle position
     ! This position gives same mu, P_phi.
     ! similar energy if B is similar.
     ! for exact calculation, 1st order correction is needed.
     psi=psi_interpol_gpu(new_phase(1),new_phase(2),0,0)
     r=old_phase(1)    ! radius of old position
     !     if((psi <= sml_inpsi .or. sml_bounce==2 .and. psi >= sml_outpsi) .and. &
     !        r < bnc_max_r .and. r > bnc_min_r .and. psi < eq_x_psi ) then                
     if(( (sml_bounce==1 .or. sml_bounce==2) .and. psi < sml_inpsi .or. (sml_bounce==2 .or. sml_bounce==3) .and. psi > sml_outpsi) .and. &
          r < bnc_max_r .and. r > bnc_min_r .and. psi < eq_x_psi ) then                
        rtn=1 ! a big jump will be happen
        if(is_rgn1_gpu(r,z,psi)) then
           inner=1
        else
           inner=0
        endif
        psi=psi_interpol_gpu(old_phase(1),old_phase(2),0,0)  !  psi is old position value
        !        if(psi > (sml_inpsi-psitol) .and. psi < (sml_outpsi+psitol)) then  !old psi is inside simulation range
        if(psi >= (sml_inpsi-psitol) .and. psi <=(sml_outpsi+psitol)) then  !old psi is inside simulation range
           z=old_phase(2)
           z_pmin=z_psi_min_gpu(r)

           ! bug fixed. sign should be -1 
           !                 if(r>eq_axis_r) then
           !                    sign=1
           !                 else
           sign=-1
           !                 endif

           ! set boundary of binary search
           z1=z_pmin              
           if(z>z_pmin) then                                   
              !             z2=max(z_pmin - (z-z_pmin+0.01)*2. , sml_bd_min_z) ! heuristic formula
              if( z_pmin- (z-z_pmin+0.01)*2. > sml_bd_min_z ) then
                 z2= z_pmin- (z-z_pmin+0.01)*2.
              else
                 z2=sml_bd_min_z
              endif
           else
              !                z2=min(z_pmin - (z-z_pmin+0.01/sml_norm_r)*2.  , sml_bd_max_z) ! heuristic formula
              if(z_pmin - (z-z_pmin+0.01)*2. < sml_bd_max_z) then
                 z2=z_pmin - (z-z_pmin+0.01)*2
              else
                 z2=sml_bd_max_z
              endif
           endif

           ! find z value using binary search.-------------------------
           do while (abs(z1-z2) > ZTOL)
              mid=0.5_work_p*(z1+z2)
#ifdef XGC_BOUNCE_DEBUG
              if( .not. (mid > sml_bd_min_z .and. mid < sml_bd_max_z) ) then
                !print *, 'invaild z in psi_interpol(1)',mid, sml_bd_min_z, sml_bd_max_z
              endif
              if( .not. (r > sml_bd_min_r .and. r < sml_bd_max_r) ) then
                  !print *, 'invaild r in psi_interpol(1)',mid, sml_bd_min_r, sml_bd_max_r
                endif
#endif
              if(sign *(psi_interpol_gpu(r,mid,0,0)-psi) < 0 ) then
                 z2=mid
              else
                 z1=mid
              endif
              !                    write(400,*) mid,z1,z2,z1-z2,psi,psi_interpol(r,mid,0,0)
           enddo

           !           z1=0.5D0*(z1+z2)
           if(inner==1) then  ! z1 gives larger psi for inner bd.
              z1=z2
           endif
           !------------------------------------------------------------

           new_phase(1:ptl_nphase)=old_phase(1:ptl_nphase)
           new_phase(2)=z1              

           ! 1st order correction  
           !           delz=(psi-new_phase(9))/psi_interpol(new_phase(1),z1,0,1)
           !           new_phase(2)=new_phase(2)+delz
           !           new_phase(9)=psi_interpol(new_phase(1),new_phase(2),0,0)

           ! (dR,dZ) correction for attaining position with same B and psi values as old position  
           oldB=B_interpol_sym_gpu(old_phase(1),old_phase(2))
           oldP=psi_interpol_gpu(old_phase(1),old_phase(2),0,0)
           deltab2=oldB-B_interpol_sym_gpu(new_phase(1),new_phase(2))
#ifdef XGC_BOUNCE_DEBUG
           if( .not. (z1 > sml_bd_min_z .and. z1 < sml_bd_max_z) ) then
             !print *, 'invaild z in psi_interpol(2)',z1, sml_bd_min_z, sml_bd_max_z
           endif
           if( .not. (new_phase(1) > sml_bd_min_r .and. new_phase(1) < sml_bd_max_r) ) then
             !print *, 'invaild r in psi_interpol(2)',new_phase(1), sml_bd_min_r, sml_bd_max_r
           endif
#endif
           deltap2=oldP-psi_interpol_gpu(new_phase(1),z1,0,0)

           !init
           if( .not. (z > sml_bd_min_z .and. z < sml_bd_max_z) ) then
              rtn=-1
              return
           endif
           use_current=.false.
           diverge=.false.
           do j=1, JMAX
              ! Newton-Raphson procedure is iterated 
#ifdef XGC_BOUNCE_DEBUG
              if( .not. (new_phase(2) > sml_bd_min_z .and. new_phase(2) < sml_bd_max_z) ) then
                !print *, 'invaild z in psi_interpol(3)',new_phase(2), sml_bd_min_z, sml_bd_max_z
              endif
              if( .not. (new_phase(1) > sml_bd_min_r .and. new_phase(1) < sml_bd_max_r) ) then
                !print *, 'invaild r in psi_interpol(3)',new_phase(1), sml_bd_min_r, sml_bd_max_r
              endif
#endif
              psi=psi_interpol_gpu(new_phase(1),new_phase(2),0,0)
              dpsi_dr=psi_interpol_gpu(new_phase(1),new_phase(2),1,0)
              dpsi_dz=psi_interpol_gpu(new_phase(1),new_phase(2),0,1)
              RBt=I_interpol_gpu(psi,0,1)
              newB=sqrt(RBt**2+dpsi_dr**2+dpsi_dz**2)/new_phase(1)
              deltab=oldB-newB
              deltap=oldP-psi
              if(((dabs(deltab)/oldB)<BTOL).and.((dabs(deltap)/oldP)<PTOL)) then 
                 use_current=.true.
                 exit
              endif
              
              d2psi_d2r=psi_interpol_gpu(new_phase(1),new_phase(2),2,0)
              d2psi_drdz=psi_interpol_gpu(new_phase(1),new_phase(2),1,1)
              d2psi_d2z=psi_interpol_gpu(new_phase(1),new_phase(2),0,2)
              
              db_dr=-newB/new_phase(1)+1D0/(newB*(new_phase(1)**2))*(dpsi_dr*d2psi_d2r+dpsi_dz*d2psi_drdz)
              db_dz=1D0/(newB*(new_phase(1)**2))*(dpsi_dr*d2psi_drdz+dpsi_dz*d2psi_d2z)
              
              denom=dpsi_dr*db_dz-dpsi_dz*db_dr
              deltar=(db_dz*deltap-dpsi_dz*deltab)/denom
              deltaz=(-db_dr*deltap+dpsi_dr*deltab)/denom
              
              ! move to new (r,z) 
              new_phase(1)=new_phase(1)+deltar
              new_phase(2)=new_phase(2)+deltaz
              
              ! check diverge
              if(new_phase(1) < sml_bd_min_r .or. new_phase(1) > sml_bd_max_r &
                   .or. new_phase(2) < sml_bd_min_z .or. new_phase(2) > sml_bd_max_z ) then
                 diverge=.true.
                 use_current=.false. ! for safety
                 exit
              endif
           enddo
           
           if(.NOT. use_current .and. .NOT. diverge ) then
              ! loop ended not satisfying the TOL codition, nor diverge
              psi=psi_interpol_gpu(new_phase(1),new_phase(2),0,0)
              dpsi_dr=psi_interpol_gpu(new_phase(1),new_phase(2),1,0)
              dpsi_dz=psi_interpol_gpu(new_phase(1),new_phase(2),0,1)
              RBt=I_interpol_gpu(psi,0,1)
              newB=sqrt(RBt**2+dpsi_dr**2+dpsi_dz**2)/new_phase(1)
              deltab=oldB-newB
              deltap=oldP-psi
              ! use original (binary search position)
              if(((deltab2/oldB)**2+(deltap2/oldP)**2)<((deltab/oldB)**2+(deltap/oldP)**2)) then
                 use_current=.false.
              endif
           endif
           
           if(.NOT. use_current) then  
              new_phase(1:ptl_nphase)=old_phase(1:ptl_nphase)
              new_phase(2)=z1
           endif
           
           psi=psi_interpol_gpu(new_phase(1),new_phase(2),0,0)
           ! end of second (dR,dZ) correction for same B

           if(psi < (sml_inpsi-psitol) .or. psi > (sml_outpsi+psitol)) then
              !           if(new_phase(9) <= sml_inpsi-psitol .or. new_phase(9) >= sml_outpsi+psitol) then
              !print *, 'Fail finding proper psi', psi/eq_x_psi, oldP/eq_x_psi, psi/eq_x_psi
              !print *, 'oldB, newB', oldB, newB
              rtn=-1
           endif
           if(sml_bounce_zero_weight==1) new_phase(piw1:piw2)=0D0           
        else ! old position was outside of boundary
           rtn=-1
           !           print *, 'ptl_elimination', i, sml_mype
        endif
     endif
  endif

end subroutine bounce_gpu
