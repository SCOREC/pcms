!!particle (ion) pushing routine using R-K 2nd order + RK4 hybrid
subroutine push(istep,ipc,grid,psn,sp,phase0,ptl,diag_on)
  use sml_module
  use ptl_module
  use fld_module
  use grid_class
  use psn_class
  use omp_module , only : split_indices
  use perf_monitor
  use eq_module
  implicit none
  integer, intent(in) :: istep, ipc !! RK4 index
  type(grid_type) :: grid
  type(psn_type) :: psn
  type(species_type) :: sp
  real (kind=8),intent(inout) :: phase0(ptl_nphase,sp%maxnum)   ! sp%phase0 -- for speed
  type(ptl_type), intent(inout) ::  ptl(sp%maxnum)    ! sp%phase
  logical, intent(in) :: diag_on
  !
  type(fld_type) :: fld
  real (kind=8) :: dt_now,dt,time_now,new_phase(ptl_nphase),old_phase(ptl_nphase)
  integer :: i,rtn,j
  integer :: ith, i_beg(sml_nthreads), i_end(sml_nthreads)
  real (kind=8) , external :: psi_interpol
  character (len=5) :: err_str(2)
  logical, parameter :: USE_SEARCH_TR2 = .true.
  err_str(1)='ion'
  err_str(2)='elec'
  

  if(sp%num==0) return  ! nothing to push

  call split_indices(sp%num, sml_nthreads, i_beg, i_end)

   
  if(sp%type/=0) then
     ! for ion
     dt=sml_dt              
  else
     ! for electron
     dt=sml_dt*0.5D0/real(sml_ncycle_half)
  endif
  
 if(sml_nrk==4)then
#ifdef PURE_RK4
  select case(ipc)
  case(1)
     dt_now=0.5D0*sml_dt
     !$OMP PARALLEL DO &
     !$OMP PRIVATE( ITH, I )
     do ith=1,sml_nthreads
        do i=i_beg(ith), i_end(ith)
           phase0(:,i)=ptl(i)%ph
        enddo
     enddo
     !do i=1, sp%num
     !   ptli=>sp%ptl(i)
     !   ptli%ph0=ptli%ph
     !enddo
  case(2)
     dt_now=0.5D0*sml_dt
  case(3)
     dt_now=sml_dt
  case(4)
     dt_now=sml_dt/6D0
  end select
#endif
 else 
  !save phase0 information
  select case(ipc)
  case(1)
     dt_now=0.5D0*dt
     !$OMP PARALLEL DO &
     !$OMP PRIVATE( ITH, I )
     do ith=1,sml_nthreads
        do i=i_beg(ith), i_end(ith)
           phase0(:,i)=ptl(i)%ph
        enddo
     enddo
  case(2)
     dt_now=dt
  end select
 endif   
  
  !$OMP PARALLEL DO &
  !$OMP PRIVATE( ITH, I, NEW_PHASE, &
  !$OMP          OLD_PHASE, RTN )
  do ith=1, sml_nthreads
     call t_startf("PUSH_LOOP")
     do i=i_beg(ith), i_end(ith)
        ! for alive particles only
        if(sp%ptl(i)%gid>0) then
           
           !******************************************************
           ! actual particle push           
           !******************************************************
           call push_single(grid,psn,sp,i,phase0(:,i),new_phase,dt_now,ith,diag_on,ipc)
           
           ! check r-z boundary validity and update psi variables
           if(new_phase(1)<eq_min_r .or. new_phase(1)>eq_max_r .or. new_phase(2)<eq_min_z .or. new_phase(2)>eq_max_z)then
              call remove_particle(sp,i,-1,ith)
!              print *, 'particle eliminated due to rz_outside :', i, sml_mype, sp%type, sp%ptl(i)%gid, new_phase(1),new_phase(2)
           else                            
              ! bounce 
              if(ipc==sml_nrk .and. sml_bounce/=0) then
                 old_phase(:)=phase0(:,i)
                 call bounce(new_phase,old_phase,rtn)
                 if(rtn<0)  then
                    call remove_particle(sp,i,-2,ith)
                 endif
              endif
              
              !******************************************************
              ! time advance one step
              !******************************************************
              ptl(i)%ph= new_phase(:)
              
           endif
           
        endif
     enddo
     call t_stopf("PUSH_LOOP")
  enddo
  
end subroutine push


! single particle push -- get new_phase using rk4 with initial E-field
subroutine push_single(grid,psn,sp,i,y,new_phase,dt,ith,diag_on,ipc)
  use sml_module
  use ptl_module
  use fld_module
  use grid_class
  use psn_class
  use perf_monitor
  implicit none
  type(grid_type), intent(in) :: grid
  type(psn_type), intent(in) :: psn
  type(species_type), intent(inout) :: sp
  real (8), intent(in) :: y(ptl_nphase)
  real (8), intent(out) :: new_phase(ptl_nphase)
  real (8), intent(in) :: dt
  integer, intent(in) :: ith
  logical, intent(in) :: diag_on
  integer, intent(in) :: ipc !! RK4 index
  !
  type(ptl_type) :: ptli
  type(fld_type) :: fld
  real (8) :: dy(ptl_nphase),dyt(ptl_nphase),dym(ptl_nphase)
  integer :: i,rtn,j,i1
  logical :: rz_outside
  real (8) :: vf_diag(sml_n_vf_diag)
  real (kind=8) , external :: psi_interpol
  character (len=5) :: err_str(0:1)
  real (8) :: hdt, time, th  ! -- time should be input 
  real (8) :: bp,B, E_mag(3)
  err_str(1)='ion'
  err_str(0)='elec'

  
  time=sml_time ! need to update -- rk4 time ###

  hdt=dt*0.5D0
  th=time + hdt

  if(sp%ptl(i)%gid>0) then
     
     !set ptli -- ptli%ct does not change
     ptli%ct=sp%ptl(i)%ct
     ptli%ph=sp%ptl(i)%ph
     ptli%gid=sp%ptl(i)%gid 
     !get derivs with updating E-field information - assigned on fld
     !diag-ports are called when diag_on is .true.
     call derivs_single(grid,psn,sp,ptli,dy,i,time,fld,ith,diag_on)


     
#ifdef PURE_RK2
     ! Only simple RK2
     
     new_phase = y + dt * dy     
     !call restrict_weight(new_phase(piw1:piw2))
#else
#ifdef PURE_RK4
     ! update phase variable -- RK4
     ! replace derivs by dy
     ! replace dt_now by dt
     ! replace ptli%ph= by new_phase=
     ! replace ptli%phi0 by y
     select case(ipc)
     case(1)
        new_phase= y + dy*dt ! y0 + k1*h/2
        !ptli%ph = y + dy*dt
        !ptli%dph = dy
        sp%ptl(i)%dph=dy     ! k1
     case(2)
        new_phase= y + dy*dt ! y0 + k2*h/2
        !ptli%ph = y + dy*dt
        !ptli%dpht= dy
        sp%ptl(i)%dpht=dy    ! k2
     case(3)
        new_phase= y + dy*dt ! y0 + k3*h
        !ptli%ph = y + dy*dt
        !ptli%dpht= dy+ptli%dpht
        sp%ptl(i)%dpht=dy+sp%ptl(i)%dpht !k4 +k3
     case(4)
        new_phase= y + dt*(sp%ptl(i)%dph + dy + 2D0*sp%ptl(i)%dpht) !y(t+h)=y(t)+h(k1+k4+2*(k3+k4))/6
     end select
#else

     ! RK2 - RK4 hybrid -- time advance with RK4 with time-constant E-field
     
     ! get E-field in magnetic field
     bp=sqrt(fld%br**2 + fld%bz**2)     
     B=sqrt(bp**2 + fld%Bphi**2 )
     E_mag(2)=(fld%Er*fld%Br + fld%Ez*fld%Bz)/bp
     E_mag(3)=(E_mag(2)*bp   + fld%Ephi*fld%Bphi)/B    ! parallel field
     E_mag(1)=(fld%Er*fld%dpsidr + fld%Ez*fld%dpsidz )/(fld%dpsidr**2 + fld%dpsidz**2)
     
     
     ! get derivs with existing E-field
     ptli%ph=y
     call derivs_single_with_e(sp,ptli,dy ,i,time    ,fld,E_mag,ith)
     
     ptli%ph = y + hdt * dy
     call derivs_single_with_e(sp,ptli,dyt,i,th      ,fld,E_mag,ith)
     
     ptli%ph = y + hdt * dyt
     call derivs_single_with_e(sp,ptli,dym,i,th      ,fld,E_mag,ith)
     
     ptli%ph = y + dt * dym
     dym = dyt + dym
     call derivs_single_with_e(sp,ptli,dyt,i,time+dt ,fld,E_mag,ith)
     
     ! Obtain new_phase
     new_phase = y + dt/6D0 * ( dy + dyt + 2D0*dym )
     !call restrict_weight(new_phase(piw1:piw2))
#endif
#endif
  endif
  
  contains
    ! set minimum weight to prevent weight to be smaller than 0
    subroutine restrict_weight(w)
      implicit none
      real (8) :: w(2)
      real (8), parameter :: weight_min = -50D0
      real (8), parameter :: weight_max = 0.999D0

      w(1)=max(w(1),weight_min)
      w(2)=max(w(2),weight_min)
      
      w(1)=min(w(1),weight_max)
      w(2)=min(w(2),weight_max)


    end subroutine restrict_weight
end subroutine push_single

!obtain derivatives of phase variable : actual calculation is done in derivs_sp. 
! prepare E-field and B-field
subroutine derivs_single(grid,psn,sp,ptli,dy,i,time,fld,ith,diag_on)
  use sml_module
  use grid_class
  use psn_class
  use fld_module
  use ptl_module
  implicit none
  type(grid_type), intent(in) :: grid
  type(psn_type), intent(in) :: psn
  type(species_type), intent(in) :: sp
  type(ptl_type),intent(in) :: ptli
  real (8), intent(out) :: dy(ptl_nphase)
  integer, intent(in) :: i
  real (8), intent(in) :: time
  type(fld_type), intent(out) :: fld
  integer, intent(in) :: ith
  logical, intent(in) :: diag_on
  !
  logical rz_outside
  real (8) :: vf_diag(sml_n_vf_diag)   ! variables for diagnosis 

  ! Save space information
  fld%r=ptli%ph(1)
  fld%z=ptli%ph(2)
  fld%phi=ptli%ph(3)

  ! obtain B-field information 
  call field(fld,time,rz_outside)
  if(.not. rz_outside) then
     ! obtain gyro-averaged E-field for each particle : use position information from 'charge'
     call efield(grid,psn,sp,i,fld,time)
!     if(sml_deltaf) call f0_info(fld)
     
     ! get time derivatives
     call derivs_sp(fld,time,ptli,sp%type,dy,diag_on,vf_diag)
     
  else
     call remove_particle(sp,i,-1,ith)
!     print *, 'particle eliminated due to rz_outside', i, sml_mype, sp%type, sp%ptl(i)%gid
  endif
     
!call  check_point('before port1')
  if(diag_on) call diag_1d_port1(sp%ptl(i),dy,sp%type,vf_diag,ith)
!call  check_point('after port1')
end subroutine derivs_single

subroutine derivs_single_with_e(sp,ptli,dy,i,time,fld,E_mag,ith)
  use sml_module
  use fld_module
  use ptl_module
  implicit none
  type(species_type), intent(in) :: sp
  type(ptl_type),intent(in) :: ptli
  real (8) :: dy(ptl_nphase)
  integer, intent(in) :: i
  real (8), intent(in) :: time 
  type(fld_type), intent(inout) :: fld
  real (8), intent(in) :: E_mag(3)
  integer, intent(in) :: ith
  !
  logical :: rz_outside
  real (8) :: vf_diag(sml_n_vf_diag)
  logical, parameter :: diag_on=.false.
  !
  real (8) :: dpsi(2), E(3), bp, dtheta_norm(2), B

  ! Save space information
  fld%r=ptli%ph(1)
  fld%z=ptli%ph(2)
  fld%phi=ptli%ph(3)

  ! obtain B-field information 
  call field(fld,time,rz_outside)
  if(.not. rz_outside) then
     ! obtain gyro-averaged E-field from previous E-field

     ! same E-field considering B-field curvature
     dpsi(1:2)=(/ fld%dpsidr, fld%dpsidz /)
     E(1:2)=E_mag(1)*dpsi

     bp=sqrt(fld%br**2+fld%bz**2)
     dtheta_norm(1:2)=(/ fld%br, fld%bz /)/bp
     E(1:2)=E(1:2) + E_mag(2)*dtheta_norm

     B=sqrt(fld%br**2+fld%bz**2+fld%bphi**2)
     E(3)=(E_mag(3)*B- E(1)*fld%Br - E(2)*fld%Bz)/fld%Bphi

     fld%Er=E(1)
     fld%Ez=E(2)
     fld%Ephi=E(3)
     
!     if(sml_deltaf) call f0_info(fld)
     
     ! get time derivatives
     call derivs_sp(fld,time,ptli,sp%type,dy,diag_on,vf_diag)
     
  else
     call remove_particle(sp,i,-1,ith)
!     print *, 'particle eliminated due to rz_outside', i, sml_mype, sp%type, sp%ptl(i)%gid
  end if
end subroutine derivs_single_with_e


subroutine shift_sp(grid,psn,sp)
  use sml_module, only : sml_electron_on, sml_pol_decomp
  use ptl_module
  use grid_class
  use psn_class
  use pol_decomp_module
  use perf_monitor
  use col_module
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  integer :: isp
  type(species_type) :: sp
  integer :: i
  logical :: finished

  if(.NOT. sml_pol_decomp) then
     call t_startf("SHIFT_IE")
     call shift_ie(grid,psn,sp,sp%phase0,sp%ptl,sp%shift_opt)
     call t_stopf("SHIFT_IE")
#ifdef VPIC_COL
     if(col_mode .eq. 3) then
        call t_startf("COL_3_DECOMP")
        call col_3_decomp(grid,sp)  ! collision 3 without f0_grid
        call t_stopf("COL_3_DECOMP")
     endif
#endif
  else
     call t_startf("SHIFT_IE")
     call shift_ie(grid,psn,sp,sp%phase0,sp%ptl,sp%shift_opt)
     call t_stopf("SHIFT_IE")
  endif
  
end subroutine shift_sp
! Code from gtc
subroutine shift_ie_tonly(grid,sp)
  use sml_module
  use ptl_module
  use grid_class
  use diag_module  
  use perf_monitor
  implicit none
  include 'mpif.h'
  type(grid_type) :: grid
  type(species_type) :: sp
  !
  integer i,m,msendleft(2),msendright(2),mrecvleft(2),mrecvright(2),mtop,&
       m0,msend,mrecv,idest,isource,isendtag,irecvtag,&
       isendcount,irecvcount,istatus(MPI_STATUS_SIZE),ierror,iteration
  integer,dimension(:),allocatable :: ihole
  real (kind=8),dimension(:,:),allocatable :: recvleft,recvright,sendleft,sendright
  integer (kind=8),dimension(:),allocatable :: recvleftid,recvrightid,sendleftid,sendrightid
  real (kind=8) :: zetaright,zetaleft,zetamid,pi_inv
  character(len=8) cdum
  integer :: st,last
  integer :: ptl_shiftmax
  integer, pointer :: ptl_rshift(:), ptl_lshift(:),ptl_hole(:)
  integer :: sendsize,arr_start, arr_end, pc_index
  integer :: iphs, iphe, icts, icte, igid, iph0s, iph0e
#ifdef PURE_RK4
  integer :: idphs, idphe, idphts, idphte
#endif
  logical :: init_alloc=.false.
  save ptl_rshift,ptl_lshift,ptl_hole, init_alloc

  !init
  call t_startf("SHIFT_IET_INIT")
  if(.not. init_alloc) then
     init_alloc=.true.
     ptl_shiftmax=sp%maxnum 
     allocate(ptl_rshift(ptl_shiftmax),ptl_lshift(ptl_shiftmax),ptl_hole(2*ptl_shiftmax))
  endif
      
  iphs=1
  iphe=ptl_nphase
  
  icts=iphe+1
  icte=iphe+ptl_nconst  
  
  iph0s=icte+1
  iph0e=icte+ptl_nphase  
  sendsize=iph0e

#ifdef PURE_RK4  
  idphs=iph0e+1
  idphe=iph0e+ptl_nphase
  
  idphts=idphe+1
  idphte=idphe+ptl_nphase
  sendsize=idphte
#endif
 
  !igid=idphte+1

  ! Modulo operation
  do i=1, sp%num
     if(sp%ptl(i)%ph(3)>= sml_2pi_wedge_n .or. sp%ptl(i)%ph(3)< 0D0 ) then
        sp%ptl(i)%ph(3)=modulo(sp%ptl(i)%ph(3),sml_2pi_wedge_n)
     endif
  enddo

  ! Shift operation
  pi_inv=1.D0/sml_pi
  m0=1
  iteration=0
  mrecv=1
  call t_stopf("SHIFT_IET_INIT")
  
  INFLOOP: do while (1==1)
     
     iteration=iteration+1
     if(iteration>sml_totalpe)then
        print *,'endless particle sorting loop at PE=',sml_mype
        stop
     endif
     
     msend=0
     msendleft=0
     msendright=0
     !        print *, 'iter=',iteration,sml_mype
     
     !************* Phase 1 *******************************
     ! Finding particle to be shifted
     !*****************************************************
     call t_startf("SHIFT_IET_PHASE1")
     if(m0 <= sp%num)then
        
        do m=m0,sp%num
           zetaright=min(sml_2pi_wedge_n,sp%ptl(m)%ph(3))-grid%phimax
           zetaleft=sp%ptl(m)%ph(3)-grid%phimin
           
           if( zetaright*zetaleft > 0 )then
              !                 print *, '[',sml_mype,']',msend,phase(m,3)*pi_inv,grid%phimax*pi_inv, grid%phimin*pi_inv
              msend=msend+1
              ptl_hole(msend)=m
              zetaright=zetaright*0.5D0*pi_inv
              zetaright=zetaright-real(floor(zetaright))
              
              if( zetaright < 0.5 )then
                 ! # of particle to move right
                 msendright(1)=msendright(1)+1
                 ptl_rshift(msendright(1))=m
                 ! # of particle to move left
              else
                 msendleft(1)=msendleft(1)+1
                 ptl_lshift(msendleft(1))=m
              endif
           endif
        enddo
     endif
     call t_stopf("SHIFT_IET_PHASE1") 
    
     ! total # of particles to be shifted for whole CPUs
     call t_startf("SHIFT_IET_RED")
     mrecv=0
     msend=msendright(1)+msendleft(1)
     call MPI_ALLREDUCE(msend,mrecv,1,MPI_INTEGER,MPI_SUM,sml_comm,ierror)
     call t_stopf("SHIFT_IET_RED")

     !        print *, '[',sml_mype,']','msend,right,left,recv',msend,msendright(1),msendleft(1),mrecv
     ! no particle to be shifted, free memory and return
     
     if ( mrecv == 0 ) then
        if(msend/=0)  then
           print *,'Error in shfit', msendright(1), msendleft(1), sml_mype
           stop
        endif
        !           print *, sml_mype,'interation # :', iteration
        !     if(sml_mype==5) print *, m0,ptl_num
        !           print *, pnum,sml_mype
        
        call shift_check(grid,sp) ! debug only
        
        exit INFLOOP
     endif
     
     call t_startf("SHIFT_IET_PACK")
     ! an extra space to prevent zero size when msendright(1)=msendleft(1)=0
     !        print *, 'mem allocation',sml_mype !debug
     !allocate memory for particle sending
     sendsize=ptl_nphase2 + ptl_nconst
#ifdef PURE_RK4
     sendsize=sendsize + ptl_nphase2  ! for RK4 
     !Julien 1. I checked and this computation of sendsize is equivalent to previous calculation, so why doing it. 
     !       2. I don't know why RK4 was always on.
#endif
     allocate(sendright(sendsize,max(msendright(1),1)),sendleft(sendsize,max(msendleft(1),1)))
     allocate(sendrightid(max(msendright(1),1)),sendleftid(max(msendleft(1),1)))
     
     ! pack particle to move right
     do m=1, msendright(1)
        sendright(iphs:iphe,m)=sp%ptl(ptl_rshift(m))%ph
        sendright(icts:icte,m)=sp%ptl(ptl_rshift(m))%ct
        sendrightid(m)=sp%ptl(ptl_rshift(m))%gid
        sendright(iph0s:iph0e,m)=sp%phase0(:,ptl_rshift(m))
#ifdef PURE_RK4
        sendright(idphs:idphe,m)=sp%ptl(ptl_rshift(m))%dph
        sendright(idphts:idphte,m)=sp%ptl(ptl_rshift(m))%dpht
#endif
     enddo
     do m=1, msendleft(1)
        sendleft(iphs:iphe,m)=sp%ptl(ptl_lshift(m))%ph
        sendleft(icts:icte,m)=sp%ptl(ptl_lshift(m))%ct
        sendleftid(m)=sp%ptl(ptl_lshift(m))%gid
        sendleft(iph0s:iph0e,m)=sp%phase0(:,ptl_lshift(m))           
#ifdef PURE_RK4
        sendleft(idphs:idphe,m)=sp%ptl(ptl_lshift(m))%dph
        sendleft(idphts:idphte,m)=sp%ptl(ptl_lshift(m))%dpht
#endif
     enddo
     !        print *, '[',sml_mype,']','sendright',pnum,mtop,last
     
     mtop=sp%num
     ! # of particles remain on local PE
     sp%num=sp%num-msendleft(1)-msendright(1)
     ! fill the hole
     last=msend
     !        print *, '[',sml_mype,']','pnum,mtop,last',pnum,mtop,last
     FILLHOLE : do i=1, msend
        m=ptl_hole(i)
        if( m > sp%num )  exit FILLHOLE
        !when empty space in the end - possible only for i=1
        do while( mtop == ptl_hole(last) )
           mtop=mtop-1
           last=last-1
        enddo
        sp%ptl(m)  = sp%ptl(mtop)
        sp%phase0(:,m) = sp%phase0(:,mtop)
#ifdef PURE_RK4
        sp%ptl(m)%dph=sp%ptl(mtop)%dph
        sp%ptl(m)%dpht=sp%ptl(mtop)%dpht
#endif
        mtop=mtop-1
        if(mtop == sp%num) exit FILLHOLE
     enddo FILLHOLE
     call t_stopf("SHIFT_IET_PACK")
     
     !*************** To Right ****************************
     !*****************************************************
     !        print *, 'send_num right',sml_mype
     ! send # of particle to move right
     call t_startf("SHIFT_IET_SR_RL")
     mrecvleft=0
     !        idest=mod(sml_mype+1,sml_totalpe)
     !        isource=mod(sml_mype-1+sml_totalpe,sml_totalpe)
     idest=modulo(sml_intpl_mype+1,sml_intpl_totalpe)
     isource=modulo(sml_intpl_mype-1,sml_intpl_totalpe)
     isendtag=sml_intpl_mype
     irecvtag=isource
     call MPI_SENDRECV(msendright,2,MPI_INTEGER,idest,isendtag,&
          mrecvleft,2,MPI_INTEGER,isource,irecvtag,sml_intpl_comm,istatus,ierror)
     
     allocate(recvleft(sendsize,max(mrecvleft(1),1)))
     recvleft=0D0 !debug
     allocate(recvleftid(max(mrecvleft(1),1)))
     recvleftid=0 !debug
     
     !        print *, 'send_data right',sml_mype,mrecvleft(1),msendright(1)
     ! send particle to right and receive from left
     recvleft=0D0
     isendcount=msendright(1)*sendsize
     irecvcount=mrecvleft(1)*sendsize
     call MPI_SENDRECV(sendright,isendcount,MPI_REAL8,idest,&
          isendtag,recvleft,irecvcount,MPI_REAL8,isource,&
          irecvtag,sml_intpl_comm,istatus,ierror)
     call MPI_SENDRECV(sendrightid,msendright(1),MPI_INTEGER8,idest,&
          isendtag,recvleftid,mrecvleft(1),MPI_INTEGER8,isource,&
          irecvtag,sml_intpl_comm,istatus,ierror)
     
     
     !        call MPI_BARRIER(sml_comm,ierror)
     
     
     
     !************** To Left *****************************
     !****************************************************
     !        print *, 'send_num left',sml_mype
     ! send # of particle to move left
     mrecvright=0
     !        idest=mod(sml_mype-1+sml_totalpe,sml_totalpe)
     !        isource=mod(sml_mype+1,sml_totalpe)
     idest=modulo(sml_intpl_mype-1,sml_intpl_totalpe)
     isource=modulo(sml_intpl_mype+1,sml_intpl_totalpe)
     isendtag=sml_intpl_mype
     irecvtag=isource
     call MPI_SENDRECV(msendleft,2,MPI_INTEGER,idest,isendtag,&
          mrecvright,2,MPI_INTEGER,isource,irecvtag,sml_intpl_comm,istatus,ierror)
     
     allocate(recvright(sendsize,max(mrecvright(1),1)))
     recvright=0D0 !debug
     allocate(recvrightid(max(mrecvright(1),1)))
     recvrightid=0 !debug
     
     ! send particle to left and receive from right
     recvright=0.0
     isendcount=msendleft(1)*sendsize
     irecvcount=mrecvright(1)*sendsize
     call MPI_SENDRECV(sendleft,isendcount,MPI_REAL8,idest,&
          isendtag,recvright,irecvcount,MPI_REAL8,isource,&
          irecvtag,sml_intpl_comm,istatus,ierror)
     call MPI_SENDRECV(sendleftid,msendleft(1),MPI_INTEGER8,idest,&
          isendtag,recvrightid,mrecvright(1),MPI_INTEGER8,isource,&
          irecvtag,sml_intpl_comm,istatus,ierror)
     
     
     call t_stopf("SHIFT_IET_SR_RL")
     
     !        call mpi_barrier(sml_comm,ierror) !debug
     !        print *, '[',sml_mype,']','mrecvl,mrevr',mrecvleft(1),mrecvright(1)

     
     
     ! need extra particle array
     if(sp%num+mrecvleft(1)+mrecvright(1) > sp%maxnum)then
        print *, 'not enough ptl array size', sp%num,sp%num+mrecvleft(1)+mrecvright(1),sp%maxnum,sml_mype, sp%type
        call MPI_ABORT(sml_comm,1,ierror)
     endif
     
     !**************** Unpack *********************88
     !***********************************************
     call t_startf("SHIFT_IET_UNPACK")
     
!        print *, 'particle  unpack', sml_mype
     ! unpack particle, particle moved from left
     
     do m=1,mrecvleft(1)
        sp%ptl(m+sp%num)%ph=recvleft(iphs:iphe,m)
        sp%ptl(m+sp%num)%ct=recvleft(icts:icte,m)
        sp%ptl(m+sp%num)%gid=recvleftid(m)
        sp%phase0(:,m+sp%num)=recvleft(iph0s:iph0e,m)
#ifdef PURE_RK4
        sp%ptl(m+sp%num)%dph=recvleft(idphs:idphe,m)
        sp%ptl(m+sp%num)%dpht=recvleft(idphts:idphte,m)
#endif
     enddo
     sp%num=sp%num+mrecvleft(1)
     
     ! particle moved from right
     do m=1,mrecvright(1)
        sp%ptl(m+sp%num)%ph =recvright(iphs:iphe,m)
        sp%ptl(m+sp%num)%ct =recvright(icts:icte,m)
        sp%ptl(m+sp%num)%gid=recvrightid(m)
        sp%phase0(:,m+sp%num)=recvright(iph0s:iph0e,m)
#ifdef PURE_RK4
        sp%ptl(m+sp%num)%dph=recvright(idphs:idphe,m)
        sp%ptl(m+sp%num)%dpht=recvright(idphts:idphte,m)
#endif
     enddo
     sp%num=sp%num+mrecvright(1)
     
     
     deallocate(sendleft,sendright,recvleft,recvright)
     deallocate(sendleftid,sendrightid,recvleftid,recvrightid)
     m0=sp%num-mrecvright(1)-mrecvleft(1)+1
!        call MPI_BARRIER(sml_comm,ierror)
     call t_stopf("SHIFT_IET_UNPACK")
     
  enddo INFLOOP

end subroutine shift_ie_tonly

subroutine shift_check(grid,sp)
  use ptl_module
  use grid_class
  use sml_module
  implicit none
  integer :: i
  type(grid_type):: grid
  type(species_type) :: sp
  
  do i= 1, sp%num
     if(grid%phimin > sp%ptl(i)%ph(3) .or. grid%phimax < sp%ptl(i)%ph(3) ) then
        print *, 'err in shift check', sp%ptl(i)%ph(3),grid%phimin, grid%phimax, sml_mype,i,sp%type
     endif
  enddo

end subroutine shift_check



subroutine shift_check2(grid,sp,message,flag)
  use ptl_module
  use grid_class
  use sml_module
  implicit none  
  character (len=30) :: message
  integer :: i,flag
  type(grid_type):: grid
  type(species_type) sp
  integer :: count
  

  count=0
  do i= 1, sp%num
     if(grid%phimin > sp%ptl(i)%ph(3) .or. grid%phimax < sp%ptl(i)%ph(3) ) then
        count=count+1
     endif
  enddo
  if(count>0) then
     print *, 'shift check :: ',message, count, sml_mype,sp%type,flag
  endif

end subroutine shift_check2




!!****************************************************************************
!!> derivatives of electron phase
!!
!!<***************************************************************************
subroutine derivs_sp(fld,t, ptli, sp_type, yprime, diag_on, vf_diag)
  use fld_module
  use sml_module
  use ptl_module
  use eq_module
  implicit none
  type(fld_type), intent(in) :: fld   ! field variable
  real (8),       intent(in) :: t     ! time
  type(ptl_type), intent(in) :: ptli  ! particle info
  integer,  intent(in)  :: sp_type      ! particle species type (ion/elec)
  real (8), intent(out) :: yprime(ptl_nphase)  !
  logical,  intent(in)  :: diag_on
  real (8), intent(out) :: vf_diag(sml_n_vf_diag)   ! variables for diagnosis 
  !
  real (8) :: mass, charge,c_m  ! charge and mass 
  real (8) :: r, z, phi, rho, mu,inv_r
  real (8) :: B, B2, over_B, over_B2
  real (8) :: D, nb_curl_nb, curl_B(3)
  real (8) :: dbdr, dbdz, dbdphi
  real (8) :: cmrho2, cmrho, murho2b, murho2b_c, vp
  real (8) :: fr, fp, fz
  real (8) :: fr_exb, fp_exb, fz_exb, yp_exb(3)

  ! for weight calculation
  real (8) :: energy, vmag, pitch, dvpdt, denergy_dt, dpitch_dt
  real (8) :: psi, den, dden, temp, dtemp, dfdp, tmp, envelop, df0(5)
  real (8) :: one_m_w
  real (8) :: total_ddpotdt    ! for electron -- from adiabatic reponse
  !
  real (8) :: bp, f0_1_lt

  ! some global parameter to be
  real (kind=8) :: sml_wdot_energy_max, sml_f0_psi_c, sml_f0_1_psi_w

  !
  real (8) :: i_factor, ep_r, ep_p, ep_z

  
  ! these are constant for deltaf method -- maybe should be moved to sml_module
  sml_wdot_energy_max=10D0
  sml_f0_psi_c=0.5*(sqrt(sml_inpsi)+sqrt(sml_outpsi))
  sml_f0_1_psi_w=1D0/( 0.4*(sqrt(sml_outpsi)-sqrt(sml_inpsi)) )

  ! prepare constant
  mass=ptl_mass(sp_type) !
  charge=ptl_charge(sp_type) !-1.D0/ptl_charge
  c_m=charge/mass

  r=ptli%ph(pir)
  z=ptli%ph(piz)
  phi=ptli%ph(pip)
  rho=ptli%ph(pirho)
  mu=ptli%ct(pim)
  if (sml_cylindrical) then
    inv_r=1D0/eq_axis_r
  else
    inv_r=1D0/r
  endif

  B = sqrt( fld%br**2 + fld%bz**2 + fld%bphi **2 )
  b2=b*b
  over_B=1/B
  over_B2=over_B*over_B
  
  ! normalized b dot curl of normalized b
  if (sml_cylindrical) then
    ! Cylindrical limit
    nb_curl_nb= -1.D0*over_B2 * ( fld%dbpdz*fld%br + (fld%dbzdr-fld%dbrdz)*fld%bphi - fld%dbpdr*fld%bz )
  else
    nb_curl_nb= -1.D0*over_B2 * ( fld%dbpdz*fld%br + (fld%dbzdr-fld%dbrdz)*fld%bphi - (fld%bphi/r  + fld%dbpdr)*fld%bz )
  endif
!  D=1.D0/ ( 1.D0 + rho * nb_curl_nb )
!#warning 'D=1.0'
  D=1.D0
      
  dbdr = ( fld%br*fld%dbrdr + fld%bphi*fld%dbpdr + fld%bz*fld%dbzdr) *over_B
  dbdz = ( fld%br*fld%dbrdz + fld%bphi*fld%dbpdz + fld%bz*fld%dbzdz) *over_B
  dbdphi=0D0  ! no B perturbation

  ! R
  curl_B(1)  = fld%dbzdp*inv_r - fld%dbpdz
  ! Z
  if (sml_cylindrical) then
    curl_B(2)  = fld%dbpdr-fld%dbrdp*inv_r
  else
    curl_B(2)  = fld%bphi*inv_r + fld%dbpdr-fld%dbrdp*inv_r
  endif
  ! phi
  curl_B(3)  = fld%dbrdz - fld%dbzdr

  ! bug fix -- rho**2 B grad B term added 2002/1/15
  !! optimization variables
  
  cmrho2=c_m*rho**2
  cmrho =c_m*rho
  murho2b=(mu+charge*cmrho2 *B)
  murho2b_c=murho2b/charge
  vp=cmrho*B


  ! F is grad H/q
  fr = fld%Er - (murho2B_c) * dbdr ! fr is force over q . -gradient of  Hamiltonian
  fp = fld%Ephi - (murho2B_c) * dbdphi/r  ! modified by shlee 5/30/2001 -- /r added
  fz = fld%Ez - (murho2B_c) * dbdz
  
  ! ignore all drifts when it close to divertor
  ! ion & electron --> pushe for electron
  ! now drift and parallel E-field are ignored.
  if(sml_ignore_drift_near_wall ) then 
     call ignore_factor(r,z,i_factor)
     !ep_r = (fld%Er  *fld%br  )*fld%br  /b2
     !ep_p = (fld%Ephi*fld%bphi)*fld%bphi/b2 
     !ep_z = (fld%Ez  *fld%bz  )*fld%bz  /b2 
     
     fr = fr*i_factor !+ ep_r*(1D0-i_factor)
     fp = fp*i_factor !+ ep_p*(1D0-i_factor)
     fz = fz*i_factor !+ ep_z*(1D0-i_factor)     
  endif
  
  yprime(1)= D*( (fld%bz*Fp - fld%Bphi * Fz) * over_B2         &
       +  cmrho * fld%br                       &
       +  cmrho2 * curl_B(1) )
  yprime(piz)= D*( (fld%bphi * fr - fld%br * fp ) * over_B2      &
       +  cmrho * fld%bz                       &
       +  cmrho2 * curl_B(2) )
  yprime(pip)= D*( (fld%br * fz - fld%bz * fr) * over_B2         &
       +  cmrho * fld%bphi                     &
       +  cmrho2 * curl_B(3) ) * inv_r


  fr_exb = fld%Er  ! fr is force over q . -gradient of  Hamiltonian
  fp_exb = fld%Ephi! modified by shlee 5/30/2001 -- /r added
  fz_exb = fld%Ez 
  
  yp_exb(1)= D*(fld%bz   *fp_exb - fld%Bphi * fz_exb) * over_B2  
  yp_exb(2)= D*(fld%bphi *fr_exb - fld%br   * fp_exb) * over_B2      
  yp_exb(3)= D*(fld%br   *fz_exb - fld%bz   * fr_exb) * over_B2 *inv_r

  yprime(pirho)=D*over_B2 *( &
       fld%br*fr + fld%bz*fz + fld%bphi*fp &
       + rho*( fr*curl_B(1) + fz*curl_B(2) + fp*curl_B(3)) )
  
  if( ptl_deltaf_sp(sp_type) .and. sml_extra_dwdt ) then
     energy=(0.5D0*mass*vp*vp + mu*b)
     vmag=sqrt(2D0*energy/mass)
     pitch = vp/vmag
     dvpdt=c_m*yprime(pirho)*B + (dbdr*yprime(1)+dbdz*yprime(2))*vp/B  !
     denergy_dt=mu*(dbdr*yprime(1)+dbdz*yprime(2)+dbdphi*yprime(3))+(mass*vp*dvpdt)
     dpitch_dt=dvpdt/vmag - 0.5D0*pitch*denergy_dt/energy


          ! dlog(f0) ------ need to be separate routine later     
     ! temp routine -- finish it later
     ! no space gradient
     psi=fld%psi

     den=eq_ftn(psi,r,z,eq_den)
     dden=eq_dftn(psi,r,z,eq_den)
     if(sp_type==0) then !electron
        temp=eq_ftn(psi,r,z,eq_tempe)*sml_e_charge
        dtemp=eq_dftn(psi,r,z,eq_tempe)*sml_e_charge
     else
        temp=eq_ftn(psi,r,z,eq_tempi)*sml_e_charge
        dtemp=eq_dftn(psi,r,z,eq_tempi)*sml_e_charge
     endif
     
     if(sml_deltaf_f0_mode==-2) then  !consistent grad f with real profile
        dfdp=dden/den + dtemp/temp*(energy/temp - 1.5D0)
     elseif(sml_deltaf_f0_mode==-1) then
        tmp=1D0/sqrt(fld%dpsidr**2+fld%dpsidz**2)  ! drdpsi 
        !envelop= exp( - ((sqrt(psi) - sml_f0_psi_c)*sml_f0_1_psi_w )**8 )
        !envelop= exp( - ((sqrt(psi) - sqrt(0.387))*100D0 )**8 )
        envelop=cosh((sqrt(psi) - sqrt(0.387))/0.0407)**(-2); 
        !envelop=1D0
        if(sp_type==0) then
           f0_1_Lt=sml_f0_1_Lt_e
        else
           f0_1_Lt=sml_f0_1_Lt
        endif
        dfdp=-envelop*(sml_f0_1_Ln + f0_1_Lt*(energy/temp - 1.5D0))*tmp
     endif
     df0(1)=dfdp *fld%dpsidr !/(2D0*q+1D-99)
     df0(2)=dfdp *fld%dpsidz !/(2D0*q+1D-99)
     df0(3)=0D0  * r     ! r : yprime(3) is phi dot --> dimension matching
#ifdef DELTAF_MODE2
     df0(4)=-1D0/temp ! energy deriv
#else
     df0(4)=0D0 ! zero for extra dwdt
#endif
     df0(5)=0D0              ! pitch deriv

     if(sml_dwdt_fix_bg)then 
        one_m_w=1D0
     else
        one_m_w=1D0 - ptli%ph(piw2)
     end if

     if( .not. sml_dwdt_exb_only ) then
        yprime(piw1)= -one_m_w* (&
             df0(1)*yprime(1) + df0(2)*yprime(2) + df0(3)*yprime(3) + df0(4)*denergy_dt + df0(5)*dpitch_dt &
             )
#ifndef DELTAF_MODE2
     else if(sml_deltaf_f0_mode == -1) then
#else
     else if(sml_deltaf_f0_mode == -1 .or. sml_deltaf_f0_mode == -2) then  ! else-if below will be ignored.
#endif
        yprime(piw1)= -one_m_w* (&
             df0(1)*yp_exb(1) + df0(2)*yp_exb(2) + df0(3)*yp_exb(3) + df0(4)*denergy_dt + df0(5)*dpitch_dt &
             )

     else if(sml_deltaf_f0_mode == -2) then  ! remove weight evolution from v_grad B  - total will be updated later
        yprime(piw1)= -one_m_w* (&
             df0(1)*(yp_exb(1)-yprime(1)) + df0(2)*(yp_exb(2)-yprime(2)) )

     endif
     yprime(piw1) = (1D0-sml_krook_rate)*yprime(piw1)
     ! electron -- minus adibatic response
     if(sp_type==0) then
        total_ddpotdt=fld%ddpotdt - yprime(1)*(fld%Er-fld%Er00) - yprime(2)*(fld%Ez-fld%Ez00)  -r*yprime(3)*fld%Ephi

        yprime(piw1) = yprime(piw1) - one_m_w*total_ddpotdt/temp*sml_e_charge
     endif
     yprime(piw2)=yprime(piw1)
  else
     yprime(piw1:piw2)=0D0
  endif
  if(diag_on) then
     vf_diag(1)=fld%psi
     vf_diag(2)=B
     vf_diag(3)=fld%Bphi
     vf_diag(4)=yp_exb(1)*fld%dpsidr + yp_exb(2)*fld%dpsidz   ! V_ExB dot grad psi
     vf_diag(5)=yprime(1)*fld%dpsidr + yprime(2)*fld%dpsidz   ! V_d dot grad psi
     vf_diag(6)=fld%dpsidr*fld%dpsidr + fld%dpsidz*fld%dpsidz  ! |grad psi|^2
     bp=sqrt( fld%br*fld%br + fld%bz*fld%bz ) 
     vf_diag(7)=(yprime(1)*fld%br + yprime(2)*fld%bz)/bp
     vf_diag(8)=(yp_exb(1)*fld%br + yp_exb(2)*fld%bz)/bp
  endif

end subroutine derivs_sp

    ! return factor=1D0 for normal 
    ! return factor=0D0 for near sheath
subroutine ignore_factor(x,y,factor)      
  use sml_module
  implicit none
  real (8), intent(in) :: x, y
  real (8), intent(out) :: factor
  real (8) :: y2
  real (8) :: dx
  
    
  dx=x-sml_ignore_drift_r0

  if(dx<0D0) then
     y2= sml_ignore_drift_z0 + sml_ignore_drift_slope1*dx
  else
     y2= sml_ignore_drift_z0 + sml_ignore_drift_slope2*dx
  endif
  
  if(y>y2) then
     factor=1D0
  else
     factor=0D0
  endif
  
end subroutine ignore_factor

#ifdef USE_BICUB_MOD
subroutine psi_der_all(r,z,ret)
use bicub_mod
implicit none
  real (8) :: r,z,ret(6)
integer, parameter :: i00 = 1
integer, parameter :: i10 = 2
integer, parameter :: i01 = 3
integer, parameter :: i20 = 4
integer, parameter :: i11 = 5
integer, parameter :: i02 = 6

call bicub_interpol(psi_bicub,r,z,                                           &
  ret(i00),ret(i10),ret(i01),                                          &
  ret(i11),ret(i20),ret(i02))
end subroutine psi_der_all
#else
subroutine psi_der_all(r,z,ret)
  implicit none
  real (8) :: r,z,ret(6)
  real (8) :: psi_interpol


  ret(1) = psi_interpol(r,z,0,0)
  ret(2) = psi_interpol(r,z,1,0)
  ret(3) = psi_interpol(r,z,0,1)
  ret(4) = psi_interpol(r,z,2,0)
  ret(5) = psi_interpol(r,z,1,1)
  ret(6) = psi_interpol(r,z,0,2)

end subroutine psi_der_all
#endif

!!****************************************************************************
!!> calculate field of a given point using interpolation funtions
!!  adopted from xorbit
!!
!!  first created : 2000/10/19
!!  last modified : 2006/02/01
!!  B->-B routine added (commented out)
!!  time dependance of field is added (2002/6/19)
!!  2006/02/01 fld module update
!!  2002/09/10 code optimization for speed
!!  2002/11/18 code modification for gxc
!!****************************************************************************
subroutine field(fld,t,rz_outside)
    use fld_module
    use sml_module, only :sml_bp_sign, sml_time, sml_cylindrical
    use eq_module, only : eq_min_r,eq_max_r,eq_min_z,eq_max_z, eq_x_psi, is_rgn12, eq_axis_r
    implicit none
    type(fld_type) :: fld  !! Field information
    real (kind=8), intent(in) :: t  !! time
    logical , intent(out) :: rz_outside
    real (kind=8) :: psi, dpsi_dr, dpsi_dz, d2psi_d2r, d2psi_drdz, d2psi_d2z
    real (kind=8) :: r,z, ret(6)
    real (kind=8) , external :: psi_interpol
    real (kind=8) , external :: I_interpol
    real (kind=8) :: r2, over_r,over_r2 !! variables for opimization
    real (8) :: cos_rip, sin_rip, ripp, dripp_dr, dripp_dz, &
         rippbphi, drippbphi_dr, drippbphi_dz, drippbphi_dphi

    r=fld%r
    z=fld%z
    r2=r**2
    if (sml_cylindrical) then
      over_r = 1.D0/eq_axis_r
    else
      over_r=1/r
    endif
    over_r2=over_r**2
    if(r<eq_min_r) then
	r=eq_min_r
	rz_outside=.true.
    else if (r>eq_max_r)then
	r=eq_max_r
	rz_outside=.true.
    endif
    if(z<eq_min_z) then
	z=eq_min_z
	rz_outside=.true.
    else if (z>eq_max_z)then
       z=eq_max_z
       rz_outside=.true.
    else
       rz_outside=.false.
    endif

    call psi_der_all(r,z,ret)
    psi        =ret(1)
    dpsi_dr    =ret(2)
    dpsi_dz    =ret(3)
    d2psi_d2r  =ret(4)
    d2psi_drdz =ret(5)
    d2psi_d2z  =ret(6)

    fld%psi=psi
    fld%dpsidr=dpsi_dr
    fld%dpsidz=dpsi_dz
    ! added 2001/06/01  - lower bound of psi
    if(psi<0D0) then
	psi=0D0
    endif

    !fld_q=q_interpol(psi,0)  --> no need
    
!    if(psi<eq_x_psi .AND. z<eq_x_z) then
    if(.not. is_rgn12(r,z,psi) ) then
	fld%I=I_interpol(psi,0,3)
	fld%dIdpsi = I_interpol(psi,1,3)
    else
        fld%I=I_interpol(psi,0,1)
	fld%dIdpsi = I_interpol(psi,1,1)
    endif

    
    fld%br=- dpsi_dz *over_r * sml_bp_sign
    fld%bz= dpsi_dr *over_r  * sml_bp_sign
    fld%bphi=fld%I *over_r   
	
    !derivativs
    fld%dbrdr=(dpsi_dz *over_r2 - d2psi_drdz *over_r) * sml_bp_sign
    fld%dbrdz=- d2psi_d2z *over_r                   * sml_bp_sign
    fld%dbrdp=0D0                                   * sml_bp_sign
    
    fld%dbzdr= (- dpsi_dr * over_r2 + d2psi_d2r *over_r) * sml_bp_sign
    fld%dbzdz= d2psi_drdz *over_r                      * sml_bp_sign
    fld%dbzdp=0D0                                      * sml_bp_sign
           
    fld%dbpdr= dpsi_dr * fld%dIdpsi *over_r - fld%I *over_r2
    fld%dbpdz= fld%dIdpsi * dpsi_dz *over_r
    fld%dbpdp=0D0
!    call efield(r,z,fld_psi,t)  ! set fld_er, ez, ephi value for a give r,z value  -> move to derivs

end subroutine


subroutine modulo_ie(grid, spall)
  use sml_module
  use ptl_module
  use grid_class
  use diag_module
  implicit none
  include 'mpif.h'
  type(grid_type) :: grid
  type(species_type) :: spall(0:ptl_nsp_max)
  integer i,m
  integer :: isp

  real (8) :: phi
  do isp=ptl_isp, ptl_nsp

  ! Modulo operation
     do i=1, spall(isp)%num
        phi=spall(isp)%ptl(i)%ph(3)
        if( phi >= sml_2pi_wedge_n .or. phi < 0D0 ) then
           spall(isp)%ptl(i)%ph(3)=modulo(phi,sml_2pi_wedge_n)
        endif
     enddo
  enddo
end subroutine modulo_ie


subroutine remove_particle(sp,i,flag,ith)
  use ptl_module
  use sml_module
  use neu_module, only : neu_weight_sum_lost
  implicit none
  type(species_type) :: sp
  integer, intent(in) :: i,flag,ith
  !error diagnosis
  real (kind=8) :: b_interpol, b,ekin,pitch
  
  if(sp%ptl(i)%gid <= 0 ) then
!     print *, 'minus gid particle in remove_particle'
     return
  endif
  sp%ptl(i)%gid=-sp%ptl(i)%gid

  !### debug
  if(sp%ptl(i)%gid > 0 ) then
     print *, 'something wrong in remove_particle'
  endif

!!!$OMP CRITICAL (REMOVE_PARTICLE)
!!  sp%lost_num=sp%lost_num + 1
!!  if(sp%lost_num <= sp%lost_nummax) then
!!     sp%lost_index(sp%lost_num)=i
!!  endif  
!!!$OMP END CRITICAL (REMOVE_PARTICLE)

  if(flag==-2 .and. sml_mype==0) then
     b=b_interpol(sp%ptl(i)%ph(1),sp%ptl(i)%ph(2),0D0)
     call rho_mu_to_ev_pitch2(sp%ptl(i)%ph(4),sp%ptl(i)%ct(pim),b,ekin,pitch,sp%type)
     write(400+sp%type,*) sp%ptl(i)%ph(1:2),&
          ekin, pitch,&          
          sp%phase0(1,i),sp%phase0(2,i) 
  endif

  if(sml_neutral .and. flag==-1.and.(sp%ptl(i)%gid<0).and.(sp%type==1))then
     neu_weight_sum_lost(ith)=neu_weight_sum_lost(ith) + sp%ptl(i)%ct(piw0)
  endif

  if(flag==-1 .and. sml_mype==0) then
!     write(450+sp%type,*) sp%phase(1,i), sp%phase(2,i),&
!          ekin, pitch,&          
!          sp%phase0(1,i),sp%phase0(2,i) 
  endif

!  if(flag==-3) then
!     print *,  'unexpected search fail in search_pid', sp%phase(1,i), sp%phase(2,i)
!  endif
end subroutine remove_particle

subroutine memory_cleaning(spall)
  use ptl_module
  use sml_module, only :sml_mype, sml_electron_on
  implicit none
  type(species_type) :: spall(0:ptl_nsp_max)
  integer :: isp

!$OMP PARALLEL DO &
!$OMP PRIVATE( ISP )
  do isp=ptl_isp, ptl_nsp
     call memory_cleaning_one_sp(spall(isp))
  enddo

contains
  subroutine memory_cleaning_one_sp(sp)
    implicit none
    type(species_type) :: sp
    integer :: i,j,ptmp


    
    if(sp%lost_num<=sp%lost_nummax) then
       ! algorithm 
       ! BASIC : for every empty space, move last particle to this 
       ! EXCEPTION : last particle can be lost particle
       !           : If so, eliminate the particle and find next non-empty space
       !           : When the index indicate erased space, just ignore it
       do j=1, sp%lost_num
          !find non-empty last particle
          if(sp%num>0) then
             do while( sp%ptl(sp%num)%gid<=0 ) 
                ptmp=sp%num
                call eliminate_one(ptmp,sp%ptl(ptmp),isp)
                sp%num=sp%num-1
                if(sp%num < 1) exit
             enddo
          endif
          i=sp%lost_index(j)
          
          ! ignore erased space 
          ! i cannot be sp%num because sp%num has positive gid value
          
          if( i < sp%num ) then
             call eliminate_one(i,sp%ptl(i),isp)
             !move last particle to i-th memory
             call copy_ptl(sp,sp%num,i)
             sp%num=sp%num-1
          endif
       enddo
    else
       !sequencial search whole particle
       i=1
       do while ( i < sp%num)
          if(sp%ptl(i)%gid <= 0) then
             !calculate weight summation of lost particle
             call eliminate_one(i,sp%ptl(i),isp)
             !eliminate i-th particle and move last particle to i-th memory
             call copy_ptl(sp,sp%num,i)
             sp%num=sp%num-1
          else
             ! next particle
             i=i+1        
          endif
       enddo
       
       !last particle exception
       if(sp%ptl(sp%num)%gid <= 0 ) then
          call eliminate_one(sp%num,sp%ptl(sp%num),isp)
          sp%num=sp%num-1
       endif
     endif
     
     ! debug - check gid again
!!$     do i=1, pnum
!!$        if(gid(i)<0) then
!!$           print *, 'error in memory cleaning',i,pnum, sml_mype
!!$           stop
!!$        endif
!!$     enddo
     sp%lost_num=0
   end subroutine memory_cleaning_one_sp
end subroutine memory_cleaning

subroutine copy_ptl(sp, src, dest)
  use ptl_module
  implicit none
  type(species_type) :: sp
  integer, intent(in) :: src, dest

  sp%ptl(dest)=sp%ptl(src)
  sp%phase0(:,dest)=sp%phase0(:,src)
  sp%tr_save(dest)=sp%tr_save(src)
  sp%p_save(:,dest)=sp%p_save(:,src)
  if(sp%type/=0) then
     sp%rhoi(dest)=sp%rhoi(src)
  endif

end subroutine copy_ptl

subroutine memory_cleaning_simple(spall)
  use ptl_module
  use sml_module, only :sml_mype, sml_electron_on
  implicit none
  type(species_type) :: spall(0:ptl_nsp_max)
  integer :: isp

  do isp=ptl_isp, ptl_nsp
     call one_sp(spall(isp))
  enddo

contains
  subroutine one_sp(sp)
    implicit none
    type(species_type) :: sp
    integer :: i,j,ptmp


    !sequencial search whole particle
    i=1
    do while ( i < sp%num)
       if(sp%ptl(i)%gid <= 0) then
          !calculate weight summation of lost particle
          call eliminate_one(i,sp%ptl(i),isp)
          !eliminate i-th particle and move last particle to i-th memory
          call copy_ptl(sp,sp%num,i)
          sp%num=sp%num-1
       else
          ! next particle
          i=i+1        
       endif
    enddo

    !last particle exception
    if(sp%num>0) then
       if(sp%ptl(sp%num)%gid <= 0) then
          call eliminate_one(sp%num,sp%ptl(sp%num),isp)
          sp%num=sp%num-1
       endif
    endif
  end subroutine one_sp

end subroutine memory_cleaning_simple


! Particle elimination - side effect calculation
subroutine eliminate_one(i,ptli,sp_type)
  use ptl_module
  implicit none
  integer,intent(in) :: i,sp_type
  type(ptl_type), intent(in) :: ptli


end subroutine eliminate_one

subroutine shift_pc_index(index,n)
  implicit none
  integer, intent(in) :: n
  integer :: index(n)
  integer :: tmp(n)

  tmp(2:n) = index(1:n-1)
  tmp(1)=index(n)
  
  index=tmp

end subroutine shift_pc_index


! save electron phase variables when ipc==1 and ihybrid==1
! It restores electron phase variabes otherwise.
! Saved variables are :
!     sp%num
!     sp%ptl(i)
! following variables are not required :
!     sp%phase0, sp%derivs, sp%dph, sp%dpht, tr_save, p_save
subroutine save_or_load_electron_phase(sp,ipc,ihybrid)
  use ptl_module 
  implicit none
  type(species_type) :: sp
  integer, intent(in) :: ipc, ihybrid
  integer :: i

  ! save electron info
  if(ipc==1 .and. ihybrid==1) then
!     if(diag_tracer_sp==0) ptl_tracer_n_save=diag_tracer_n     
     ptl_enum_save=sp%num
!$OMP PARALLEL DO &
!$OMP PRIVATE( I )
     do i=1, sp%num
        ptl_ephase_save(i)=sp%ptl(i)
     enddo
     

  else ! restore saved electron info
!     if(diag_tracer_sp==0) diag_tracer_n=ptl_tracer_n_save
     sp%num=ptl_enum_save
!$OMP PARALLEL DO &
!$OMP PRIVATE( I )
     do i=1, sp%num
        sp%ptl(i)=ptl_ephase_save(i)
     enddo
  endif
end subroutine save_or_load_electron_phase


! update time derivative of potential
! ddpotdt is  d(delta phi) / dt 
! obtained from saved potential (previous time step) and predicted potential
! ihybrid is the iteration number
! saved potential has nrk*nhybrid data sets.
subroutine ddpotdt_update(grid,psn,ipc,ihybrid)
  use grid_class
  use psn_class
  use sml_module
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  integer, intent(in) :: ipc, ihybrid
  integer :: i,iphi, j
  real (8) :: inv_dt
  !
  real (8) :: tmp(grid%nnode,0:1)
  
  inv_dt=1D0/sml_dt

  j=ipc + sml_nrk*(ihybrid-1)

  do iphi=0, 1
!$OMP PARALLEL DO &
!$OMP PRIVATE( I )
     do i=1, grid%nnode
        psn%ddpotdt(i,iphi)=(psn%dpot(i,iphi)-psn%dpotsave(i,iphi,j))*inv_dt
        psn%dpotsave(i,iphi,j)=psn%dpot(i,iphi)
     enddo
  enddo

  ! psn%dpot is not field following coord
  ! ddpotdt is converted to field following coord.
  ! for speed up this algorithm need to be optimized.
  !### simple upgrade : tmp = dpot -dpotsave , cnv_grid(..,.., tmp, ddpotdt)
  call cnvt_grid_real2ff(grid,psn%ff_hdp_tr,psn%ff_hdp_p,psn%ddpotdt,tmp)
  psn%ddpotdt=tmp

end subroutine ddpotdt_update




!!$subroutine single_phase_derivs_with_E(ptli,yprime,sp,i,time,E_mag)
!!$  use ptl_module
!!$  use fld_module
!!$  use sml_module
!!$  implicit none
!!$  type(ptl_type) :: ptli
!!$  type(species_type) :: sp
!!$  real (kind=8) ::  yprime(ptl_nphase),yp_exb(ptl_nphase),time
!!$  integer :: i
!!$  type(fld_type) :: fld
!!$  real (8) :: dpsi(2), E(3), dtheta_norm(2),B, E_mag(3), bp
!!$  logical :: rz_outside
!!$
!!$  fld%r=ptli%ph(1)
!!$  fld%z=ptli%ph(2)
!!$  fld%phi=ptli%ph(3)
!!$  call field(fld,time,rz_outside)
!!$  if(.not. rz_outside) then
!!$
!!$     dpsi(1:2)=(/ fld%dpsidr, fld%dpsidz /)
!!$     E(1:2)=E_mag(1)*dpsi
!!$
!!$     bp=sqrt(fld%br**2+fld%bz**2)
!!$     dtheta_norm(1:2)=(/ fld%br, fld%bz /)/bp
!!$     E(1:2)=E(1:2) + E_mag(2)*dtheta_norm
!!$
!!$     B=sqrt(fld%br**2+fld%bz**2+fld%bphi**2)
!!$     E(3)=(E_mag(3)*B- E(1)*fld%Br - E(2)*fld%Bz)/fld%Bphi
!!$
!!$     fld%Er=E(1)
!!$     fld%Ez=E(2)
!!$     fld%Ephi=E(3)
!!$     if(sml_deltaf) then 
!!$        call f0_info(fld) ! some mode is not working -triangle related ftn - cannot use this
!!$     endif
!!$     call derivs_sp(fld,time,ptli,yprime,yp_exb,sp%type)
!!$  else
!!$     call remove_particle(sp,i,-1)
!!$     return
!!$  endif
!!$  
!!$end subroutine single_phase_derivs_with_E

