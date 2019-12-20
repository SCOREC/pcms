#ifdef USE_CALC_GRADIENT

#include "calc_Er_Ez.F90"
#include "calc_E_para.F90"
#include "calc_E_phi_ff.F90"

#endif 

subroutine pushe(istep,ihybrid,ncycle,grid,psn,sp,diag_on)
  use grid_class
  use psn_class
  use ptl_module
  use sml_module
  use omp_module , only : split_indices
  use perf_monitor
  implicit none
  integer, intent(in) :: istep, ihybrid
  integer, intent(in) :: ncycle
  type(grid_type) :: grid
  type(psn_type) :: psn
  type(species_type) :: sp
  logical :: diag_on
  !
  integer :: icycle
  integer :: epc, i
  integer :: ith, i_beg(sml_nthreads), i_end(sml_nthreads)
  
  logical, parameter :: use_sort_particles = .true.

  if (use_sort_particles) then
     call t_startf("pushe_sort_particles")
     call sort_particles(grid,sp)
     call t_stopf("pushe_sort_particles")
  endif

  do icycle=1, ncycle
     !                        
     do epc=1,sml_nrk
        call t_startf("ELECTRON_LOOP")

        sml_epc=epc

        ! search electron position        
        call t_startf("PUSHE_SEARCH_INDEX")
        call chargee_search_index(grid,psn,sp)
        call t_stopf("PUSHE_SEARCH_INDEX")
               
        !call determine_diag_on(istep,ipc,diag_on) ! redundant - but for safety
        if(ihybrid>1 .or. icycle >1 .or. epc > 1) diag_on=.false.

        call t_startf("PUSHE_1step")
#ifdef RK4_ACCURATE_EFIELD
        call pushe_1step2(istep,epc,grid,psn,sp,sp%phase0,sp%ptl,diag_on)
#else
        call pushe_1step(istep,epc,grid,psn,sp,sp%phase0,sp%ptl,diag_on)
#endif
        call t_stopf("PUSHE_1step")

        ! Modulo operation
        call split_indices(sp%num, sml_nthreads, i_beg, i_end)
  !$OMP PARALLEL DO &
  !$OMP PRIVATE( ITH, I )
        do ith=1, sml_nthreads
           call t_startf("PUSHE_MOD_LOOP") 
           do i=i_beg(ith), i_end(ith)
              if(sp%ptl(i)%ph(3)>= sml_2pi_wedge_n .or. sp%ptl(i)%ph(3)< 0D0 ) then
                 sp%ptl(i)%ph(3)=modulo(sp%ptl(i)%ph(3),sml_2pi_wedge_n)
              endif
           enddo
           call t_stopf("PUSHE_MOD_LOOP") 
        enddo
        call t_stopf("ELECTRON_LOOP")        
     enddo
  enddo

  call t_startf("PUSHE_SEARCH_INDEX")
  call chargee_search_index(grid,psn,sp)
  call t_stopf("PUSHE_SEARCH_INDEX")

  sml_epc=1  ! for diagnostic purpose -- sheath_calculatoin
end subroutine pushe

subroutine gather_field_info(grid,psn) 
  use grid_class
  use psn_class
  use sml_module
  use perf_monitor
  implicit none
  type(grid_type):: grid
  type(psn_type) :: psn
  include 'mpif.h'
  
  integer :: iphi, ipe, sz, ierr
#ifdef USE_CALC_GRADIENT
  integer, parameter :: idebug = 0
  integer, parameter :: iodev = 31
  character(len=80) :: filename
  integer :: i
  real(8), dimension(0:(sml_nphi_total-1)) :: maxerr_E_r
  real(8), dimension(0:(sml_nphi_total-1)) :: maxerr_E_z
  real(8), dimension(0:(sml_nphi_total-1)) :: maxerr_E_para
  real(8) :: global_maxerr_E_r, err_E_r, maxval_maxerr_E_r 
  real(8) :: global_maxerr_E_z, err_E_z, maxval_maxerr_E_z
  real(8) :: global_maxerr_E_para, err_E_para, maxval_maxerr_E_para 
  real(8) :: E_phi_ff(1:3,0:1), psn_E_phi_ff(1:3,0:1)
  real(8) :: global_list(0:(sml_plane_totalpe+sml_intpl_totalpe))
  integer :: isize, comm
  logical :: is_ionode

  integer :: pid, iplane, iplane_max, iphi_max
  real(8) :: dval
  real(8) :: global_err(0:(sml_plane_totalpe-1),0:(sml_intpl_totalpe-1))
  integer :: root, tag, irequest
  integer :: array_of_request( size(global_err) )
  integer :: array_of_status(mpi_status_size, size(array_of_request))
#endif

#ifdef USE_OLD_GATHER_FIELD

  do iphi=0, sml_nphi_total - 1

     call t_startf("GATHER_FIELD_INFO_COPY")
     ipe = sml_pe_per_plane * iphi
     if(sml_mype==ipe) then
#ifdef USE_CALC_GRADIENT
#else
        psn%E_phi_ff(:,0,:,iphi)=psn%E_rho_ff(:,0,0,:)
        psn%E_phi_ff(:,1,:,iphi)=psn%E_rho_ff(:,1,0,:)
#endif

        if(sml_extra_dwdt) then
           psn%pot_phi_ff(0,:,iphi)=psn%pot_rho_ff(0,0,:)
           psn%pot_phi_ff(1,:,iphi)=psn%pot_rho_ff(1,0,:)

           psn%ddpotdt_phi(:,0,iphi)=psn%ddpotdt(:,0)
           psn%ddpotdt_phi(:,1,iphi)=psn%ddpotdt(:,1)
        endif
     endif
     call t_stopf("GATHER_FIELD_INFO_COPY")

     call t_startf("GATHER_FIELD_INFO_BCAST")
     sz=grid%nnode * 2 * 3
#ifdef USE_CACL_GRADIENT
#else
     call mpi_bcast(psn%E_phi_ff(:,:,:,iphi),sz,mpi_real8,ipe,sml_comm,ierr)
#endif

     if(sml_extra_dwdt) then
        sz=grid%nnode * 2
        call mpi_bcast(psn%pot_phi_ff(:,:,iphi),sz,mpi_real8,ipe,sml_comm,ierr)

        sz=grid%nnode * 2
        call mpi_bcast(psn%ddpotdt_phi(:,:,iphi),sz,mpi_real8,ipe,sml_comm,ierr)
     endif
     call t_stopf("GATHER_FIELD_INFO_BCAST")
  enddo

#else

  call t_startf("GATHER_FIELD_INFO_COPY")
  do iphi=0, sml_nphi_total - 1
     if(sml_intpl_mype==iphi) then


#ifdef USE_CALC_GRADIENT
#else
        psn%E_phi_ff(:,0,:,iphi)=psn%E_rho_ff(:,0,0,:)
        psn%E_phi_ff(:,1,:,iphi)=psn%E_rho_ff(:,1,0,:)
#endif

        if(sml_extra_dwdt) then
           psn%pot_phi_ff(0,:,iphi)=psn%pot_rho_ff(0,0,:)
           psn%pot_phi_ff(1,:,iphi)=psn%pot_rho_ff(1,0,:)

           psn%ddpotdt_phi(:,0,iphi)=psn%ddpotdt(:,0)
           psn%ddpotdt_phi(:,1,iphi)=psn%ddpotdt(:,1)


           !###### E00_ff required for sml_extra_dwdt
           !##### maybe 
           

        endif
     endif
  enddo
  call t_stopf("GATHER_FIELD_INFO_COPY")

  call t_startf("GATHER_FIELD_INFO_GATHER")
 
  sz=grid%nnode * 2 * 3
#ifdef USE_CACL_GRADIENT
#else
  call mpi_allgather(MPI_IN_PLACE,sz,MPI_REAL8,psn%E_phi_ff(:,:,:,:),sz,MPI_REAL8,sml_intpl_comm,ierr)
#endif

  if(sml_extra_dwdt) then
     sz=grid%nnode * 2
     call mpi_allgather(MPI_IN_PLACE,sz,MPI_REAL8,psn%pot_phi_ff(:,:,:),sz,MPI_REAL8,sml_intpl_comm,ierr)     
     sz=grid%nnode * 2 
     call mpi_allgather(MPI_IN_PLACE,sz,MPI_REAL8,psn%ddpotdt_phi(:,:,:),sz,MPI_REAL8,sml_intpl_comm,ierr)
  endif

  call t_stopf("GATHER_FIELD_INFO_GATHER")

#endif

#ifdef USE_CALC_GRADIENT
!  ----------------------------------------
!  double check values of gradient in
!  E_rho(1:3,0:1,1:grid%nnode,0)
!  compared to computed using 
!  psn%dot_phi_real(1:grid%nnode,0:nphim1)
!  ----------------------------------------

 if (idebug >= 1) then

  is_ionode = (sml_mype .eq. 0)

  maxval_maxerr_E_r  = 0
  maxval_maxerr_E_z  = 0
  maxval_maxerr_E_para  = 0

  global_maxerr_E_r  = 0
  global_maxerr_E_z  = 0
  global_maxerr_E_para  = 0

  if (idebug >= 2) then
    filename =  ' '
    write(filename,8010) sml_mype
 8010 format('pE_para.',i3.3,'.txt')
    open(unit=iodev,file=trim(filename),access='sequential', &
         form='formatted',status='replace')
    rewind(iodev)
    write(iodev,8012) sml_plane_mype, sml_intpl_mype
 8012 format('% sml_plane_mype ',i9,' sml_intpl_mype ',i9 )
    write(iodev,*) '% i,calc_E_para(0:1), psn_E_para(0:1) '
  endif


  do iphi=0, sml_nphi_total-1
  if (sml_intpl_mype ==  iphi) then

  maxerr_E_r(iphi) = 0
  maxerr_E_z(iphi) = 0
  maxerr_E_para(iphi) = 0

  do i=1,grid%nnode
    call calc_E_phi_ff(grid,psn,    i, iphi,E_phi_ff )
!    -------------------------------------------------------------------
!    err_E_r = maxval(abs(E_phi_ff(1,0:1) - psn%E_phi_ff(1,0:1,i,iphi)))
!    err_E_z = maxval(abs(E_phi_ff(2,0:1) - psn%E_phi_ff(2,0:1,i,iphi)))
!    err_E_para = maxval(abs(E_phi_ff(3,0:1) - psn%E_phi_ff(3,0:1,i,iphi)))
!    -------------------------------------------------------------------
    psn_E_phi_ff(1:3,0:1) = psn%E_rho_ff(1:3,0:1,0,i)

    err_E_r = max(abs(E_phi_ff(1,0) - psn_E_phi_ff(1,0)),     &
                  abs(E_phi_ff(1,1) - psn_E_phi_ff(1,1)) )
    err_E_z = max(abs(E_phi_ff(2,0) - psn_E_phi_ff(2,0)),     &
                  abs(E_phi_ff(2,1) - psn_E_phi_ff(2,1)) )
    err_E_para = max(abs(E_phi_ff(3,0) - psn_E_phi_ff(3,0)),     &
                     abs(E_phi_ff(3,1) - psn_E_phi_ff(3,1)) )

    maxerr_E_r(iphi) = max(maxerr_E_r(iphi),err_E_r)
    maxerr_E_z(iphi) = max(maxerr_E_z(iphi),err_E_z)
    maxerr_E_para(iphi) = max(maxerr_E_para(iphi),err_E_para)

    if (idebug >= 2) then
      write(iodev,8020)  i, E_phi_ff(3,0:1), psn_E_phi_ff(3,0:1)
 8020 format(1x,i9,4(1x,1pe16.6))
    endif

  enddo

  if (idebug >= 2) then
    close(iodev)
  endif


  maxval_maxerr_E_r = max(maxval_maxerr_E_r, maxerr_E_r(iphi))
  maxval_maxerr_E_z = max(maxval_maxerr_E_z, maxerr_E_z(iphi))
  maxval_maxerr_E_para = max(maxval_maxerr_E_para, maxerr_E_para(iphi))

  endif
  enddo

! -------------------------------------
! compute maximum across all processors
! -------------------------------------

       

! ---------------------------
! global error in same plane 
! ---------------------------

   
   isize = 1
   comm = sml_plane_comm


   call MPI_Allreduce( maxval_maxerr_E_r, global_maxerr_E_r, isize, &
       MPI_REAL8, MPI_MAX, comm, ierr )

   call MPI_Allreduce( maxval_maxerr_E_z, global_maxerr_E_z, isize, &
       MPI_REAL8, MPI_MAX, comm, ierr )

   call MPI_Allreduce( maxval_maxerr_E_para, global_maxerr_E_para, isize, &
       MPI_REAL8, MPI_MAX, comm, ierr )

   if (is_ionode) then
     write(*,*) 'maxerr_E_r_plane ', global_maxerr_E_r
     write(*,*) 'maxerr_E_z_plane ', global_maxerr_E_z
     write(*,*) 'maxerr_E_para_plane ', global_maxerr_E_para
   endif



!  -------------------
!  global inter plane
!  -------------------
   
   isize = 1
   comm = sml_intpl_comm

   call MPI_Allreduce( maxval_maxerr_E_r, global_maxerr_E_r, isize, &
       MPI_REAL8, MPI_MAX, comm, ierr )

   call MPI_Allreduce( maxval_maxerr_E_z, global_maxerr_E_z, isize, &
       MPI_REAL8, MPI_MAX, comm, ierr )

   call MPI_Allreduce( maxval_maxerr_E_para, global_maxerr_E_para, isize, &
       MPI_REAL8, MPI_MAX, comm, ierr )

   if (is_ionode) then
     write(*,*) 'maxerr_E_r_intpl ', global_maxerr_E_r
     write(*,*) 'maxerr_E_z_intpl ', global_maxerr_E_z
     write(*,*) 'maxerr_E_para_intpl ', global_maxerr_E_para
   endif




!  -------------------
!  global across all
!  -------------------
   
   isize = 1

   comm = sml_comm

   call MPI_Allreduce( maxval_maxerr_E_r, global_maxerr_E_r, isize, &
       MPI_REAL8, MPI_MAX, comm, ierr )

   call MPI_Allreduce( maxval_maxerr_E_z, global_maxerr_E_z, isize, &
       MPI_REAL8, MPI_MAX, comm, ierr )

   call MPI_Allreduce( maxval_maxerr_E_para, global_maxerr_E_para, isize, &
       MPI_REAL8, MPI_MAX, comm, ierr )

   if (is_ionode) then
     write(*,*) 'maxerr_E_r_all ', global_maxerr_E_r
     write(*,*) 'maxerr_E_z_all ', global_maxerr_E_z
     write(*,*) 'maxerr_E_para_all ', global_maxerr_E_para
   endif


  endif


#endif
end subroutine gather_field_info



!!particle electron pushing routine using R-K 2nd order + RK4 hybrid
subroutine pushe_1step(istep,ipc,grid,psn,sp,phase0,ptl,diag_on)
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
  real (kind=8) :: phase0(ptl_nphase,sp%maxnum)   ! sp%phase0 -- for speed
  type(ptl_type) ::  ptl(sp%maxnum)    ! sp%phase
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

!   
!  if(sp%type/=0) then
!     ! for ion
!     dt=sml_dt              
!  else
!     ! for electron
!     dt=sml_dt*0.5D0/real(sml_ncycle_half)
!  endif

!** electron only
  dt=sml_dt*0.5D0/real(sml_ncycle_half)
  
  
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
    
  
  !$OMP PARALLEL DO &
  !$OMP PRIVATE( ITH, I, NEW_PHASE, &
  !$OMP          OLD_PHASE, RTN )
  do ith=1, sml_nthreads
     call t_startf("PUSHE_LOOP") 
     do i=i_beg(ith), i_end(ith)
        ! for alive particles only
        if(sp%ptl(i)%gid>0) then
           
           !******************************************************
           ! actual particle push           
           !******************************************************
           call pushe_single(grid,psn,sp,i,phase0(:,i),new_phase,dt_now,ith,diag_on)
           
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
     call t_stopf("PUSHE_LOOP") 
  enddo
  
end subroutine pushe_1step

!!particle electron pushing routine using pure RK4
subroutine pushe_1step2(istep,ipc,grid,psn,sp,phase0,ptl,diag_on)
  use sml_module
  use ptl_module
  use fld_module
  use grid_class
  use psn_class
  use omp_module , only : split_indices
  use perf_monitor
  use eq_module
  implicit none
  integer, intent(in) :: istep, ipc !
  type(grid_type) :: grid
  type(psn_type) :: psn
  type(species_type) :: sp
  real (kind=8) :: phase0(ptl_nphase,sp%maxnum)   ! sp%phase0 -- for speed
  type(ptl_type) ::  ptl(sp%maxnum)    ! sp%phase
  logical, intent(in) :: diag_on
  !
  type(fld_type) :: fld
  real (kind=8) :: dt_now,time_now,new_phase(ptl_nphase),old_phase(ptl_nphase)
  integer :: i,rtn,j
  integer :: ith, i_beg(sml_nthreads), i_end(sml_nthreads)
  real (kind=8) , external :: psi_interpol
  character (len=5) :: err_str(2)
  logical, parameter :: USE_SEARCH_TR2 = .true.
  err_str(1)='ion'
  err_str(2)='elec'
  

  if(sp%num==0) return  ! nothing to push

  call split_indices(sp%num, sml_nthreads, i_beg, i_end)

  dt_now=sml_dt*0.5D0/real(sml_ncycle_half)
  
  


!$OMP PARALLEL DO &
!$OMP PRIVATE( ITH, I, NEW_PHASE, &
!$OMP          OLD_PHASE, RTN )
  do ith=1, sml_nthreads
     call t_startf("PUSHE_LOOP") 
     do i=i_beg(ith), i_end(ith)

        !save phase0 information
        phase0(:,i)=ptl(i)%ph


        ! for alive particles only
        if(sp%ptl(i)%gid>0) then
           
           !******************************************************
           ! actual particle push           
           !******************************************************
           call pushe_single(grid,psn,sp,i,phase0(:,i),new_phase,dt_now,ith,diag_on)
           
           ! check r-z boundary validity and update psi variables
           if(new_phase(1)<eq_min_r .or. new_phase(1)>eq_max_r .or. new_phase(2)<eq_min_z .or. new_phase(2)>eq_max_z)then
              call remove_particle(sp,i,-1,ith)
!              print *, 'particle eliminated due to rz_outside :', i, sml_mype, sp%type, sp%ptl(i)%gid, new_phase(1),new_phase(2)
           else                            
              ! bounce 
              !if(ipc==sml_nrk .and. sml_bounce/=0) then
              if(sml_bounce/=0) then
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
     call t_stopf("PUSHE_LOOP") 
  enddo
  
end subroutine pushe_1step2



! single particle push -- get new_phase using rk4 with initial E-field
subroutine pushe_single(grid,psn,sp,i,y,new_phase,dt,ith,diag_on)
  use sml_module
  use ptl_module
  use fld_module
  use grid_class
  use psn_class
  use perf_monitor
  implicit none
  type(grid_type), intent(in) :: grid
  type(psn_type), intent(in) :: psn
  type(species_type), intent(in) :: sp
  real (8), intent(in) :: y(ptl_nphase)
  real (8), intent(out) :: new_phase(ptl_nphase)
  real (8), intent(in) :: dt
  integer, intent(in) :: ith
  logical, intent(in) :: diag_on
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
     call derivs_single_elec(grid,psn,sp,ptli,dy,i,time,fld,ith,diag_on)


     
#ifdef PURE_RK2
     ! Only simple RK2
     
     new_phase = y + dt * dy     
     !call restrict_weight(new_phase(piw1:piw2))
     
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
     call derivs_single_with_e_elec(grid,psn,sp,ptli,dy ,i,time    ,fld,E_mag,ith)
     
     ptli%ph = y + hdt * dy
     call derivs_single_with_e_elec(grid,psn,sp,ptli,dyt,i,th      ,fld,E_mag,ith)
     
     ptli%ph = y + hdt * dyt
     call derivs_single_with_e_elec(grid,psn,sp,ptli,dym,i,th      ,fld,E_mag,ith)
     
     ptli%ph = y + dt * dym
     dym = dyt + dym
     call derivs_single_with_e_elec(grid,psn,sp,ptli,dyt,i,time+dt ,fld,E_mag,ith)
     
     ! Obtain new_phase
     new_phase = y + dt/6D0 * ( dy + dyt + 2D0*dym )
     !call restrict_weight(new_phase(piw1:piw2))
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
end subroutine pushe_single

!obtain derivatives of phase variable : actual calculation is done in derivs_sp. 
! prepare E-field and B-field
subroutine derivs_single_elec(grid,psn,sp,ptli,dy,i,time,fld,ith,diag_on)
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
     call efield_elec(grid,psn,sp,i,fld,time)
!     if(sml_deltaf) call f0_info(fld)
     
     ! get time derivatives
     call derivs_elec(fld,time,ptli,sp%type,dy,diag_on,vf_diag)
     
  else
     call remove_particle(sp,i,-1,ith)
!     print *, 'particle eliminated due to rz_outside', i, sml_mype, sp%type, sp%ptl(i)%gid
  endif
     
!call  check_point('before port1')
  if(diag_on) call diag_1d_port1(sp%ptl(i),dy,sp%type,vf_diag,ith)
!call  check_point('after port1')
end subroutine derivs_single_elec

subroutine derivs_single_with_e_elec(grid,psn,sp,ptli,dy,i,time,fld,E_mag,ith)
  use grid_class
  use psn_class
  use sml_module
  use fld_module
  use ptl_module
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
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
  real (8) :: p(3), x(2), phi, phi_mid, xff(2)
  integer :: itr

  ! Save space information
  fld%r=ptli%ph(1)
  fld%z=ptli%ph(2)
  fld%phi=ptli%ph(3)

  ! obtain B-field information 
  call field(fld,time,rz_outside)
  if(.not. rz_outside) then
     ! obtain gyro-averaged E-field from previous E-field

#ifdef RK4_ACCURATE_EFIELD
     !rh test performance implications of getting the actual field ad every time step
     x(1)=fld%r
     x(2)=fld%z

     if(fld%phi>= sml_2pi_wedge_n .or. fld%phi< 0D0 ) then
        fld%phi=modulo(fld%phi,sml_2pi_wedge_n)
     endif
     phi =fld%phi
     phi_mid=(floor(phi/grid%delta_phi) + 0.5D0) * grid%delta_phi

     call field_following_pos2(x,phi,phi_mid,xff)


     call search_tr2(grid,xff,itr,p)

     if (itr .gt. 0) then
       call efield_gk_elec2(grid,psn,sp,i,fld,itr,p)
     else
#endif
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
#ifdef RK4_ACCURATE_EFIELD
     endif
#endif

!     if(sml_deltaf) call f0_info(fld)
     
     ! get time derivatives
     call derivs_elec(fld,time,ptli,sp%type,dy,diag_on,vf_diag)
     
  else
     call remove_particle(sp,i,-1,ith)
!     print *, 'particle eliminated due to rz_outside', i, sml_mype, sp%type, sp%ptl(i)%gid
  end if
end subroutine derivs_single_with_e_elec




!!****************************************************************************
!!> derivatives of electron phase
!!
!!<***************************************************************************
subroutine derivs_elec(fld,t, ptli, sp_type, yprime, diag_on, vf_diag)
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


  if(sml_deltaf_f0_mode==-1) then
     ! these are constant for deltaf method -- maybe should be moved to sml_module
     sml_wdot_energy_max=10D0
     sml_f0_psi_c=0.5*(sqrt(sml_inpsi)+sqrt(sml_outpsi))
     sml_f0_1_psi_w=1D0/( 0.4*(sqrt(sml_outpsi)-sqrt(sml_inpsi)) )
  endif
     
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
  D=1.D0/ ( 1.D0 + rho * nb_curl_nb )
      
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
        envelop= exp( - ((sqrt(psi) - sml_f0_psi_c)*sml_f0_1_psi_w )**8 )
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

end subroutine derivs_elec

subroutine efield_elec(grid,psn,sp,i,fld,time)
  use sml_module
  use grid_class
  use psn_class
  use fld_module
  use eq_module
  use ptl_module
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  type(fld_type) :: fld
  integer, intent(in) :: i
  type(species_type) :: sp
  real (kind=8) :: time

  
  if(sml_00_efield .or. sml_turb_efield) then
     !Gyrokinetic E
     call efield_gk_elec(grid,psn,sp,i,fld)
  else
       fld%Er=0D0
       fld%Ez=0D0
       fld%Ephi=0D0
       fld%Epot=0D0
       fld%Er00=0D0
       fld%Ez00=0D0
       fld%ddpotdt=0D0
  endif

end subroutine efield_elec


subroutine efield_gk_elec(grid,psn,sp,i,fld) !require fld%(r,z,phi,br,bz,bphi)
  use sml_module
  use grid_class
  use psn_class
  use fld_module
  use ptl_module
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  type(species_type) :: sp
  integer, intent(in) :: i
  type(fld_type), intent(inout) :: fld

  integer :: itr,ip,node,irho, iphi
  real (kind=8) :: wphi(0:1), wp, rho, rhon, wrho(2)
  real (kind=8) :: pot, E(3),E00(2), B , ddpotdt !, D(3)
  real (kind=8) :: E_phi_ff(3,0:1)

  iphi=floor(fld%phi/grid%delta_phi)
  !wphi(1)=(fld%phi/grid%delta_phi)  - grid%iphi_offset
  wphi(1)=(fld%phi/grid%delta_phi) - iphi
  wphi(0)=1D0 - wphi(1)



  pot=0D0
  E=0D0
  E00=0D0
  !D=0D0 !debug
  ddpotdt=0D0

  !get E
  itr=sp%tr_save(i)
  if(itr>0) then        
     do ip=1, 3
        node=grid%nd(ip,itr)
        wp=sp%p_save(ip,i)
        
        ! for electron -- rho=0 case, but iphi=/=0.           

#ifdef USE_CALC_GRADIENT
        call calc_E_phi_ff(grid,psn,node,iphi,E_phi_ff )
        E = E + wp*wphi(0)*E_phi_ff(:,0)
        E = E + wp*wphi(1)*E_phi_ff(:,1)
#else

        E   = E   + wp*wphi(0)*psn%E_phi_ff(:,0,node,iphi)
        E   = E   + wp*wphi(1)*psn%E_phi_ff(:,1,node,iphi)
#endif
        
        if(sml_extra_dwdt) then           
           E00 = E00 + wp*wphi(0)*psn%E00_ff(:,0,node)
           E00 = E00 + wp*wphi(1)*psn%E00_ff(:,1,node)
           ddpotdt = ddpotdt + wp*wphi(0)*psn%ddpotdt_phi(node,0,iphi)
           ddpotdt = ddpotdt + wp*wphi(1)*psn%ddpotdt_phi(node,1,iphi)
           pot = pot + wp*wphi(0)*psn%pot_phi_ff(0,node,iphi) 
           pot = pot + wp*wphi(1)*psn%pot_phi_ff(1,node,iphi)
        endif

     enddo
  else
     print *, 'E-field ion invalid tr : (i,itr,mype,gid)=', i,itr,sml_mype,sp%ptl(i)%gid
  endif
  
  !E(3) was para E-field and becomes Ephi
  if(sml_turb_efield) then
     B=sqrt(fld%br**2 + fld%bz**2 + fld%bphi**2)
     E(3)=(E(3)*B- E(1)*fld%Br - E(2)*fld%Bz)/fld%Bphi
     
  else
     E(3)=0D0
  endif
  
  fld%Er=E(1)
  fld%Ez=E(2)
  fld%Ephi=E(3)
  if(sml_extra_dwdt) then
     fld%Epot=pot
     fld%Er00=E00(1)
     fld%Ez00=E00(2) 
     fld%ddpotdt=ddpotdt
  endif
end subroutine efield_gk_elec

#ifdef RK4_ACCURATE_EFIELD
subroutine efield_gk_elec2(grid,psn,sp,i,fld,itr,p) !require fld%(r,z,phi,br,bz,bphi)
  use sml_module
  use grid_class
  use psn_class
  use fld_module
  use ptl_module
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  type(species_type) :: sp
  integer, intent(in) :: i, itr
  real (kind=8), intent(in) :: p(3)
  type(fld_type), intent(inout) :: fld

  integer :: ip,node,irho, iphi
  real (kind=8) :: wphi(0:1), wp, rho, rhon, wrho(2)
  real (kind=8) :: pot, E(3),E00(2), B , ddpotdt !, D(3)
  real (kind=8) :: E_phi_ff(3,0:1)

  iphi=floor(fld%phi/grid%delta_phi)
  !wphi(1)=(fld%phi/grid%delta_phi)  - grid%iphi_offset
  wphi(1)=(fld%phi/grid%delta_phi) - iphi
  wphi(0)=1D0 - wphi(1)



  pot=0D0
  E=0D0
  E00=0D0
  !D=0D0 !debug
  ddpotdt=0D0

  !get E
!  itr=sp%tr_save(i)
  if(itr>0) then        
     do ip=1, 3
        node=grid%nd(ip,itr)
        !wp=sp%p_save(ip,i)
         wp=p(ip)

        ! for electron -- rho=0 case, but iphi=/=0.           

#ifdef USE_CALC_GRADIENT
        call calc_E_phi_ff(grid,psn,node,iphi,E_phi_ff )
        E = E + wp*wphi(0)*E_phi_ff(:,0)
        E = E + wp*wphi(1)*E_phi_ff(:,1)
#else

        E   = E   + wp*wphi(0)*psn%E_phi_ff(:,0,node,iphi)
        E   = E   + wp*wphi(1)*psn%E_phi_ff(:,1,node,iphi)
#endif
        
        if(sml_extra_dwdt) then           
           E00 = E00 + wp*wphi(0)*psn%E00_ff(:,0,node)
           E00 = E00 + wp*wphi(1)*psn%E00_ff(:,1,node)
           ddpotdt = ddpotdt + wp*wphi(0)*psn%ddpotdt_phi(node,0,iphi)
           ddpotdt = ddpotdt + wp*wphi(1)*psn%ddpotdt_phi(node,1,iphi)
           pot = pot + wp*wphi(0)*psn%pot_phi_ff(0,node,iphi) 
           pot = pot + wp*wphi(1)*psn%pot_phi_ff(1,node,iphi)
        endif

     enddo
  else
     print *, 'E-field ion invalid tr : (i,itr,mype,gid)=', i,itr,sml_mype,sp%ptl(i)%gid
  endif
  
  !E(3) was para E-field and becomes Ephi
  if(sml_turb_efield) then
     B=sqrt(fld%br**2 + fld%bz**2 + fld%bphi**2)
     E(3)=(E(3)*B- E(1)*fld%Br - E(2)*fld%Bz)/fld%Bphi
     
  else
     E(3)=0D0
  endif
  
  fld%Er=E(1)
  fld%Ez=E(2)
  fld%Ephi=E(3)
  if(sml_extra_dwdt) then
     fld%Epot=pot
     fld%Er00=E00(1)
     fld%Ez00=E00(2) 
     fld%ddpotdt=ddpotdt
  endif
end subroutine efield_gk_elec2
#endif

