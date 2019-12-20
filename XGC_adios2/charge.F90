subroutine chargee(grid,psn,sp)
  use grid_class
  use psn_class
  use sml_module
  use eq_module, only : is_rgn12, eq_x_psi,eq_x_z
  use ptl_module
  use smooth_module
  use random_xgc
  use pol_decomp_module
  use omp_module, only: split_indices
  use perf_monitor
  implicit none
  include 'mpif.h'
  type(grid_type):: grid
  type(psn_type) :: psn
  type(species_type) :: sp
  real (8), parameter :: weight_min = -1000D0 ! minimum weight (w1) in chargee

  call t_startf("CHARGEE_SEARCH_INDEX")
  call chargee_search_index(grid,psn,sp)
  call t_stopf("CHARGEE_SEARCH_INDEX")

#ifndef DELTAF_MODE2
  call t_startf("UPDATE_ELEC_WEIGHT")
  call update_elec_weight      ! update electron weight
  call t_stopf("UPDATE_ELEC_WEIGHT")
#endif
  call t_startf("CHARGEE_SCATTER")
  call chargee_scatter         ! sub-subroutine
  call t_stopf("CHARGEE_SCATTER")

  call t_startf("CHARGEE_MPISUM")
  call chargee_mpisum(grid,psn,sp)   ! corresponding chargei_gyro_average
#ifdef F0_CHARGE_N0
  if(sml_f0_grid) then
    call chargee_mpisum_f0(grid,psn)
  endif
#endif
  call t_stopf("CHARGEE_MPISUM")

#ifdef XGC1_EM
!  if(sml_gstep>=1)call remove_chargei(psn)
  if(sml_mode_select_on==1)call mode_selection_comb(sml_mode_select_n,grid,psn,psn%idensity(:,1))
#endif
contains
  subroutine chargee_scatter
    implicit none
    real (8),allocatable :: den(:,:,:), den00(:,:)
    real :: inv_delta_phi ! for fast calculation -- const
    integer :: ith, i_beg(sml_nthreads), i_end(sml_nthreads)
    integer :: nn_beg(sml_nthreads), nn_end(sml_nthreads)
    ! in loop
    integer ::  i,irho, j, node
    real (8) :: phi, wphi(0:1), particle_weight
    real (8) :: wp
    ! for 1D
    integer :: larmor, ip
    real (8) :: x(2), psi, pn
    !
    integer :: jth, iphi
    real (8), external :: psi_interpol

    !allocate memory and deallocate it at the end of subroutine
    allocate(den(0:1,grid%nnode,sml_nthreads))  ! indexing is different for cache performance
    allocate(den00(grid%npsi00,sml_nthreads))

    inv_delta_phi=1D0/grid%delta_phi

    ! for reduction, divide particles among OpenMP threads
    ! and calculate contributions to nodes for each subset
    ! of particles independently. This will introduce a
    ! roundoff difference in the results for different numbers
    ! of threads.
    call split_indices(sp%num, sml_nthreads, i_beg, i_end)

!$OMP PARALLEL DO &
!$OMP PRIVATE( ITH, I, IRHO, J, NODE, &
!$OMP          PHI, WPHI, PARTICLE_WEIGHT, WP, &
!$OMP          IP, X, PSI, PN )
    do ith=1,sml_nthreads
       call t_startf("CHARGEE_SCATTER_LOOP1")

       den(:,:,ith) =0D0
       den00(:,ith) =0D0

       do i=i_beg(ith),i_end(ith)
          if(sp%ptl(i)%gid<=0) cycle

          if((sp%tr_save(i)) > 0) then

             x=sp%ptl(i)%ph(1:2)  ! for 1D potential


             ! phi weight
             phi=sp%ptl(i)%ph(3)
             wphi(1) = (phi*inv_delta_phi) - grid%iphi_offset ! larger index weight
             wphi(0) = 1D0 - wphi(1)  ! smaller index weight


             ! particle weight
             if(sml_deltaf_elec) then  !
                particle_weight=max(sp%ptl(i)%ph(piw1),weight_min)*sp%ptl(i)%ct(piw0)
             else
                particle_weight=sp%ptl(i)%ct(piw0) ! for full f simulation only
             endif

             do j=1, 3
                node=grid%nd(j,sp%tr_save(i))
                wp=sp%p_save(j,i)

                den(0:1,node,ith) = den(0:1,node,ith) +  wp*particle_weight*wphi(0:1)
!                den(0:1,node,ith) = den(0:1,node,ith) +  wp*particle_weight*wphi(0:1)
             enddo


             !1D ---

             psi = psi_interpol(x(1),x(2),0,0)
             pn=(psi-grid%psi00min)/grid%dpsi00
             ip=floor(pn)+1
!             if(0 < ip .and. ip < grid%npsi00 .and. (x(2) > eq_x_z .or. psi > eq_x_psi)) then
             if(0 < ip .and. ip < grid%npsi00 .and. is_rgn12(x(1),x(2),psi) .and. (.not. sml_00_xz_up .or. x(2)>eq_x_z)) then
                wp=1D0 - ( pn - real(ip-1,8) )  ! recylce of wp -- differenet meaning from above : w of p(1:3) vs w of psi
                den00(ip  ,ith)=den00(ip  ,ith) + particle_weight* wp
                den00(ip+1,ith)=den00(ip+1,ith) + particle_weight*(1D0-wp)
             endif


          else

             !eliminate particle
             call remove_particle(sp,i,-1,ith) ! no need due to sheath_calculation

          endif

       enddo

       call t_stopf("CHARGEE_SCATTER_LOOP1")
    enddo  ! end of particle-thread loop




    ! open mp sum

    ! combine results
    call split_indices(grid%nnode, sml_nthreads, nn_beg, nn_end)

!$OMP PARALLEL DO &
!$OMP PRIVATE( JTH, ITH, NODE, IRHO, IPHI )
    do jth=1,sml_nthreads
       call t_startf("CHARGEE_SCATTER_LOOP2")

       do ith=2,sml_nthreads
          do node=nn_beg(jth), nn_end(jth)
             do iphi=0, 1
                den(iphi,node,1)=den(iphi,node,1)+ den(iphi,node,ith)
             enddo
          enddo
          ! compare performance
          ! with     den(:,:,:,1)=den(:,:,:,1)+den(:,:,:,ith)
       enddo

       ! indexing order change
       do iphi=0, 1
          do node=nn_beg(jth), nn_end(jth)
             psn%edensity(node,iphi) = den(iphi,node,1)
          enddo
       enddo
       call t_stopf("CHARGEE_SCATTER_LOOP2")
    enddo

    ! 1D
    do ith=2, sml_nthreads
       den00(:,1)=den00(:,1)+den00(:,ith)
    enddo
    psn%eden00_1d=den00(:,1)

    deallocate(den,den00)

  end subroutine chargee_scatter

  ! update electron weight from f0 change
  subroutine update_elec_weight
    implicit none
    integer :: i
    integer :: ith, i_beg(sml_nthreads), i_end(sml_nthreads)
    real (8) :: r,z,B,psi,En,new_f0,w2,dw
    real (8), external :: get_f0_elec
    logical, external :: is_nan
    call split_indices(sp%num, sml_nthreads, i_beg, i_end)

!$OMP PARALLEL DO &
!$OMP PRIVATE( ITH, I, NEW_F0, W2, DW )
    do ith=1, sml_nthreads
       call t_startf("UPD_ELEC_WT_LOOP")
       do i=i_beg(ith), i_end(ith)
          if(sp%ptl(i)%gid>0) then

             new_f0 = get_f0_elec(grid,psn,sp%ptl(i),sp%tr_save(i),sp%p_save(:,i))

             w2= 1D0 - new_f0/sp%ptl(i)%ct(pif0)

             dw = w2 - sp%ptl(i)%ph(piw2)
             sp%ptl(i)%ph(piw1)= sp%ptl(i)%ph(piw1) + dw
             sp%ptl(i)%ph(piw2)= w2


             !call restrict_weight(sp%ptl(i)%ph(piw1:piw2), new_f0)

             !if(sml_ipc==2) sp%ptl(i)%ct(pif0)=new_f0
             !sp%ptl(i)%ct(pif0)=new_f0  ! update f0 always unlike ion --> f0 is saved on ptl_ephase_save

             !NaN test
             if(is_nan(w2) .or. is_nan(dw)) then
                print *, 'NaN found in update_elec_weight : new_f0, w2, dw', new_f0, w2, dw
                stop
             endif

          endif
       enddo
       call t_stopf("UPD_ELEC_WT_LOOP")
    enddo
  end subroutine update_elec_weight


end subroutine chargee

! set minimum weight to prevent weight to be smaller than 0
subroutine restrict_weight(w,f0)
  implicit none
  real (8) :: w(2),f0, tmp1,tmp2
  real (8), parameter :: weight_min = -100D0
  real (8), parameter :: weight_max = 0.99999D0

  tmp1=f0/(1D0-w(1))
  tmp2=f0/(1D0-w(2))

  w(1)=max(w(1),weight_min)
  w(2)=max(w(2),weight_min)

  w(1)=min(w(1),weight_max)
  w(2)=min(w(2),weight_max)


  f0=(1D0-w(1))*tmp1
  f0=(1D0-w(2))*tmp2

end subroutine restrict_weight


subroutine chargee_search_index(grid,psn,sp)
  use grid_class
  use psn_class
  use ptl_module
  use sml_module
  use omp_module, only: split_indices
  use perf_monitor
  implicit none
  type(grid_type):: grid
  type(psn_type) :: psn
  type(species_type) :: sp
  !
  integer :: ith, i
  integer :: i_beg(sml_nthreads), i_end(sml_nthreads)
  real (8) :: phi_mid, x(2), phi, mu, rho, xff(2)
  real (8) :: p(3)
  integer :: itr, ip
  real (kind=8), external :: gyro_radius
  logical, parameter :: USE_SEARCH_TR2 = .true.

  !phi_mid=0.5D0*(grid%phimin+grid%phimax)

  call split_indices(sp%num, sml_nthreads, i_beg, i_end)

!$OMP PARALLEL DO &
!$OMP PRIVATE( ITH, I, PHI, PHI_MID, &
!$OMP          X, MU, RHO, XFF,  &
!$OMP          P, ITR, IP )
  do ith=1,sml_nthreads
     call t_startf("CHARGEE_SRCHLOOP")
     do i=i_beg(ith),i_end(ith)
        if(sp%ptl(i)%gid>0)  then

           ! get proper toroidal angle index and weight
           x=sp%ptl(i)%ph(1:2)
           phi=sp%ptl(i)%ph(3)
           phi_mid=(floor(phi/grid%delta_phi) + 0.5D0) * grid%delta_phi
!           mu=sp%ptl(i)%ct(pim)
!           rho=gyro_radius(x,mu)  !gyro radius
!           sp%rhoi(i)=rho

             ! get field following posision at 1/2 angle
           call field_following_pos2(x,phi,phi_mid,xff)

           ! find triangle
           if (USE_SEARCH_TR2) then
              call search_tr2(grid,xff,itr,p)
           else
              call search_tr(grid,xff,itr,p)
           endif


           !remove particle or sheath calculation
           if(itr<0) then
              if(sml_sheath_mode==0 .or. sml_gstep <= 0 ) then
                 call remove_particle(sp,i,-1,ith)
              else
!$omp critical
                 call sheath_calculation(grid,psn,sp,i,0,itr,p,ith)
!$omp end critical
              endif
           endif

           sp%tr_save(i)=itr
           sp%p_save(:,i)=p
           ! compare performance
           !do ip=1,3
           !   sp%p_save(ip,i)=p(ip)
           !enddo
        endif
     enddo
     call t_stopf("CHARGEE_SRCHLOOP")
  enddo
end subroutine chargee_search_index

#ifdef XGC1_EM
subroutine ijpar_mpisum(grid,psn,sp)
  use grid_class
  use psn_class
  use sml_module
  use ptl_module
  use perf_monitor
  implicit none
  type(grid_type):: grid
  type(psn_type) :: psn
  type(species_type) :: sp
  include 'mpif.h'

  integer :: nn, ierr
  integer :: iphi, ipe
  integer :: idest, isendtag, isource, irecvtag
  integer, dimension(MPI_STATUS_SIZE) :: istatus
  real (8) :: tmp(grid%nnode,0:1)
  real (8) :: inv_nphi_total

  nn=grid%nnode

  call mpi_reduce(psn%ijpar_ff,tmp,nn*2, &
       MPI_REAL8,mpi_sum,0,sml_plane_comm,ierr)

  ! coordinate conversion
  if(sml_plane_mype==0) then
    tmp=tmp/grid%node_vol_ff

  ! coordinate conversion
    call cnvt_grid_ff2real(grid,psn%ff_hdp_tr,psn%ff_hdp_p,tmp,psn%ijpar_ff)
  endif

  ! send 0 plane receive 1 plane
  if(sml_plane_mype==0) then
    idest=mod(sml_intpl_mype-1+sml_intpl_totalpe,sml_intpl_totalpe)
    isendtag=sml_intpl_mype

    isource=mod(sml_intpl_mype+1,sml_intpl_totalpe)
    irecvtag=isource
    call mpi_sendrecv(psn%ijpar_ff(:,0),nn,MPI_REAL8,  idest,isendtag, &
         tmp(:,0),nn,MPI_REAL8,isource,irecvtag, &
         sml_intpl_comm,istatus,ierr)
    psn%ijpar_ff(:,1) = 0.5D0*(psn%ijpar_ff(:,1) + tmp(:,0)) ! tmp(:,0) is ijpar received

  endif

  call mpi_bcast(psn%ijpar_ff(:,1),nn,mpi_real8,0,sml_plane_comm,ierr)
end subroutine ijpar_mpisum

subroutine chargee_hyb_mpisum(grid,psn)
  use grid_class
  use psn_class
  use sml_module
!  use eq_module, only : eq_x_z
  use ptl_module
  use perf_monitor
  use smooth_module
  implicit none
  type(grid_type):: grid
  type(psn_type) :: psn
  include 'mpif.h'
  !
  integer :: nn, ierr
  integer :: irho, iphi, ipe
  integer :: idest, isendtag, isource, irecvtag
  integer, dimension(MPI_STATUS_SIZE) :: istatus
  real (8) :: tmp(grid%nnode)
  real (8) :: dum00(grid%npsi00)
  real (8) :: inv_nphi_total
  !
  character (len=30) :: filename
  integer :: i

  nn=grid%nnode
!  call continutity_eq(psn, grid)

  ! idensity0 -- sum-up
  call t_startf("CHARGEE_SUM_RED")
  call mpi_allreduce(psn%eden_hyb,tmp,nn, &
       MPI_REAL8,mpi_sum,sml_intpl_comm,ierr)
  call t_stopf("CHARGEE_SUM_RED")

  inv_nphi_total=1D0/real(sml_nphi_total)
  psn%edensity0=tmp*inv_nphi_total


end subroutine chargee_hyb_mpisum

subroutine remove_chargei(psn)
  use psn_class
  use sml_module
  implicit none
  type(psn_type) :: psn
  if(sml_electron_hyb) then
     psn%idensity0=0D0
     psn%idensity=0D0
     psn%iden00_1d=0D0
     psn%cden00_1d=0D0
     psn%ijpar_ff=0D0
  endif

end subroutine remove_chargei
#endif

subroutine chargee_mpisum(grid,psn,sp)
  use grid_class
  use psn_class
  use sml_module
!  use eq_module, only : eq_x_z
  use ptl_module
  use perf_monitor
  use smooth_module
  use f0_module
  implicit none
  type(grid_type):: grid
  type(psn_type) :: psn
  type(species_type) :: sp
  include 'mpif.h'
  !
  integer :: nn, ierr
  integer :: irho, iphi, ipe
  integer :: idest, isendtag, isource, irecvtag
  integer, dimension(MPI_STATUS_SIZE) :: istatus
  real (8) :: tmp(grid%nnode,0:1)
  real (8) :: dum00(grid%npsi00)
  real (8) :: inv_nphi_total
  !
  character (len=30) :: filename
  integer :: i

  nn=grid%nnode

#ifndef F0_CHARGE_N0
  if(sml_f0_grid) then ! omp?  -- 0.5 is multiplied already in eden_f0
#ifndef F0_TOR_LINEAR
    psn%edensity(:,0)=psn%edensity(:,0) + psn%eden_f0*grid%node_vol_ff(:,0)
    psn%edensity(:,1)=psn%edensity(:,1) + psn%eden_f0*grid%node_vol_ff(:,1)
#else
    psn%edensity(:,0)=psn%edensity(:,0) + psn%eden_f0(:,0)*grid%node_vol_ff(:,0)
    psn%edensity(:,1)=psn%edensity(:,1) + psn%eden_f0(:,1)*grid%node_vol_ff(:,1)
#endif
  endif
#endif

  ! send 0-rho to proc 0 of plane_comm
  call t_startf("CHARGEE_SUM_RED")
  call mpi_reduce(psn%edensity,tmp,nn*2, &
       MPI_REAL8,mpi_sum,0,sml_plane_comm,ierr)
  call t_stopf("CHARGEE_SUM_RED")

  ! coordinate conversion
  if(sml_plane_mype==0) then
     call t_startf("CHARGEE_SUM_CNVRT")
     !make it density
     ! omp?
     tmp=tmp/grid%node_vol_ff


     ! coordinate conversion
     call cnvt_grid_ff2real(grid,psn%ff_hdp_tr,psn%ff_hdp_p,tmp,psn%edensity)
     call t_stopf("CHARGEE_SUM_CNVRT")
  endif


  ! send 0 plane receive 1 plane
  if(sml_plane_mype==0) then

     call t_startf("CHARGEE_SUM_SR_RED")
     idest=mod(sml_intpl_mype-1+sml_intpl_totalpe,sml_intpl_totalpe)
     isendtag=sml_intpl_mype

     isource=mod(sml_intpl_mype+1,sml_intpl_totalpe)
     irecvtag=isource
     call mpi_sendrecv(psn%edensity(:,0),nn,MPI_REAL8,  idest,isendtag, &
          tmp(:,0),nn,MPI_REAL8,isource,irecvtag, &
          sml_intpl_comm,istatus,ierr)

     psn%edensity(:,1) = 0.5D0*(psn%edensity(:,1) + tmp(:,0)) ! tmp(:,0) is density received

     ! idensity0 -- sum-up
     call mpi_allreduce(psn%edensity(:,1),psn%edensity0,nn, mpi_real8, mpi_sum, sml_intpl_comm,ierr)
     call t_stopf("CHARGEE_SUM_SR_RED")

     inv_nphi_total=1D0/real(sml_nphi_total)
     psn%edensity0=psn%edensity0*inv_nphi_total

     ! smoothing
!pw  call t_startf("CHARGEE_SUM_SMOOTH_R")
     call smooth_r(psn%edensity0(:),tmp(:,0),smooth_r1,grid)
!pw  call t_stopf("CHARGEE_SUM_SMOOTH_R")
     psn%edensity0=tmp(:,0)
  endif


  ! 00 mode
  call t_startf("CHARGEE_SUM_RED_BCAST")

  !if(.false.) then
  !call my_mpi_allreduce(psn%eden00_1d,dum00,grid%npsi00)
  !psn%eden00_1d=dum00/psn%vol00
  !endif

  call mpi_bcast(psn%edensity(:,1),nn,mpi_real8,0,sml_plane_comm,ierr)
  call mpi_bcast(psn%edensity0    ,nn,mpi_real8,0,sml_plane_comm,ierr)
  call t_stopf("CHARGEE_SUM_RED_BCAST")

  if(.true.) then ! Use edensity0 to get eden00_1d
    call convert_grid_2_001d(grid,psn%edensity0,psn%eden00_1d)
  endif

#ifdef F0_CHARGE_N0

  f0_den0_ptl(:,0) = psn%edensity0

#endif


#ifndef F0_CHARGE_N0
  if(sml_f0_grid) then
  !prevent total density becomes below 0  -- following formula is not exact for electron
   do i=1, grid%nnode
      if(psn%edensity(i,1) < -f0_den_global(i)) then
          psn%edensity(i,1)=-f0_den_global(i)
      endif
      if(psn%edensity0(i) < -f0_den_global(i)) then
        psn%edensity0(i)=-f0_den_global(i)
      endif
   enddo
  endif
#endif

end subroutine chargee_mpisum

!only required with F0_CHARGE_N0
subroutine chargee_mpisum_f0(grid,psn)
  use grid_class
  use psn_class
  use sml_module
!  use eq_module, only : eq_x_z
  use ptl_module
  use perf_monitor
  use smooth_module
  use f0_module
  implicit none
  type(grid_type):: grid
  type(psn_type) :: psn
  include 'mpif.h'
  !
  integer :: nn, ierr
  integer :: irho, iphi, ipe
  integer :: idest, isendtag, isource, irecvtag
  integer, dimension(MPI_STATUS_SIZE) :: istatus
  real (8) :: tmp(grid%nnode,0:1)
  real (8) :: inv_nphi_total
  !
  character (len=30) :: filename
  integer :: i

  nn=grid%nnode

  ! send 0-rho to proc 0 of plane_comm
  call t_startf("CHARGEE_SUM_RED")
  call mpi_reduce(psn%eden_f0,tmp(:,0),nn, &
       MPI_REAL8,mpi_sum,0,sml_plane_comm,ierr)
  call t_stopf("CHARGEE_SUM_RED")

  ! coordinate conversion
  if(sml_plane_mype==0) then
     call t_startf("CHARGEE_SUM_CNVRT")

     tmp(:,1)=tmp(:,0)

     ! coordinate conversion
     call cnvt_grid_ff2real(grid,psn%ff_hdp_tr,psn%ff_hdp_p,tmp,psn%edensity_f0)
     call t_stopf("CHARGEE_SUM_CNVRT")
  endif


  ! send 0 plane receive 1 plane
  if(sml_plane_mype==0) then

     call t_startf("CHARGEE_SUM_SR_RED")
     idest=mod(sml_intpl_mype-1+sml_intpl_totalpe,sml_intpl_totalpe)
     isendtag=sml_intpl_mype

     isource=mod(sml_intpl_mype+1,sml_intpl_totalpe)
     irecvtag=isource
     call mpi_sendrecv(psn%edensity(:,0),nn,MPI_REAL8,  idest,isendtag, &
          tmp(:,0),nn,MPI_REAL8,isource,irecvtag, &
          sml_intpl_comm,istatus,ierr)

     psn%edensity_f0(:,1) = 0.5D0*(psn%edensity_f0(:,1) + tmp(:,0)) ! tmp(:,0) is density received

     ! extract n=0
     call mpi_allreduce(psn%edensity_f0(:,1),tmp(:,0),nn, mpi_real8, mpi_sum, sml_intpl_comm,ierr)
     call t_stopf("CHARGEE_SUM_SR_RED")

     inv_nphi_total=1D0/real(sml_nphi_total)
     tmp(:,0)=tmp(:,0)*inv_nphi_total

     psn%edensity_f0(:,1)=psn%edensity_f0(:,1) - tmp(:,0) + f0_density_n0_add(:,0) ! last 0 in f0_density_n0_add is electron index
  endif

  call t_startf("CHARGEE_SUM_RED_BCAST")
  call mpi_bcast(psn%edensity_f0(:,1),nn,mpi_real8,0,sml_plane_comm,ierr)
  call t_stopf("CHARGEE_SUM_RED_BCAST")

  psn%edensity(:,1) = psn%edensity(:,1) + psn%edensity_f0(:,1)
  psn%edensity0     = psn%edensity0     + f0_density_n0_add(:,0)


  !prevent total density becomes below 0  -- following formula is not exact for electron
  do i=1, grid%nnode
    if(psn%edensity(i,1) < -f0_den_global(i)) then
        psn%edensity(i,1)=-f0_den_global(i)
    endif
    if(psn%edensity0(i) < -f0_den_global(i)) then
        psn%edensity0(i)=-f0_den_global(i)
    endif
  enddo

end subroutine chargee_mpisum_f0

#ifdef OLD_INIT_FF

subroutine init_bfollow(grid,psn)
  use sml_module
  use grid_class
  use psn_class
  use eq_module, only: eq_axis_r
  use comm_mod, only : allgatherv_2d_ring, allgatherv_2i_ring
  implicit none
  include 'mpif.h'
  integer :: i,j
  type(psn_type) :: psn
  type(grid_type) :: grid
  real (kind=8), allocatable :: xp(:,:),xm(:,:)
  integer , allocatable :: nold(:)
  integer :: itr,dir,min_dist_node,nd,node_new,node_old,dum
  real (kind=8) :: min_dist,p(3),x(2),xtmp(4),dist,dphi,phi
  integer :: iseg_offset,segnum, segend, segpe, ierror

  logical, parameter :: USE_SEARCH_TR2 = .true.
  logical, parameter :: use_bcast = .true.  ! .false. does not work on edison
  integer :: n1, n2

  !allocate memeory
  allocate(nold(grid%nnode), xp(2,grid%nnode), xm(2,grid%nnode) )
  nold=0 !safty
  xp=0D0 !safty
  xm=0D0 !safty

  if(sml_bfollow==1) then
     ! Generates bfollow
     nold=0
     phi=0D0 ! Cannot specifiy toroidal angle - symmetry is assumed.
     do i=1, grid%nnode
        dir=1
        if(sml_bt_sign<0D0) dir=-dir

        dphi=real(dir)*grid%delta_phi
        ! field following direction
        call field_following_pos2(grid%x(:,i),phi,phi+dphi,xp(:,i))

        ! opposite direction to B
        dphi=-dphi
        call field_following_pos2(grid%x(:,i),phi,phi+dphi,xm(:,i))

     enddo
  else
     xp(:,:)=grid%x(:,:)
     xm(:,:)=grid%x(:,:)
     nold(:)=0
  endif

  !  print *, 'end b-following file'

  do i=1, grid%nnode
     if(nold(i)==-1) then
        xp(:,i)=grid%x(:,i)
        xm(:,i)=grid%x(:,i)
     endif
  enddo

  !Parallizing below routine
  ! algorithm - each pe find each segment

  segnum = grid%nnode/sml_totalpe
  if(grid%nnode > sml_totalpe*segnum ) then
     segnum=segnum+1
  endif
  iseg_offset = sml_mype*segnum
  segend = min(grid%nnode, iseg_offset + segnum)



  do i=1+iseg_offset, segend
     !     print *, i
     if (sml_cylindrical) then
       psn%bfollow_1_dx(i)=1D0/sqrt( (xp(1,i)-xm(1,i))**2 + (xp(2,i)-xm(2,i))**2 + (eq_axis_r*grid%delta_phi*2D0)**2)
     else
       psn%bfollow_1_dx(i)=1D0/sqrt( (xp(1,i)-xm(1,i))**2 + (xp(2,i)-xm(2,i))**2 + (grid%x(1,i)*grid%delta_phi*2D0)**2)
     endif
     ! What is exact formula for curved derivative of curved system ??
     ! Here we assumed that it is the same with XYZ derivative of XYZ coordinate

     ! for two directions
     do dir=1,2
        if(dir==1) then
           x=xp(:,i)
        else
           x=xm(:,i)
        endif
        !        print *, dir , 'search tr',x

        if (USE_SEARCH_TR2) then
           call search_tr2(grid,x,itr,p)
        else
           call search_tr(grid,x,itr,p)
        endif

        !find  node and p
        !        print *, dir, 'search tr end'
        if(itr<=0) then
           !find nearest node point
           !           print *, 'find_nearest',dir
           min_dist=1D50
           do nd=1, grid%nnode
              dist=(grid%x(1,nd)-x(1))**2 + (grid%x(2,nd)-x(2))**2
              if(min_dist > dist) then
                 min_dist=dist
                 min_dist_node=nd
              endif
           enddo
           !           print *, 'min',min_dist, min_dist_node
           !choose one triangle that has this node point
           itr=grid%tr_node(1,min_dist_node)
           do j=1,3
              if(min_dist_node==grid%nd(j,itr)) then
                 p(j)=1D0
              else
                 p(j)=0D0
              endif
           enddo
        endif
        if(itr<=0 .or. itr> grid%ntriangle) then
           print *, 'Wrong itr', itr,grid%ntriangle
           call err_count
        endif
        psn%bfollow_tr(dir,i) = itr
        psn%bfollow_p(:,dir,i) = p
     enddo
  enddo

  !propagate result

  if (use_bcast) then
  do segpe=0, sml_totalpe-1
     iseg_offset = segpe*segnum
     segend = min(grid%nnode, iseg_offset + segnum)
     if(segend-iseg_offset >0 ) then
         call MPI_BCAST(psn%bfollow_tr(:,iseg_offset+1:segend) , (segend-iseg_offset)*2  ,&
              MPI_INTEGER, segpe, sml_comm, ierror)
         call MPI_BCAST(psn%bfollow_p(:,:,iseg_offset+1:segend) , (segend-iseg_offset)*6  ,&
              MPI_REAL8, segpe, sml_comm, ierror)
         call MPI_BCAST(psn%bfollow_1_dx(iseg_offset+1:segend), (segend-iseg_offset) ,&
              MPI_REAL8, segpe, sml_comm, ierror)
     endif
  enddo
  else

    n1 = size(psn%bfollow_tr,1)
    n2 = min(grid%nnode, size( psn%bfollow_tr, 2))
    call allgatherv_2i_ring( psn%bfollow_tr, n1,n2, segnum, sml_comm )

    n1 = size( psn%bfollow_p, 1)
    n2 = min(grid%nnode, size(psn%bfollow_p,2))
    call allgatherv_2d_ring( psn%bfollow_p, n1, n2, segnum, sml_comm )

    n1 = 1
    n2 = min(grid%nnode, size(psn%bfollow_1_dx,1))
    call allgatherv_2d_ring( psn%bfollow_1_dx, n1, n2, segnum, sml_comm)
 

  endif

  !deallocate memeory
  deallocate(nold,xp,xm)

!1000 format (2I8,1x,4(e19.13,1x))

end subroutine init_bfollow

subroutine init_ff(grid,psn,half)
  use sml_module
  use grid_class
  use psn_class

  use comm_mod, only : allgatherv_2d_ring, allgatherv_2i_ring
  implicit none
  include 'mpif.h'
  integer :: i,j
  type(psn_type) :: psn
  type(grid_type) :: grid
  logical, intent(in) :: half
!  real (kind=8), allocatable :: xp(:,:),xm(:,:)
  integer :: itr,dir,min_dist_node,nd,node_new,node_old,dum
  real (kind=8) :: min_dist,p(3),x(2),xtmp(4),dist,dphi,phi
  integer :: iseg_offset,segnum, segend, segpe, ierror
  real (8) :: delta_phi
  logical, parameter :: USE_SEARCH_TR2 = .true.

  logical, parameter :: use_bcast = .true.  ! .false. does not work on edison
  integer :: n1, n2

  if(half) then
     delta_phi = 0.5D0* grid%delta_phi
  else
     delta_phi = grid%delta_phi
  end if


  ! Generates bfollow
  phi=0D0 ! Cannot specifiy toroidal angle - symmetry is assumed.





  do dir=0,1

     !Parallizing below routine
     ! algorithm - each pe find each segment
     segnum = grid%nnode/sml_totalpe
     if(grid%nnode > sml_totalpe*segnum ) then
        segnum=segnum+1
     endif
     iseg_offset = sml_mype*segnum
     segend = min(grid%nnode, iseg_offset + segnum)


     do i=1+iseg_offset, segend

        ! find field following position
        if(dir==1) then
           dphi=delta_phi
           ! phi direction
           call field_following_pos2(grid%x(:,i),phi,phi+dphi,x)
        else
           ! -phi direction
           dphi=-delta_phi
           call field_following_pos2(grid%x(:,i),phi,phi+dphi,x)
        endif

        ! triangle search
        if (USE_SEARCH_TR2) then
           call search_tr2(grid,x,itr,p)
        else
           call search_tr(grid,x,itr,p)
        endif


        !find  node and p
        !        print *, dir, 'search tr end'
        if(itr<=0) then
           !find nearest node point
           !           print *, 'find_nearest',dir
           min_dist=1D50
           do nd=1, grid%nnode
              dist=(grid%x(1,nd)-x(1))**2 + (grid%x(2,nd)-x(2))**2
              if(min_dist > dist) then
                 min_dist=dist
                 min_dist_node=nd
              endif
           enddo
           !           print *, 'min',min_dist, min_dist_node
           !choose one triangle that has this node point
           itr=grid%tr_node(1,min_dist_node)
           do j=1,3
              if(min_dist_node==grid%nd(j,itr)) then
                 p(j)=1D0
              else
                 p(j)=0D0
              endif
           enddo
        endif

        if(itr<=0 .or. itr> grid%ntriangle) then
           print *, 'Wrong itr', itr,grid%ntriangle
           call err_count
        endif

        !store it
        if(half) then
           psn%ff_hdp_tr(i,dir)= itr
           psn%ff_hdp_p(:,i,dir)=p
        else
           psn%ff_1dp_tr(i,dir) = itr
           psn%ff_1dp_p(:,i,dir) = p
        endif
     enddo

     !propagate result
     if (use_bcast) then

     do segpe=0, sml_totalpe-1
        iseg_offset = segpe*segnum
        segend = min(grid%nnode, iseg_offset + segnum)
        if(segend-iseg_offset >0 ) then
           if(half) then
              call MPI_BCAST(psn%ff_hdp_tr(iseg_offset+1:segend,dir) , (segend-iseg_offset)  ,&
                   MPI_INTEGER, segpe, sml_comm, ierror)
              call MPI_BCAST(psn%ff_hdp_p(:,iseg_offset+1:segend,dir) , (segend-iseg_offset)*3  ,&
                   MPI_REAL8, segpe, sml_comm, ierror)
           else
              call MPI_BCAST(psn%ff_1dp_tr(iseg_offset+1:segend,dir) , (segend-iseg_offset)  ,&
                   MPI_INTEGER, segpe, sml_comm, ierror)
              call MPI_BCAST(psn%ff_1dp_p(:,iseg_offset+1:segend,dir) , (segend-iseg_offset)*3  ,&
                   MPI_REAL8, segpe, sml_comm, ierror)
           endif
        endif
     enddo

    else

     if (half) then
        n1 = 1
        n2 = min(grid%nnode,size(psn%ff_hdp_tr,1))
        call allgatherv_2i_ring( psn%ff_hdp_tr(:,dir), n1,n2, segnum, sml_comm)

        n1 = size(psn%ff_hdp_p,1)
        n2 = min(grid%nnode,size(psn%ff_hdp_p,2))
        call allgatherv_2d_ring( psn%ff_hdp_p(:,:,dir), n1,n2, segnum,sml_comm)
      else

        n1 = 1
        n2 = min(grid%nnode, size(psn%ff_1dp_tr,1))
        call allgatherv_2i_ring( psn%ff_1dp_tr(:,dir), n1,n2,segnum, sml_comm)

        n1 = size(psn%ff_1dp_p,1)
        n2 = min(grid%nnode, size(psn%ff_1dp_p,2))
        call allgatherv_2d_ring( psn%ff_1dp_p(:,:,dir), n1,n2,segnum, sml_comm)
      endif
    endif

!  debug
!     if(sml_mype==0) then
!        do i=1, grid%nnode
!           if(psn%ff_hdp_tr(i,dir)==0) then
!              print *, 'zero itr', i, dir, psn%ff_hdp_tr(i,dir)
!           endif
!        enddo
!     endif
  enddo


!1000 format (2I8,1x,4(e19.13,1x))

end subroutine init_ff

#else
! not(OLD_INIT_FF)

subroutine init_ff(grid,psn)
  use sml_module
  use grid_class
  use psn_class
  use eq_module, only: eq_axis_r
  implicit none
  include 'mpif.h'
  integer :: i,j
  type(psn_type) :: psn
  type(grid_type) :: grid
  integer :: half
!  real (kind=8), allocatable :: xp(:,:),xm(:,:)
  integer :: itr,dir,min_dist_node,nd,node_new,node_old,dum
  real (kind=8) :: min_dist,p(3),x(2,0:2),dist,dphi,phi,lens(0:2)
  integer :: iseg_offset,segnum, segend, segpe, ierror
  real (8) :: delta_phi
  logical, parameter :: USE_SEARCH_TR2 = .true.

! first half, then another half

  delta_phi = 0.5D0*grid%delta_phi
 
  ! Generates bfollow
  phi=0D0 ! Cannot specifiy toroidal angle - symmetry is assumed.

  do dir=0,1

     !Parallizing below routine
     ! algorithm - each pe find each segment
     segnum = grid%nnode/sml_totalpe
     if(grid%nnode > sml_totalpe*segnum ) then
        segnum=segnum+1
     endif
     iseg_offset = sml_mype*segnum
     segend = min(grid%nnode, iseg_offset + segnum)


     do i=1+iseg_offset, segend
        do half=0,1
          if(half==0) then
              x(:,0)=grid%x(:,i)
          endif
          ! find field following position

          if(dir==1) then
             dphi=+delta_phi
          else
             dphi=-delta_phi
          endif

          ! rotational symmetry assumed
          call field_following_pos2(x(:,half),phi,phi+dphi,x(:,half+1))

          ! triangle search
          if (USE_SEARCH_TR2) then
            call search_tr2(grid,x(:,half+1),itr,p)
          else
            call search_tr(grid,x(:,half+1),itr,p)
          endif

          !find  node and p
          !        print *, dir, 'search tr end'
          if(itr<=0) then
            !find nearest node point
            !           print *, 'find_nearest',dir
            min_dist=1D50
            do nd=1, grid%nnode
              dist=sum((grid%x(:,nd)-x(:,half+1))**2)
              if(min_dist > dist) then
                min_dist=dist
                min_dist_node=nd
              endif
            enddo
            !           print *, 'min',min_dist, min_dist_node
            !choose one triangle that has this node point
            itr=grid%tr_node(1,min_dist_node)
            do j=1,3
              if(min_dist_node==grid%nd(j,itr)) then
                p(j)=1D0
              else
                p(j)=0D0
              endif
            enddo
          endif

          if(itr<=0 .or. itr> grid%ntriangle) then
            print *, 'Wrong itr', itr,grid%ntriangle
            call err_count
          endif

          !store it
          if(half==0) then
            psn%ff_hdp_tr(i,dir)= itr
            psn%ff_hdp_p(:,i,dir)=p
            ! Cartesian length -- good enough for short steps
            !   mid-point \rho for good measure
            if (sml_cylindrical) then
              psn%ff_hdp_dx(i,dir)=sqrt(sum((x(:,1)-x(:,0))**2)+(eq_axis_r*delta_phi)**2)
            else
              psn%ff_hdp_dx(i,dir)=sqrt(sum((x(:,1)-x(:,0))**2)+0.25d0*((x(1,0)+x(1,1))*delta_phi)**2)
            endif
          else
            psn%ff_1dp_tr(i,dir) = itr
            psn%ff_1dp_p(:,i,dir) = p
            ! Simpson's rule for parabolic segment in cylinder coords
            lens(0)=sum((x(:,1)-x(:,0))**2)
            lens(1)=sum((x(:,2)-x(:,0))**2)
            lens(2)=sum((x(:,2)-x(:,1))**2)
            if (sml_cylindrical) then
              dist=sqrt(lens(1)+(eq_axis_r*2d0*delta_phi)**2)
              dist=dist+sqrt(0.75d0*lens(0)+0.25d0*lens(2)-.1875d0*lens(1)+0.25d0*(eq_axis_r*delta_phi)**2)
              dist=dist+sqrt(0.25d0*lens(0)+0.75d0*lens(2)-.1875d0*lens(1)+0.25d0*(eq_axis_r*delta_phi)**2)
              psn%ff_1dp_dx(i,dir)=2d0/3d0*dist
            else
              dist=sqrt(lens(1)+(x(1,1)*2d0*delta_phi)**2)
              dist=dist+sqrt(0.75d0*lens(0)+0.25d0*lens(2)-.1875d0*lens(1)+0.25d0*(x(1,0)*delta_phi)**2)
              dist=dist+sqrt(0.25d0*lens(0)+0.75d0*lens(2)-.1875d0*lens(1)+0.25d0*(x(1,2)*delta_phi)**2)
              psn%ff_1dp_dx(i,dir)=2d0/3d0*dist
            endif
          endif
       enddo   ! do half=0,1
     enddo

     !propagate result
     do segpe=0, sml_totalpe-1
        iseg_offset = segpe*segnum
        segend = min(grid%nnode, iseg_offset + segnum)
        if(segend-iseg_offset >0 ) then
              call MPI_BCAST(psn%ff_hdp_tr(iseg_offset+1:segend,dir) , (segend-iseg_offset)  ,&
                   MPI_INTEGER, segpe, sml_comm, ierror)
              call MPI_BCAST(psn%ff_hdp_p(:,iseg_offset+1:segend,dir) , (segend-iseg_offset)*3  ,&
                   MPI_REAL8, segpe, sml_comm, ierror)
              call MPI_BCAST(psn%ff_hdp_dx(iseg_offset+1:segend,dir) , (segend-iseg_offset)  ,&
                   MPI_REAL8, segpe, sml_comm, ierror)
              call MPI_BCAST(psn%ff_1dp_tr(iseg_offset+1:segend,dir) , (segend-iseg_offset)  ,&
                   MPI_INTEGER, segpe, sml_comm, ierror)
              call MPI_BCAST(psn%ff_1dp_p(:,iseg_offset+1:segend,dir) , (segend-iseg_offset)*3  ,&
                   MPI_REAL8, segpe, sml_comm, ierror)
              call MPI_BCAST(psn%ff_1dp_dx(iseg_offset+1:segend,dir) , (segend-iseg_offset)  ,&
                   MPI_REAL8, segpe, sml_comm, ierror)
        endif
     enddo

!  debug
!     if(sml_mype==0) then
!        do i=1, grid%nnode
!           if(psn%ff_hdp_tr(i,dir)==0) then
!              print *, 'zero itr', i, dir, psn%ff_hdp_tr(i,dir)
!           endif
!        enddo
!     endif
  enddo
! 1000 format (2I8,1x,4(e19.13,1x))
end subroutine init_ff

#endif
! not(OLD_INIT_FF)

! coordinate conversion -- field following coord to real coord
subroutine cnvt_grid_ff2real(grid,tr,p,var_ff,var_real)
  use grid_class
  use psn_class
  implicit none
  type(grid_type), intent(in) :: grid
  integer, intent(in) :: tr(grid%nnode,0:1)
  real (8), intent(in) :: p(3,grid%nnode,0:1)
  real (8), intent(in) :: var_ff(grid%nnode,0:1)
  real (8), intent(out) :: var_real(grid%nnode,0:1)
  real (8), allocatable :: vol(:,:)
  integer :: dir, i, j, cdir, nd
  real (8) :: wt
  !debug
  integer :: count
  count=0
  allocate(vol(grid%nnode,0:1))

  ! intialize
  var_real=0D0
  vol=0D0

  do dir=0, 1
     ! for all node point of ff
     do i=1, grid%nnode

        !### debug
        !if(tr(i,dir)<=0) print *, 'wrong tr', tr(i,dir),i,dir

        do j=1, 3
           ! node and weight
           nd=grid%nd(j,tr(i,dir))
           wt=p(j,i,dir)*grid%node_vol_ff(nd,dir)

           ! accumulate values
           var_real(nd,dir) =  var_real(nd,dir) +  wt*var_ff(i,dir)
           vol(nd,dir) = vol(nd,dir) + wt
        enddo
     enddo
  enddo

  !var_real=var_real/vol
  where(vol /= 0D0 ) 
    var_real=var_real/vol
  elsewhere
    var_real=0D0
  endwhere

  !fill empty area --- replace nan values
  do dir=0, 1
     cdir=1-dir

     do i=1, grid%nnode
        if(vol(i,dir)==0D0) then
           count=count+1  !debug
           var_real(i,dir)=0D0
           do j=1, 3
              nd=grid%nd(j,tr(i,cdir))
              wt=p(j,i,cdir)

              ! obtain interpolated values
              var_real(i,dir)=var_real(i,dir) + wt*var_ff(nd,dir)

           enddo
        endif
     enddo
  enddo

  deallocate(vol)
  !debug
  !print *, 'cnvt_grid_ff2real : empty count=',count

end subroutine cnvt_grid_ff2real


! get data from interpolation
subroutine cnvt_grid_real2ff(grid,tr,p,var_real,var_ff)
  use grid_class
  use psn_class
  implicit none
  type(grid_type), intent(in) :: grid
  integer, intent(in) :: tr(grid%nnode,0:1)
  real (8), intent(in) :: p(3,grid%nnode,0:1)
  real (8), intent(in) :: var_real(grid%nnode,0:1)
  real (8), intent(out) :: var_ff(grid%nnode,0:1)
  integer :: dir, i, j, cdir, nd
  real (8) :: wt

  var_ff=0D0

  do dir=0, 1
     cdir=1-dir

     do i=1, grid%nnode
        do j=1, 3
           nd=grid%nd(j,tr(i,dir))
           wt=p(j,i,dir)

           ! obtain interpolated values
           var_ff(i,dir)=var_ff(i,dir) + wt*var_real(nd,dir)

        enddo
     enddo
  enddo

end subroutine cnvt_grid_real2ff



subroutine bfollow_test(xp,xm,grid)
  use grid_class
  use sml_module
  implicit none
  type(grid_type) :: grid
  real (kind=8) :: xp(2,grid%nnode), xm(2,grid%nnode)
  integer :: i
  real (kind=8) :: dx(3), r,z,phi,b(3),cost1, cost2, cost3

  if(sml_mype==0) then

     do i=1, grid%nnode
        r=grid%x(1,i)
        z=grid%x(2,i)
        dx(1:2)= xp(:,i)-xm(:,i)
        dx(3)=- 2D0*r *   grid%delta_phi ! r * delta phi
        call bvec_interpol(r,z,phi,b(1),b(2),b(3))

        cost1 = (dx(1)*b(1)+dx(2)*b(2)+dx(3)*b(3))/sqrt( (dx(1)**2+dx(2)**2+dx(3)**2) * (b(1)**2+b(2)**2+b(3)**2) )
        cost2 = (-dx(1)*b(1)-dx(2)*b(2)+dx(3)*b(3))/sqrt( (dx(1)**2+dx(2)**2+dx(3)**2) * (b(1)**2+b(2)**2+b(3)**2) )
        cost3 = (dx(3)*b(3))/sqrt( (dx(3)**2) * (b(1)**2+b(2)**2+b(3)**2) )
        write(1357,1000) i, r,z,1D0-cost1, 1D0-cost2, 1D0-cost3
        write(1358,1000) i, r,z, xp(1,i), xp(2,i), xm(1,i), xm(2,i)
     enddo
     close(1357)
  endif

1000 format (I6,1x,6(e19.13,1x))

end subroutine

subroutine chargee_background(grid,psn,spall)
  use grid_class
  use psn_class
  use sml_module
  use ptl_module
  use smooth_module
  implicit none
  type(grid_type):: grid
  type(psn_type) :: psn
  type(species_type) :: spall(0:ptl_nsp_max)
  integer :: i
  real (kind=8) :: psi
  real (8) :: tmp(grid%npsi00)

  if(sml_deltaf) then
     psn%edensity0=0D0
#ifdef ADIOS
  elseif(sml_restart) then
     ! load files
     call background_edensity0_read(grid,psn)
#endif
  else if(.true.) then
     call chargei(grid,psn,spall(1))
     psn%edensity0(:)=psn%idensity0(:)
     psn%eden00_1d(:)=psn%iden00_1d(:)
     if(sml_flat_electron_density==1) then
        !call smooth_pol0(grid, psn%edensity0, smooth00)  !???? how? call it or not
        call convert_grid_2_001d(grid,psn%edensity0,tmp)
        call convert_001d_2_grid(grid,tmp,psn%edensity0)
     endif
     if(sml_mype==0) then
#ifdef ADIOS
        call background_edensity0_output(grid,psn)
#endif
     endif
  else
     call chargee_background_monte(grid,psn)
  endif

end subroutine chargee_background

! this subroutine obtains gyro averaged ion guiding center density
! it is set to electron background density.
! does not work with canonical distribution
subroutine chargee_background_monte(grid,psn)
  use grid_class
  use psn_class
  use sml_module
  use ptl_module
  use eq_module
  use random_xgc

  implicit none
  type(grid_type):: grid
  type(psn_type) :: psn
  integer :: i, j, larmor, num
  real (8) :: rdim, roffset, zdim, zoffset
  real (8) :: r,z, angle, x(2), x_ring(2), dx_unit(2,sml_nlarmor), dx_unit2(2)
  real (8) :: ti,ti_ev, ni, b, zmax, zdum, v, mu, psi,rho, cosa, sina
  real (8), external ::  b_interpol, gyro_radius, psi_interpol
  integer :: nodes(3), itr, init, count
  real (8) :: p(3)
  logical, parameter :: use_search_tr2=.true.

  rdim=sml_bd_max_r - sml_bd_min_r
  roffset=sml_bd_min_r
  zdim=sml_bd_max_z - sml_bd_min_z
  zoffset=sml_bd_min_z

  do larmor=1, sml_nlarmor
     angle=sml_2pi/real(sml_nlarmor)*real(larmor-1)
     dx_unit(:,larmor)=(/cos(angle),sin(angle)/)
  enddo

  psn%edensity0=0D0
  num=5*sml_monte_num

  do i=1, num
     ! uniform distribution
     r=sqrt(roffset**2 + rdim*ranx()*(2D0*roffset + rdim) )
     z=zdim*ranx() + zoffset
     psi=psi_interpol(r,z,0,0)

     ! get ion temperature and density, B
     ni=eq_ftn(psi,r,z,eq_den)
     ti_ev=eq_ftn(psi,r,z, eq_tempi)
     b=b_interpol(r,z,0D0)

     ! generate mu
     ti=ti_ev*sml_ev2j
     zmax=1.d0 - dexp(-7.d0)
     zdum=zmax*ranx()
     v= sqrt(-2.d0/ptl_mass(1)*dlog(1.d0-zdum)*ti)
     mu=0.5D0*ptl_mass(1)*v**2/b

     x(1)=r
     x(2)=z
     rho=gyro_radius(x,mu)

     angle=sml_2pi*ranx()  ! base angle for 4-point average
     cosa=cos(angle)
     sina=sin(angle)

     ! search 4-point
     do larmor=1, sml_nlarmor
        dx_unit2(1)= dx_unit(1,larmor)*cosa+dx_unit(2,larmor)*sina
        dx_unit2(2)=-dx_unit(1,larmor)*sina+dx_unit(2,larmor)*cosa


        x_ring = x + rho* dx_unit2(:)

        ! find position for
        if (USE_SEARCH_TR2) then
           call search_tr2(grid,x_ring,itr,p)
        else
           if(larmor==1) then
              call search_tr(grid,x_ring,itr,p)
              init=itr
           else
              call search_tr_with_guess(grid,x_ring,init,itr,p,count)
           endif
        endif

        ! add weight onto edensity0
        if(itr>0 .and. itr <=grid%ntriangle) then
           nodes=grid%nd(:,itr)
           do j=1, 3
              psn%edensity0(nodes(j)) = psn%edensity0(nodes(j)) + &
                   p(j)*ni
           enddo
        endif
     enddo

  enddo
  ! all reduce
  call my_mpi_allreduce(psn%edensity0,grid%rtmp1,grid%nnode)

  psn%edensity0=grid%rtmp1*rdim*zdim*sml_2pi_wedge_n*(roffset+0.5D0*rdim)*grid%inv_node_vol &
       /(real(sml_nlarmor*sml_totalpe,8)*real(num,8)*real(sml_nphi_total,8))

  if(sml_mype==0) then
#ifdef ADIOS
     call background_edensity0_output(grid,psn)
#endif
  endif
end subroutine chargee_background_monte

subroutine field_following_pos(x_org,phi,dir,w,delta_phi,x_dest)
  implicit none
  real (kind=8), intent(in) :: x_org(2), w, delta_phi,phi
  real (kind=8), intent(out) :: x_dest(2)
  real (kind=8) :: dphi
  integer :: dir

  if(dir==2) then   ! dir==2 -> larger index. To larger phi plane. dphi >0
     dphi=delta_phi*(1D0-w)
  else              ! dir==1 -> smaller index. To smaller phi plane. dphi <0
     dphi=delta_phi*(w-1D0)
  endif

  call field_following_pos2(x_org,phi,phi+dphi,x_dest)

end subroutine field_following_pos



#ifndef POS2_ORG
subroutine field_following_pos2(x_org,phi_org,phi_dest,x_dest)
  use sml_module, only : sml_nphi_total, sml_cylindrical
  implicit none
  real (kind=8), intent(in) :: x_org(2),phi_org, phi_dest
  real (kind=8), intent(out) :: x_dest(2)
  real (kind=8) :: phi, x(2),dphi
  real (kind=8) :: hh,h6,dx1(2),dx2(2),dx3(2),x_tmp(2)
  integer :: sml_ff_step, order
  integer :: i

#ifdef POS2_ORG_HACK
   x_dest=x_org
   return
#endif

  sml_ff_step=1 ! default
  order=2
  if(sml_nphi_total<=4) then
    order=4
    sml_ff_step=3
  else if(sml_nphi_total<=8) then
    order=4
    sml_ff_step=2
  else if(sml_nphi_total<=16) then
    order=4
  endif
  ! otherwise default value step=1, order=2


  dphi=(phi_dest-phi_org)/real(sml_ff_step)

  ! 0 step
  phi=phi_org
  x=x_org

  do i=1, sml_ff_step
     ! get first derivative
     call derivs(x,phi,dx1)

     if( order==1 ) then ! first order calculation
        x = x + dx1*dphi
     else if( order==2 ) then  ! second order calculation - rk2
        ! obtain mid point
        hh=dphi*0.5D0

        x_tmp = x + dx1*hh

        ! get new derivative
        call derivs(x_tmp,phi+hh,dx2)

        ! advance one step using mid-point derivative
        x = x + dx2*dphi
     else if( order==4 ) then
        ! 4th Order Calculation - rk4
        !
        hh=dphi*0.5D0
        h6=dphi/6D0
        ! derivative 1 (from x)- obtained already
        x_tmp=x + hh*dx1             ! yt=y+hh*dydx      : yt -> x_tmp
        !
        ! derivative 2 (from x_tmp)
        call derivs(x_tmp,phi+hh,dx2)! dyt from yt       : dyt -> dx2
        x_tmp=x + hh*dx2             ! yt=y+hh*dyt       : yt -> x_tmp
        !
        ! derivative 3 (from x_tmp)
        call derivs(x_tmp,phi+hh,dx3)! dym from yt       : dym -> dx3
        x_tmp=x + dphi*dx3           ! yt=y + h*dym      : yt -> x_tmp
        dx3 = dx2 + dx3              ! dym = dyt + dym   : dym -> dx3 , dyt -> dx2
        !
        ! derivative 4 (from x_tmp)
        call derivs(x_tmp,phi+dphi,dx2)    ! dyt from yt       : dyt -> dx2, yt -> x_tmp
        x = x + h6 * (dx1 + dx2 + 2D0*dx3) ! yout = y + h6* (dydx+dyt+2D0*dym)
     else
      print *, 'Wrong Order in field_following_pos2.  order=', order
     endif
     phi=phi+dphi
  enddo
  x_dest=x
contains
  subroutine derivs(x,phi,dx)
    use sml_module
    use eq_module, only: eq_axis_r
    implicit none
    real (8), intent(in)  :: x(2), phi
    real (8), intent(out) :: dx(2)
    real (8) :: b(2), bphi, r,z

    r=min(max(x(1),sml_bd_min_r),sml_bd_max_r)
    z=min(max(x(2),sml_bd_min_z),sml_bd_max_z)


    call bvec_interpol(r,z,phi,b(1),b(2),bphi)
    if (sml_cylindrical) then
      dx = b/bphi * eq_axis_r
    else
      dx = b/bphi * x(1)
    endif

  end subroutine derivs
end subroutine field_following_pos2

#else

subroutine field_following_pos2(x_org,phi0,phi,x_dest)
  use sml_module, only: sml_cylindrical
  use eq_module, only: eq_axis_r
  implicit none
  real (kind=8), intent(in) :: x_org(2), phi,phi0
  real (kind=8), intent(out) :: x_dest(2)
  real (kind=8) :: B(2),Bphi, x_mid(2), dphi
  real (kind=8) :: hh,h6,dx1(2),dx2(2),dx3(2),x_tmp(2)


  dphi=phi-phi0

  call bvec_interpol(x_org(1),x_org(2),phi0,b(1),b(2),bphi)

  if( .false. ) then ! first order calculation
     if (sml_cylindrical) then
       x_dest(:)=x_org(:) + b(:)/bphi*(eq_axis_r*dphi)
     else
       x_dest(:)=x_org(:) + b(:)/bphi*(x_org(1)*dphi)
     endif
  else if( .true. ) then
     ! second order calculation - rk2
     ! obtain mid point
     hh=dphi*0.5D0
     x_mid(:)=x_org(:) + b(:)/bphi*(x_org(1)*hh)

     ! get new derivative
     call bvec_interpol(x_mid(1),x_mid(2),phi0+hh,b(1),b(2),bphi)

     ! advance one step using mid-point derivative
     if (sml_cylindrical) then
       x_dest(:)=x_org(:) + b(:)/bphi*(eq_axis_r*dphi)
     else
       x_dest(:)=x_org(:) + b(:)/bphi*(x_mid(1)*dphi)
     endif
  else
     ! 4th Order Calculation - rk4
     !
     hh=dphi*0.5D0
     h6=dphi/6D0
     ! derivative 1 (from x_org)- obtained already
     if (sml_cylindrical) then
       dx1=b/bphi*eq_axis_r
     else
       dx1=b/bphi*x_org(1)     ! dydx from y       : y-> x_org , dx1 -> dydx
     endif
     x_tmp=x_org + hh*dx1    ! yt=y+hh*dydx      : yt -> x_tmp
     !
     ! derivative 2 (from x_tmp)
     call bvec_interpol(x_tmp(1),x_tmp(2),phi0+hh,b(1),b(2),bphi)
     if (sml_cylindrical) then
       dx2=b/bphi*eq_axis_r
     else
       dx2=b/bphi*x_tmp(1)     ! dyt from yt       : dyt -> dx2
     endif
     x_tmp=x_org + hh*dx2    ! yt=y+hh*dyt       : yt -> x_tmp
     !
     ! derivative 3 (from x_tmp)
     call bvec_interpol(x_tmp(1),x_tmp(2),phi0+hh,b(1),b(2),bphi)
     if (sml_cylindrical) then
       dx3=b/bphi*eq_axis_r
     else
       dx3=b/bphi*x_tmp(1)     ! dym from yt       : dym -> dx3
     endif
     x_tmp=x_org + dphi*dx3     ! yt=y + h*dym      : yt -> x_tmp
     dx3 = dx2 + dx3         ! dym = dyt + dym   : dym -> dx3 , dyt -> dx2
     !
     ! derivative 4 (from x_tmp)
     call bvec_interpol(x_tmp(1),x_tmp(2),phi0+dphi,b(1),b(2),bphi)
     if (sml_cylindrical) then
       dx2=b/bphi*eq_axis_r
     else
       dx2=b/bphi*x_tmp(1)    ! dyt from yt       : dyt -> dx2, yt -> x_tmp
     endif
     x_dest=x_org+h6* ( dx1 + dx2 + 2D0*dx3)     ! yout = y + h6* (dydx+dyt+2D0*dym)  : yout -> x_dest, dydx -> dx1, dyt-> dx2 , dym -> dx3
  endif

  !x_dest(:)=x_org(:)
end subroutine field_following_pos2
#endif


subroutine smooth_r_init2(smooth_r,grid)
  use smooth_module
  use grid_class
  implicit none
  type(smooth_r_type) :: smooth_r
  type(grid_type) :: grid
  real (kind=8) :: vring(smooth_r%n), fring(0:smooth_r%n)
  real (kind=8) :: weight,p(3),dpsi(2),r_norm(2),psi_interpol,x0(2),x(2),tmp
  integer :: i,j,kr,dir,itr,nodes(3)

  if(smooth_r%n<=0) then
     return ! no smoothing
  else if(smooth_r%n==1) then
     vring(1)=1.22474D0
     fring(0)=2D0/3D0
     fring(1)=1D0/6D0
  else if(smooth_r%n==2) then
     vring(1)= 0.958572D0
     vring(2)= 2.02018D0
     fring(0)= 8D0/15D0
     fring(1)= (7D0+2D0*sqrt(10D0))/60
     fring(2)= (7D0-2D0*sqrt(10D0))/60
  else
     !simple linear gauss smoothing
     do i=1, smooth_r%n
        vring(i)= real(2.2*i)/real(smooth_r%n)
        fring(i)= exp( - vring(i)**2 )
     enddo
     fring(0)=1D0
     tmp=1D0 + 2D0 * sum(fring(1:smooth_r%n))
     fring=fring/tmp
     !print *, 'error : invalid smooth_r%n, not implimented yet', smooth_r%n
     !stop
  endif

  do i=1, grid%nnode

     ! 1st point is the original grid point
     call set_value(smooth_r%mat,i,i,fring(0),1)

     ! position of grid points
     x0=grid%x(:,i)

     ! B-field vector
     dpsi(1)=psi_interpol(x0(1),x0(2),1,0)
     dpsi(2)=psi_interpol(x0(1),x0(2),0,1)
     r_norm=dpsi / sqrt(dpsi(1)**2+dpsi(2)**2)

     do kr=1,smooth_r%n
        do dir=-1,1,2
           x=x0+r_norm*real(dir)*smooth_r%d0
           call search_tr2(grid,x,itr,p)
           if(itr>0) then
              nodes=grid%nd(:,itr)

              do j=1,3
                 weight=fring(kr)*p(j)
                 if(weight<0.00001) then
                    call set_value(smooth_r%mat,i,i,weight,1)  ! add to diagonal
                 else
                    call set_value(smooth_r%mat,i,nodes(j),weight,1)
                 endif
              enddo !end of 3-point interpolation loop
           else ! cannot find smoothing grid
              call set_value(smooth_r%mat,i,i,fring(kr),1)
           endif
        enddo !end of direction loop
     enddo ! end of n ring loop
  end do ! end of node loop


end subroutine smooth_r_init2


subroutine smooth_r(in,out,smoothr,grid)
  use grid_class
  use sml_module
  use smooth_module
  use omp_module, only: split_indices
  use perf_monitor
  implicit none
  type(smooth_r_type) :: smoothr
  type(grid_type) :: grid
  real (kind=8) , dimension(grid%nnode) :: in,out
  integer :: i,j,ierr
  integer :: ith, i_beg(sml_nthreads), i_end(sml_nthreads)

  call split_indices(grid%nnode, sml_nthreads, i_beg, i_end)
stop
  if(smoothr%n<=0) then
!$OMP PARALLEL DO &
!$OMP PRIVATE( ITH, I )
     do ith=1,sml_nthreads
        do i=i_beg(ith),i_end(ith)
!pw        do i=1, grid%nnode
           out(i) = in(i)
        enddo
     enddo
     return
  end if

!$OMP PARALLEL DO &
!$OMP PRIVATE( ITH, I, J )
  do ith=1,sml_nthreads
     do i=i_beg(ith),i_end(ith)
!pw     do i=1, grid%nnode
        out(i)=0D0
        do j=1, smoothr%mat%nelement(i)
#ifndef OPTIM_GYRO_AVG_MAT
           out(i)=out(i)+smoothr%mat%value(j,i)*in(smoothr%mat%eindex(j,i))
#else
           out(i)=out(i)+smoothr%mat%value(smoothr%mat%a(i)+j)*in(smoothr%mat%eindex(smoothr%mat%a(i)+j))
#endif
        enddo
     enddo
  enddo

end subroutine smooth_r

!!$subroutine smooth_nearx(in,snx,grid)
!!$  use grid_class
!!$  use smooth_module
!!$  implicit none
!!$  type(smooth_nearx_type) :: snx
!!$  type(grid_type) :: grid
!!$  real (kind=8) , dimension(grid%nnode) :: in
!!$!
!!$  integer :: i,j,ierr, loop
!!$
!!$
!!$  if(snx%n<=0) then
!!$     return
!!$  end if
!!$
!!$  do loop=1, snx%nsmth
!!$     ! get smoothed value
!!$     do i=1, snx%n
!!$        snx%tmp(i)=0D0
!!$        do j=1, snx%mat%nelement(i)
!!$           snx%tmp(i)=snx%tmp(i)+snx%mat%value(j,i)*in(snx%mat%eindex(j,i))
!!$        enddo
!!$     enddo
!!$
!!$     ! save to in with space conversion
!!$     do i=1, snx%n
!!$        in(snx%nodes(i))=snx%tmp(i)
!!$     enddo
!!$  enddo
!!$end subroutine smooth_nearx
!!$
!!$subroutine smooth_nearx_init(snx,grid,d0,dr,dz, nsmth)
!!$  use eq_module
!!$  use sml_module
!!$  use grid_class
!!$  use smooth_module
!!$  implicit none
!!$  type(smooth_nearx_type) :: snx
!!$  type(grid_type) :: grid
!!$  integer :: nsmth
!!$  real (kind=8) :: d0, dr, dz
!!$  !
!!$  integer, parameter :: nlarmor=8
!!$  integer :: ct, larmor, i,j, n
!!$  real (8) :: x0(2), angle, x(2), dx_unit(2,nlarmor)
!!$  integer :: itr, nodes(3)
!!$  real (8) :: p(3), w0, w1, weight
!!$  snx%d0=d0
!!$  snx%dr=dr
!!$  snx%dz=dz
!!$  snx%nsmth=nsmth
!!$
!!$  ! no smoothing
!!$  if(nsmth<=0) then
!!$     snx%n=0
!!$     return
!!$  endif
!!$
!!$  ! get number of nearx nodes
!!$  call search_nearx(.false.,ct)
!!$  snx%n=ct
!!$
!!$  ! exception handle for saftey
!!$  if(ct<=0) return
!!$
!!$  ! memory allocation
!!$  allocate(snx%nodes(snx%n),snx%tmp(snx%n))
!!$
!!$  ! search nearx nodes
!!$  call search_nearx(.true.,ct)
!!$
!!$  !matrix initialize
!!$  call new_mat(snx%mat,snx%n,(8*snx%n+1)*3)
!!$
!!$  ! prepare parameters
!!$  do larmor=1, nlarmor
!!$     angle=sml_2pi/real(sml_nlarmor)*real(larmor-1)
!!$     dx_unit(:,larmor)=(/cos(angle),sin(angle)/)
!!$  enddo
!!$  w0=0.5D0
!!$  w1=(1D0 - w0)/real(nlarmor)
!!$
!!$  do i=1, snx%n
!!$     ! 1st point is the original grid point
!!$     n=snx%nodes(i)
!!$     call set_value(snx%mat,i,n,w0,1)
!!$
!!$     ! position of grid points
!!$     x0=grid%x(:,n)
!!$
!!$     ! larmor loop
!!$     do larmor=1,nlarmor
!!$        x = x0 + dx_unit(:,larmor)*snx%d0
!!$        call search_tr(grid,x,itr,p)
!!$        if(itr>0) then
!!$           nodes=grid%nd(:,itr)
!!$           do j=1,3
!!$              weight=p(j)*w1
!!$              call set_value(snx%mat,i,nodes(j),weight,1)
!!$           enddo !end of 3-point interpolation loop
!!$        else ! cannot find smoothing grid
!!$           call set_value(snx%mat,i,n,weight,1)
!!$        endif
!!$     enddo !end of nlarmor loop
!!$
!!$  end do ! end of node loop
!!$
!!$  contains
!!$    subroutine search_nearx(store,count)
!!$      implicit none
!!$      integer :: count
!!$      logical :: store
!!$      !
!!$      integer :: i
!!$      real (8) :: r,z, rmin, rmax, zmin, zmax
!!$      rmin=eq_x_r - dr
!!$      rmax=eq_x_r + dr
!!$      zmin=eq_x_z - dz
!!$      zmax=eq_x_z + dz
!!$
!!$      count=0
!!$
!!$      do i=1, grid%nnode
!!$         r=grid%x(1,i)
!!$         z=grid%x(2,i)
!!$         if( rmin < r .and. r<rmax .and. zmin < z .and. z < zmax ) then
!!$            count=count+1
!!$            if(store) then
!!$               snx%nodes(count)=i
!!$            endif
!!$         endif
!!$      enddo
!!$
!!$    end subroutine search_nearx
!!$end subroutine smooth_nearx_init


!! Obtain charge density from gyro-center particle
subroutine chargei(grid,psn,sp)
  use grid_class
  use psn_class
  use sml_module
  use eq_module, only : is_rgn12, eq_x_psi, eq_x_z
  use ptl_module
  use smooth_module
  use random_xgc
  use pol_decomp_module
  use omp_module, only: split_indices
  use perf_monitor
#ifdef XGC_COUPLING_CORE_EDGE
!  use coupling_core_edge
  use new_coupling
#endif
  implicit none
  include 'mpif.h'
  type(grid_type):: grid
  type(psn_type) :: psn
  type(species_type) :: sp

  logical, save :: first_flag = .true.
  real*8 :: inv_nphi_total
  integer :: nn,ierr

  call t_startf("CHARGEI_SEARCH_INDEX")
  call chargei_search_index(grid,psn,sp) ! sub-subroutine
  call t_stopf("CHARGEI_SEARCH_INDEX")

#ifndef DELTAF_MODE2
  call t_startf("UPDATE_ION_WEIGHT")
  call update_ion_weight      ! sub-subroutine
  call t_stopf("UPDATE_ION_WEIGHT")
#endif
  call t_startf("CHARGEI_SCATTER")
  call chargei_scatter         ! sub-subroutine
  call t_stopf("CHARGEI_SCATTER")

#ifdef XGC1_EM
  call ijpar_mpisum(grid,psn,sp)
#endif

  call t_startf("CHARGEI_GYRO_AVG")
  call chargei_gyro_average(grid,psn,sp)
#ifdef F0_CHARGE_N0
  if(sml_f0_grid) then
    call chargei_gyro_average_f0(grid,psn,sp)  !need sp?
  endif
#endif
  call t_stopf("CHARGEI_GYRO_AVG")

#ifdef XGC1_EM
  if(sml_gstep>=1)call remove_chargei(psn)
  if(sml_mode_select_on==1)call mode_selection_comb(sml_mode_select_n,grid,psn,psn%idensity(:,1))
#else 
!Test Julien
  !if(sml_mode_select_on==1)call mode_selection_comb(sml_mode_select_n,grid,psn,psn%idensity(:,1))
#endif

#ifdef XGC_COUPLING_CORE_EDGE
  !Seung-Hoe quick patch for BCs: apply charge boundary condition -- additional
  call set_boundary2_values(psn%idensity(:,1),0D0,psn%cbd0_tmp)

  nn=grid%nnode

#ifdef SC17DEMO
  !! jyc: temporary fix for SC17 demo
  if(first_flag)then
    first_flag=.false.
    call cce_initialize()
  endif
#endif
  if(sml_plane_mype==0)then
#ifndef SC17DEMO
    !! jyc: temporary fix for SC17 demo
    if(first_flag)then
      first_flag=.false.
      call cce_initialize()
    endif
#endif
    call cce_send_density(psn%idensity(:,1))
    call cce_receive_density()
    call cce_process_density(psn%idensity(:,1))

!  endif
     ! idensity0 -- sum-up
     call mpi_allreduce(psn%idensity(:,1),psn%idensity0,nn, mpi_real8, mpi_sum, sml_intpl_comm,ierr)
     call t_stopf("CHARGEI_GA_SR_RED")
     inv_nphi_total=1D0/real(sml_nphi_total)
     psn%idensity0=psn%idensity0*inv_nphi_total

  endif

#ifdef F0_CHARGE_N0

  f0_den0_ptl(:,1) = psn%idensity0

#endif

  ! 00 mode
  call t_startf("CHARGEI_GA_RED_BCAST")

  call mpi_bcast(psn%idensity(:,1),nn,mpi_real8,0,sml_plane_comm,ierr)
  call mpi_bcast(psn%idensity0    ,nn,mpi_real8,0,sml_plane_comm,ierr)
  call t_stopf("CHARGEI_GA_RED_BCAST")

  if(.true.) then ! Use idensity0 to get iden00_1d
    call convert_grid_2_001d(grid,psn%idensity0,psn%iden00_1d)
  endif

#endif

contains


  subroutine update_ion_weight
    implicit none
    integer :: i
    integer :: ith, i_beg(sml_nthreads), i_end(sml_nthreads)
    real (8) :: r,z,B,psi,En,new_f0,w2,dw
    real (8), external :: get_f0_ion
    logical, external :: is_nan
    call split_indices(sp%num, sml_nthreads, i_beg, i_end)

!$OMP PARALLEL DO &
!$OMP PRIVATE( ITH, I, NEW_F0, W2, DW )
    do ith=1, sml_nthreads
       call t_startf("UPD_ION_WT_LOOP")
       do i=i_beg(ith), i_end(ith)
          if(sp%ptl(i)%gid>0) then

             new_f0 = get_f0_ion(grid,sp%ptl(i),sp%tr_save(i),sp%p_save(:,i),sp%type)

             w2= 1D0 - new_f0/sp%ptl(i)%ct(pif0)

             dw = w2 - sp%ptl(i)%ph(piw2)
             sp%ptl(i)%ph(piw1)= sp%ptl(i)%ph(piw1) + dw
             sp%ptl(i)%ph(piw2)= w2

             !call restrict_weight(sp%ptl(i)%ph(piw1:piw2), new_f0)

             !if(sml_ipc==2) sp%ptl(i)%ct(pif0)=new_f0

             !NaN test
             !if(is_nan(w2) .or. is_nan(dw)) then
             !   print *, 'NaN found in update_ion_weight : new_f0, w2, dw', new_f0, w2, dw
             !   stop
             !endif
          endif
       enddo
       call t_stopf("UPD_ION_WT_LOOP")
    enddo
  end subroutine update_ion_weight

  subroutine chargei_scatter
    implicit none
    real (8),allocatable :: iden(:,:,:,:), iden00(:,:)
#ifdef XGC1_EM
    real (8), allocatable :: ijpar(:,:,:)
#endif
    real :: inv_delta_phi, inv_drho, dx_unit(2,sml_nlarmor)  ! for fast calculation -- const
    integer :: ith, i_beg(sml_nthreads), i_end(sml_nthreads)
    integer :: nn_beg(sml_nthreads), nn_end(sml_nthreads)
    ! in loop
    integer ::  i,irho, j, node
    real (8) :: phi, wphi(0:1), particle_weight, rho
    real (8) :: rhon, wrho(2), wp, upar
    ! for 1D
    integer :: larmor, ip
    real (8) :: x(2),angle, cosa, sina, dx_unit2(2), x_ring(2), psi_ring, pn
    !
    integer :: jth, iphi
    real (8), external :: psi_interpol
    real (8) :: dpdr, dpdz, dp
    !allocate memory and deallocate it at the end subroutine
    allocate(iden(0:1,0:grid%nrho,grid%nnode, sml_nthreads))  ! indexing is different for cache performance
    allocate(iden00(grid%npsi00,sml_nthreads))
#ifdef XGC1_EM
    allocate(ijpar(0:1, grid%nnode, sml_nthreads))
#endif

    inv_delta_phi=1D0/grid%delta_phi
    inv_drho = 1D0/ grid%drho

    do larmor=1, sml_nlarmor
       angle=sml_2pi/real(sml_nlarmor)*real(larmor-1)
       dx_unit(:,larmor)=(/cos(angle),sin(angle)/)
    enddo


    ! for reduction, divide particles among OpenMP threads
    ! and calculate contributions to nodes for each subset
    ! of particles independently. This will introduce a
    ! roundoff difference in the results for different numbers
    ! of threads.
    call split_indices(sp%num, sml_nthreads, i_beg, i_end)

!$OMP PARALLEL DO &
!$OMP PRIVATE( ITH, I, IRHO, J, NODE, &
!$OMP          PHI, WPHI, PARTICLE_WEIGHT, RHO, RHON, WRHO, WP, &
!$OMP          LARMOR, IP, X, ANGLE, COSA, SINA,DX_UNIT2, X_RING, PSI_RING, PN )
    do ith=1,sml_nthreads
       call t_startf("CHARGEI_SCATTER_LOOP1")

       iden(:,:,:,ith) = 0.0D0
       iden00(:,ith) =0D0
#ifdef XGC1_EM
       ijpar(:,:,ith) = 0D0
#endif
       do i=i_beg(ith),i_end(ith)
          if(sp%ptl(i)%gid<=0) cycle

          if((sp%tr_save(i)) > 0) then ! particle is inside of grid

             x=sp%ptl(i)%ph(1:2)  ! for 1D potential


             ! phi weight
             phi=sp%ptl(i)%ph(3)
             wphi(1) = (phi*inv_delta_phi) - grid%iphi_offset ! larger index weight
             wphi(0) = 1D0 - wphi(1)  ! smaller index weight


             ! particle weight
             if(sml_deltaf) then
                particle_weight=sp%ptl(i)%ph(piw1)*sp%ptl(i)%ct(piw0)
             else
                particle_weight=sp%ptl(i)%ct(piw0) ! for full f simulation only
             endif
             ! rho index
             rho=sp%rhoi(i)
             rhon=min(rho,grid%rhomax)*inv_drho
             irho=min(floor(rhon),grid%nrho-1)
             wrho(2)=rhon - real(irho)
             wrho(1)=1D0-wrho(2)
             do j=1, 3
                node=grid%nd(j,sp%tr_save(i))
                wp=sp%p_save(j,i)
                iden(0:1,irho  ,node,ith) = iden(0:1,irho  ,node,ith) +  wp*particle_weight*wphi(0:1)*wrho(1)
                iden(0:1,irho+1,node,ith) = iden(0:1,irho+1,node,ith) +  wp*particle_weight*wphi(0:1)*wrho(2)
#ifdef XGC1_EM
                ! particle parallel velocity: for parallel ion current calculation
                upar = sp%ptl(i)%ph(4)*ptl_charge(1)*grid%bfield(4,node)/ptl_mass(1)
                ijpar(0:1,node,ith) = ijpar(0:1,node,ith) + wp*particle_weight*wphi(0:1)*ptl_charge(1)*upar
#endif

             enddo

             !1D ---

#ifdef GYRO_AVG_RADIAL
             dpdr=psi_interpol(x(1),x(2),1,0)
             dpdz=psi_interpol(x(1),x(2),0,1)
             dp=sqrt(dpdr**2 + dpdz**2)
             cosa=dpdr/dp
             sina=dpdz/dp
#else
             ! random angle
             angle=sml_2pi*ranx()
             cosa=cos(angle)
             sina=sin(angle)
#endif

             do larmor=1, sml_nlarmor
                ! gyro-ring position
                dx_unit2(1)= dx_unit(1,larmor)*cosa+dx_unit(2,larmor)*sina
                dx_unit2(2)=-dx_unit(1,larmor)*sina+dx_unit(2,larmor)*cosa
                x_ring = x + rho* dx_unit2(:)   ! 1D does not require field following


                psi_ring = psi_interpol(x_ring(1),x_ring(2),0,0)
                pn=(psi_ring-grid%psi00min)/grid%dpsi00
                ip=floor(pn)+1
!                if(0 < ip .and. ip < grid%npsi00 .and. (x_ring(2) > eq_x_z .or. psi_ring > eq_x_psi)) then
                if(0 < ip .and. ip < grid%npsi00 .and. is_rgn12(x_ring(1),x_ring(2),psi_ring) &
                     &.and. (.not. sml_00_xz_up .or. x_ring(2)>eq_x_z) ) then
                   wp=1D0 - ( pn - real(ip-1,8) )  ! recylce of wp -- differenet meaning from above : w of p(1:3) vs w of psi
                   iden00(ip  ,ith)=iden00(ip  ,ith) + particle_weight* wp
                   iden00(ip+1,ith)=iden00(ip+1,ith) + particle_weight*(1D0-wp)
                endif
             enddo

          else

             !eliminate particle
             call remove_particle(sp,i,-1,ith) ! no need due to sheath_calculation

          endif

       enddo

       call t_stopf("CHARGEI_SCATTER_LOOP1")
    enddo  ! end of particle-thread loop

    

    ! open mp sum

    ! combine results
    call split_indices(grid%nnode, sml_nthreads, nn_beg, nn_end)


!$OMP PARALLEL DO &
!$OMP PRIVATE( JTH, ITH, NODE, IRHO, IPHI )
    do jth=1,sml_nthreads
       call t_startf("CHARGEI_SCATTER_LOOP2")

       do ith=2,sml_nthreads
          do node=nn_beg(jth), nn_end(jth)
             do irho=0,grid%nrho
                do iphi=0, 1
                   iden(iphi,irho,node,1)=iden(iphi,irho,node,1)+ iden(iphi,irho,node,ith)
#ifdef XGC1_EM
                   ijpar(iphi,node,1)=ijpar(iphi,node,1) + ijpar(iphi,node,ith)
#endif
                enddo
             enddo
          enddo
          ! compare performance
          ! with     iden(:,:,:,1)=iden(:,:,:,1)+iden(:,:,:,ith)
       enddo

       ! indexing order change
       do irho=0,grid%nrho
          do iphi=0, 1
             do node=nn_beg(jth), nn_end(jth)
                psn%iden_rho_ff(node,iphi,irho) = iden(iphi,irho,node,1)
             enddo
          enddo
       enddo
       call t_stopf("CHARGEI_SCATTER_LOOP2")
    enddo

    ! 1D
    do ith=2, sml_nthreads
       iden00(:,1)=iden00(:,1)+iden00(:,ith)
    enddo
    psn%iden00_1d=iden00(:,1)

    deallocate(iden,iden00)
#ifdef XGC1_EM
    deallocate(ijpar)
#endif
  end subroutine chargei_scatter
end subroutine chargei



subroutine chargei_search_index(grid,psn,sp)
  use grid_class
  use psn_class
  use ptl_module
  use sml_module
  use omp_module, only: split_indices
  use perf_monitor
  implicit none
  type(grid_type):: grid
  type(psn_type) :: psn
  type(species_type) :: sp
  !
  integer :: ith, i
  integer :: i_beg(sml_nthreads), i_end(sml_nthreads)
  real (8) :: phi_mid, x(2), phi, mu, rho, xff(2)
  real (8) :: p(3)
  integer :: itr, ip
  real (kind=8), external :: gyro_radius
  logical, parameter :: USE_SEARCH_TR2 = .true.

  phi_mid=0.5D0*(grid%phimin+grid%phimax)

  call split_indices(sp%num, sml_nthreads, i_beg, i_end)

!$OMP PARALLEL DO &
!$OMP PRIVATE( ITH, I, &
!$OMP          X, PHI, MU, RHO, XFF,  &
!$OMP          P, ITR, IP )
  do ith=1,sml_nthreads
     call t_startf("CHARGEI_SRCHLOOP")
     do i=i_beg(ith),i_end(ith)
        if(sp%ptl(i)%gid>0)  then

           ! get proper toroidal angle index and weight
           x=sp%ptl(i)%ph(1:2)
           phi=sp%ptl(i)%ph(3)
           mu=sp%ptl(i)%ct(pim)
           rho=gyro_radius(x,mu)  !gyro radius
           sp%rhoi(i)=rho

           ! get field following posision at 1/2 angle
           call field_following_pos2(x,phi,phi_mid,xff)

           ! find triangle
           if (USE_SEARCH_TR2) then
              call search_tr2(grid,xff,itr,p)
           else
              call search_tr(grid,xff,itr,p)
           endif

           !remove particle or sheath calculation
           if(itr<0) then
              if(sml_sheath_mode==0 .or. sml_gstep <= 0 ) then
                 call remove_particle(sp,i,-1,ith)
              else
!$omp critical
                 call sheath_calculation(grid,psn,sp,i,1,itr,p,ith)
!$omp end critical
              endif
           endif

           sp%tr_save(i)=itr
           sp%p_save(:,i)=p
           ! compare performance
           !do ip=1,3
           !   sp%p_save(ip,i)=p(ip)
           !enddo

        endif
     enddo
     call t_stopf("CHARGEI_SRCHLOOP")
  enddo
end subroutine chargei_search_index



! get charge density using gyro_average
! sml_pe_per_plane >= sml_nrho * 2   -- important condition
subroutine chargei_gyro_average(grid,psn,sp)
  use grid_class
  use psn_class
  use sml_module
  use f0_module
!  use eq_module, only : eq_x_z
  use ptl_module
  use perf_monitor
  use smooth_module
#ifdef XGC_COUPLING_CORE_EDGE_VARPI2
!  use coupling_core_edge, only : cce_varpi_grid
  use new_coupling, only : cce_varpi_grid
#endif
  implicit none
  type(grid_type):: grid
  type(psn_type) :: psn
  type(species_type) :: sp
  include 'mpif.h'
  !
  integer :: nn, ierr
  integer :: irho, iphi, ipe
  integer :: idest, isendtag, isource, irecvtag
  integer, dimension(MPI_STATUS_SIZE) :: istatus
  real (8) :: tmp(grid%nnode,0:1)
  real (8) :: dum00(grid%npsi00)
  real (8) :: inv_nphi_total
  !
  character (len=30) :: filename
  integer :: i


  nn=grid%nnode

#ifndef F0_CHARGE_N0
  !### check cofficient is correct
  !### optimization -- f0_inode1:f0_inode2 may be enough
  if(sml_f0_grid) then
     do irho=0, grid%nrho
        do iphi=0,1
#ifndef F0_TOR_LINEAR
           psn%iden_rho_ff(:,iphi,irho)=psn%iden_rho_f0(:,irho)*grid%node_vol_ff(:,iphi)&
                & + psn%iden_rho_ff(:,iphi,irho)
#else
           psn%iden_rho_ff(:,iphi,irho)=psn%iden_rho_f0(:,iphi,irho)*grid%node_vol_ff(:,iphi)&
                & + psn%iden_rho_ff(:,iphi,irho)
#endif
        enddo
     enddo
  endif
#endif


#ifdef XGC_COUPLING_CORE_EDGE_VARPI2
  do irho=0, grid%nrho
     do iphi=0,1
       call cce_varpi_grid(psn%iden_rho_ff(:,iphi,irho))
     enddo
  enddo
#endif

  call t_startf("CHARGEI_GA_RED")
  ! send 0-rho to proc 0 of plane_comm
  call mpi_reduce(psn%iden_rho_ff(:,0:1,0),tmp,nn*2, &
       MPI_REAL8,mpi_sum,0,sml_plane_comm,ierr)
  if(sml_plane_mype==0) then
     !psn%iden_rho_ff(:,0:1,0)=tmp
     psn%iden_rho_ff(:,0:1,0)=tmp/grid%node_vol_ff
  endif

  ! (irho-1) * 2 + iphi ==> working plane pe index

  do irho=1,grid%nrho
     do iphi=0,1
        ipe=(irho-1)*2 + iphi
        call mpi_reduce(psn%iden_rho_ff(:,iphi,irho),tmp(:,0),nn, &
             MPI_REAL8,mpi_sum,ipe,sml_plane_comm,ierr)
     enddo
  enddo
  call t_stopf("CHARGEI_GA_RED")

  ! only following proc calculates
  if(sml_plane_mype<grid%nrho*2) then
     ! Do gyro average

     call t_startf("CHARGEI_GA_MAT")
     tmp(:,0)=tmp(:,0)/grid%node_vol_ff(:,mod(sml_plane_mype,2))
     call mat_mult(psn%gyro_avg_mat,tmp(:,0),tmp(:,1))
     !call mat_transpose_mult(psn%gyro_avg_mat,tmp(:,0),tmp(:,1))
     !tmp(:,1)=tmp(:,0) ! debug
     call t_stopf("CHARGEI_GA_MAT")

#ifdef IDEN_DEBUG
     !for debuging
     if(sml_istep==200 .and. sml_ipc==1) then
        do irho=1, grid%nrho
           do iphi=0,1
              if(sml_mype==(irho-1)*2 + iphi) then
                 write(filename,'("debug.",i2.2,".",i1.1,".txt")') irho,iphi
                 open(unit=333,file=filename,status='replace')
                 do i=1, nn
                    write(333,1000) tmp(i,0), tmp(i,1)
                 enddo
              endif
           enddo
        enddo
     endif
1000 format(18(e19.13,1x))
#endif

     call t_startf("CHARGEI_GA_SR")
     if(sml_plane_mype/=0) then
        ! send data to 0 proc
        idest=0
        isendtag=sml_plane_mype
        call mpi_send(tmp(:,1),nn,MPI_REAL8,0,isendtag,sml_plane_comm,ierr)
        !
     else  ! plane_mype==0
        ! from 0 proc to 0 proc -- isource=0 case
        psn%iden_rho_ff(:,0,1)=tmp(:,1)

        ! receive data from other procs
        do isource=1,grid%nrho*2-1
           iphi=mod(isource,2)
           irho=isource/2  + 1
           call mpi_recv(psn%iden_rho_ff(:,iphi,irho),nn,MPI_REAL8,isource,isource,sml_plane_comm, istatus, ierr)
        enddo

        ! sum-up to irho = 0
        ! omp ??
        do irho=1, grid%nrho
           psn%iden_rho_ff(:,:,0)=psn%iden_rho_ff(:,:,0) + psn%iden_rho_ff(:,:,irho)
        enddo
     endif
     call t_stopf("CHARGEI_GA_SR")
  endif


  !prevent total density becomes below 0
  if(sml_f0_grid) then
    do iphi=0, 1
      do i=1, grid%nnode
        if(psn%iden_rho_ff(i,iphi,0) < -f0_den_global(i)) then
            psn%iden_rho_ff(i,iphi,0)=-f0_den_global(i)
        endif
      enddo
    enddo
  endif
  ! coordinate conversion
  if(sml_plane_mype==0) then
     call t_startf("CHARGEI_GA_CNVRT")
     !make it density
     ! omp?
     !psn%iden_rho_ff(:,:,0)=psn%iden_rho_ff(:,:,0)/grid%node_vol_ff

     ! coordinate conversion
     call cnvt_grid_ff2real(grid,psn%ff_hdp_tr,psn%ff_hdp_p,psn%iden_rho_ff(:,:,0),psn%idensity)
     call t_stopf("CHARGEI_GA_CNVRT")
  endif


  ! send 0 plane receive 1 plane
  if(sml_plane_mype==0) then

     call t_startf("CHARGEI_GA_SR_RED")
     idest=mod(sml_intpl_mype-1+sml_intpl_totalpe,sml_intpl_totalpe)
     isendtag=sml_intpl_mype

     isource=mod(sml_intpl_mype+1,sml_intpl_totalpe)
     irecvtag=isource
     call mpi_sendrecv(psn%idensity(:,0),nn,MPI_REAL8,  idest,isendtag, &
                                tmp(:,0),nn,MPI_REAL8,isource,irecvtag, &
                                sml_intpl_comm,istatus,ierr)

     psn%idensity(:,1) = 0.5D0*(psn%idensity(:,1) + tmp(:,0)) ! tmp(:,0) is density received

#ifdef XGC_COUPLING_CORE_EDGE
  endif
#else
     ! idensity0 -- sum-up
     call mpi_allreduce(psn%idensity(:,1),psn%idensity0,nn, mpi_real8, mpi_sum, sml_intpl_comm,ierr)
     call t_stopf("CHARGEI_GA_SR_RED")

     inv_nphi_total=1D0/real(sml_nphi_total)
     psn%idensity0=psn%idensity0*inv_nphi_total

     ! smoothing
!pw     call t_startf("CHARGEI_GA_SMOOTH_R")
!     call smooth_r(psn%idensity0(:),tmp(:,0),smooth_r1,grid)
!pw     call t_stopf("CHARGEI_GA_SMOOTH_R")
!     psn%idensity0=tmp(:,0)
  endif

#ifdef F0_CHARGE_N0

  f0_den0_ptl(:,1) = psn%idensity0

#endif

  ! 00 mode
  call t_startf("CHARGEI_GA_RED_BCAST")

  if(.false.) then !This part does not have iden_rho_ff
  call my_mpi_allreduce(psn%iden00_1d,dum00,grid%npsi00)
  psn%iden00_1d=dum00/(psn%vol00*real(sml_nlarmor))
  endif

  call mpi_bcast(psn%idensity(:,1),nn,mpi_real8,0,sml_plane_comm,ierr)
  call mpi_bcast(psn%idensity0    ,nn,mpi_real8,0,sml_plane_comm,ierr)
  call t_stopf("CHARGEI_GA_RED_BCAST")

  if(.true.) then ! Use idensity0 to get iden00_1d
    call convert_grid_2_001d(grid,psn%idensity0,psn%iden00_1d)
  endif
#endif
end subroutine chargei_gyro_average

#ifdef F0_CHARGE_N0
! get charge density using gyro_average of psn%iden_rho_f0
! sml_pe_per_plane >= sml_nrho * 2   -- important condition
subroutine chargei_gyro_average_f0(grid,psn,sp)
  use grid_class
  use psn_class
  use sml_module
  use f0_module
!  use eq_module, only : eq_x_z
  use ptl_module
  use perf_monitor
  use smooth_module
  implicit none
  type(grid_type):: grid
  type(psn_type) :: psn
  type(species_type) :: sp
  include 'mpif.h'
  !
  integer :: nn, ierr
  integer :: irho, iphi, ipe
  integer :: idest, isendtag, isource, irecvtag
  integer, dimension(MPI_STATUS_SIZE) :: istatus
  real (8) :: tmp(grid%nnode)
  real (8), allocatable :: iden_f0_ff(:,:)
  real (8) :: inv_nphi_total
  !
  character (len=30) :: filename
  integer :: i


  nn=grid%nnode


  allocate(iden_f0_ff(nn,0:1))


  call t_startf("CHARGEI_GA_RED")
  ! send 0-rho to proc 0 of plane_comm
  call mpi_reduce(psn%iden_rho_f0(:,0),tmp,nn,MPI_REAL8,mpi_sum,0,sml_plane_comm,ierr)
  if(sml_plane_mype==0) then
    iden_f0_ff(:,0)=tmp
    iden_f0_ff(:,1)=tmp
  endif

  ! (irho-1) * 2 + iphi ==> working plane pe index
  do irho=1,grid%nrho
     do iphi=0,1
        ipe=(irho-1)*2 + iphi
        call mpi_reduce(psn%iden_rho_f0(:,irho),psn%iden_rho_f0(:,0),nn, &
             MPI_REAL8,mpi_sum,ipe,sml_plane_comm,ierr)   ! psn%iden_rho_f0(:,0) as temp variable
     enddo
  enddo
  call t_stopf("CHARGEI_GA_RED")

  ! only following proc calculates
  if(sml_plane_mype<grid%nrho*2) then
     ! Do gyro average

     call t_startf("CHARGEI_GA_MAT")
     call mat_mult(psn%gyro_avg_mat,psn%iden_rho_f0(:,0),tmp)
     call t_stopf("CHARGEI_GA_MAT")

     call t_startf("CHARGEI_GA_SR")
     if(sml_plane_mype/=0) then
        ! send data to 0 proc
        idest=0
        isendtag=sml_plane_mype
        call mpi_send(tmp,nn,MPI_REAL8,0,isendtag,sml_plane_comm,ierr)
        !
     else  ! plane_mype==0
        ! from 0 proc to 0 proc -- isource=0 case
        iden_f0_ff(:,0)=iden_f0_ff(:,0) + tmp

        ! receive data from other procs
        do isource=1,grid%nrho*2-1
           iphi=mod(isource,2)
           irho=isource/2  + 1
           call mpi_recv(tmp,nn,MPI_REAL8,isource,isource,sml_plane_comm, istatus, ierr)

           iden_f0_ff(:,iphi)=iden_f0_ff(:,iphi) + tmp
        enddo

     endif
     call t_stopf("CHARGEI_GA_SR")
  endif

  ! coordinate conversion
  if(sml_plane_mype==0) then
     call t_startf("CHARGEI_GA_CNVRT")
     !make it density
     ! omp?
     !psn%iden_rho_ff(:,:,0)=psn%iden_rho_ff(:,:,0)/grid%node_vol_ff

     ! coordinate conversion
     call cnvt_grid_ff2real(grid,psn%ff_hdp_tr,psn%ff_hdp_p,iden_f0_ff,psn%idensity_f0)
     call t_stopf("CHARGEI_GA_CNVRT")
  endif


  ! send 0 plane receive 1 plane
  if(sml_plane_mype==0) then

     call t_startf("CHARGEI_GA_SR_RED")
     idest=mod(sml_intpl_mype-1+sml_intpl_totalpe,sml_intpl_totalpe)
     isendtag=sml_intpl_mype

     isource=mod(sml_intpl_mype+1,sml_intpl_totalpe)
     irecvtag=isource
     call mpi_sendrecv(psn%idensity_f0(:,0),nn,MPI_REAL8,  idest,isendtag, &
                                tmp,nn,MPI_REAL8,isource,irecvtag, &
                                sml_intpl_comm,istatus,ierr)

     psn%idensity_f0(:,1) = 0.5D0*(psn%idensity_f0(:,1) + tmp) ! tmp is density received

     ! idensity_f0 n=0 -- sum-up
     call mpi_allreduce(psn%idensity_f0(:,1),tmp,nn, mpi_real8, mpi_sum, sml_intpl_comm,ierr)
     call t_stopf("CHARGEI_GA_SR_RED")

     inv_nphi_total=1D0/real(sml_nphi_total)
     tmp=tmp*inv_nphi_total

     psn%idensity_f0(:,1) = psn%idensity_f0(:,1) - tmp + f0_density_n0_add(:,1)  ! last index 1 is for ion
  endif

  call t_startf("CHARGEI_GA_RED_BCAST")
  call mpi_bcast(psn%idensity_f0(:,1),nn,mpi_real8,0,sml_plane_comm,ierr)
  call t_stopf("CHARGEI_GA_RED_BCAST")

  psn%idensity(:,1) = psn%idensity(:,1) + psn%idensity_f0(:,1)
  psn%idensity0     = psn%idensity0     + f0_density_n0_add(:,1)

  deallocate(iden_f0_ff)

end subroutine chargei_gyro_average_f0
#endif

real (8) function get_f0_ion(grid,ptli,itr,p,sp_type)
  use grid_class
  use eq_module
  use ptl_module
  use f0_module
  use sml_module
  implicit none
  type(grid_type) :: grid
  type(ptl_type) :: ptli
  integer, intent(in) :: itr, sp_type
  real(8), intent(in) :: p(3)
  real (8) :: r,z,B, psi, En, den, temp_ev
  real (8) :: mu_n, vp_n, phi
  real (8) :: f0_g
  real (8), external :: b_interpol, psi_interpol
  real (8) :: enp
  logical :: err
  logical, external :: is_nan
  ! 1. f0_analytic - Maxwellian
  r=ptli%ph(pir)
  z=ptli%ph(piz)
  phi=ptli%ph(pip)
  B=b_interpol(r,z,phi)
  psi=psi_interpol(r,z,0,0)

  en=ptli%ct(pim)*B + 0.5D0*(ptl_charge(sp_type)*ptli%ph(pirho)*B)**2/ptl_mass(sp_type)

  den=eq_ftn(psi,r,z,eq_den)
  temp_ev=eq_ftn(psi,r,z,eq_tempi)
  !den=f0_iden(psi,r,z)
  !temp_ev=f0_tempi(psi,r,z)

!rh #ifndef MU_LINEAR
!rh   enp=ptli%ct(pim)*eq_axis_b*sml_j2ev/temp_ev
!rh   ! The factor sqrt(enp) can lead to numerical problems in the weight update for particles with mu-->0
!rh   !get_f0_ion = den/temp_ev*exp(-en*sml_j2ev/temp_ev)*2D0*sqrt(enp)
!rh   ! Use this instead -->
!rh   get_f0_ion = den*sqrt(1D0/temp_ev**3)*exp(-en*sml_j2ev/temp_ev)
!rh #else
  ! f is defined in sqrt(m/(2 pi e^3)) B dmu dvp space
  get_f0_ion = den*sqrt(1D0/temp_ev**3)*exp(-en*sml_j2ev/temp_ev)
!rh #endif
  ! 2. f0_grid
  if(sml_f0_grid) then
#ifndef V_PERP
     mu_n=ptli%ct(pim)*(2.*eq_axis_b)/(temp_ev*sml_ev2j)  ! mu_n here is mu2B_n in f0_get_f0g
#else
     ! mu_n becomes (v_perp/v_th)^2 here --->
     mu_n=ptli%ct(pim)*(2.*b/(temp_ev*sml_ev2j))
#endif
     vp_n=ptl_c_m(sp_type)*ptli%ph(pirho)*B/sqrt(temp_ev*sml_ev2j/ptl_mass(sp_type))

     call f0_get_f0g(grid,sp_type,itr,p,phi,mu_n,vp_n,f0_g,err)
     if(err) then
        call for_debugging
     endif
     get_f0_ion=get_f0_ion + f0_g
  endif


  if(is_nan(get_f0_ion)) then
     print *, 'NaN found in get_f0_ion (en, psi, den, temp_ev,mu_n, vp_n, f0_g)', en,psi, den, temp_ev, mu_n, vp_n, f0_g
     print *, 'sml_mype,f0_inode1,f0_inode2',sml_mype, f0_inode1, f0_inode2
     stop
  endif

contains
subroutine for_debugging
  implicit none
  real (8) :: x(2), xff(2), phi_mid, p2(3)
  integer :: itr2
  logical, parameter :: use_search_tr2=.true.

  x(1)=r
  x(2)=z

  phi_mid=(floor(phi/grid%delta_phi) + 0.5D0) * grid%delta_phi

  ! get field following posision at 1/2 angle
  call field_following_pos2(x,phi,phi_mid,xff)

  ! find triangle again
  if (USE_SEARCH_TR2) then
      call search_tr2(grid,xff,itr2,p2)
  else
      call search_tr(grid,xff,itr2,p2)
  endif

  print *, 'r,z,phi,sp_type,itr,p=',r,z,phi,sp_type,itr,p
  print *, 'r_ff,z_ff,phi_mid=', xff(1), xff(2), phi_mid
  print *, 'itr2, p2',itr2, p2

end subroutine
end function get_f0_ion

real (8) function get_f0_elec(grid,psn,ptli,itr,p)
  use grid_class
  use psn_class
  use eq_module
  use ptl_module
  use f0_module
  use sml_module
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  type(ptl_type) :: ptli
  integer, intent(in) :: itr
  real(8), intent(in) :: p(3)
  integer, parameter :: sp_type=0
  real (8) :: r,z,B, psi, En, den, temp_ev
  real (8) :: mu_n, vp_n, phi
  real (8) :: f0_g
  real (8), external :: b_interpol, psi_interpol
  real (8) :: dpot
  real (8) :: tmp, enp
  logical :: err
  logical, external :: is_nan
  ! 1. f0_analytic - Maxwellian
  r=ptli%ph(pir)
  z=ptli%ph(piz)
  phi=ptli%ph(pip)
  B=b_interpol(r,z,phi)
  psi=psi_interpol(r,z,0,0)
  en=ptli%ct(pim)*B + 0.5D0*(ptl_charge(sp_type)*ptli%ph(pirho)*B)**2/ptl_mass(sp_type)

  den=eq_ftn(psi,r,z,eq_den)
  temp_ev=eq_ftn(psi,r,z,eq_tempe)
  !den=f0_eden(psi,r,z)
  !temp_ev=f0_tempe(psi,r,z)

  call get_dpot(dpot)


!  get_f0_elec = den*sqrt(1D0/temp_ev**3)*exp(-en*sml_j2ev/temp_ev)*exp(dpot/temp_ev)
  !maybe (1D0+dpot/temp_ev) is better than exp(dpot/temp_ev) to be consistent with poisson equation --> but f0_a can be minus when dpot<-temp_ev


!  get_f0_elec = den*sqrt(1D0/temp_ev**3)*exp(-en*sml_j2ev/temp_ev)*(1D0 + dpot / temp_ev)
! did not work well



  tmp=dpot/temp_ev
  tmp=min(tmp, 1D0)
!rh #ifndef MU_LINEAR
!rh   enp=ptli%ct(pim)*eq_axis_b*sml_j2ev/temp_ev
!rh   ! The factor sqrt(enp) can lead to numerical problems in the weight update for particles with mu-->0
!rh   !get_f0_elec = den/temp_ev*exp(-en*sml_j2ev/temp_ev)*exp(tmp)*2D0*sqrt(enp)
!rh   ! Use this instead -->
!rh   get_f0_elec = den*sqrt(1D0/temp_ev**3)*exp(-en*sml_j2ev/temp_ev)*exp(tmp)
!rh #else
  get_f0_elec = den*sqrt(1D0/temp_ev**3)*exp(-en*sml_j2ev/temp_ev)*exp(tmp)
!rh #endif



  ! 2. f0_grid

  if(sml_f0_grid) then
#ifndef V_PERP
     mu_n=ptli%ct(pim)*(2.*eq_axis_b)/(temp_ev*sml_ev2j)
#else
     ! mu_n becomes (v_perp/v_th)^2 here --->
     mu_n=ptli%ct(pim)*(2*b/(temp_ev*sml_ev2j))
#endif
     vp_n=ptl_c_m(sp_type)*ptli%ph(pirho)*B/sqrt(temp_ev*sml_ev2j/ptl_mass(sp_type))

     call f0_get_f0g(grid,sp_type,itr,p,phi,mu_n,vp_n,f0_g,err)
     if(err) then
        print *, 'r,z,phi,sp_type,itr,p=',r,z,phi,sp_type,itr,p
     endif
     get_f0_elec=get_f0_elec + f0_g
  endif

  if(is_nan(get_f0_elec)) then
    print *, 'NaN found in get_f0_elec (en, psi, den, temp_ev,mu_n, vp_n, f0_g)', en,psi, den, temp_ev, mu_n, vp_n, f0_g
    print *, 'itr, dpot, sml_mype,f0_inode1,f0_inode2',itr, dpot, sml_mype, f0_inode1, f0_inode2
    stop
  endif

contains
  subroutine get_dpot(dpot)
    implicit none
    real (8) , intent(out) :: dpot
    integer :: ip, node
    real (8) :: wp, wphi(0:1)


    wphi(1)=(phi/grid%delta_phi)  - grid%iphi_offset
    wphi(0)=1D0 - wphi(1)

    dpot=0D0

    if(sml_turb_efield) then
       if(itr>0) then
          do ip=1, 3
             node=grid%nd(ip,itr)
             wp=p(ip)

             ! for electron -- rho=0 case, optimized.
             dpot = dpot + wp*wphi(0)*psn%dpot_ff(node,0)
             dpot = dpot + wp*wphi(1)*psn%dpot_ff(node,1)

          enddo
       endif
    endif
  end subroutine get_dpot
end function get_f0_elec


subroutine sheath_calculation(grid,psn,sp,iptl,type,itrout,pout,ith)
  use grid_class
  use psn_class
  use sml_module
  use ptl_module
  use neu_module, only : neu_weight_sum_lost
  use diag_module
  implicit none
  type(grid_type):: grid
  type(psn_type) :: psn
  type(species_type) :: sp
  integer, intent(in) :: type, iptl
  integer :: itrout
  real (kind=8) :: pout(3)
  integer, intent(in) :: ith
  !
  integer :: i,l, node,itr
  real (kind=8) :: p(3),psave(3)
  integer, parameter :: rgn_wall=100
  real (kind=8), parameter :: minus_val=-1D50
  real (kind=8) :: rho,b,en_para, x(2),phi, phi_mid, xff(2)
  real (8) :: time_now, xn(2), dist_sqr, dist_min
  integer :: node_min
  real (kind=8) , external :: b_interpol
  logical, parameter :: USE_SEARCH_TR2 = .true.
  real (kind=8) :: new_phase(ptl_nphase)
  integer :: widx
  real (8) :: w1_change, en_perp
  ! find nearest wall point

  new_phase=sp%phase0(:,iptl)
  x = new_phase(1:2)
  phi=new_phase(3)
  phi_mid=(floor(phi/grid%delta_phi) + 0.5D0) * grid%delta_phi

  ! get field following posision at 1/2 angle
  call field_following_pos2(x,phi,phi_mid,xff)

  ! find position of previous time step
  !if (USE_SEARCH_TR2) then
  call search_tr2(grid,xff,itr,p)
  !else
  !   call search_tr(grid,xff,itr,p)
  !endif
  psave=p

  ! if old position is also outside of grid --> remove particle and return
  if(itr<0) then
     call remove_particle(sp,iptl,-1,ith)
     itrout=itr
     pout=p
     return
  endif

  ! search three nodes of the triangle and check if it is wall nodes
  do i=1, 3
     l=maxloc(p,1)
     node = grid%nd(l,itr)

     if(grid%rgn(node)==rgn_wall) then
        exit
     else
        p(l)=minus_val
     endif
  enddo

  !if the triangle does not have a wall node
  ! find nearest one using primitive scan
  if(grid%rgn(node)/=rgn_wall) then
     dist_min = 1D99
     do i=1, psn%nwall ! for all wall nodes
        ! check distance
        xn=grid%x(:,psn%wall_nodes(i))
        dist_sqr=(xn(1)-xff(1))**2 + (xn(2)-xff(2))**2
        ! check minimum
        if(dist_min > dist_sqr) then
           dist_min=dist_sqr
           node_min = i
        endif
     enddo
     node=psn%wall_nodes(node_min)
  endif


  !
  !check potential and energy

  rho=new_phase(4)
  b=b_interpol(x(1),x(2),0D0)

  en_para=ptl_c2_2m(sp%type) * (rho*b)**2

  ! reflection
  new_phase(pirho)=-new_phase(pirho)


  widx=psn%node_to_wall(node)
  if(en_para < psn%sheath_pot(widx) *sml_ev2j .and. sp%type==0) then
     ! do nothing -- just reflection due to sheath potential
  else
     w1_change = 1D0 - new_phase(piw2) + new_phase(piw1)

     ! sum-up to wall node
     if(sml_ipc==2 .and. (sml_epc==1 .or. sp%type/=0 ) ) then
       psn%sheath_lost(widx,ith)=psn%sheath_lost(widx,ith) + (w1_change)*sp%charge*sp%ptl(iptl)%ct(piw0)
       if(sp%type/=0) then
          psn%sheath_ilost(widx,ith)=psn%sheath_ilost(widx,ith) + (w1_change)*sp%charge*sp%ptl(iptl)%ct(piw0)
       endif
     endif

     ! --- reflect with weight change
     new_phase(piw1)=-1D0+new_phase(piw2)
     ! w2 does not change

     ! for neutral: 
     ! for ion electron simulations .and. .not. sml_neutral_use_ion_loss --> count electrons only
     ! for ion only simulations .or. sml_neutral_use_ion_loss --> count ions
     ! Two separate if statements for better code reading --> sinlge if would be better for code performance
     if (sml_neutral) then
        if(sml_ipc==2) then
          if(sp%type/=0 .and. (sml_neutral_use_ion_loss .or. .not. sml_electron_on)) then  !count ion
              neu_weight_sum_lost(ith) = neu_weight_sum_lost(ith) + (w1_change)*sp%ptl(iptl)%ct(piw0)
          elseif(sp%type==0 .and. .not. sml_neutral_use_ion_loss .and. sml_epc==1) then    ! count electron
              neu_weight_sum_lost(ith) = neu_weight_sum_lost(ith) + (w1_change)*sp%ptl(iptl)%ct(piw0)
          endif
        endif
     endif

     ! for heat load diagnosie
     if(sml_ipc==2 .and. diag_heat_on .and. (sp%type/=0 .or. sml_epc==1)) then

        en_perp = sp%ptl(iptl)%ct(pim)*b
        ! arguments:
        ! 1. weight_change, 2. potential, 3. en_para, 4. en_perp, 5. ct, 6. phase(old), 7. phase(new), 8. stype
        ! new_phase is old position
        call diag_heat_port(w1_change,psn%sheath_pot(widx),en_para, en_perp, sp%ptl(iptl)%ct, &
             & new_phase, sp%ptl(iptl)%ph,grid%delta_phi, type, ith)
     endif
  endif

  sp%ptl(iptl)%ph(:)=new_phase

  if(sp%ptl(iptl)%gid>0) then
     itrout=itr
     pout=psave
  endif

end subroutine sheath_calculation

subroutine sheath_adjust(grid,psn)
  use sml_module
  use grid_class
  use psn_class
  use eq_module
  use diag_module
  implicit none
  include 'mpif.h'
  type(grid_type) :: grid
  type(psn_type) :: psn
  integer :: i, ierr
  real(8) :: lost_sum(psn%nwall)
  real(8) :: factor

  !openmp reduction
  do i=2, sml_nthreads
    psn%sheath_ilost(:,1)=psn%sheath_ilost(:,1)+psn%sheath_ilost(:,i)
    psn%sheath_lost(:,1) =psn%sheath_lost(:,1) +psn%sheath_lost(:,i)
  enddo


  !mpi reduce of sheath_lost & sheath_ilost -- assuming n=0 sheath potential
  lost_sum=0D0 !safety
  call my_mpi_reduce(psn%sheath_ilost(:,1),lost_sum,psn%nwall)
  psn%sheath_ilost(:,1)=lost_sum
  ! sheath_ilost is for diagnosis only

  lost_sum=0D0
  call my_mpi_reduce(psn%sheath_lost(:,1),lost_sum,psn%nwall)


  ! T / (n * A * l)
  ! T= 10 eV  , n = eq_den_v2, A= 2pi R rho_i, l = rho_i
  !  sml_e_charge * T / (eq_den_v2 * sml_2pi * eq_axis_r * 0.005 * 0.005)
  factor=sml_sheath_adjust_factor* 10. / ( eq_den_v2 * sml_2pi * eq_axis_r * 0.005 * 0.005 ) ! charge_e is already mutliplied. check sign###
  !rh factor seems too small, try to divide it by sml_e_charge
  factor=factor/sml_e_charge

  if(sml_mype==0) then
     print *,maxval(abs(lost_sum)),factor
     do i=1, psn%nwall
        psn%sheath_pot(i) = psn%sheath_pot(i) - factor * lost_sum(i)
        !rh          it has to be "-" here ---^
     enddo
  endif

  !diagnosis
  psn%sheath_lost(:,1)=lost_sum
  if(sml_mype==0 .and. mod(sml_gstep,diag_1d_period)==0 ) call diag_sheath(grid,psn)


  !rh Reset psn%sheath_lost
  psn%sheath_lost=0D0
  psn%sheath_ilost=0D0
  call mpi_bcast(psn%sheath_pot, psn%nwall, mpi_real8, 0 , sml_comm,ierr)


end subroutine sheath_adjust

!call after update_f0_sp
subroutine chargei_f0(grid,psn)
  use sml_module
  use grid_class
  use psn_class
  use f0_module
  use eq_module
  use ptl_module
  use omp_module
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn

  integer,parameter :: isp=1
  real (8),allocatable, save :: iden(:,:,:,:)   ! one dummy dim without F0_TOR_LINEAR
  integer :: ith,jth, i_beg(sml_nthreads), i_end(sml_nthreads)
  integer ::  nd, imu, irho
  real (8) :: inv_drho,mu, rho, smu_n, particle_weight
  real (8) :: rhon, wrho(2), wp
  real (8) :: tmp

#ifndef F0_TOR_LINEAR
  allocate(iden(1,0:grid%nrho,f0_inode1:f0_inode2, sml_nthreads))  ! indexing is different for cache performance
#else
  allocate(iden(0:1,0:grid%nrho,f0_inode1:f0_inode2, sml_nthreads))  ! indexing is different for
#endif

  inv_drho = 1D0/ grid%drho
  call split_indices(f0_inode2-f0_inode1+1, sml_nthreads, i_beg, i_end)

  iden=0D0

!$OMP PARALLEL DO &
!$OMP PRIVATE( ith, nd, imu, &
!$OMP          PARTICLE_WEIGHT, smu_n, mu, RHO, RHON, irho, WRHO)
  do ith=1,sml_nthreads
     do nd=i_beg(ith)+f0_inode1-1,i_end(ith)+f0_inode1-1
        do imu=f0_imu1, f0_imu2

           smu_n=real(imu)*f0_dsmu     ! sqrt(mu) normalized
           rho=smu_n * sqrt(ptl_mass(1)*psn%tempi_ev(nd)*sml_ev2j)/(sml_e_charge*grid%bfield(4,nd))

           rhon=min(rho,grid%rhomax)*inv_drho
           irho=min(floor(rhon),grid%nrho-1)
           wrho(2)=rhon - real(irho)
           wrho(1)=1D0-wrho(2)

#ifndef F0_TOR_LINEAR
           particle_weight=sum(f0_f0g(:,nd,imu,isp))*f0_grid_vol(nd,isp)
           iden(1,irho  ,nd,ith) = iden(1,irho  ,nd,ith) +  particle_weight*wrho(1)
           iden(1,irho+1,nd,ith) = iden(1,irho+1,nd,ith) +  particle_weight*wrho(2)
#else
           tmp=sum(f0_f0g(:,nd,imu,0,isp))*f0_grid_vol_vonly(nd,isp)
           particle_weight=tmp*grid%node_vol_ff(nd,0)
           iden(0,irho  ,nd,ith) = iden(0,irho  ,nd,ith) +  particle_weight*wrho(1)
           iden(0,irho+1,nd,ith) = iden(0,irho+1,nd,ith) +  particle_weight*wrho(2)

           tmp=sum(f0_f0g(:,nd,imu,1,isp))*f0_grid_vol_vonly(nd,isp)
           particle_weight=tmp*grid%node_vol_ff(nd,1)
           iden(1,irho  ,nd,ith) = iden(1,irho  ,nd,ith) +  particle_weight*wrho(1)
           iden(1,irho+1,nd,ith) = iden(1,irho+1,nd,ith) +  particle_weight*wrho(2)
#endif
        enddo
     enddo
  enddo
  do ith=2,sml_nthreads
     iden(:,:,:,1)=iden(:,:,:,1)+iden(:,:,:,ith)
  enddo

  !rh We have to set psn%iden_rho_f0 to zero, because f0_inode1 and f0_inode2
  !rh may change during the simulation
  psn%iden_rho_f0=0D0

  ! make it density with 0.5 factor
  ! indexing order changed
  do irho=0,grid%nrho
!$OMP PARALLEL DO &
!$OMP PRIVATE( nd)
     do nd=f0_inode1, f0_inode2
#ifndef F0_TOR_LINEAR
        psn%iden_rho_f0(nd,irho) = 0.5D0*iden(1,irho,nd,1)/grid%node_vol_nearest(nd)
#else
        psn%iden_rho_f0(nd,0:1,irho) = iden(0:1,irho,nd,1)/grid%node_vol_ff(nd,0:1)
#endif
     enddo
  enddo

  deallocate(iden)

end subroutine chargei_f0

!call after update_f0_sp
subroutine chargee_f0(grid,psn)
  use sml_module
  use grid_class
  use psn_class
  use f0_module
  use eq_module
  use ptl_module
  use omp_module
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn

  integer,parameter :: isp=0
  real (8),allocatable, save :: eden(:,:,:)
  integer :: ith,jth, i_beg(sml_nthreads), i_end(sml_nthreads)
  integer ::  nd, imu, iphi
  real (8) :: mu, particle_weight
  real (8) :: wp

#ifndef F0_TOR_LINEAR
  allocate(eden(1,f0_inode1:f0_inode2, sml_nthreads))  ! indexing is different for cache performance
#else
  allocate(eden(0:1,f0_inode1:f0_inode2, sml_nthreads))  ! indexing is different for cache performance
#endif

  call split_indices(f0_inode2-f0_inode1+1, sml_nthreads, i_beg, i_end)
  eden=0D0

!$OMP PARALLEL DO &
!$OMP PRIVATE( ith, nd, imu, &
!$OMP          PARTICLE_WEIGHT, mu )
  do ith=1,sml_nthreads
     do nd=i_beg(ith)+f0_inode1-1,i_end(ith)+f0_inode1-1
        particle_weight=0D0

#ifndef F0_TOR_LINEAR
        do imu=f0_imu1, f0_imu2
           particle_weight=particle_weight + sum(f0_f0g(:,nd,imu,isp))
        enddo  ! make sum(sum(f0_f0g )) ???
        particle_weight=particle_weight*f0_grid_vol(nd,isp)

        eden(1,nd,ith) = eden(1,nd,ith) +  particle_weight
#else
        do iphi=0,1
          do imu=f0_imu1, f0_imu2
           particle_weight=particle_weight + sum(f0_f0g(:,nd,imu,iphi,isp))
          enddo  ! make sum(sum(f0_f0g )) ???
          particle_weight=particle_weight*f0_grid_vol_vonly(nd,isp)*grid%node_vol_ff(nd,iphi)

          eden(iphi,nd,ith) = eden(iphi,nd,ith) +  particle_weight
        enddo


#endif

     enddo
  enddo

  do ith=2, sml_nthreads
     eden(:,:,1)=eden(:,:,1)+eden(:,:,ith)
  enddo

  !rh We have to set psn%eden_f0 to zero, because f0_inode1 and f0_inode2
  !rh may change during the simulation
  psn%eden_f0=0D0

  ! make it density with 0.5 factor
!$OMP PARALLEL DO &
!$OMP PRIVATE( nd)
  do nd=f0_inode1, f0_inode2
#ifndef F0_TOR_LINEAR
     psn%eden_f0(nd)=0.5D0*eden(1,nd,1)/grid%node_vol_nearest(nd)
#else
     psn%eden_f0(nd,0:1)=eden(0:1,nd,1)/grid%node_vol_ff(nd,0:1)
#endif
  enddo


  !mpi all reduce will be done together with particle density
  deallocate(eden)

end subroutine chargee_f0


!******************************************************************************************************
!********  Funtions for debugging ********************************************************************
!******************************************************************************************************
subroutine enforce_modeled_charge(grid,psn)
  use psn_class
  use grid_class
  use sml_module
  implicit none
  type(grid_type):: grid
  type(psn_type) :: psn



  psn%idensity=1D15
  psn%idensity0=1D15


end subroutine enforce_modeled_charge

