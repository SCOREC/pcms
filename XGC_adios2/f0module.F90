module f0_module
  use eq_module
  use mat_class
  use pol_decomp_module

  !!PRIVATE MEMBER FUNCTIONS:
  private ceil2    
  private pair     
  private partition_intersect

  integer :: f0_nmu  !!=31  now set in input file/setup
  integer :: f0_nvp  !!=15
#ifndef F0_TOR_LINEAR
  real (8), allocatable :: f0_f0g(:,:,:,:)
#else
  real (8), allocatable :: f0_f0g(:,:,:,:,:)
  real (8), allocatable :: f0_df0g3(:,:,:,:,:)
#endif
  real (8), allocatable :: f0_f(:,:,:,:)
  real (8), allocatable :: f0_df0g(:,:,:,:)  !! stores the result of the collision operator
  real (8), allocatable :: f0_df0g2(:,:,:,:)  !! for giving df back to particles
#ifndef F_USE_MARKER_DEN2
  real (8), allocatable :: f0_n(:,:,:,:) !! Normalization for mesh --> particle interpolation
#else
  ! F0_TOR_LINEAR should be defined together
  real (8), allocatable :: f0_n(:,:,:,:,:) !! Normalization for mesh --> particle interpolation
#endif
  
#ifdef DIAG_NOISE
  real (8), allocatable :: f0_diag_n(:,:,:,:), f0_diag_w0(:,:,:,:), f0_diag_w0s(:,:,:,:), f0_diag_w(:,:,:,:), f0_diag_ws(:,:,:,:)
#endif



  real (8) :: f0_smu_max       !
  real (8) :: f0_vp_max  !!=3  ! normalized



  real (8) :: f0_dsmu, f0_dvp
!  real (8) :: f0_mu_norm, f0_vp_norm
!moved  integer :: f0_inode1,f0_inode2
  integer :: f0_imu1,f0_imu2
  real (8), allocatable :: f0_grid_vol(:,:), f0_grid_vol_vonly(:,:)
  real (8), allocatable :: f0_B_B0(:), f0_n_Ta(:,:), f0_den(:), f0_t_ev(:,:)
  real (8), allocatable :: f0_den_global(:)

  logical :: f0_col_change_weight !! if true, transfer of collisional f0_f0g to the particle weights
  logical :: f0_f_correction

  real (8), parameter :: f0_mu0_factor=3D0  !! Set value of lowest mu in grid --> 1/f0_mu0_factor

  type(eq_ftn_type) :: f0_edge_envelope

  real (8), allocatable :: f0_den0_ptl(:,:)
  real (8), allocatable :: f0_density_n0_add(:,:)
contains
  subroutine f0_mem_allocation
    use ptl_module
    implicit none
    integer :: alloc_stat

    if (allocated(f0_f0g)) deallocate(f0_f0g)
#ifndef F0_TOR_LINEAR
    allocate(f0_f0g(-f0_nvp:f0_nvp, &
         f0_inode1:f0_inode2,&
         f0_imu1:f0_imu2,&
         ptl_isp:ptl_nsp),stat=alloc_stat)
#else
    allocate(f0_f0g(-f0_nvp:f0_nvp, &
         f0_inode1:f0_inode2,&
         f0_imu1:f0_imu2,0:1,&
         ptl_isp:ptl_nsp),stat=alloc_stat)
#endif
     call assert(alloc_stat .eq. 0, &
                'alloc(f0_f0g) return istat=',alloc_stat)

    if (allocated(f0_f)) deallocate(f0_f)
    allocate(f0_f(-f0_nvp:f0_nvp, &
         f0_inode1:f0_inode2,&
         f0_imu1:f0_imu2,&
         ptl_isp:ptl_nsp),stat=alloc_stat)
     call assert(alloc_stat .eq. 0, &
                'alloc(f0_f) return istat=',alloc_stat)

    if (allocated(f0_df0g)) deallocate(f0_df0g)
    allocate(f0_df0g(-f0_nvp:f0_nvp, &
         f0_inode1:f0_inode2,&
         f0_imu1:f0_imu2,&
         ptl_isp:ptl_nsp),stat=alloc_stat)
     call assert(alloc_stat .eq. 0, &
                'alloc(f0_df0g) return istat=',alloc_stat)

    if (allocated(f0_n)) deallocate(f0_n)
#ifndef F_USE_MARKER_DEN2
    allocate(f0_n(-f0_nvp:f0_nvp, &
         f0_inode1:f0_inode2,&
         f0_imu1:f0_imu2,&
         ptl_isp:ptl_nsp),stat=alloc_stat)
#else
    allocate(f0_n(-f0_nvp:f0_nvp, &
         f0_inode1:f0_inode2,&
         f0_imu1:f0_imu2,0:1,&
         ptl_isp:ptl_nsp),stat=alloc_stat)
#endif
     call assert(alloc_stat .eq. 0, &
                'alloc(f0_n) return istat=',alloc_stat)

    if (f0_col_change_weight) then
      if (allocated(f0_df0g2)) deallocate(f0_df0g2)
      allocate(f0_df0g2(-f0_nvp:f0_nvp, &
           f0_inode1:f0_inode2,&
           f0_imu1:f0_imu2,&
           ptl_isp:ptl_nsp),stat=alloc_stat)
       call assert(alloc_stat .eq. 0, &
                   'alloc(f0_df0g2) return istat=',alloc_stat)
    endif
#ifdef F0_TOR_LINEAR
      if (allocated(f0_df0g3)) deallocate(f0_df0g3)
      allocate(f0_df0g3(-f0_nvp:f0_nvp, &
           f0_imu1:f0_imu2, &
           f0_inode1:f0_inode2,0:1,&
           ptl_isp:ptl_nsp),stat=alloc_stat)
       call assert(alloc_stat .eq. 0, &
                   'alloc(f0_df0g3) return istat=',alloc_stat)
#endif
    if (allocated(f0_node_cost)) deallocate(f0_node_cost)
    allocate(f0_node_cost(f0_inode1:f0_inode2), stat=alloc_stat)
     call assert(alloc_stat .eq. 0, &
                 'alloc(f0_node_cost) return istat=',alloc_stat)
    f0_node_cost = 0.D0

#ifdef DIAG_NOISE
    if (allocated(f0_diag_n)) deallocate(f0_diag_n)
    allocate(f0_diag_n(-f0_nvp:f0_nvp, &
         f0_inode1:f0_inode2,&
         f0_imu1:f0_imu2,&
         ptl_isp:ptl_nsp),stat=alloc_stat)
     call assert(alloc_stat .eq. 0, &
                 'alloc(f0_diag_n) return istat=',alloc_stat)

    if (allocated(f0_diag_w0)) deallocate(f0_diag_w0)
    allocate(f0_diag_w0(-f0_nvp:f0_nvp, &
         f0_inode1:f0_inode2,&
         f0_imu1:f0_imu2,&
         ptl_isp:ptl_nsp),stat=alloc_stat)
     call assert(alloc_stat .eq. 0, &
                 'alloc(f0_diag_w0) return istat=',alloc_stat)

    if (allocated(f0_diag_w0s)) deallocate(f0_diag_w0s)
    allocate(f0_diag_w0s(-f0_nvp:f0_nvp, &
         f0_inode1:f0_inode2,&
         f0_imu1:f0_imu2,&
         ptl_isp:ptl_nsp),stat=alloc_stat)
     call assert(alloc_stat .eq. 0, &
                 'alloc(f0_diag_w0s) return istat=',alloc_stat)

    if (allocated(f0_diag_w)) deallocate(f0_diag_w)
    allocate(f0_diag_w(-f0_nvp:f0_nvp, &
         f0_inode1:f0_inode2,&
         f0_imu1:f0_imu2,&
         ptl_isp:ptl_nsp),stat=alloc_stat)
     call assert(alloc_stat .eq. 0, &
                 'alloc(f0_diag_w) return istat=',alloc_stat)

    if (allocated(f0_diag_ws)) deallocate(f0_diag_ws)
    allocate(f0_diag_ws(-f0_nvp:f0_nvp, &
         f0_inode1:f0_inode2,&
         f0_imu1:f0_imu2,&
         ptl_isp:ptl_nsp),stat=alloc_stat)
     call assert(alloc_stat .eq. 0, &
                 'alloc(f0_diag_ws) return istat=',alloc_stat)

    f0_diag_n=0D0
    f0_diag_w0=0D0
    f0_diag_w0s=0D0
    f0_diag_w=0D0
    f0_diag_ws=0D0
#endif

  end subroutine f0_mem_allocation

  ! memory allocation, volume calculation, etc
  subroutine f0_initialize(grid,inode1,inode2,imu1,imu2)
    use sml_module
    use grid_class
    implicit none
    type(grid_type) :: grid
    integer, intent(in) :: inode1, inode2, imu1, imu2


    f0_inode1=inode1
    f0_inode2=inode2
    f0_imu1=imu1
    f0_imu2=imu2

    call f0_mem_allocation

 !   allocate(f0_f0a(-f0_nvp:f0_nvp, &
 !        f0_inode1:f0_inode2,&
 !        f0_imu1:f0_imu2,&
 !        ptl_isp:ptl_nsp))


    ! set zero initial value -- read in at restart_read
    f0_f0g=0D0
    f0_f =0D0
    f0_df0g=0.D0
    f0_n=0D0

    if (f0_col_change_weight) then
      f0_df0g2=0.D0
    endif

    call f0_init_rest(grid)
  end subroutine f0_initialize


  ! After node weights have been updated, f0 has to be redistributed.
  ! This routine sends the grid points a process lost and receives the grid points
  ! a process gained during the update of the poloidal domain decomposition.
  ! To do: make general by using f0_imu1 and f0_imu2!!!
  subroutine f0_redistribute(grid,imu1,imu2)
    use sml_module
    use ptl_module, only: ptl_isp,ptl_nsp
    use grid_class
    use perf_monitor
    implicit none
    include 'mpif.h'
    type(grid_type) :: grid
    integer, intent(in) :: imu1, imu2

    ! local variables
    integer :: alloc_stat, ldim1, ldim2, ldim3
    integer :: is, ie, i, j, k, isendtag, istatus(MPI_STATUS_SIZE), dum1(1), dum
    integer :: msize, nsp, ierr
    integer :: my_nodes(2), my_nodes_old(2), my_nodes_cnt, offset
    integer :: nfields
    integer :: pid, ssignal, rsignal, send_counts, recv_counts
    integer :: keep_min, keep_max
    integer :: send_left_min, send_left_max, send_right_min, send_right_max
    integer :: recv_left_min, recv_left_max, recv_right_min, recv_right_max
    integer :: send_count(0:sml_pe_per_plane-1), recv_count(0:sml_pe_per_plane-1)
    integer :: send_index(0:sml_pe_per_plane-1), recv_index(0:sml_pe_per_plane-1)
    integer :: rrequest(0:sml_pe_per_plane-1)
    integer, parameter :: F0COST  = 1
#ifndef F0_TOR_LINEAR
    integer, parameter :: F0F0G   = 2
    integer, parameter :: F0F     = 3
    integer, parameter :: F0DF0G  = 4
    integer, parameter :: F0DF0G2 = 5
#else
    integer, parameter :: F0F0G0  = 2
    integer, parameter :: F0F0G1  = 3
    integer, parameter :: F0F     = 4
    integer, parameter :: F0DF0G  = 5
    integer, parameter :: F0DF0G2 = 6
#endif
    integer, parameter :: MAXPREPOSTS = 20   ! nothing magical to this number - just be reasonable
    real (kind=8), allocatable, dimension(:,:,:,:,:) :: sendarr
    real (kind=8), allocatable, dimension(:,:,:,:,:) :: temp

    call t_startf("F0_REDIST_CHKPT0")
    call check_point('f0_redistribute: start')
    call t_stopf("F0_REDIST_CHKPT0")

    nsp=ptl_nsp-ptl_isp+1

    ! send/receive only nodes that need to be shifted
    my_nodes(1)=gvid0_pid(sml_plane_mype)
    my_nodes(2)=gvid0_pid(sml_plane_mype+1)-1
    offset = my_nodes(1) - 1
    my_nodes_cnt = my_nodes(2) - offset

    my_nodes_old(1)=gvid0_pid_old(sml_plane_mype)
    my_nodes_old(2)=gvid0_pid_old(sml_plane_mype+1)-1

#ifndef F0_TOR_LINEAR
    nfields = 4
#else
    nfields = 5
#endif
    if (allocated(f0_df0g2)) then
      nfields = nfields + 1
    endif

    allocate(temp(-f0_nvp:f0_nvp,0:f0_nmu,ptl_isp:ptl_nsp,nfields,my_nodes_cnt), stat=alloc_stat)
     call assert(alloc_stat .eq. 0, &
                'alloc(temp) in f0_redistribute return istat=',alloc_stat)
    ldim1 = lbound(temp,1)
    ldim2 = lbound(temp,2)
    ldim3 = lbound(temp,3)

    ! Nodes retained are {max(old_min,new_min), min(old_max,new_max)}
    keep_min = max(my_nodes_old(1), my_nodes(1))
    keep_max = min(my_nodes_old(2), my_nodes(2))

    ! Nodes sent to the left are {old_min, min(old_max, new_min-1)}
    send_left_min = my_nodes_old(1)
    send_left_max = min(my_nodes_old(2),my_nodes(1)-1)
    
    ! Nodes sent to the right are {max(old_min,new_max+1), old_max}
    send_right_min = max(my_nodes_old(1),my_nodes(2)+1)
    send_right_max = my_nodes_old(2)

    ! Determine intersection of left/right send intervals with new partition, to determine
    ! what to send to neighboring processes in the same plane
    call partition_intersect(send_left_min, send_left_max, send_right_min, send_right_max, &
                             grid%nnode, sml_pe_per_plane, gvid0_pid, send_index, send_count, &
                             send_counts, ierr)
    if (ierr .ne. 0) then
      write(6,*) "error in left/right send intervals:", sml_mype, sml_plane_mype, &
        send_left_min, send_right_max, send_right_min, send_right_max
      stop
    endif

    ! Nodes received from left are {new_min, min(old_min-1,new_max)}
    recv_left_min = my_nodes(1)
    recv_left_max = min(my_nodes_old(1)-1,my_nodes(2))

    ! Nodes received from right are {max(old_max+1,new_min), new_max}
    recv_right_min = max(my_nodes_old(2)+1,my_nodes(1))
    recv_right_max = my_nodes(2)

    ! Determine intersection of left/right receive intervals with old partition, to determine
    ! what to receive from neighboring processes in the same plane
    call partition_intersect(recv_left_min, recv_left_max, recv_right_min, recv_right_max, &
                             grid%nnode, sml_pe_per_plane, gvid0_pid_old, recv_index, recv_count, &
                             recv_counts, ierr)

    if (ierr .ne. 0) then
      write(6,*) "error in left/right receive intervals:", sml_mype, sml_plane_mype, &
        recv_left_min, recv_right_max, recv_right_min, recv_right_max
      stop
    endif

    ! Determine MPI message size
    msize=(2*f0_nvp+1) * (f0_nmu+1) * nsp * nfields

    ! allocate send buffer
    dum1=maxval(send_count)
    dum=max(1,dum1(1))
    allocate(sendarr(-f0_nvp:f0_nvp,0:f0_nmu,ptl_isp:ptl_nsp,nfields,0:dum-1), stat=alloc_stat)
     call assert(alloc_stat .eq. 0, &
                'alloc(sendarr) in f0_redistribute return istat=',alloc_stat)
     sendarr = 0.D0

    ! If not too many, prepost all receive requests
    ! (using XOR ordering and flow control)
    if (recv_counts <= MAXPREPOSTS) then
      do i=1,ceil2(sml_pe_per_plane)-1
        pid = pair(sml_pe_per_plane,i,sml_plane_mype)
        if (pid >= 0) then
          if (recv_count(pid) > 0) then
            call mpi_irecv(temp(-f0_nvp,0,ptl_isp,1,(recv_index(pid)-offset)), &
                           recv_count(pid)*msize, MPI_REAL8, pid, pid, &
                           sml_plane_comm, rrequest(pid), ierr)
            call mpi_send (rsignal, 1, MPI_INTEGER, pid, pid, &
                           sml_plane_comm, ierr)
          endif
        endif
      enddo
    endif

    ! Copy information that stays on the same process
    do i=keep_min,keep_max
      temp(ldim1,ldim2,ldim3,F0COST,i-offset) = f0_node_cost(i)
#ifndef F0_TOR_LINEAR
      temp(:,:,:,F0F0G,i-offset)  =  f0_f0g(:,i,:,:)
#else
      temp(:,:,:,F0F0G0,i-offset)  =  f0_f0g(:,i,:,0,:)
      temp(:,:,:,F0F0G1,i-offset)  =  f0_f0g(:,i,:,1,:)
#endif
      temp(:,:,:,F0F,i-offset)    =    f0_f(:,i,:,:)
      temp(:,:,:,F0DF0G,i-offset) = f0_df0g(:,i,:,:)

      if (allocated(f0_df0g2)) then
        temp(:,:,:,F0DF0G2,i-offset)=f0_df0g2(:,i,:,:)
      endif

    enddo

    ! send/recv data using XOR ordering and flow control
    do i=1,ceil2(sml_pe_per_plane)-1
      pid = pair(sml_pe_per_plane,i,sml_plane_mype)
      if (pid >= 0) then

        if (recv_count(pid) > 0) then

          if (recv_counts > MAXPREPOSTS) then
            ! Post receive request
            call mpi_irecv(temp(-f0_nvp,0,ptl_isp,1,(recv_index(pid)-offset)), &
                           recv_count(pid)*msize, MPI_REAL8, pid, pid, &
                           sml_plane_comm, rrequest(pid), ierr)
            call mpi_send (rsignal, 1, MPI_INTEGER, pid, pid, &
                           sml_plane_comm, ierr)

          endif

        endif

        if (send_count(pid) > 0) then

          ! Fill send buffer
          do j=0,send_count(pid)-1
            sendarr(ldim1,ldim2,ldim3,F0COST,j) = f0_node_cost(send_index(pid)+j)
          enddo
          do j=0,send_count(pid)-1
#ifndef F0_TOR_LINEAR
            sendarr(:,:,:,F0F0G,j)  =  f0_f0g(:,send_index(pid)+j,:,:)
#else
            sendarr(:,:,:,F0F0G0,j)  =  f0_f0g(:,send_index(pid)+j,:,0,:)
            sendarr(:,:,:,F0F0G1,j)  =  f0_f0g(:,send_index(pid)+j,:,1,:)
#endif
          enddo
          do j=0,send_count(pid)-1
            sendarr(:,:,:,F0F,j)    =    f0_f(:,send_index(pid)+j,:,:)
          enddo
          do j=0,send_count(pid)-1
            sendarr(:,:,:,F0DF0G,j) = f0_df0g(:,send_index(pid)+j,:,:)
          enddo

          if (allocated(f0_df0g2)) then
            do j=0,send_count(pid)-1
              sendarr(:,:,:,F0DF0G2,j) = f0_df0g2(:,send_index(pid)+j,:,:)
            enddo
          endif

          ! Wait for signal, and then send
          call mpi_recv (ssignal, 1, MPI_INTEGER, pid, sml_plane_mype, &
                         sml_plane_comm, istatus, ierr)

          call mpi_rsend(sendarr, send_count(pid)*msize, &
                         MPI_REAL8, pid, sml_plane_mype, &
                         sml_plane_comm, ierr)

        endif

        if (recv_count(pid) > 0) then
          ! wait for message
          ! Note: since using XOR ordering and flow control, no reason to delay
          !       the wait
          call mpi_wait(rrequest(pid), istatus, ierr)
        endif

      endif

    enddo


    ! Deallocate old f0_f and f0_f0g and re-initialize
    deallocate(f0_f0g,f0_f,f0_grid_vol,f0_grid_vol_vonly,f0_B_B0,f0_n_Ta, f0_den,f0_T_ev)
    deallocate(f0_df0g)
    ! Note: f0_n does not have to be redistributed, only re-allocated
    deallocate(f0_n)
    if (allocated(f0_df0g2)) then
      deallocate(f0_df0g2)
    endif

    call f0_initialize(grid,my_nodes(1),my_nodes(2),imu1,imu2)

    !call check_point('f0_redistribute: re-allocated f0 memory')

    do i=my_nodes(1),my_nodes(2)
      f0_node_cost(i) = temp(ldim1,ldim2,ldim3,F0COST,i-offset)
#ifndef F0_TOR_LINEAR
      f0_f0g(:,i,:,:) = temp(:,:,:,F0F0G,i-offset)
#else
      f0_f0g(:,i,:,0,:) = temp(:,:,:,F0F0G0,i-offset)
      f0_f0g(:,i,:,1,:) = temp(:,:,:,F0F0G1,i-offset)
#endif
      f0_f(:,i,:,:)   = temp(:,:,:,F0F,i-offset)
      f0_df0g(:,i,:,:)= temp(:,:,:,F0DF0G,i-offset)
      if (allocated(f0_df0g2)) then
        f0_df0g2(:,i,:,:) = temp(:,:,:,F0DF0G2,i-offset)
      endif
    enddo

    deallocate(temp,sendarr)

    call t_startf("F0_REDIST_CHKPT1")
    call check_point('f0_redistribute: finished')
    call t_stopf("F0_REDIST_CHKPT1")

#ifdef DIAG_NOISE
    print *, 'Not implemented f0_redistribute, yet'

#endif
  end subroutine f0_redistribute


  subroutine f0_init_rest(grid)
    use sml_module
    use grid_class
    use ptl_module
    use eq_module
    implicit none
    type(grid_type) :: grid

    integer :: alloc_stat
    integer :: i, isp
    real (8) :: r,z,b,psi,T_ev
    real (8), external :: b_interpol
    real (8) :: den

    f0_dsmu = f0_smu_max / real(f0_nmu,8)
    f0_dvp = f0_vp_max / real(f0_nvp,8)

    if(allocated(f0_grid_vol)) deallocate(f0_grid_vol)
    allocate(f0_grid_vol(f0_inode1:f0_inode2,ptl_isp:ptl_nsp), stat=alloc_stat)
     call assert(alloc_stat .eq. 0, &
                 'alloc(f0_grid_vol) return istat=',alloc_stat)

    if(allocated(f0_grid_vol_vonly)) deallocate(f0_grid_vol_vonly)
    allocate(f0_grid_vol_vonly(f0_inode1:f0_inode2,ptl_isp:ptl_nsp), stat=alloc_stat)
     call assert(alloc_stat .eq. 0, &
                 'alloc(f0_grid_vol_vonly) return istat=',alloc_stat)

    if(allocated(f0_B_B0)) deallocate(f0_B_B0,f0_n_Ta,f0_den)
    allocate(f0_B_B0(f0_inode1:f0_inode2), &
             f0_n_Ta(f0_inode1:f0_inode2,ptl_isp:ptl_nsp), &
             f0_den(f0_inode1:f0_inode2), stat=alloc_stat)
     call assert(alloc_stat .eq. 0, &
                 'alloc(f0_B_B0,f0_n_Ta,F0_den) return istat=',alloc_stat)

    if(allocated(f0_T_ev)) deallocate(f0_T_ev)
    allocate(f0_T_ev(f0_inode1:f0_inode2,ptl_isp:ptl_nsp), stat=alloc_stat)
     call assert(alloc_stat .eq. 0, &
                 'alloc(f0_T_ev) return istat=',alloc_stat)

    if (allocated(f0_node_cost)) deallocate(f0_node_cost)
    allocate(f0_node_cost(f0_inode1:f0_inode2), stat=alloc_stat)
     call assert(alloc_stat .eq. 0, &
                 'alloc(f0_node_cost) return istat=',alloc_stat)
    f0_node_cost = 0.D0

    !setup grid_vol, B and n/T^(3/2)
    do i=f0_inode1, f0_inode2
       r=grid%x(1,i)
       z=grid%x(2,i)
       B=B_interpol(r,z,0D0) ! should be B* to be exact
       psi=grid%psi(i)


       f0_B_B0(i)=B/eq_axis_b

       den=eq_ftn(psi,r,z,eq_den)
       f0_den(i)=den

       do isp=ptl_isp, ptl_nsp
          if(isp==0) then
             T_ev=eq_ftn(psi,r,z,eq_tempe)
          else
             T_ev=eq_ftn(psi,r,z,eq_tempi)
          endif
          f0_T_ev(i,isp)=T_ev


          ! When MU_LINEAR is defined -- original code
          ! Phase-Space volume of f
          ! f is defined in sqrt(m/(2pi e^3)) B dmu dv|| space
          ! dmu = T_ev * e / B_0 dmu_n
          ! dvp = sqrt(T_ev*e/m) dvp_n
          ! _n mean normalized quantity
          ! sqrt(T_ev^3/2pi) B/B_0 delta mu_n delta vp_n is the volme of velocity grid.
          ! space volume is obtained from grid%node_vol_nearest
          ! why factor of 2 in front ??

          ! When MU_LINEAR is not defined
          ! f is defined in sqrt(m/(2pi e^2)) B/sqrt(B0) dsmu dv|| space
          ! In general f is not conserved in this space.
          ! Since mu is constant in time, f is conserved following path.
          ! With normalization, the volume element becomes sqrt(1/2pi) B/B0 T_ev dsmu_n dvp_n

!          f0_grid_vol(i,isp)=grid%node_vol_nearest(i)*sqrt(T_ev**3/sml_2pi)*B/eq_axis_b*f0_dsmu*f0_dvp
          f0_grid_vol_vonly(i,isp)=T_ev*sqrt(1D0/sml_2pi)*f0_dsmu*f0_dvp  ! v-space term only
          f0_grid_vol(i,isp)=grid%node_vol_nearest(i)*f0_grid_vol_vonly(i,isp)


          f0_n_Ta(i,isp)=den/T_ev

       enddo
    enddo


    ! Allocate f0_den_global and initialize it
    if(.not. allocated(f0_den_global)) then
      allocate(f0_den_global(grid%nnode))

      ! following do loop can be opitmized using mpi_reduce
      do i=1, grid%nnode
        r=grid%x(1,i)
        z=grid%x(2,i)
        psi=grid%psi(i)
        den=eq_ftn(psi,r,z,eq_den)
        f0_den_global(i)=den
      enddo
    endif

#ifdef  F0_CHARGE_N0
    if(.not. allocated(f0_den0_ptl)) then
       allocate(f0_den0_ptl(grid%nnode,ptl_isp:ptl_nsp))
       allocate(f0_density_n0_add(grid%nnode,ptl_isp:ptl_nsp))
       f0_den0_ptl=0D0
       f0_density_n0_add=0D0
    endif
#endif


  end subroutine f0_init_rest

  !
  subroutine f0_get_f0g(grid,isp,itr,p,phi,mu_n_in,vp_n,f0g,err)
    use grid_class
    implicit none
    type(grid_type) :: grid
    integer,  intent(in) :: isp,itr
    real (8), intent(in) :: p(3)
    real (8), intent(in) :: phi, mu_n_in, vp_n
    real (8), intent(out) :: f0g
    logical, intent(out) :: err
    !
    real (8) :: mu_n   ! v_perp^2
    integer :: i_mu, i_vp, ip, node, ml(1)
    real (8) ::  wmu(0:1), wvp(0:1), smu, smu0, smu1
    logical, external :: is_nan
    real (8) :: wphi(0:1) ! for F0_TOR_LINEAR

    !Set default value / initial value of f0g
    f0g = 0D0

    !exclude out-of-grid particle -- 1st
    if(itr<=0) then
       return ! zero f0g
    endif

    !find nearest node
    ml=maxloc(p)
    node= grid%nd(ml(1),itr)
    if(node < f0_inode1 .or. node > f0_inode2) then
       print *, 'wrong particle found, isp=',isp
       print *, 'f0_inode1, f0_inode2 node', f0_inode1, f0_inode2, node
       print *, 'grid%nd',grid%nd(:,itr), itr
       print *, 'p', p
       print *, 'x from nd and p=', grid%x(:,grid%nd(1,itr))*p(1) + grid%x(:,grid%nd(2,itr))*p(2) + grid%x(:,grid%nd(3,itr))*p(3)
       print *, 'x1=', grid%x(:,grid%nd(1,itr)), p(1)
       print *, 'x2=', grid%x(:,grid%nd(2,itr)), p(2)
       print *, 'x3=', grid%x(:,grid%nd(3,itr)), p(3)
       err=.true.
       !stop
       return
    else
       err=.false.
    endif

    !get mu_n (It is actually v_perp^2/vth^2)
    !when V_PERP is not defined, mu_n is 2*B(grid)*mu / T(particle)
#ifndef V_PERP
    mu_n = mu_n_in*f0_B_B0(node)
#else
    mu_n = mu_n_in
#endif


#ifdef F0_TOR_LINEAR
    wphi(1)= phi/grid%delta_phi - grid%iphi_offset
    wphi(0)= 1D0 - wphi(1)
#endif

    smu=sqrt(mu_n)       ! square root mu

    wmu(0)=smu/f0_dsmu   ! temp. variable - normalized sqrt(mu) with delta
    i_mu = floor(wmu(0)) ! get index for lower grid
    wmu(1)=wmu(0)-real(i_mu,8) ! wmu(1) is now weight for upper grid
    wmu(0)=1D0 - wmu(1)        ! wmu(0) is now weight for lower grid

    if(i_mu==0) then
       smu0=f0_dsmu/f0_mu0_factor
       smu1=f0_dsmu
    elseif(i_mu==f0_nmu-1) then
       smu0=f0_dsmu*real(i_mu,8)
       smu1=f0_dsmu*(real(f0_nmu,8)-1D0/f0_mu0_factor)
    else
       smu0=f0_dsmu*real(i_mu,8)
       smu1=f0_dsmu*real(i_mu+1,8)
    endif


    wvp(0)=vp_n/f0_dvp
    i_vp= floor(wvp(0))
    wvp(1)=wvp(0)-real(i_vp,8)
    wvp(0)=1D0-wvp(1)

    !exclude out-of-grid particle - 2nd
    if(i_mu >= f0_nmu .or. i_vp >= f0_nvp .or. i_vp < -f0_nvp) then
       return ! zero f0g
    endif

    !sum-up

    ! This is for v_perp-v_parallel grid --->
#ifndef F0_TOR_LINEAR
    f0g=f0g+ wmu(0)/(smu0)*wvp(0)*f0_f0g(i_vp+0,node,i_mu+0,isp)
    f0g=f0g+ wmu(0)/(smu0)*wvp(1)*f0_f0g(i_vp+1,node,i_mu+0,isp)
    f0g=f0g+ wmu(1)/(smu1)*wvp(0)*f0_f0g(i_vp+0,node,i_mu+1,isp)
    f0g=f0g+ wmu(1)/(smu1)*wvp(1)*f0_f0g(i_vp+1,node,i_mu+1,isp)
#else
    f0g=f0g+ wmu(0)/(smu0)*wvp(0)*wphi(0)*f0_f0g(i_vp+0,node,i_mu+0,0,isp)
    f0g=f0g+ wmu(0)/(smu0)*wvp(1)*wphi(0)*f0_f0g(i_vp+1,node,i_mu+0,0,isp)
    f0g=f0g+ wmu(1)/(smu1)*wvp(0)*wphi(0)*f0_f0g(i_vp+0,node,i_mu+1,0,isp)
    f0g=f0g+ wmu(1)/(smu1)*wvp(1)*wphi(0)*f0_f0g(i_vp+1,node,i_mu+1,0,isp)

    f0g=f0g+ wmu(0)/(smu0)*wvp(0)*wphi(1)*f0_f0g(i_vp+0,node,i_mu+0,1,isp)
    f0g=f0g+ wmu(0)/(smu0)*wvp(1)*wphi(1)*f0_f0g(i_vp+1,node,i_mu+0,1,isp)
    f0g=f0g+ wmu(1)/(smu1)*wvp(0)*wphi(1)*f0_f0g(i_vp+0,node,i_mu+1,1,isp)
    f0g=f0g+ wmu(1)/(smu1)*wvp(1)*wphi(1)*f0_f0g(i_vp+1,node,i_mu+1,1,isp)
#endif
    f0g=f0g/sqrt(f0_T_ev(node,isp))

    if(is_nan(f0g)) then
      print *, 'isp,itr,p,phi,mu_n,vp_n =', isp,itr,p,phi,mu_n,vp_n
      print *, 'wmu(0:1), wvp(0:1), smu0, smu1=',wmu(0:1), wvp(0:1), smu0, smu1
#ifndef F0_TOR_LINEAR
      print *, 'f0_f0g(i_vp+0~1,node,i_mu+0~1,isp)', f0_f0g(i_vp+0,node,i_mu+0,isp),  &
        f0_f0g(i_vp+1,node,i_mu+0,isp), &
        f0_f0g(i_vp+0,node,i_mu+1,isp), &
        f0_f0g(i_vp+1,node,i_mu+1,isp)

      ! conflict of memory access should not be harmfull - zeroing out
      call set_zero_for_nan(f0_f0g(i_vp+0,node,i_mu+0,isp))
      call set_zero_for_nan(f0_f0g(i_vp+1,node,i_mu+0,isp))
      call set_zero_for_nan(f0_f0g(i_vp+0,node,i_mu+1,isp))
      call set_zero_for_nan(f0_f0g(i_vp+1,node,i_mu+1,isp))
      f0g=0D0

#else

  ! not implimented.
#endif
    endif

  end subroutine f0_get_f0g

  subroutine set_zero_for_nan(var)
    implicit none
    real (8), intent(inout) :: var
    logical, external :: is_nan

    if(is_nan(var)) var=0D0

  end subroutine

#ifndef DIAG_NOISE
  subroutine f0_update_f0g(grid,isp,itr,p,phi,mu_n,vp_n,df0g,df2,s1g, s2g, s1f, s2f,iflag)
#else
  subroutine f0_update_f0g(grid,isp,itr,p,phi,mu_n,vp_n,df0g,df2,s1g, s2g, s1f, s2f,iflag, w0, w1w0)
#endif
    use grid_class
    use sml_module
    implicit none
    type(grid_type) :: grid
    integer,  intent(in) :: isp,itr
    real (8), intent(in) :: p(3)
    real (8), intent(in) :: phi, mu_n, vp_n
    real (8), intent(in) :: df0g, df2
    real (8), dimension(-f0_nvp:f0_nvp,f0_inode1:f0_inode2), intent(inout) :: s1g, s2g, s1f, s2f
    integer, intent(out) :: iflag
    integer :: i_mu, i_vp, ip, node, ml(1)
    real (8) :: wmu(0:1), wvp(0:1), wp(0:1), wp2(0:1), smu, mu_vol(0:1)
    real (8) :: t00, t01, t10, t11, smu0, smu1
    real (8) :: w0, w1w0
    real (8) :: wphi(0:1), wp_p0(0:1), wp_p1(0:1) ! for F0_TOR_LINEAR

#ifdef F0_TOR_LINEAR
    wphi(1)= phi/grid%delta_phi - grid%iphi_offset
    wphi(0)= 1D0 - wphi(1)
#endif


    smu=sqrt(mu_n)       ! square root mu


    wmu(0)=smu/f0_dsmu   ! temp. variable - normalized sqrt(mu) with delta
    i_mu = floor(wmu(0)) ! get index for lower grid
    wmu(1)=wmu(0)-real(i_mu,8) ! wmu(1) is now weight for upper grid
    wmu(0)=1D0 - wmu(1)        ! wmu(0) is now weight for lower grid


    wvp(0)=vp_n/f0_dvp
    i_vp= floor(wvp(0))
    wvp(1)=wvp(0)-real(i_vp,8)
    wvp(0)=1D0-wvp(1)

    !exclude out-of-grid particle
    if(itr<=0 .or. i_mu >= f0_nmu .or. i_vp >= f0_nvp .or. i_vp < -f0_nvp) then
       iflag=-1
       return
    endif

    ml=maxloc(p)
    node= grid%nd(ml(1),itr)

    mu_vol=1D0 ! 0 and 1 index
    if(i_mu==0) then
       mu_vol(0)=0.5D0
       smu0=f0_dsmu/f0_mu0_factor
       smu1=f0_dsmu
    elseif (i_mu==f0_nmu-1) then
       mu_vol(1)=0.5D0
       smu0=f0_dsmu*real(i_mu,8)
       smu1=f0_dsmu*(real(f0_nmu,8)-1D0/f0_mu0_factor)
    else
       smu0=f0_dsmu*real(i_mu,8)
       smu1=f0_dsmu*real(i_mu+1,8)
    endif

    !rh For f0g: phase space volume element
#ifndef F0_TOR_LINEAR
    wp(0:1)=1D0/(f0_grid_vol(node,isp)*mu_vol(0:1))
#else
#ifndef F_USE_MARKER_DEN2
    wp_p0(0:1)=1D0/(f0_grid_vol_vonly(node,isp)*grid%node_vol_ff(node,0)*mu_vol(0:1))
    wp_p1(0:1)=1D0/(f0_grid_vol_vonly(node,isp)*grid%node_vol_ff(node,1)*mu_vol(0:1))
#else
    wp_p0(0)=smu0*sqrt(f0_T_ev(node,isp))
    wp_p0(1)=smu1*sqrt(f0_T_ev(node,isp))
    wp_p1=wp_p0
#endif
#endif
    !rh For f0_f: conversion from mu-vpar to sqrt(mu)-vpar
#ifdef F_USE_MARKER_DEN
    wp2(0)=smu0*sqrt(f0_T_ev(node,isp))
    wp2(1)=smu1*sqrt(f0_T_ev(node,isp))
#else
    wp(0:1)=1D0/(f0_grid_vol(node,isp)*mu_vol(0:1))
    wp2=wp
#endif

    !temp variable for optimization
    t00=wmu(0)*wvp(0)
    t01=wmu(0)*wvp(1)
    t10=wmu(1)*wvp(0)
    t11=wmu(1)*wvp(1)

    ! for f0_grid
#ifndef F0_TOR_LINEAR
    f0_f0g(i_vp+0,node,i_mu+0,isp)=f0_f0g(i_vp+0,node,i_mu+0,isp)+df0g*t00*wp(0)
    f0_f0g(i_vp+1,node,i_mu+0,isp)=f0_f0g(i_vp+1,node,i_mu+0,isp)+df0g*t01*wp(0)
    f0_f0g(i_vp+0,node,i_mu+1,isp)=f0_f0g(i_vp+0,node,i_mu+1,isp)+df0g*t10*wp(1)
    f0_f0g(i_vp+1,node,i_mu+1,isp)=f0_f0g(i_vp+1,node,i_mu+1,isp)+df0g*t11*wp(1)
#else
    f0_df0g3(i_vp+0,i_mu+0,node,0,isp)=f0_df0g3(i_vp+0,i_mu+0,node,0,isp)+df0g*t00*wp_p0(0)*wphi(0)
    f0_df0g3(i_vp+1,i_mu+0,node,0,isp)=f0_df0g3(i_vp+1,i_mu+0,node,0,isp)+df0g*t01*wp_p0(0)*wphi(0)
    f0_df0g3(i_vp+0,i_mu+1,node,0,isp)=f0_df0g3(i_vp+0,i_mu+1,node,0,isp)+df0g*t10*wp_p0(1)*wphi(0)
    f0_df0g3(i_vp+1,i_mu+1,node,0,isp)=f0_df0g3(i_vp+1,i_mu+1,node,0,isp)+df0g*t11*wp_p0(1)*wphi(0)

    f0_df0g3(i_vp+0,i_mu+0,node,1,isp)=f0_df0g3(i_vp+0,i_mu+0,node,1,isp)+df0g*t00*wp_p1(0)*wphi(1)
    f0_df0g3(i_vp+1,i_mu+0,node,1,isp)=f0_df0g3(i_vp+1,i_mu+0,node,1,isp)+df0g*t01*wp_p1(0)*wphi(1)
    f0_df0g3(i_vp+0,i_mu+1,node,1,isp)=f0_df0g3(i_vp+0,i_mu+1,node,1,isp)+df0g*t10*wp_p1(1)*wphi(1)
    f0_df0g3(i_vp+1,i_mu+1,node,1,isp)=f0_df0g3(i_vp+1,i_mu+1,node,1,isp)+df0g*t11*wp_p1(1)*wphi(1)
#endif

#ifdef F_USE_MARKER_DEN2
  !F0_LINEAR is defined together
       f0_n(i_vp+0,node,i_mu+0,0,isp)=f0_n(i_vp+0,node,i_mu+0,0,isp)+t00*wphi(0)
       f0_n(i_vp+1,node,i_mu+0,0,isp)=f0_n(i_vp+1,node,i_mu+0,0,isp)+t01*wphi(0)
       f0_n(i_vp+0,node,i_mu+1,0,isp)=f0_n(i_vp+0,node,i_mu+1,0,isp)+t10*wphi(0)
       f0_n(i_vp+1,node,i_mu+1,0,isp)=f0_n(i_vp+1,node,i_mu+1,0,isp)+t11*wphi(0)

       f0_n(i_vp+0,node,i_mu+0,1,isp)=f0_n(i_vp+0,node,i_mu+0,1,isp)+t00*wphi(1)
       f0_n(i_vp+1,node,i_mu+0,1,isp)=f0_n(i_vp+1,node,i_mu+0,1,isp)+t01*wphi(1)
       f0_n(i_vp+0,node,i_mu+1,1,isp)=f0_n(i_vp+0,node,i_mu+1,1,isp)+t10*wphi(1)
       f0_n(i_vp+1,node,i_mu+1,1,isp)=f0_n(i_vp+1,node,i_mu+1,1,isp)+t11*wphi(1)
#endif

#ifdef COL_F_NAN_CHECK
    if(.not. (f0_f(2,node,1,isp) > 1D0 .or. f0_f(2, node, 1, isp) < 2D0) ) then
       print *, sml_mype, '] NAN FOUND f0_f0g (f0_update_f0g)', f0_f0g(2, node,1, isp), i_vp, i_mu, df0g
       stop
    endif
#endif
    ! for f_grid
    if(.not. sml_no_fp_in_f .and. mod(sml_istep,sml_f_source_period)==0) then
       f0_f(i_vp+0,node,i_mu+0,isp)=f0_f(i_vp+0,node,i_mu+0,isp)+df2*t00*wp2(0)
       f0_f(i_vp+1,node,i_mu+0,isp)=f0_f(i_vp+1,node,i_mu+0,isp)+df2*t01*wp2(0)
       f0_f(i_vp+0,node,i_mu+1,isp)=f0_f(i_vp+0,node,i_mu+1,isp)+df2*t10*wp2(1)
       f0_f(i_vp+1,node,i_mu+1,isp)=f0_f(i_vp+1,node,i_mu+1,isp)+df2*t11*wp2(1)

       ! Normalization for reverse interpolation --> This needs to be done
       ! regardless of the value of df2!!!
#ifndef F_USE_MARKER_DEN2
       f0_n(i_vp+0,node,i_mu+0,isp)=f0_n(i_vp+0,node,i_mu+0,isp)+t00
       f0_n(i_vp+1,node,i_mu+0,isp)=f0_n(i_vp+1,node,i_mu+0,isp)+t01
       f0_n(i_vp+0,node,i_mu+1,isp)=f0_n(i_vp+0,node,i_mu+1,isp)+t10
       f0_n(i_vp+1,node,i_mu+1,isp)=f0_n(i_vp+1,node,i_mu+1,isp)+t11
#else
       !done already above
#endif
    endif


! for mu-parallelization. based on old code. Need updates to use
!!$    ! for boundary operation to send to other processor
!!$    if(i_mu/=0 .and. i_mu==f0_imu1) then
!!$       s1g(i_vp+0,node)=s1g(i_vp+0,node)+df0g*t00
!!$       s1g(i_vp+1,node)=s1g(i_vp+1,node)+df0g*t01
!!$       if(.not. sml_no_fp_in_f) then
!!$          s1f(i_vp+0,node)=s1f(i_vp+0,node)+df2 *t00
!!$          s1f(i_vp+1,node)=s1f(i_vp+1,node)+df2 *t01
!!$       endif
!!$    endif
!!$
!!$    if(i_mu/=f0_nmu .and. i_mu==f0_imu2) then
!!$       s2g(i_vp+0,node)=s2g(i_vp+0,node)+df0g*t10
!!$       s2g(i_vp+1,node)=s2g(i_vp+1,node)+df0g*t11
!!$       if(.not. sml_no_fp_in_f) then
!!$          s2f(i_vp+0,node)=s2f(i_vp+0,node)+df2 *t10
!!$          s2f(i_vp+1,node)=s2f(i_vp+1,node)+df2 *t11
!!$       endif
!!$    endif


!    enddo

#ifdef DIAG_NOISE
    wp(0:1)=1D0/(f0_grid_vol(node,isp)*mu_vol(0:1))
    f0_diag_n(i_vp+0,node,i_mu+0,isp)=f0_diag_n(i_vp+0,node,i_mu+0,isp)+1D0*t00/wp(0)
    f0_diag_n(i_vp+1,node,i_mu+0,isp)=f0_diag_n(i_vp+1,node,i_mu+0,isp)+1D0*t01/wp(0)
    f0_diag_n(i_vp+0,node,i_mu+1,isp)=f0_diag_n(i_vp+0,node,i_mu+1,isp)+1D0*t10/wp(1)
    f0_diag_n(i_vp+1,node,i_mu+1,isp)=f0_diag_n(i_vp+1,node,i_mu+1,isp)+1D0*t11/wp(1)


    f0_diag_w0(i_vp+0,node,i_mu+0,isp)=f0_diag_w0(i_vp+0,node,i_mu+0,isp)+w0*t00/wp(0)
    f0_diag_w0(i_vp+1,node,i_mu+0,isp)=f0_diag_w0(i_vp+1,node,i_mu+0,isp)+w0*t01/wp(0)
    f0_diag_w0(i_vp+0,node,i_mu+1,isp)=f0_diag_w0(i_vp+0,node,i_mu+1,isp)+w0*t10/wp(1)
    f0_diag_w0(i_vp+1,node,i_mu+1,isp)=f0_diag_w0(i_vp+1,node,i_mu+1,isp)+w0*t11/wp(1)
    

    f0_diag_w0s(i_vp+0,node,i_mu+0,isp)=f0_diag_w0s(i_vp+0,node,i_mu+0,isp)+w0*w0*t00/wp(0)
    f0_diag_w0s(i_vp+1,node,i_mu+0,isp)=f0_diag_w0s(i_vp+1,node,i_mu+0,isp)+w0*w0*t01/wp(0)
    f0_diag_w0s(i_vp+0,node,i_mu+1,isp)=f0_diag_w0s(i_vp+0,node,i_mu+1,isp)+w0*w0*t10/wp(1)
    f0_diag_w0s(i_vp+1,node,i_mu+1,isp)=f0_diag_w0s(i_vp+1,node,i_mu+1,isp)+w0*w0*t11/wp(1)

    f0_diag_w(i_vp+0,node,i_mu+0,isp)=f0_diag_w(i_vp+0,node,i_mu+0,isp)+w1w0*t00/wp(0)
    f0_diag_w(i_vp+1,node,i_mu+0,isp)=f0_diag_w(i_vp+1,node,i_mu+0,isp)+w1w0*t01/wp(0)
    f0_diag_w(i_vp+0,node,i_mu+1,isp)=f0_diag_w(i_vp+0,node,i_mu+1,isp)+w1w0*t10/wp(1)
    f0_diag_w(i_vp+1,node,i_mu+1,isp)=f0_diag_w(i_vp+1,node,i_mu+1,isp)+w1w0*t11/wp(1)

    f0_diag_ws(i_vp+0,node,i_mu+0,isp)=f0_diag_ws(i_vp+0,node,i_mu+0,isp)+w1w0*w1w0*t00/wp(0)
    f0_diag_ws(i_vp+1,node,i_mu+0,isp)=f0_diag_ws(i_vp+1,node,i_mu+0,isp)+w1w0*w1w0*t01/wp(0)
    f0_diag_ws(i_vp+0,node,i_mu+1,isp)=f0_diag_ws(i_vp+0,node,i_mu+1,isp)+w1w0*w1w0*t10/wp(1)
    f0_diag_ws(i_vp+1,node,i_mu+1,isp)=f0_diag_ws(i_vp+1,node,i_mu+1,isp)+w1w0*w1w0*t11/wp(1)

#endif



    iflag=1

  end subroutine f0_update_f0g

  subroutine reset_f0_f
    implicit none

    f0_f=0D0
    f0_n=0D0

#ifdef DIAG_NOISE
    f0_diag_n=0D0
    f0_diag_w0=0D0
    f0_diag_w0s=0D0
    f0_diag_w=0D0
    f0_diag_ws=0D0
#endif

  end subroutine reset_f0_f

  ! add
!  subroutine updata_f_grid_part  ! should be called before update_f0_sp -- otherwise, double counted.
!    implicit none

!    f0_f = f0_f + f0_f0g

!  end subroutine updata_f_grid_part

  ! Determine the intersection between node subsets {l1,r1} and {l2,r2} with a 
  ! psize ordered (1D) partition of {1,nnodes}. Return index of first node and total 
  ! number of nodes in intersection for each partition element.
  ! Note: Because of the construction of intervals in calling routine, r1 and l2 can 
  ! not be in the same partition element (there is always one intervening partition 
  ! element that neither are in) so do not have to worry about that case.
  subroutine partition_intersect(l1, r1, l2, r2, nnodes, psize, partition, index, count, counts, ierror)
  integer, intent(in) :: l1, r1, l2, r2, nnodes, psize
  integer, intent(in) :: partition(psize+1)
  integer, intent(out) :: index(psize)
  integer, intent(out) :: count(psize)
  integer, intent(out) :: counts
  integer, intent(out) :: ierror

  ! local variables
  integer :: a, b, c, d, i

  index(:) = nnodes
  count(:) = 0
  counts   = 0
  ierror   = 0

  ! determine the intersection of the first interval with the partition.
  if (r1 >= l1) then

    if (l1 < 1) then
      write(6,*) "error in partition_intersect (l1 < 1):", l1
      ierror = 1
      return
    endif

    if (l1 > nnodes) then
      write(6,*) "error in partition_intersect (l1 > nnodes):", l1, nnodes
      ierror = 1
      return
    endif

    do i=1,psize
      if (l1 < partition(i+1)) then
        a = i
        exit
      endif
    enddo  

    if (r1 > nnodes) then
      write(6,*) "error in partition_intersect (r1 > nnodes):", r1, nnodes
      ierror = 1
      return
    endif

    do i=a,psize
      if (r1 < partition(i+1)) then
        b = i
        exit
      endif
    enddo  

  else

    a = psize+1
    b = 1

  endif

  ! determine the intersection of the second interval with the partition
  if (r2 >= l2) then

    if (l2 < 1) then
      write(6,*) "error in partition_intersect (l2 > 1):", l2
      ierror = 1
      return
    endif

    if (l2 > nnodes) then
      write(6,*) "error in partition_intersect (l2 > nnodes):", l2, nnodes
      ierror = 1
      return
    endif

    do i=b,psize
      if (l2 < partition(i+1)) then
        c = i
        exit
      endif
    enddo  

    if (r2 > nnodes) then
      write(6,*) "error in partition_intersect (r2 > nnodes):", r2, nnodes
      ierror = 1
      return
    endif

    do i=c,psize
      if (r2 < partition(i+1)) then
        d = i
        exit
      endif
    enddo  

  else

    c = psize+1
    d = 1

  endif

  ! Verify that r1 and l2 are not in the same partition.
  if (b == c) then
    write(6,*) "error in partition_intersect (b == c)"
    ierror = 1
    return
  endif

  ! calculate index and count contributions from the first interval
  if (a <= b) then
    if (a == b) then
      index(a) = l1
      count(a) = r1 - l1 +1
      counts   = counts + 1
    else
      index(a) = l1
      count(a) = partition(a+1) - l1
      counts   = counts + 1

      do i=a+1,b-1
        count(i) = partition(i+1) - partition(i)
        index(i) = partition(i)
        counts   = counts + 1
      enddo

      index(b) = partition(b)
      count(b) = r1 - partition(b) + 1
      counts   = counts + 1
    endif
  endif

  ! calculate index and count contributions from the second interval
  if (c <= d) then
    if (c == d) then
      index(c) = l2
      count(c) = r2 - l2 +1
      counts   = counts + 1
    else
      index(c) = l2
      count(c) = partition(c+1) - l2
      counts   = counts + 1

      do i=c+1,d-1
        index(i) = partition(i)
        count(i) = partition(i+1) - partition(i)
        counts   = counts + 1
      enddo

      index(d) = partition(d)
      count(d) = r2 - partition(d) + 1
      counts   = counts + 1
    endif
  endif

  end subroutine partition_intersect

  !------------------------------------------------------------------------------
  !
  !------------------------------------------------------------------------------
  !BOP
  ! !ROUTINE: ceil2
  !
  ! !INTERFACE:
  integer function ceil2(n)
    !
    ! !INPUT PARAMETERS:
    integer :: n
    ! !DESCRIPTION:
    !
    !     Smallest power of 2 greater than or equal to the argument
    !
    ! !REVISION HISTORY: 
    !    2008.08.21   Worley         Imported from spmdutils in
    !                                Community Atmosphere Model
    !
    !EOP
    !------------------------------------------------------------------------------
    !BOC
    !
    ! !LOCAL VARIABLES:
    integer p
    
    p=1
    do while ( p < n )
       p=p*2
    enddo
    ceil2=p
    
    return
    !EOC
  end function ceil2

  !------------------------------------------------------------------------------
  !
  !------------------------------------------------------------------------------
  !BOP
  ! !ROUTINE: pair 
  !
  ! !INTERFACE:
  integer function pair(np,p,k)
    !
    ! !INPUT PARAMETERS:
    integer :: np
    integer :: p
    integer :: k
    ! !DESCRIPTION:
    !
    !     Bitwise XOR of arguments p and k, if less than upper bound np
    !
    ! !REVISION HISTORY: 
    !    2008.08.21   Worley         Imported from spmdutils in 
    !                                Community Atmosphere Model
    !
    !EOP
    !------------------------------------------------------------------------------
    !BOC
    !
    ! !LOCAL VARIABLES:
    integer q
    !
    q = ieor(p,k)
    if ( q > np-1 ) then
       pair = -1
    else
       pair = q
    endif
    
    return
    
    !EOC
  end function pair

end module f0_module


subroutine update_f0_sp(grid,sp)
  use sml_module
  use grid_class
  use ptl_module
  use eq_module
  use f0_module
  use perf_monitor
  implicit none
  type(grid_type), intent(in) :: grid
  type(species_type) :: sp
  type(eq_ftn_type) :: temp_ftn
  integer :: i
  real (8) :: r, z, phi, b, psi, temp, mu_n, vp_n, w, alpha, dw, df, w0, df2, temp_g, psi_g
  integer :: iflag
  real (8) , external :: psi_interpol, b_interpol
  real (8), allocatable :: send1_g(:,:), send2_g(:,:), send1_f(:,:), send2_f(:,:)
  integer :: node, ml(1)

  !copying structure will not take much time
  if(sp%type==0) then
     temp_ftn=eq_tempe
  else
     temp_ftn=eq_tempi
  endif

! The following arrays are not used when the commented if-statemens
! are ture, but using those in subroutine arguments without allocation
! may cause error messages
!
!  if(f0_imu1>0) then ! send lower part
     allocate(send1_g(-f0_nvp:f0_nvp,f0_inode1:f0_inode2))
     allocate(send1_f(-f0_nvp:f0_nvp,f0_inode1:f0_inode2))
     send1_g=0D0
     send1_f=0D0
!  endif

!  if(f0_imu2<f0_nmu) then
     allocate(send2_g(-f0_nvp:f0_nvp,f0_inode1:f0_inode2))
     allocate(send2_f(-f0_nvp:f0_nvp,f0_inode1:f0_inode2))
     send2_g=0D0
     send2_f=0D0
!  endif

#ifdef F0_TOR_LINEAR
    f0_df0g3(:,:,:,:,sp%type)=0D0
#ifdef F_USE_MARKER_DEN2
    f0_n(:,:,:,:,sp%type)=0D0
#else
    f0_n(:,:,:,sp%type)=0D0
#endif
#endif

  call t_startf("UPDATE_F0_SP_LOOP")
  ! omp need?
  ! how to use OpenMP ?

  do i=1, sp%num
     if(sp%ptl(i)%gid<=0) cycle

     if( sp%tr_save(i)>0 ) then

        r=sp%ptl(i)%ph(pir)
        z=sp%ptl(i)%ph(piz)
        phi=sp%ptl(i)%ph(pip)
        B=b_interpol(r,z,phi)
        ! psi
        psi=psi_interpol(r,z,0,0)

        ! T
        temp=eq_ftn(psi,r,z,temp_ftn)*sml_ev2j

        ! get normalized quantity, mu_n and vp_n, from T


#ifndef V_PERP
        ml=maxloc(sp%p_save(:,i))
        node= grid%nd(ml(1),sp%tr_save(i))
        mu_n=sp%ptl(i)%ct(pim)*(2.*eq_axis_b*f0_B_B0(node)/temp)
#else
        ! mu_n becomes (v_perp/v_th)^2 here --->
        mu_n=sp%ptl(i)%ct(pim)*(2.*b/temp)
#endif
        vp_n=ptl_c_m(sp%type)*sp%ptl(i)%ph(pirho)*B/sqrt(temp/ptl_mass(sp%type))


        if(mu_n <= f0_smu_max*f0_smu_max .and. vp_n <= f0_vp_max .and. vp_n >= -f0_vp_max) then

           ! determine alpha from weight (+ position and velocity)
           w=sp%ptl(i)%ph(piw1)
           w0=sp%ptl(i)%ct(piw0)

           !rh Calculate dynamic alpha by taking into account the 
           !rh particles toroidal transit time --->
           !if (sml_gstep .ge. sml_f0_grid_alpha_start) then
             !vpara = abs(ptl_c_m(sp%type)*sp%ptl(i)%ph(pirho)*B)
             !if (sml_inpsi .le. psi .and. psi .lt. sml_outpsi .and. psi .lt. eq_x_psi) then
             !  ! alpha = transit_time/sml_dt
             !  ! determine local safety factor and R_major
             !  dum=(psi-sml_inpsi)*real(sml_00_npsi-1,8)/(sml_outpsi-sml_inpsi)
             !  dum=max(0D0,min(real(sml_00_npsi-1,8),dum))
             !  j=min(floor(dum)+1,sml_00_npsi)
             !  bb=dum-real(j-1,8)
             !  aa=1D0-bb
             !  qloc = aa*grid%qsafety(j) + bb*grid%qsafety(j+1)
             !  Rmaj = aa*grid%epspar(j,2) + bb*grid%epspar(j+1,2)
             !  alpha = (vpara*sml_dt)/(sml_2pi*qloc*Rmaj)
             !else
             !  alpha = (vpara*sml_dt)/(sml_2pi*eq_axis_r)
             !endif
             ! Limit alpha to something reasonable
             ! Actually, if this truncation is necessary, the time step
             ! may be too large anyway...
             !alpha = min(alpha,sml_f0_grid_alpha)
             !alpha = sml_dt * sqrt(temp*sml_j2ev/200D0)
             !alpha = sml_f0_grid_alpha ! + sml_f0_alpha2*(1D0 - exp(- w*w)) ! hard coded -- change later
           !else
           !  alpha=0D0
           !endif

           alpha = sml_f0_grid_alpha * min(real(sml_gstep,8),real(sml_f0_grid_alpha_start,8)) &
                /real(sml_f0_grid_alpha_start,8)
           !alpha = sml_f0_grid_alpha ! + sml_f0_alpha2*(1D0 - exp(- w*w)) ! hard coded -- change later

           !
           dw= alpha*w
           dw = min(max(dw,-alpha),alpha) ! if |w1|>1 , make |dw|=alpha
           sp%ptl(i)%ph(piw1)=w-dw   !! account for weight transferred to the grid
#ifdef F_USE_MARKER_DEN2
           df=dw*sp%ptl(i)%ct(pif0)
#else
           df=dw*w0
#endif
           if(mod(sml_istep, sml_f_source_period)==0) then
#ifdef F_USE_MARKER_DEN
              df2=w*sp%ptl(i)%ct(pif0)
#else
              df2=w*w0
#endif
           else
              df2=0D0
           endif

           !if(df/=0D0 .or. df2/=0D0) then  !for performance
#ifndef DIAG_NOISE
              call f0_update_f0g(grid,sp%type,sp%tr_save(i),sp%p_save(:,i),&
                   &phi,mu_n,vp_n,df,df2,send1_g, send2_g, send1_f, send2_f,iflag)
#else
              call f0_update_f0g(grid,sp%type,sp%tr_save(i),sp%p_save(:,i),&
                   &phi,mu_n,vp_n,df,df2,send1_g, send2_g, send1_f, send2_f,iflag, w0, w*w0)
#endif
              if(iflag==-1) then
                 print *, 'out of range mu_n or vp_p', mu_n, vp_n  !something wrong
              else
                 ! Does not need.
                 ! time variation of f will be captured in next charge operation
                 !sp%ptl(i)%ph(piw1)=sp%ptl(i)%ph(piw1) - dw
                 !sp%ptl(i)%ph(piw2)=sp%ptl(i)%ph(piw2) - dw
              endif
           !endif
           
        endif

     endif
  enddo
  call t_stopf("UPDATE_F0_SP_LOOP")


#ifdef F0_TOR_LINEAR

  ! ### this is simple form. 
  ! ## to be precise, f0_df0g3 need to be sent to next planes and ff transformed to be added.
  do node=f0_inode1, f0_inode2
#ifndef F_USE_MARKER_DEN2
    f0_f0g(:,node,:,:,sp%type) = f0_f0g(:,node,:,:,sp%type) + f0_df0g3(:,:,node,:,sp%type)
#else
    f0_f0g(:,node,:,:,sp%type) = f0_f0g(:,node,:,:,sp%type) + f0_df0g3(:,:,node,:,sp%type)/(f0_n(:,node,:,:,sp%type)+1D-10)
#endif
  enddo


#endif




! following code block is not used -- disabled for now.
  ! send-receive mu boundary data
!  if(f0_imu1 > 0) then
!    call t_startf("UPDATE_F0_SP_SR1")
!    call send_receive_mu_boundary(f0_f0g(:,:,:,sp%type),send1_g,send2_g)
!    call t_stopf("UPDATE_F0_SP_SR1")
!  endif
!  if(f0_imu2 < f0_nmu) then
!    call t_startf("UPDATE_F0_SP_SR2")
!    call send_receive_mu_boundary(f0_f(:,:,:,sp%type),send1_f,send2_f)
!    call t_stopf("UPDATE_F0_SP_SR2")
!  endif

  if(allocated(send1_g)) deallocate(send1_g, send1_f)
  if(allocated(send2_g)) deallocate(send2_g, send2_f)

#ifdef F0_CHARGE_N0
  ! charge update
  !#### THIS HAS TIME DEALY OF CHARGE DENSITY 
  !#### CHARGE DENSITY is calculated from previous particle wieghts and position
  !#### assuming that n=0 mode does not have much time variation
  !-- alpha should be constant over particle

  f0_density_n0_add(:,sp%type) = f0_density_n0_add(:,sp%type) + alpha * f0_den0_ptl(:,sp%type)

#endif

contains
  subroutine send_receive_mu_boundary(array, send1, send2)
    implicit none
    real (8) :: array(-f0_nvp:f0_nvp, f0_inode1:f0_inode2, f0_imu1:f0_imu2)
    real (8), intent(in) :: send1(-f0_nvp:f0_nvp,f0_inode1:f0_inode2)
    real (8), intent(in) :: send2(-f0_nvp:f0_nvp,f0_inode1:f0_inode2)
    real (8), allocatable :: recv1(:,:),recv2(:,:)
    integer :: rreq1, rreq2, sreq1, sreq2
    integer :: isp, icount
    integer :: src1, dest1, recv_tag1, send_tag1
    integer :: src2, dest2, recv_tag2, send_tag2
    include 'mpif.h'
    integer :: istatus(MPI_STATUS_SIZE), ierr

    isp=sp%type
    icount=(2*f0_nvp+1)*(f0_inode2-f0_inode1+1)

    !prepare memory and  receive first
    if(f0_imu1>0) then !recv lower part
       allocate(recv1(-f0_nvp:f0_nvp,f0_inode1:f0_inode2))
       ! receive
       src1=sml_mype - sml_f0_nmu_decomp
       recv_tag1=src1
       call mpi_irecv(recv1,icount, MPI_REAL8, src1, recv_tag1,&
            sml_comm,rreq1,ierr)
    endif


    if(f0_imu2<f0_nmu) then
       allocate(recv2(-f0_nvp:f0_nvp,f0_inode1:f0_inode2))
       !receive
       src2=sml_mype + sml_f0_nmu_decomp
       recv_tag2=src2
       call mpi_irecv(recv2,icount, MPI_REAL8, src2, recv_tag2,&
            sml_comm,rreq2,ierr)
    endif


    ! send
    if(f0_imu1>0) then

       !send1 = array(:,:,f0_imu1)
       dest1 = sml_mype - sml_f0_nmu_decomp
       send_tag1=sml_mype

       call mpi_isend(send1,icount,MPI_REAL8,dest1,send_tag1,&
            sml_comm,sreq1,ierr)
    endif

    if(f0_imu2<f0_nmu) then

       !send2 = array(:,:,f0_imu2)
       dest2 = sml_mype + sml_f0_nmu_decomp
       send_tag1=sml_mype

       call mpi_isend(send2,icount,MPI_REAL8,dest2,send_tag2,&
            sml_comm,sreq2,ierr)
    endif


    if(f0_imu1>0) then
       call mpi_wait(rreq1,istatus,ierr)
       array(:,:,f0_imu1)= array(:,:,f0_imu1) + recv1
    endif
    if(f0_imu2<f0_nmu) then
       call mpi_wait(rreq2,istatus,ierr)
       array(:,:,f0_imu2)= array(:,:,f0_imu2) + recv2
    endif

    if(f0_imu1>0)  call mpi_wait(sreq1,istatus,ierr)
    if(f0_imu2<f0_nmu) call mpi_wait(sreq2,istatus,ierr)

    if(allocated(recv1)) deallocate(recv1)
    if(allocated(recv2)) deallocate(recv2)

  end subroutine send_receive_mu_boundary

end subroutine update_f0_sp

subroutine add_f0_analytic(stype)
  use f0_module
  implicit none
  integer, intent(in) :: stype
  integer :: imu, node, ivp
  real (8) :: smu,mu, vp, en


  do imu=f0_imu1, f0_imu2
     do node=f0_inode1, f0_inode2
        do ivp=-f0_nvp, f0_nvp
           ! normalized
           smu=(imu*f0_dsmu) !  v_perp (normalized)
           mu=smu*smu         ! now mu is v_perp^2

           vp=ivp*f0_dvp

           if(imu==0) smu=f0_dsmu/f0_mu0_factor
           en=0.5D0 * (vp*vp + mu)
           ! 'smu' is v_perp/v_th
           f0_f(ivp,node,imu,stype) = f0_f(ivp,node,imu,stype)  + f0_n_Ta(node,stype)*exp(-en)*smu

        enddo
     enddo
  enddo

end subroutine add_f0_analytic

!electron f0 - adiabatic response added
subroutine add_f0_analytic_elec(grid,psn)
  use f0_module
  use grid_class
  use psn_class
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
!  integer, intent(in) :: stype
  integer, parameter :: stype=0
  integer :: imu, node, ivp
  real (8) :: smu, mu, vp, en, pot


  do imu=f0_imu1, f0_imu2
     do node=f0_inode1, f0_inode2
        !normalized potential
        pot=0.5D0*(psn%dpot_ff(node,0)+psn%dpot_ff(node,1))/f0_t_ev(node,stype)
        do ivp=-f0_nvp, f0_nvp
           ! normalized
           smu=(imu*f0_dsmu) ! v_perp

           mu=smu*smu         ! now mu is (v_perp/v_th)^2

           vp=ivp*f0_dvp

           if(imu==0) smu=f0_dsmu/f0_mu0_factor


           en=0.5D0 * (vp*vp +mu)
           f0_f(ivp,node,imu,stype) = f0_f(ivp,node,imu,stype)  + f0_n_Ta(node,stype)*exp(-en+pot)*smu
        enddo
     enddo
  enddo

end subroutine add_f0_analytic_elec



! non hamiltonian term of kinetic equation
subroutine f_source(grid, psn, spall)
  use sml_module
  use grid_class
  use psn_class
  use ptl_module
  use f0_module
  use perf_monitor
  use col_module, only : col_mode
  use src_module, only : src_narea, src_narea_e
  implicit none
  type(grid_type), intent(in) :: grid
  type(psn_type), intent(in) :: psn
  type(species_type), intent(inout) :: spall(0:ptl_nsp_max)
  !debug
  integer :: node, st
  integer :: index_diag_f0_df

   interface
     subroutine diag_f0(istep,grid,psn,flag)
       use sml_module
       use diag_module
       use grid_class
       use psn_class
       use f0_module
       implicit none
       integer, intent(in) :: istep
       type(grid_type) :: grid
       type(psn_type) :: psn
       integer, intent(in), optional :: flag ! flag is for debug purpose
     end subroutine diag_f0
   end interface

  call t_startf("F_SOURCE")

#ifdef COL_F_NAN_CHECK
  call t_startf("F_NAN_CHECK")
  st=1
  do node=f0_inode1, f0_inode2
     if(.not. (f0_f(2,node,1,st) > 1D0 .or. f0_f(2, node, 1, st) < 2D0) ) then
        print *, sml_mype, '] NAN FOUND f0_f (f_source)', f0_f(2, node,1, st), node
        stop
     endif
  enddo
  call t_stopf("F_NAN_CHECK")
#endif

#ifndef DELTAF_MODE2
  call t_startf("F_UPD_W_PTL")
  ! update w1 and w2 - which is from particle motion
  call update_w_ion(spall(1),1)
  if(sml_electron_on) call update_w_elec(spall(0),1)
  call t_stopf("F_UPD_W_PTL")
#endif

  ! Update f0_f0g and f0_f
  call t_startf("F_UPD_F0_SP")
  call update_f0_sp(grid,spall(1))
  if (sml_electron_on) call update_f0_sp(grid,spall(0))
  call t_stopf("F_UPD_F0_SP")


#ifdef COL_F_NAN_CHECK
  call t_startf("F_NAN_CHECK")
  do node=f0_inode1, f0_inode2
     if(.not. (f0_f(2,node,1,st) > 1D0 .or. f0_f(2, node, 1, st) < 2D0) ) then
        print *, sml_mype, '] NAN FOUND f0_f (after update_w_ion)', f0_f(2, node,1, st), node
        stop
     endif
  enddo
  call t_stopf("F_NAN_CHECK")
#endif

  ! This statement was in main.F90 previously
  ! but had to be shifted to this position
  ! for the gradual update of f0_f0g.
  if(mod(sml_istep,sml_f_source_period)==0) then
     call t_startf("F_SOURCE_FIRST_PART")
     ! Reset f0_df0g
     f0_df0g=0.D0

#ifdef F_USE_MARKER_DEN
     f0_n=f0_n+1D-10  ! to avoid nan
#ifndef F_USE_MARKER_DEN2
     f0_f=f0_f/f0_n   ! average operation 
#else
     f0_f=f0_f/(f0_n(:,:,:,0,:)+f0_n(:,:,:,1,:))   ! average operation
#endif
#endif
     
     !rh Add f0_f0g to f0_f to make f0_f the full distribution function
#ifndef F0_TOR_LINEAR
     f0_f=f0_f+f0_f0g
#else
     f0_f=f0_f+(f0_f0g(:,:,:,0,:)+f0_f0g(:,:,:,1,:))*0.5D0
#endif
     call t_stopf("F_SOURCE_FIRST_PART")


     ! add maxwellian (f0_a) to make whole f
     call t_startf("ADD_F0_ANALYTIC")
     if(sml_electron_on) call add_f0_analytic_elec(grid,psn) ! sp%type==0 is assumed in the subroutine
     call add_f0_analytic(spall(1)%type)
     call t_stopf("ADD_F0_ANALYTIC")

#ifdef COL_F_NAN_CHECK
     call t_startf("F_NAN_CHECK")
     do node=f0_inode1, f0_inode2
        if(.not. (f0_f(2,node,1,st) > 1D0 .or. f0_f(2, node, 1, st) < 2D0) ) then
           print *, sml_mype, '] NAN FOUND f0_f (after f_heat_torque)', f0_f(2, node,1, st), node
           stop
        endif
     enddo
     call t_stopf("F_NAN_CHECK")
#endif

     call t_startf("DIAG_F0")
     call diag_f0(sml_istep,grid,psn,1) ! 1 for safety. to avoid compiler bug.
     call t_stopf("DIAG_F0")
     call t_startf("DIAG_3D_F0_F")
     call diag_3d_f0_f(grid,psn)
     call t_stopf("DIAG_3D_F0_F")

     !Make symmetric f0_f
     if(sml_symmetric_f) then
        call t_startf("SYMMETRIC_F")
        call symmetric_f   !symmetic f0_n??
        call t_stopf("SYMMETRIC_F")
     endif


     ! collision
     ! Collision comes first since f_collision will have additional parallelism when f0_f is symmetric. --> df0g=0
     if(col_mode==4) then
        call t_startf("F_COLLISON")
        call f_collision(grid)
        call t_stopf("F_COLLISON")

        call t_startf("DIAG_F0_DF_PORT1_1")
        !diagnosis of collisional change of density, energy and momentum of each species
        index_diag_f0_df=1
        call diag_f0_df_port1(index_diag_f0_df,grid,psn)
        call t_stopf("DIAG_F0_DF_PORT1_1")


     endif
#ifdef COL_F_NAN_CHECK
     call t_startf("F_DF0G_NAN_CHECK")
     call df0g_nan_check('After_collision')
     call t_stopf("F_DF0G_NAN_CHECK")
#endif

     ! source
     if(sml_source) then
        call t_startf("F0_HEAT_TORQUE")
        if(sml_electron_on) then
#ifdef HEAT_TORQUE_MAXWELL
           call f_heat_torque_maxwellian(spall(0), psn%tempe_ev,src_narea_e)
#else
           call f_heat_torque_org(spall(0), psn%tempe_ev,src_narea_e)
#endif
        endif
#ifdef HEAT_TORQUE_MAXWELL
        call f_heat_torque_maxwellian(spall(1), psn%tempi_ev, src_narea)
#else
        call f_heat_torque_org(spall(1), psn%tempi_ev, src_narea)
#endif
        call t_stopf("F0_HEAT_TORQUE")

        call t_startf("DIAG_F0_DF_PORT1_2")
        index_diag_f0_df=2
        call diag_f0_df_port1(index_diag_f0_df,grid,psn)
        call t_stopf("DIAG_F0_DF_PORT1_2")
     endif

#ifdef COL_F_NAN_CHECK
     call t_startf("F_NAN_CHECK")
     do node=f0_inode1, f0_inode2
        if(.not. (f0_f(2,node,1,st) > 1D0 .or. f0_f(2, node, 1, st) < 2D0) ) then
           print *, sml_mype, '] NAN FOUND f0_f (after f_heat_torque', f0_f(2, node,1, st), node
           stop
        endif
     enddo
     call t_stopf("F_NAN_CHECK")
#endif

     
     ! neutral
     if(sml_neutral .and. sml_gstep > sml_neutral_start_step) then
      ! call check_point('neutral')
        call t_startf("F_NEUTRAL3")
        call neutral3(sml_istep,grid,psn,spall(0))
        call t_stopf("F_NEUTRAL3")

        call t_startf("DIAG_F0_DF_PORT1_3")
        index_diag_f0_df=3
        call diag_f0_df_port1(index_diag_f0_df,grid,psn)
        call t_stopf("DIAG_F0_DF_PORT1_3")
     endif

#ifdef COL_F_NAN_CHECK
     call t_startf("F_DF0G_NAN_CHECK")
     call df0g_nan_check('After_neutral')
     call t_stopf("F_DF0G_NAN_CHECK")
#endif

     !radiation colloing of electron
     if(sml_electron_on .and. sml_radiation) then
        call t_startf("F_RADIATION")
        call f_radiation(grid)
        call t_stopf("F_RADIATION")

        call t_startf("DIAG_F0_DF_PORT1_4")
        index_diag_f0_df=4
        call diag_f0_df_port1(index_diag_f0_df,grid,psn)
        call t_stopf("DIAG_F0_DF_PORT1_4")
     endif

#ifdef COL_F_NAN_CHECK
     call t_startf("F_DF0G_NAN_CHECK")
     call df0g_nan_check('After_radiation')
     call t_stopf("F_DF0G_NAN_CHECK")
#endif

     ! NaN remove
     call t_startf("F_DF0G_NAN_REMOVE")
     call df0g_nan_remove(f0_df0g)
     if(f0_col_change_weight .and. col_mode==4) call df0g_nan_remove(f0_df0g2)
     call t_stopf("F_DF0G_NAN_REMOVE")

     
     ! set zero df0g for near wall nodes
     call t_startf("DF0G_NEAR_WALL_REMOVE")
     call df0g_near_wall_remove
     call t_stopf("DF0G_NEAR_WALL_REMOVE")

     !******************************************************
     !*****    Apply change of f                      ******
     !******************************************************

     ! apply df0g on particle -- for all f0g modification
     if(f0_col_change_weight) then
        f0_df0g2=f0_df0g
        f0_df0g=0D0   ! f0_df0g - f0_df0g2  = 0
     endif

     ! on particle
    if (f0_col_change_weight) then
      ! Change particle weights for every phase space cell with a particle in it
      ! instead of leaving the collisional information on the grid (f0_df0g)
      call t_startf("distribute_f0g_i")
      call distribute_f0g(grid,spall(1))
      call t_stopf("distribute_f0g_i")
      if (sml_electron_on) then
         call t_startf("distribute_f0g_e")
         call distribute_f0g(grid,spall(0))
         call t_stopf("distribute_f0g_e")
      endif
    endif

    ! Update f0_f0g abruptly
    call t_startf("UPDATE_F0_F0G")
#ifndef F0_TOR_LINEAR
    f0_f0g=f0_f0g+f0_df0g
#else
    !when f0_col_change_weight is true,  f0_df0g is zero.
    !if(f0_col_change_weight==.false.) then
    !    print *, 'Not implemented with F0_TOR_LINEAR'
    !    stop
    !endif

    !rather simple way
    f0_f0g(:,:,:,0,:)=f0_f0g(:,:,:,0,:)+f0_df0g
    f0_f0g(:,:,:,1,:)=f0_f0g(:,:,:,1,:)+f0_df0g

#endif
    call t_stopf("UPDATE_F0_F0G")


    ! write diagnosis file of f_source
    call t_startf("DIAG_F0_DF")
    call diag_f0_df(sml_istep,grid,psn)
    call t_stopf("DIAG_F0_DF")

  endif


#ifndef DELTAF_MODE2
  ! update w2 - which reflect the change of f0
  call t_startf("UPDATE_W_ION_ELEC2")
  call update_w_ion(spall(1),2)
  if(sml_electron_on) call update_w_elec(spall(0),2)
  call t_stopf("UPDATE_W_ION_ELEC2")
#endif

  call t_stopf("F_SOURCE")
contains
  subroutine update_w_ion(sp,iflag)
    use omp_module
    implicit none
    type(species_type) :: sp
    integer, intent(in) :: iflag
    integer :: ith, i_beg(sml_nthreads), i_end(sml_nthreads)
    integer :: i
    real (8) ::  new_f0, w2, dw
    real (8), external :: get_f0_ion

    call split_indices(sp%num, sml_nthreads, i_beg, i_end)

!$OMP PARALLEL DO &
!$OMP PRIVATE( ITH, I, NEW_F0, W2, DW )
    do ith=1, sml_nthreads
       do i=i_beg(ith), i_end(ith)
          if(sp%ptl(i)%gid>0) then

             new_f0 = get_f0_ion(grid,sp%ptl(i),sp%tr_save(i),sp%p_save(:,i),sp%type)
             w2= 1D0 - new_f0/sp%ptl(i)%ct(pif0)
             if(iflag==1) then
                dw=w2 - sp%ptl(i)%ph(piw2)
                sp%ptl(i)%ph(piw1)= sp%ptl(i)%ph(piw1) + dw
             endif
             sp%ptl(i)%ph(piw2)= w2

          endif
       enddo
    enddo
  end subroutine update_w_ion

  subroutine update_w_elec(sp,iflag)
    use omp_module
    implicit none
    type(species_type) :: sp
    integer, intent(in) :: iflag
    integer :: ith, i_beg(sml_nthreads), i_end(sml_nthreads)
    integer :: i
    real (8) ::  new_f0, w2, dw
    real (8), external :: get_f0_elec

    call split_indices(sp%num, sml_nthreads, i_beg, i_end)

!$OMP PARALLEL DO &
!$OMP PRIVATE( ITH, I, NEW_F0, W2, DW )
    do ith=1, sml_nthreads
       do i=i_beg(ith), i_end(ith)
          if(sp%ptl(i)%gid>0) then

             new_f0 = get_f0_elec(grid,psn,sp%ptl(i),sp%tr_save(i),sp%p_save(:,i))
             w2= 1D0 - new_f0/sp%ptl(i)%ct(pif0)
             if(iflag==1) then
                dw=w2 - sp%ptl(i)%ph(piw2)
                sp%ptl(i)%ph(piw1)= sp%ptl(i)%ph(piw1) + dw
             endif
             sp%ptl(i)%ph(piw2)= w2

          endif
       enddo
    enddo
  end subroutine update_w_elec


  ! apply heat and torque  -- obsolete code -- remove later
  subroutine f_heat_torque_org(sp, t_ev,narea)
    use src_module
    use omp_module, only : split_indices
    implicit none
    type(species_type) :: sp
    real(8), intent(in) :: t_ev(grid%nnode)
    integer, intent(in) :: narea
    !
    real (8) :: c_m, mass
    real (8), dimension(narea) :: pimw, popw, pi, po, width, alpha, beta, beta2, del_en, del_mo
    real (8), dimension(narea) :: heat_power, torque
    integer, dimension(narea) :: dmu!, dvp
    integer :: idvp
    integer, parameter :: N=3
    real (8) :: l_sum(N,narea,sml_nthreads), sumall(N,narea), v2(N), vsum(N)
    integer :: imu, ivp, node, j, stype
    real (8) :: r, r_bphi_B, psi, mu, dmu2, E_perp, E_para, E, factor, vp, vth, tmp
    logical :: rgn_check
    integer :: ith
    real (8) , allocatable :: df(:,:,:)
    integer :: imu1, imu2, itmp, itmp_mu, itmp_vp
    real (8) :: epsilon1, epsilon2, f, mu_vol

    real (8), parameter :: alpha_minimum=-0.6

    !only when narea is positive
    if(narea<=0) return


    c_m=sp%c_m
    mass=sp%mass
    stype=sp%type

    if(stype==0) then
      heat_power=src_heat_power_e
      torque=src_torque_e
      pi = src_pin_e
      po = src_pout_e
      width = src_decay_width_e
    else
      heat_power=src_heat_power
      torque=src_torque
      pi = src_pin
      po = src_pout
      width = src_decay_width
    endif

    pimw = pi - width
    popw = po + width

    l_sum=1D-90 ! zero + epsilon to avoid NaN



    ! Get total E & nR for each section
!    do ith=1, sml_nthreads
    ith=1
    do node=f0_inode1,f0_inode2

      r=grid%x(1,node) ! get major radius
      rgn_check=(grid%rgn(node)==1 .or. grid%rgn(node)==2)
      !rh We have B_phi and |B| in grid%bfield, why not use it?
      !rh r_bphi_B = r*sml_bt_sign ! bphi/B assumed 1 for now. complete it later for NSTX
      r_bphi_B = r*grid%bfield(3,node)/grid%bfield(4,node)
      psi=grid%psi(node)
      ! simulation boundary check
      if(psi < sml_inpsi .or. sml_outpsi < psi) cycle


      vth=sqrt(t_ev(node)*sml_ev2j/mass)

      vsum=0D0

      do imu=f0_imu1, f0_imu2

        !!rh mu is either norm. mag. moment (sqrt(mu)-grid) or (v_perp/v_th)^2 (v_perp-grid)
        mu=(imu*f0_dsmu)**2 ! crude approximation
        if(imu==0) mu=f0_dsmu*f0_dsmu/f0_mu0_factor ! boundary effect of linear interpolation -->> wrong

        if(imu==0) then
            mu_vol=0.5D0             ! obtain from integrating of interpolation function
        elseif(imu==f0_nmu) then
            mu_vol=0.5D0         ! slightly smaller than half
        else
            mu_vol=1D0
        endif

        E_perp=0.5D0*mu*t_ev(node)*sml_ev2j   !mu is v_perp^2 - normalized

        do ivp=-f0_nvp, f0_nvp
          vp=ivp*f0_dvp  !normalized vp

          E_para=0.5D0*vp*vp* t_ev(node)*sml_ev2j
          E=E_perp+E_para

          v2(1)=E
          v2(2)=mass*r_Bphi_B
          v2(3)=mass*r_Bphi_B*vp*vth

          tmp=f0_f(ivp,node,imu,stype)*f0_grid_vol(node,stype)*mu_vol
          tmp=max(tmp,0D0)  !prevent minus f0_f



          vsum = vsum + tmp * v2

        enddo  !ivp
      enddo  !imu

      do j=1, narea
        if(pimw(j) < psi .AND. psi < popw(j) .and. rgn_check ) then
          if( pi(j) < psi .and. psi < po(j) ) then
            factor=1D0
          elseif( pimw(j) < psi .and. psi < pi(j) ) then
            factor= (psi-pimw(j))/width(j)
          else
            factor= (popw(j)-psi)/width(j)
          endif

          l_sum(:,j,ith)=l_sum(:,j,ith) + factor*vsum(:)
        endif
      enddo

    enddo  !node loop


    ! omp summation
    !do ith=2, sml_nthreads
    !   l_sum(:,:,1)=l_sum(:,:,1)+l_sum(:,:,ith)
    !enddo

    ! mpi summation
    call t_startf("SOURCE_RED")
    call my_mpi_allreduce(l_sum(:,:,1), sumall, N*narea)
    call t_stopf("SOURCE_RED")


    ! 2. find alpha and beta
    del_en(:)=heat_power(:) * sml_dt * real(sml_f_source_period)
    del_mo(:)=torque(:)     * sml_dt * real(sml_f_source_period)



    ! alpha is obtained approximately
    ! beta is more exact
    alpha=sqrt(1D0+del_en/sumall(1,:))  ! alpha is velocity enarlge factor (different from source in heat.F90)
    beta=(del_mo-(alpha-1D0)*sumall(3,:))/sumall(2,:)       ! beta is v-shift amount
    !if(sml_mype==0) print *,'del_en,alpha,beta',del_en(1),alpha,beta


    !3. calculate df (or directly modify f0_f0g

    ! send base loop

    do j=1, narea
       if(alpha(j)>1D0) then  ! heating
          dmu(j)=1
       else
          dmu(j)=-1
       endif

       !if(beta(j)>0D0) then ! positive
       !   dvp(j)=1
       !else
       !   dvp(j)=-1
       !endif
    enddo

    ! add one more grid for boundary communication
    imu1=max(f0_imu1-1, 0)
    imu2=min(f0_imu2+1, f0_nmu)

    allocate(df(-f0_nvp:f0_nvp,f0_inode1:f0_inode2,imu1:imu2))  ! how to optimize the memory access?
    df=0D0 ! initialization


    do imu=f0_imu1, f0_imu2
       do node=f0_inode1,f0_inode2
          rgn_check=(grid%rgn(node)==1 .or. grid%rgn(node)==2)
          psi=grid%psi(node)
          ! simulation boundary check
          if(psi < sml_inpsi .or. sml_outpsi < psi) cycle

          ! beta is dimensional - get beta2
          beta2=beta/sqrt(t_ev(j)*sml_ev2j/mass)


          do j=1, narea
             if(pimw(j) < psi .AND. psi < popw(j) .and. rgn_check ) then
                if( pi(j) < psi .and. psi < po(j) ) then
                   factor=1D0
                elseif( pimw(j) < psi .and. psi < pi(j) ) then
                   factor= (psi-pimw(j))/width(j)
                else
                   factor= (popw(j)-psi)/width(j)
                endif

                ! mu is sqrt(mu)
                ! dmu is delta sqrt(mu)
                if(imu==0) then
                   !rh I believe dmu2 must be 1 - 1/f0_mu0_factor
                   dmu2=1D0-1D0/f0_mu0_factor !rh 2D0/3D0
                   mu=1D0/f0_mu0_factor
                elseif(imu==1 .and. dmu(j)<0 ) then
                   dmu2=1D0-1D0/f0_mu0_factor !rh 2D0/3D0
                   mu=1D0
                else
                   dmu2=1D0
                   mu=imu
                endif



                epsilon1=mu/dmu2*factor*(alpha(j)-1D0)*real(dmu(j))
                epsilon1=min(epsilon1,1D0)                                   ! epsilon is smaller than 1

                do ivp=-f0_nvp, f0_nvp
                   vp=ivp

                   epsilon2=factor*((alpha(j)-1D0)*vp+beta2(j)/f0_dvp)
                   if(epsilon2>0D0) then
                        idvp=1
                   else
                        idvp=-1
                   endif

                   epsilon2=epsilon2*real(idvp) ! epsilon is positive always
                   epsilon2=min(epsilon2,1D0)             ! epsilon is smaller than 1. if it is large, source does not work correctly
                   !epsilon2=0D0


                   f=f0_f(ivp,node,imu,stype)
                   ! 1. shift in vp
!                   itmp=min(max(ivp+dvp(j),-f0_nvp),f0_nvp)     ! boundary of itmp - ignore shift at the vp boundary
!                   df(ivp ,node,imu)=df(ivp, node,imu) -epsilon2*f
!                   df(itmp,node,imu)=df(itmp,node,imu)+epsilon2*f

                   ! 2. shift in mu
!                   itmp=min(max(imu+dmu(j),0),f0_nmu)         ! two dimensional -- maybe use f+df
!                   df(ivp,node,imu )=df(ivp, node,imu )-epsilon1*f
!                   df(ivp,node,itmp)=df(ivp, node,itmp)+epsilon1*f




                   itmp_vp=min(max(ivp+idvp,-f0_nvp),f0_nvp)
                   itmp_mu=min(max(imu+dmu(j),0),f0_nmu)

                   df(ivp,node,imu)=df(ivp,node,imu)+ ((1D0-epsilon1)*(1D0-epsilon2)-1D0)*f
                   df(itmp_vp,node,imu)=df(itmp_vp,node,imu) + ((1D0-epsilon1)*epsilon2)*f
                   df(ivp,node,itmp_mu)=df(ivp,node,itmp_mu) + ((1D0-epsilon2)*epsilon1)*f
                   df(itmp_vp,node,itmp_mu)=df(itmp_vp,node,itmp_mu) + epsilon1*epsilon2*f


                enddo
             endif
          enddo ! j, narea
       enddo !node
    enddo ! imu



    f0_df0g(:,:,:,stype) = f0_df0g(:,:,:,stype) + df(:,:,:)

    deallocate(df)

  end subroutine f_heat_torque_org
! apply heat and torque
! the difference from f_heat_torque_org is that f_heat_Torque_maxwellian assumes f as maxwellian, and calculate df.
! This algorithm is expected to be stable when f is highly distorted. 
  subroutine f_heat_torque_maxwellian(sp, t_ev,narea)
    use src_module
    use omp_module, only : split_indices
    implicit none
    type(species_type) :: sp
    real(8), intent(in) :: t_ev(grid%nnode)
    integer :: narea
    !
    real (8), allocatable :: Tperp(:),Tpara(:),upara(:), number(:)
    !
    real (8) :: c_m, mass
    integer :: stype
    real (8), dimension(narea) :: pi, po, width, pimw, popw, heat_power, torque
    integer, parameter :: N=4
    real (8) :: l_sum(N,narea,sml_nthreads)
    integer :: ith

    integer :: node
    logical :: rgn_check
    real (8) :: r, r_bphi_B, psi, vth
    integer, parameter :: N2=4
    real (8) :: v2(N2), vsum(N2), tmp
    !
    integer :: imu, ivp
    real (8) :: mu, mu_vol, vp, E_perp,E_para
    !
    integer :: j
    real (8) :: factor, v(N)
    !
    real (8) :: sumall(N,narea)
    real (8), dimension(narea) :: alpha, beta, del_en, del_mo
    !
    real (8) , allocatable :: df(:,:,:)
    real (8) :: smu, tmp2
    real (8) :: den, T_perp1, T_para1, u_para1, T_perp2,T_para2, u_para2, E_para1, E_para2

    real (8), parameter :: alpha_minimum=-0.6


    !only when narea is positive
    if(narea<=0) return


    ! memory allocation for temperature and parallel velocity
    allocate(Tperp(f0_inode1:f0_inode2))
    allocate(Tpara(f0_inode1:f0_inode2))
    allocate(upara(f0_inode1:f0_inode2))
    allocate(number(f0_inode1:f0_inode2))



    c_m=sp%c_m
    mass=sp%mass
    stype=sp%type

    if(stype==0) then
      heat_power=src_heat_power_e
      torque=src_torque_e
      pi = src_pin_e
      po = src_pout_e
      width = src_decay_width_e
    else
      heat_power=src_heat_power
      torque=src_torque
      pi = src_pin
      po = src_pout
      width = src_decay_width
    endif

    pimw = pi - width
    popw = po + width

    l_sum=1D-90 ! zero + epsilon to avoid NaN


    ith=1  ! need to make omp loop
    do node=f0_inode1,f0_inode2
       
       r=grid%x(1,node) ! get major radius
       rgn_check=(grid%rgn(node)==1 .or. grid%rgn(node)==2)
       !rh We have B_phi and |B| in grid%bfield, why not use it?
       !rh r_bphi_B = r*sml_bt_sign ! bphi/B assumed 1 for now. complete it later for NSTX
       r_bphi_B = r*grid%bfield(3,node)/grid%bfield(4,node)
       psi=grid%psi(node)
       ! simulation boundary check
       if(psi < sml_inpsi .or. sml_outpsi < psi) cycle


       vth=sqrt(t_ev(node)*sml_ev2j/mass)

       vsum=0D0

       ! get node T_perp, T_para, u_para, particle number
       do imu=f0_imu1, f0_imu2

          mu=(imu*f0_dsmu)**2 ! crude approximation
          !rh is this a BUG ???????
          if(imu==0) mu=f0_dsmu*f0_dsmu/f0_mu0_factor


          if(imu==0) then
             mu_vol=0.5D0             ! obtain from integrating of interpolation function
          elseif(imu==f0_nmu) then
             mu_vol=0.5D0         ! slightly smaller than half
          else
             mu_vol=1D0
          endif


          E_perp=0.5D0*mu*t_ev(node)*sml_ev2j

          do ivp=-f0_nvp, f0_nvp
             vp=ivp*f0_dvp  !normalized vp


             E_para=0.5D0*vp*vp* t_ev(node)*sml_ev2j


             v2(1)=1D0    ! number
             v2(2)=E_perp ! E_perp
             v2(3)=E_para ! E_para
             v2(4)=vp*vth ! v_parallel

             ! tmp is particle number
             tmp=f0_f(ivp,node,imu,stype)*f0_grid_vol_vonly(node,stype)*mu_vol

             tmp=max(tmp,0D0)  !prevent minus f0_f

             vsum = vsum + tmp * v2

          enddo
       enddo

       number(node)=vsum(1)  ! particle number
       Tperp(node)=vsum(2)/vsum(1) ! joule
       upara(node)=vsum(4)/vsum(1)
       Tpara(node)=2D0*vsum(3)/vsum(1) - mass*upara(node)**2

       !sum with

       do j=1, narea
          if(pimw(j) < psi .AND. psi < popw(j) .and. rgn_check ) then
             if( pi(j) < psi .and. psi < po(j) ) then
                factor=1D0
             elseif( pimw(j) < psi .and. psi < pi(j) ) then
                factor= (psi-pimw(j))/width(j)
             else
                factor= (popw(j)-psi)/width(j)
             endif


            ! number(node) is density. volume multiplication is required for number
             tmp = factor * number(node) * grid%node_vol_nearest(node)
             v(1)= Tperp(node)+ Tpara(node)/2D0
             v(2)= mass * upara(node)
             v(3)= mass/2D0
             v(4)= mass*r_Bphi_B

             l_sum(:,j,ith)=l_sum(:,j,ith) + tmp * v

          endif
       enddo

    enddo !do node


    ! omp summation
    !do ith=2, sml_nthreads
    !   l_sum(:,:,1)=l_sum(:,:,1)+l_sum(:,:,ith)
    !enddo

    ! mpi summation
    call t_startf("SOURCE_RED")
    call my_mpi_allreduce(l_sum(:,:,1), sumall, N*narea)
    call t_stopf("SOURCE_RED")



    ! 2. find alpha and beta
    del_en(:)=heat_power(:) * sml_dt * real(sml_f_source_period)
    del_mo(:)=torque(:)     * sml_dt * real(sml_f_source_period)


    ! - obtain beta first
    beta = del_mo/sumall(4,:) ! sumall is rotational inertia
    alpha = (del_en - sumall(2,:)*beta - sumall(3,:)*beta*beta)/sumall(1,:)  ! Temperature loss/gain by troque

    !check minimum alpha -- if alpha < -1, is causes NaN error
    do j=1, narea
      if(alpha(j)<alpha_minimum) then
        if(sml_mype==0) then
          print *, 'Too low alpha in f_heat_torque: Cooling is too strong or something wrong. (j, alpha)=',j,alpha(j)
          print *, 'alpha(j) is enforced to', alpha_minimum
        endif
        alpha(j)=alpha_minimum
      endif

    enddo

    ! 3. calculate df

    allocate(df(-f0_nvp:f0_nvp,f0_inode1:f0_inode2,f0_imu1:f0_imu2))  ! how to optimize the memory access?
    df=0D0 ! initialization


    do node=f0_inode1,f0_inode2

       rgn_check=(grid%rgn(node)==1 .or. grid%rgn(node)==2)
       psi=grid%psi(node)
       ! simulation boundary check
       if(psi < sml_inpsi .or. sml_outpsi < psi) cycle
       vth=sqrt(t_ev(node)*sml_ev2j/mass)

       den = number(node) !bugfix: already density
       T_perp1 = Tperp(node)
       T_para1 = Tpara(node)
       u_para1 = upara(node)


       do imu=f0_imu1, f0_imu2

          smu=imu*f0_dsmu
          mu=(imu*f0_dsmu)**2 ! crude approximation
          !rh is this a BUG ?????????
          if(imu==0) mu=f0_dsmu*f0_dsmu/f0_mu0_factor

          !perp energy in joule
          E_perp=0.5D0*mu*t_ev(node)*sml_ev2j


          do ivp=-f0_nvp, f0_nvp
             vp=ivp*f0_dvp  !normalized vp


             do j=1, narea
                if(pimw(j) < psi .AND. psi < popw(j) .and. rgn_check ) then
                   if( pi(j) < psi .and. psi < po(j) ) then
                      factor=1D0
                   elseif( pimw(j) < psi .and. psi < pi(j) ) then
                      factor= (psi-pimw(j))/width(j)
                   else
                      factor= (popw(j)-psi)/width(j)
                   endif


                   T_perp2 = (1D0+alpha(j)*factor)*T_perp1
                   T_para2 = (1D0+alpha(j)*factor)*T_para1
                   u_para2 = u_para1 + beta(j)*factor
                   !Parallel energy in moving frame (joule)
                   E_para1=0.5D0*mass*(vp*vth-u_para1)**2
                   E_para2=0.5D0*mass*(vp*vth-u_para2)**2

!                   tmp = 1D0/(T_perp1*sml_j2ev)*2D0*smu*sqrt(f0_T_ev(node,stype)*sml_ev2j/T_para1)
!                   tmp2= 1D0/(T_perp2*sml_j2ev)*2D0*smu*sqrt(f0_T_ev(node,stype)*sml_ev2j/T_para2)
                   tmp = 1D0/(T_perp1*sml_j2ev)*smu*sqrt(f0_T_ev(node,stype)*sml_ev2j/T_para1)
                   tmp2= 1D0/(T_perp2*sml_j2ev)*smu*sqrt(f0_T_ev(node,stype)*sml_ev2j/T_para2)


                   df(ivp,node,imu)=df(ivp,node,imu)+  den * &
                        ( tmp2*exp(-E_perp/T_perp2-E_para2/T_para2) &
                        - tmp *exp(-E_perp/T_perp1-E_para1/T_para1) )

                endif
             enddo !j
          enddo !ivp
       enddo !imu
    enddo !node
    
    ! send - recv df at imu boundary
    f0_df0g(:,:,:,stype) = f0_df0g(:,:,:,stype) + df(:,:,:)

    deallocate(df)
    deallocate(Tperp)
    deallocate(Tpara)
    deallocate(upara)
    deallocate(number)

  end subroutine f_heat_torque_maxwellian

  subroutine distribute_f0g(grid,sp)
    ! Uses the result of the collision operator and neutral collisions (f0_df0g)
    ! to update particle weights instead of keeping it on the grid
    use sml_module
    use f0_module
    use eq_module
    use grid_class
    use ptl_module
    use omp_module
    implicit none
    include 'mpif.h'
    type(grid_type), intent(in) :: grid
    type(species_type), intent(inout) :: sp
    type(eq_ftn_type) :: temp_ftn
    real (kind=8) :: r, z, phi, B, temp, psi, smu, wmu(0:1), wvp(0:1), df0g_total1, df0g_total2
    real (kind=8) :: mu_n, vp_n, mu_vol(0:1), t00, t01, t10, t11, wp(0:1), delfrac, dum1, smu0, smu1
    integer :: i, ith, i_beg(sml_nthreads), i_end(sml_nthreads), ml(1), i_mu, i_vp
    integer :: ierr, isize
    real (kind=8), allocatable :: df0g_delete(:,:,:), dum2(:,:,:)
    real (kind=8) :: new_weight, weight, total_weight
    real (kind=8), external :: b_interpol, psi_interpol
    logical, external :: is_nan
!#ifdef NODE_TEMP
    real (kind=8) :: psi_g, temp_g
!#endif

    !copying structure will not take much time
    if(sp%type==0) then
       temp_ftn=eq_tempe
    else
       temp_ftn=eq_tempi
    endif

    allocate(df0g_delete(-f0_nvp:f0_nvp,f0_inode1:f0_inode2,0:f0_nmu))
    allocate(dum2(-f0_nvp:f0_nvp,f0_inode1:f0_inode2,0:f0_nmu))

    df0g_delete=1.D0
    df0g_total1=0D0
    df0g_total2=0D0


!rh    ! Find fraction to delete from f0_df0g if particle is found
       delfrac=1.D0

    dum1=sum(abs(f0_df0g2(:,:,:,sp%type)))

    if(sp%type==0) then
       call t_startf("DIST_F0G_ELEC_RED1")
    else
       call t_startf("DIST_F0G_ION_RED1")
    endif
    call mpi_allreduce(dum1, df0g_total1, 1, MPI_REAL8, MPI_SUM, sml_plane_comm,ierr)
    if(sp%type==0) then
       call t_stopf("DIST_F0G_ELEC_RED1")
    else
       call t_stopf("DIST_F0G_ION_RED1")
    endif

    if (sml_mype==0) then
      print '(a27,e25.10)','Initial sum of |f0_df0g2|: ', df0g_total1
    endif


    call split_indices(sp%num, sml_nthreads, i_beg, i_end)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! particle --> mesh interpolation weights gathered in f0_update_sp -> f0_n
!!!!!!! will be used for reverse interpolation here
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if(sp%type==0) then
       call t_startf("DIST_F0G_ELEC_LOOP")
    else
       call t_startf("DIST_F0G_ION_LOOP")
    endif

!$OMP PARALLEL DO &
!$OMP PRIVATE( ITH, I, PSI, R, Z, B, ML, NODE, TEMP, MU_N, VP_N,    &
!$OMP          PSI_G, TEMP_G, SMU, WMU, WVP, I_MU, I_VP,            &
!$OMP          WEIGHT, MU_VOL, WP, SMU0, SMU1, NEW_WEIGHT,          &
!$OMP          T00, T01, T10, T11)
    do ith=1, sml_nthreads
       do i=i_beg(ith), i_end(ith)
          if(sp%ptl(i)%gid>0) then

             if (sp%tr_save(i) .le. 0) then
               ! Prevent out of mesh particles
               cycle
             endif

             r=sp%ptl(i)%ph(pir)
             z=sp%ptl(i)%ph(piz)
             phi=sp%ptl(i)%ph(pip)
             B=b_interpol(r,z,phi)
             ! node
             ml=maxloc(sp%p_save(:,i))
!#ifdef TRIGRID
             node= grid%nd(ml(1),sp%tr_save(i))
!#else
!             node= sp%cell_save(ml(1),i)
!#endif
             ! psi
             psi=psi_interpol(r,z,0,0)
             ! T
             temp=eq_ftn(psi,r,z,temp_ftn)*sml_ev2j

             ! get normalized quantity, mu_n and vp_n, from T

#ifndef V_PERP
             mu_n=sp%ptl(i)%ct(pim)*(2.*eq_axis_b*f0_B_B0(node)/temp)
#else
             ! mu_n becomes (v_perp/v_th)^2 here --->
             mu_n=sp%ptl(i)%ct(pim)*(2.*b/temp)
#endif
             vp_n=ptl_c_m(sp%type)*sp%ptl(i)%ph(pirho)*B/sqrt(temp/ptl_mass(sp%type))


             smu=sqrt(mu_n)       ! v_perp/v_th else


             wmu(0)=smu/f0_dsmu   ! temp. variable - normalized sqrt(mu) with delta
             i_mu = floor(wmu(0)) ! get index for lower grid
             wmu(1)=wmu(0)-real(i_mu,8) ! wmu(1) is now weight for upper grid
             wmu(0)=1D0 - wmu(1)        ! wmu(0) is now weight for lower grid

             wvp(0)=vp_n/f0_dvp
             i_vp= floor(wvp(0))
             wvp(1)=wvp(0)-real(i_vp,8)
             wvp(0)=1D0-wvp(1)

             !exclude out-of-grid particle
             ! now we have the same conditions here and in f0_update_sp
             ! f0_n has to be larger than 0 for all particles that are
             ! not excluded by this statement
             if(i_mu >= f0_nmu .or. i_vp >= f0_nvp .or. i_vp < -f0_nvp .or.     &
                mu_n > f0_smu_max*f0_smu_max .or. vp_n > f0_vp_max .or. vp_n < -f0_vp_max) then
                cycle
             endif

             mu_vol(0:1)=1D0
             if(i_mu==0) then
                mu_vol(0)=0.5D0
                smu0=f0_dsmu/f0_mu0_factor
                smu1=f0_dsmu
             elseif (i_mu==f0_nmu-1) then
                mu_vol(1)=0.5D0
                smu0=f0_dsmu*real(i_mu,8)
                smu1=f0_dsmu*(real(f0_nmu,8)-1D0/f0_mu0_factor)
             else
                smu0=f0_dsmu*real(i_mu,8)
                smu1=f0_dsmu*real(i_mu+1,8)
             endif

             !temp variable for optimization
             t00=wmu(0)*wvp(0)
             t01=wmu(0)*wvp(1)
             t10=wmu(1)*wvp(0)
             t11=wmu(1)*wvp(1)

             !rh We do this axisymmetrically:
             !rh wsum was summed over all tor. sections ==> multiply volume by sml_nphi_total
             !sk - sml_nphi_total is not needed for XGC1
             wp(0:1)=f0_grid_vol(node,sp%type)*mu_vol(0:1) !*sml_nphi_total

             ! We found a particle in this phase space cell ---> mark df0g in that cell for deletion
             df0g_delete(i_vp+0,node,i_mu+0)=0D0 !rh 1.D0 - delfrac
             df0g_delete(i_vp+1,node,i_mu+0)=0D0 !rh 1.D0 - delfrac
             df0g_delete(i_vp+0,node,i_mu+1)=0D0 !rh 1.D0 - delfrac
             df0g_delete(i_vp+1,node,i_mu+1)=0D0 !rh 1.D0 - delfrac

#ifdef F0_DIST_WEIGHT_UPDATE
! This uses the weight evolution equation to modify the particle weight

             new_weight =  delfrac * ( t00/smu0 * f0_df0g2(i_vp,node,i_mu,sp%type)       &
                             +t01/smu0 * f0_df0g2(i_vp+1,node,i_mu,sp%type)     &
                             +t10/smu1 * f0_df0g2(i_vp,node,i_mu+1,sp%type)     &
                             +t11/smu1 * f0_df0g2(i_vp+1,node,i_mu+1,sp%type) ) &
                          / sqrt(f0_T_ev(node,sp%type))

             sp%ptl(i)%ph(piw1)=sp%ptl(i)%ph(piw1) + new_weight/sp%ptl(i)%ct(pif0)

#else
! This uses the reverse interpolation
#ifdef F_USE_MARKER_DEN
             new_weight =        f0_df0g2(i_vp+0,node,i_mu+0,sp%type) / smu0 * t00    &
                         +       f0_df0g2(i_vp+1,node,i_mu+0,sp%type) / smu0 * t01  &
                         +       f0_df0g2(i_vp+0,node,i_mu+1,sp%type) / smu1 * t10   &
                         +       f0_df0g2(i_vp+1,node,i_mu+1,sp%type) / smu1 * t11
             new_weight = new_weight/sqrt(f0_T_ev(node,sp%type))

             sp%ptl(i)%ph(piw1)=sp%ptl(i)%ph(piw1)+new_weight/sp%ptl(i)%ct(pif0)
#else
             new_weight =        f0_df0g2(i_vp+0,node,i_mu+0,sp%type) * wp(0)    &
                            * t00 / f0_n(i_vp+0,node,i_mu+0,sp%type)            &
                         +       f0_df0g2(i_vp+1,node,i_mu+0,sp%type) * wp(0)    &
                            * t01 / f0_n(i_vp+1,node,i_mu+0,sp%type)            &
                         +       f0_df0g2(i_vp+0,node,i_mu+1,sp%type) * wp(1)    &
                            * t10 / f0_n(i_vp+0,node,i_mu+1,sp%type)            &
                         +       f0_df0g2(i_vp+1,node,i_mu+1,sp%type) * wp(1)    &
                            * t11 / f0_n(i_vp+1,node,i_mu+1,sp%type)

             sp%ptl(i)%ph(piw1)=sp%ptl(i)%ph(piw1)+new_weight/sp%ptl(i)%ct(piw0)
#endif
#endif

!rh             if (is_nan(new_weight)) then
!rh                print *,sml_mype,sml_intpl_mype,node,i_mu,i_vp,     &
!rh                        f0_n(i_vp+0,node,i_mu+0,sp%type),               &
!rh                        f0_n(i_vp+1,node,i_mu+0,sp%type),               &
!rh                        f0_n(i_vp+0,node,i_mu+1,sp%type),               &
!rh                        f0_n(i_vp+1,node,i_mu+1,sp%type)
!rh                stop
!rh             endif


          endif
       enddo
    enddo

    if(sp%type==0) then
       call t_stopf("DIST_F0G_ELEC_LOOP")
       call t_startf("DIST_F0G_ELEC_RED2")
    else
       call t_stopf("DIST_F0G_ION_LOOP")
       call t_startf("DIST_F0G_ION_RED2")
    endif

    ! Gather df0g_delete information
    isize=size(df0g_delete)
    if(sml_symmetric_f) then
      call mpi_allreduce(df0g_delete,dum2,isize,MPI_REAL8,MPI_MIN,sml_intpl_comm,ierr)
      df0g_delete=dum2
    endif

    ! delete values of f0_f0g which have been used to change particle weights
    ! remainder stays on grid until a particle becomes available
    f0_df0g2(:,:,:,sp%type)=f0_df0g2(:,:,:,sp%type)*df0g_delete

    deallocate(df0g_delete,dum2)

    dum1=sum(abs(f0_df0g2(:,:,:,sp%type)))
    call mpi_allreduce(dum1, df0g_total2, 1, MPI_REAL8, MPI_SUM, sml_plane_comm,ierr)


    !f0_f0g(:,:,:,sp%type)=f0_f0g(:,:,:,sp%type)+f0_df0g2(:,:,:,sp%type)
    f0_df0g(:,:,:,sp%type)=f0_df0g2(:,:,:,sp%type)


    if (sml_mype==0) then
      print '(a54,e25.10)','Remaining sum of |f0_df0g2| (is transferred to grid): ', df0g_total2
    endif

    if(sp%type==0) then
       call t_stopf("DIST_F0G_ELEC_RED2")
    else
       call t_stopf("DIST_F0G_ION_RED2")
    endif

  end subroutine distribute_f0g

#ifdef COL_F_NAN_CHECK
  subroutine df0g_nan_check(str)
    use ptl_module
    implicit none
    character (len=*) :: str
    integer :: node, isp, ivp, imu
    logical :: is_nan
    do isp=ptl_isp, ptl_nsp
       do imu=f0_imu1, f0_imu2
          do node=f0_inode1, f0_inode2
             do ivp=-f0_nvp, f0_nvp

                if(is_nan(f0_df0g(ivp,node,imu,isp))) then
                   print *, 'NaN found in df0g (ivp,node,imu,isp)', ivp, node, imu, isp, str
                   stop
                endif
             enddo
          enddo
       enddo
    enddo
  end subroutine df0g_nan_check
#endif

  subroutine df0g_nan_remove(df0g)
    implicit none
    real (8),intent(inout) :: df0g(-f0_nvp:f0_nvp, f0_inode1:f0_inode2, f0_imu1:f0_imu2, ptl_isp:ptl_nsp)
    integer :: node, isp, ivp, imu
    logical :: is_nan
    logical :: happened
    integer :: cnt
    do isp=ptl_isp, ptl_nsp
       do imu=f0_imu1, f0_imu2
          do node=f0_inode1, f0_inode2
             happened=.false.
             cnt=0
             do ivp=-f0_nvp, f0_nvp
                if(is_nan(df0g(ivp,node,imu,isp))) then
                   df0g(ivp,node,imu,isp)=0D0
                   cnt=cnt+1
                   if(.not. happened) then
                      print *, 'NaN removed from df0g', ivp, node, imu, isp
                      happened=.true.
                   endif
                endif
             enddo
             if(happened .and. cnt>1) then
                print *, '# of NaNs,node,imu,isp',cnt, node, imu, isp
             endif
          enddo
       enddo
    enddo
  end subroutine df0g_nan_remove

  ! set zero for near wall nodes
  ! XGCa grid can have simple test
  subroutine df0g_near_wall_remove
    implicit none
    logical, save :: first=.true.
    integer, save :: minn, maxn
!pw    logical is_near_wall

    if(first) then
       !initialize
       minn=10000000  ! 10M nodes
       maxn=0

       do node=1, grid%nnode
          if(is_near_wall(node)) then
             minn=min(minn,node)
             maxn=max(maxn,node)
          endif
       enddo
              
    endif    


    do node=f0_inode1, f0_inode2
       if( minn <= node .and. node <= maxn) then
          if(is_near_wall(node)) then
             f0_df0g(:,node,:,:)=0D0
             if(f0_col_change_weight) f0_df0g2(:,node,:,:)=0D0
          endif
       endif
    enddo


  end subroutine df0g_near_wall_remove

  logical function is_near_wall(inode) !
    implicit none
    integer :: inode,j, tr, k, nd

    is_near_wall=.false.
    do j=1, grid%num_t_node(inode)
       tr=grid%tr_node(j,inode)
       do k=1,3
          nd=grid%nd(k,tr)
          if(grid%rgn(nd)==grid_rgn_wall) then
             is_near_wall=.true.
             return
          endif
       enddo
    enddo
  end function is_near_wall

end subroutine f_source


subroutine set_gvid0_pid_from_f0(nnode)
  use sml_module
  use f0_module
  implicit none
  integer, intent(in) :: nnode
  include 'mpif.h'
  integer :: sz, ierr

  gvid0_pid(sml_plane_mype)=f0_inode1

  sz=1
  call mpi_allgather(MPI_IN_PLACE,sz, mpi_integer, gvid0_pid, sz, mpi_integer, sml_plane_comm,ierr)

  gvid0_pid(sml_plane_totalpe)=nnode+1


end subroutine set_gvid0_pid_from_f0

subroutine symmetric_f0g
  use sml_module
  use f0_module
  use ptl_module
  implicit none
  include 'mpif.h'
  real (8), allocatable :: f0g_sum(:,:,:)
  integer :: size, isp, ierr

#ifndef F0_TOR_LINEAR
  allocate(f0g_sum(-f0_nvp:f0_nvp,f0_inode1:f0_inode2,f0_imu1:f0_imu2))

  size=(2*f0_nvp+1)*(f0_inode2-f0_inode1+1)*(f0_imu2-f0_imu1+1)

  do isp=ptl_isp, ptl_nsp
     call mpi_allreduce(f0_f0g(:,:,:,isp),f0g_sum,size,mpi_real8,mpi_sum,sml_intpl_comm,ierr)
     f0_f0g(:,:,:,isp)=f0g_sum/real(sml_intpl_totalpe)

  enddo
#else
  allocate(f0g_sum(-f0_nvp:f0_nvp,f0_inode1:f0_inode2,f0_imu1:f0_imu2))
  size=(2*f0_nvp+1)*(f0_inode2-f0_inode1+1)*(f0_imu2-f0_imu1+1)

  do isp=ptl_isp, ptl_nsp
     f0_f0g(:,:,:,0,isp)=f0_f0g(:,:,:,0,isp)+f0_f0g(:,:,:,1,isp)
     call mpi_allreduce(f0_f0g(:,:,:,0,isp),f0g_sum,size,mpi_real8,mpi_sum,sml_intpl_comm,ierr)
     f0_f0g(:,:,:,0,isp)=f0g_sum/real(sml_intpl_totalpe*2)
     f0_f0g(:,:,:,1,isp)=f0_f0g(:,:,:,0,isp)
  enddo
#endif
  deallocate(f0g_sum)
end subroutine symmetric_f0g


subroutine symmetric_f
  use sml_module
  use f0_module
  use ptl_module
  implicit none
  include 'mpif.h'
  real (8), allocatable :: f0g_sum(:,:,:)   ! means f_sum
  integer :: size, isp, ierr
  
  allocate(f0g_sum(-f0_nvp:f0_nvp,f0_inode1:f0_inode2,f0_imu1:f0_imu2))
  
  size=(2*f0_nvp+1)*(f0_inode2-f0_inode1+1)*(f0_imu2-f0_imu1+1)

  do isp=ptl_isp, ptl_nsp
     call mpi_allreduce(f0_f(:,:,:,isp),f0g_sum,size,mpi_real8,mpi_sum,sml_intpl_comm,ierr)
     f0_f(:,:,:,isp)=f0g_sum/real(sml_intpl_totalpe)

  enddo

  deallocate(f0g_sum)
end subroutine symmetric_f

subroutine f_radiation(grid)
  use sml_module
  use f0_module
  use grid_class
  implicit none
  type(grid_type) :: grid
  integer :: node, ivp, imu
  integer, parameter :: isp=0  ! electron only
  real (8) :: psi, mu, eperp1(f0_imu1:f0_imu2),epara(-f0_nvp:f0_nvp)
  real (8) :: f, fsum, fesum, mu_vol, eperp,b_b0, n, Tev, dE, meane, alpha
  logical :: rgn_check

  ! prepare eperp and epara
  do imu=f0_imu1,f0_imu2
     mu=imu*f0_dsmu
     mu=mu*mu !!rh mu is either norm. mag. moment (sqrt(mu)-grid) or (v_perp/v_th)^2 (v_perp-grid)
     eperp1(imu)=mu
  enddo

  do ivp=-f0_nvp, f0_nvp
     epara(ivp) = 0.5D0*(ivp*f0_dvp)**2
  enddo



  ! for each node
  do node=f0_inode1, f0_inode2
     ! check validity of f_radiation

     !rivate region check
     !rh rgn_check=(grid%rgn(node)==1 .or. grid%rgn(node)==2)
     !rh if(.not. rgn_check) cycle
     if(sml_inpsi > grid%psi(node) .or. sml_outpsi < grid%psi(node) &
        .or. (grid%x(1,node) .lt. eq_x_z .and. grid%psi(node) .lt. eq_x_psi) &
        .or. grid%rgn(node)==grid_rgn_wall ) cycle

     ! simulation boundary check
     psi = grid%psi(node)
     if(psi < sml_inpsi .or. sml_outpsi < psi) cycle

     ! obtain n & T
     !vsum=0D0
     fsum=0D0
     fesum=0D0
     b_b0=f0_b_b0(node)
     do imu=f0_imu1, f0_imu2

        ! normalized volume element
        if(imu==0 .or. imu==f0_nmu) then
           mu_vol=0.5D0
        else
           mu_vol=1D0
        endif

        eperp=0.5D0*eperp1(imu)
        do ivp=-f0_nvp,f0_nvp
           f=f0_f(ivp,node,imu,isp)
           f=max(f,0D0)
           !vsum  = vsum  + mu_vol
           fsum  = fsum  + mu_vol*f
           fesum = fesum + mu_vol*f*(eperp+epara(ivp))
        enddo
     enddo

     n=fsum*f0_grid_vol_vonly(node,isp)
     Tev=2D0/3D0*fesum/fsum*f0_T_ev(node,isp)

     !debug
     !if(sml_mype==0) print *, 'n & T', n, Tev

!!$     ! get cooling power
!!$     ! NRL radiation fomula :
!!$     ! Br : Bremsstrahlung
!!$     ! r  : recombination
!!$     ! c  : cyclotron --> ignored
!!$     ! P_Br + P_r = 1.69E-32 * Ne * sqrt(T_e_ev) * sum( Z^2 N_z ( 1 + (Ez-1)/Te ) ) watt/cm^3
!!$
!!$     ! SI
!!$     ! P = 1.69E-32 * Ne(SI) * sqrt(T_e_ev) * sum( Z^2 N_z(SI)*1E-6 ( ------------- ) )  watt/m^3
!!$     ! P = 1.69E-38 * Ne(SI) * sqrt(T_e_ev) * sum( Z^2 N_z(SI) ( ------------- ) )  watt/m^3
!!$     ! dE = P * dt /n  -- energy change per particle for given time interval (J)
!!$     dE = 1.69E-38 * sqrt(Tev) * Zsqr * n_imp * sml_dt * real(sml_f_source_period)


     ! Use formula from  ADAS GCR data for Carbon in adpak format
     ! rates as a function of electron temperature, neutral fraction
     !and residence time  ---> residence time is ignored: use infinte residence time (maximum value)

     call power_loss(dE)
     meane=3D0/2D0*Tev*sml_ev2j
     alpha = sqrt(1D0 - dE/meane) ! minus sign for alow loss  : v --> alpha v

     !if(sml_mype==0) print *, 'alpha=', alpha

     ! change of f0
     call f_shift_from_f_heat_torque


  enddo

contains
  subroutine power_loss(E_change)
    use rad_module
    use eq_module, only : eq_ftn, eq_zeff
    implicit none
    real (8) :: E_change
    real (8) :: r, z, Lz, avgZ, avgZ2, n_imp
    real (8) :: zeff, chargez, chargez2

    ! dE = P * dt / n
    ! P = n_imp n_e Lz

    r=grid%x(1,node)
    z=grid%x(2,node)
    call Rad_Emissivity2(psi,r,z,Lz,avgZ,avgZ2,TeV,n) ! r,z are for neutral density

    if(rad_use_zeff_profile) then

       zeff=eq_ftn(psi,r,z,eq_zeff)
       if(rad_use_fix_charge) then
          chargez=rad_z
          chargez2=rad_z*rad_z
       else
          chargez=avgZ
          chargez2=avgZ2
       endif

       n_imp=(zeff-1D0)/(chargez2-chargez)*n

    else
       n_imp = n*rad_impurity_fraction
    endif

    E_change = n_imp * Lz * sml_dt * real(sml_f_source_period)


  end subroutine power_loss
  subroutine f_shift_from_f_heat_torque
    implicit none
    real (8) :: dmu2, mu, epsilon1, epsilon2, vp
    integer, parameter :: dmu = -1
    real (8) :: f
    integer :: idvp,itmp_vp, itmp_mu
    real (8) :: df(-f0_nvp:f0_nvp, f0_imu1:f0_imu2)

    df=0D0

    do imu=f0_imu1, f0_imu2

    ! mu is sqrt(mu) or v_perp/v_th
    ! dmu is delta sqrt(mu) or delta v_perp/v_th

!rh #ifndef MU_LINEAR
       ! mu is sqrt(mu) or v_perp/v_th
       ! dmu is delta sqrt(mu) or delta v_perp/v_th
       if(imu==0) then
          !rh I believe dmu2 must be 1-1/f0_mu0_factor
          dmu2=1D0-1D0/f0_mu0_factor !rh 2D0/3D0
          mu=1D0/f0_mu0_factor
       elseif(imu==1 ) then
          dmu2=1D0-1D0/f0_mu0_factor !rh 2D0/3D0
          mu=1D0
       else
          dmu2=1D0
          mu=imu
       endif
!rh #else
!rh        if(imu==0) then
!rh           dmu2=2D0/3D0
!rh           mu=1D0/3D0
!rh        elseif(imu==1 ) then
!rh           dmu2=2D0/3D0
!rh           mu=1D0
!rh        else
!rh           dmu2=1D0
!rh           mu=imu
!rh        endif
!rh        ! same as without MU_LINEAR ??
!rh #endif


       epsilon1=mu/dmu2*(alpha-1D0)*real(dmu)
       epsilon1=min(epsilon1,1D0)                                   ! epsilon is smaller than 1

       do ivp=-f0_nvp, f0_nvp
          vp=ivp

          epsilon2=(alpha-1D0)*vp
          if(epsilon2>0D0) then
             idvp=1
          else
             idvp=-1
          endif

          epsilon2=epsilon2*real(idvp) ! epsilon is positive always
          epsilon2=min(epsilon2,1D0)             ! epsilon is smaller than 1. if it is large, source does not work correctly

          f=f0_f(ivp,node,imu,isp)
          f=max(f,0D0)

          itmp_vp=min(max(ivp+idvp,-f0_nvp),f0_nvp)
          itmp_mu=min(max(imu-dmu,0),f0_nmu)

          df(ivp,imu)=df(ivp,imu)+ ((1D0-epsilon1)*(1D0-epsilon2)-1D0)*f
          df(itmp_vp,imu)=df(itmp_vp,imu) + ((1D0-epsilon1)*epsilon2)*f
          df(ivp,itmp_mu)=df(ivp,itmp_mu) + ((1D0-epsilon2)*epsilon1)*f
          df(itmp_vp,itmp_mu)=df(itmp_vp,itmp_mu) + epsilon1*epsilon2*f


       enddo ! ivp
    enddo ! imu


    ! send - recv df at imu boundary

    f0_df0g(:,node,:,isp) = f0_df0g(:,node,:,isp) + df(:,f0_imu1:f0_imu2)

    !deallocate(df)


  end subroutine f_shift_from_f_heat_torque

end subroutine f_radiation


! calculate emissivity from data file
subroutine Rad_Emissivity2(psi_in,r_in,z_in,Lz_out,avgZ_out,avgZ2_out,tval_eV_in,eden) ! The input arguments are different from XGC0
  use rad_module
  use eq_module, only : eq_x_psi, eq_x_z
  implicit none


  real (kind=8), intent(in)  :: psi_in,z_in,r_in,tval_eV_in, eden
  real (kind=8), intent(out) :: Lz_out,avgZ_out,avgZ2_out
  integer :: i, j
  real (kind=8) :: tval_eV, n0_ne, delt, deln, Lz, avgZ, avgZ2, tmp1, tmp2
  real (kind=8), external ::  neutral_den

  tval_eV=tval_eV_in   !
  n0_ne=neutral_den(r_in,z_in,psi_in)/eden
!  if(z_in<eq_x_z) tval_eV=10D0

  tval_eV=max(tval_eV,rad_Te_ev(1)); tval_eV=min(tval_eV,rad_Te_ev(rad_nT))
  n0_ne=max(n0_ne,rad_n0_ne(1)); n0_ne=min(n0_ne,rad_n0_ne(rad_nn0ne))

  ! sequencial search
  do i=1, rad_nT
     if( rad_Te_ev(i) <= tval_ev .and. tval_ev <= rad_Te_ev(i+1) ) exit
  enddo

  do j=1, rad_nn0ne
     if( rad_n0_ne(j) <= n0_ne .and. n0_ne <= rad_n0_ne(j+1)) exit
  enddo
  tmp1=log(rad_Te_ev(i))
  tmp2=log(rad_Te_ev(i+1))
  delt = (log(tval_ev)-tmp1)/(tmp2-tmp1)

  tmp1=log(rad_n0_ne(j))
  tmp2=log(rad_n0_ne(j+1))
  deln = (log(n0_ne)-tmp1)/(tmp2-tmp1)

  Lz   = (1D0-delt)*(1D0-deln)*rad_Lz(i  ,j  ) &
       + (1D0-delt)* deln     *rad_Lz(i  ,j+1) &
       +    delt   *(1D0-deln)*rad_Lz(i+1,j  ) &
       +    delt   * deln     *rad_Lz(i+1,j+1)


  avgZ = (1D0-delt)*(1D0-deln)*rad_avgz(i  ,j  ) &
       + (1D0-delt)* deln     *rad_avgz(i  ,j+1) &
       +    delt   *(1D0-deln)*rad_avgz(i+1,j  ) &
       +    delt   * deln     *rad_avgz(i+1,j+1)

  avgz2= (1D0-delt)*(1D0-deln)*rad_avgz2(i  ,j  ) &
       + (1D0-delt)* deln     *rad_avgz2(i  ,j+1) &
       +    delt   *(1D0-deln)*rad_avgz2(i+1,j  ) &
       +    delt   * deln     *rad_avgz2(i+1,j+1)

  if((avgZ<2.0D0).or.(avgZ2<=avgZ)) then
     avgZ=2.0D0
     avgZ2=avgZ**2
  endif

  Lz_out=Lz
  avgZ_out=avgZ
  avgZ2_out=avgZ2
end subroutine Rad_Emissivity2


subroutine init_radiation
  use sml_module
  use rad_module
  implicit none
  include 'mpif.h'
  logical :: ok
  real (8) :: dum
  integer :: ierr, i, j


  !read radiation emissivity data for sml_radiation .and. sml_electron

  if(sml_mype==0) then
     inquire(file=rad_filename, exist=ok)
     if(.not. ok) then
        print *, 'file for radiation data is not found:',rad_filename
        stop
     endif
     open(400,file=rad_filename, status='old')
     read (400,*) rad_nT, rad_nn0ne
     allocate(rad_Te_ev(rad_nT),rad_n0_ne(rad_nn0ne))
     allocate(rad_Lz(rad_nT,rad_nn0ne), rad_avgZ(rad_nT,rad_nn0ne), rad_avgZ2(rad_nT,rad_nn0ne) )

     do j=1,rad_nn0ne
        do i=1, rad_nT
           read(400,*) rad_Te_ev(i),rad_n0_ne(j),dum,rad_Lz(i,j),rad_avgZ(i,j),rad_avgZ2(i,j)

           !debug
           if(sml_mype==0.and.i==1.and.j==1) print *, rad_Te_ev(i),rad_n0_ne(j),dum,rad_Lz(i,j),rad_avgZ(i,j),rad_avgZ2(i,j)
           if(sml_mype==0.and.i==2.and.j==1) print *, rad_Te_ev(i),rad_n0_ne(j),dum,rad_Lz(i,j),rad_avgZ(i,j),rad_avgZ2(i,j)
           if(sml_mype==0.and.i==1.and.j==2) print *, rad_Te_ev(i),rad_n0_ne(j),dum,rad_Lz(i,j),rad_avgZ(i,j),rad_avgZ2(i,j)
           if(sml_mype==0.and.i==2.and.j==2) print *, rad_Te_ev(i),rad_n0_ne(j),dum,rad_Lz(i,j),rad_avgZ(i,j),rad_avgZ2(i,j)
        enddo
     enddo
     ! rad_Lz=1D-6*rad_Lz
     rad_Lz=1D-6*rad_Lz  ! unit conversion from watt cm^3 to watt  m^3
     close(400)
  endif

  ! broadcast array size
  call mpi_bcast(rad_nT,    1, MPI_INTEGER, 0, sml_comm, ierr)
  call mpi_bcast(rad_nn0ne, 1, MPI_INTEGER, 0, sml_comm, ierr)

  ! allocate memory
  if(sml_mype/=0) then
     allocate(rad_Te_ev(rad_nT),rad_n0_ne(rad_nn0ne))
     allocate(rad_Lz(rad_nT,rad_nn0ne), rad_avgZ(rad_nT,rad_nn0ne), rad_avgZ2(rad_nT,rad_nn0ne) )
  endif

  ! broadcast data
  call mpi_bcast(rad_Te_ev, rad_nT   , MPI_REAL8, 0, sml_comm, ierr)
  call mpi_bcast(rad_n0_ne, rad_nn0ne, MPI_REAL8, 0, sml_comm, ierr)
  call mpi_bcast(rad_Lz,   rad_nT*rad_nn0ne, MPI_REAL8, 0, sml_comm, ierr)
  call mpi_bcast(rad_avgZ, rad_nT*rad_nn0ne, MPI_REAL8, 0, sml_comm, ierr)
  call mpi_bcast(rad_avgZ2,rad_nT*rad_nn0ne, MPI_REAL8, 0, sml_comm, ierr)


end subroutine init_radiation

subroutine f0_nan_check(str)
  use f0_module
  use ptl_module
  implicit none
  character (len=*) :: str
  integer :: node, isp, ivp, imu, iphi
  logical :: is_nan
  do isp=ptl_isp, ptl_nsp
     do node=f0_inode1, f0_inode2
        do imu=f0_imu1, f0_imu2
           do ivp=-f0_nvp, f0_nvp
#ifndef F0_TOR_LINEAR
              if(is_nan(f0_f0g(ivp,node,imu,isp))) then
                 print *, 'NaN found in f0g (ivp,node,imu,isp)', ivp, node, imu, isp, str
                 stop
              endif
#else
              do iphi=0, 1
                if(is_nan(f0_f0g(ivp,node,imu,iphi,isp))) then
                  print *, 'NaN found in f0g (ivp,node,imu,iphi,isp)', ivp, node, imu, iphi, isp, str
                  stop
                endif
              enddo
#endif
           enddo
        enddo
     enddo
  enddo
end subroutine f0_nan_check

logical function is_nan(a)
  implicit none
  real (8) :: a

  is_nan = .not. ( a > 1D0 .or. a < 2D0 )

end function is_nan

