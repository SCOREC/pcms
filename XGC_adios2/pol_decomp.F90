!!!!#include "assert_mod.F90"
#include "tri_pol_decomp_part.F90"

!!$subroutine decompose_triangles(grid,spall) ! by Ed
!!$  use grid_class
!!$  use ptl_module
!!$  use sml_module
!!$  use mpi
!!$  implicit none
!!$  type(grid_type) :: grid
!!$  type(species_type) :: spall(0:ptl_nsp_max)
!!$  !
!!$  integer, allocatable :: tri_weight(:), part(:)
!!$  integer :: ntr
!!$  logical :: is_ok
!!$  integer :: nparts, comm,ierror
!!$
!!$
!!$
!!$  ntr=grid%ntriangle
!!$
!!$  ! memory allocation
!!$  allocate( tri_weight(ntr), stat=ierror )
!!$  is_ok = (ierror.eq.0)
!!$  call assert( is_ok,'allocate( tri_weight ) ',ierror)        
!!$  allocate( part(ntr), stat=ierror )
!!$  is_ok = (ierror.eq.0)
!!$  call assert( is_ok,'allocate( tri_weight ) ',ierror)        
!!$
!!$  !-------------------------------
!!$  !get weight of each triangle
!!$  !-------------------------------
!!$  call get_weights(grid,spall,tri_weight)
!!$
!!$
!!$  ! ---------------------------
!!$  ! perform graph partitioning
!!$  !  ---------------------------    
!!$  nparts = sml_pe_per_plane ! NOTE: nparts is number of processors per plane
!!$  !
!!$  !	each processor independently run a "serial" copy of parmetis
!!$  !	this is the simplest and seems to give the best quality
!!$  !       -------------------------------------------------------------
!!$  comm  = MPI_COMM_SELF    
!!$  call tri_pol_decomp_part(ntr,grid%nd,nparts,tri_weight,comm,part)
!!$  !       -----------------------------------------------
!!$  !       optional: 
!!$  !       perform a broad cast to make sure all
!!$  !       processors agree on the initial partition
!!$  !       -----------------------------------------------    
!!$  comm = sml_comm
!!$  call MPI_BCAST( part, ntr, MPI_INTEGER, 0, comm,ierror)
!!$  is_ok = (ierror.eq.MPI_SUCCESS)
!!$  call assert(is_ok,'mpi_bcast(part) ',ierror)
!!$
!!$
!!$  !-----------------------------------
!!$  ! setup decomp_type data structure
!!$  !-----------------------------------
!!$  call setup_decomp_type(grid,part)
!!$
!!$
!!$  deallocate( tri_weight, part)    
!!$
!!$end subroutine decompose_triangles
!!$
!!$subroutine get_weights(grid,spall,tri_weight)
!!$  use grid_class
!!$  use ptl_module
!!$  use sml_module
!!$  use mpi
!!$  implicit none
!!$  type(grid_type),intent(in) :: grid
!!$  type(species_type) :: spall(0:ptl_nsp_max)
!!$  integer :: tri_weight(grid%ntriangle)
!!$  !
!!$  integer, allocatable :: count_ptl(:)
!!$  integer :: factor
!!$  integer :: isp, i, itr, ierror
!!$  logical :: is_ok
!!$  real (8) :: p(3)
!!$
!!$  allocate(count_ptl(grid%ntriangle))
!!$  count_ptl(:)=0
!!$
!!$  do isp=ptl_isp, ptl_nsp
!!$     if(isp==0) then
!!$        factor=5 ! electron factor
!!$     else
!!$        factor=1 ! ion factor
!!$     endif
!!$
!!$     do i=1,spall(isp)%num
!!$        if (spall(isp)%ptl(i)%gid <= 0) cycle
!!$
!!$        call search_tr2(grid,spall(isp)%ptl(i)%ph(1:2),itr,p)
!!$        if( (1 <= itr) .and. (itr <= grid%ntriangle) ) then
!!$           count_ptl(itr) = count_ptl(itr) + factor
!!$        endif
!!$     enddo
!!$  enddo
!!$
!!$  call MPI_ALLREDUCE( count_ptl, tri_weight, grid%ntriangle, MPI_INTEGER, MPI_SUM, sml_comm, ierror)
!!$  is_ok = (ierror.eq.MPI_SUCCESS)
!!$  call assert(is_ok,'mpi_allreduce(tri_weight)',ierror)
!!$
!!$  deallocate(count_ptl)
!!$
!!$end subroutine get_weights
!!$
!!$subroutine setup_decomp_type(grid,part)
!!$  use grid_class
!!$  use sml_module, only : sml_pe_per_plane, sml_plane_mype
!!$  implicit none
!!$  type(grid_type) :: grid
!!$  integer :: part(grid%ntriangle)
!!$  integer, dimension(0:sml_pe_per_plane-1) :: cnt, cnt2, cnt3
!!$  integer :: itr, pe, nd, type1
!!$ 
!!$
!!$  allocate(grid%p%gtr_2_pe_ltr(2,grid%ntriangle))  ! pe and local triangle index of a triangle(global)
!!$  allocate(grid%p%gnd_2_pe_lnd(2,grid%nnode))  !
!!$  
!!$  ! scan 'part'
!!$  cnt(:)=0  ! initialize
!!$  do itr = 1, grid%ntriangle
!!$     pe=part(itr)
!!$     cnt(pe)=cnt(pe) + 1
!!$     
!!$     grid%p%gtr_2_pe_ltr(:,itr)=(/ part(itr), cnt(pe) /)
!!$  enddo
!!$  grid%p%ntr = cnt(sml_plane_mype)
!!$
!!$  allocate(grid%p%l2g_tr(grid%p%ntr))  ! local triangle --> global triangle index
!!$  ! inverse function of gtr_2_pe_ltr
!!$  do itr=1, grid%ntriangle
!!$     if(part(itr)==sml_plane_mype) then
!!$        grid%p%l2g_tr(grid%p%gtr_2_pe_ltr(2,itr))=itr
!!$     endif
!!$  enddo
!!$  
!!$  !--- node ----
!!$  ! 1. finding nodes which are not shared with outher routine
!!$  cnt=0
!!$  cnt2=0
!!$  cnt3=0
!!$  do pe=1, sml_pe_per_plane
!!$     do nd=1, grid%nnode
!!$        ! checking type of node
!!$        call check_type(nd,pe,type1)
!!$
!!$        if(type1==1) then  ! non-sharing node
!!$           
!!$           cnt(pe)  = cnt(pe)  + 1
!!$           grid%p%gnd_2_pe_lnd(:,nd)= (/ pe, cnt(pe) /)
!!$
!!$        elseif(type1==2) then ! sharing node and put here
!!$
!!$           cnt2(pe) = cnt2(pe) + 1
!!$           grid%p%gnd_2_pe_lnd(:,nd)= (/ pe, -cnt2(pe) /)  ! put minus to distinguish --> put correct number later
!!$
!!$        elseif(type1==3) then ! sharing node and put elsewhere
!!$
!!$           cnt3(pe) = cnt3(pe) + 1
!!$
!!$        endif                        
!!$     enddo
!!$  enddo
!!$
!!$  ! set correct number for type 2 nodes (minus local node index)
!!$!  do pe=1, sml_pe_per_plane
!!$!     do nd=1, grid%nnode
!!$!        pe2  = grid%p%gnd_2_pe_lnd(1,nd)
!!$!        lnd = grid%p%gnd_2_pe_lnd(2,nd)
!!$!        if(lnd < 0) then
!!$!           grid%p%gnd_2_pe_lnd(2,nd) = -lnd + cnt(pe)
!!$!        endif
!!$!     enddo
!!$!  enddo
!!$  
!!$  ! local (=depending on processors) variables
!!$  grid%p%nn1=cnt(sml_plane_mype)
!!$  grid%p%nn2=grid%p%nn1 + cnt2(sml_plane_mype)
!!$  grid%p%nn3=grid%p%nn2 + cnt3(sml_plane_mype)
!!$  
!!$!  do nd=1, grid%nnode
!!$!     if(
!!$  
!!$
!!$
!!$  
!!$contains
!!$  subroutine check_type(nd,pe,type1)
!!$    implicit none
!!$    integer, intent(in) :: nd, pe
!!$    integer, intent(out) :: type1
!!$    integer :: i, itr, tr_pe
!!$    logical :: have_all, hit, hit_first
!!$
!!$    have_all=.true.  ! all triangle is in pe
!!$    hit=.false.      ! at least one triangle is in pe 
!!$    hit_first=.false.
!!$
!!$    ! scan all triangle that has the node and 
!!$    do i=1, grid%num_t_node(nd)
!!$       itr=grid%tr_node(i,nd)
!!$       tr_pe=grid%p%gtr_2_pe_ltr(1,itr)
!!$       if(pe==tr_pe) then
!!$          hit=.true.
!!$          if(i==1) hit_first=.true.
!!$       else
!!$          have_all=.false.
!!$       endif
!!$    enddo
!!$    
!!$    if(have_all) then
!!$       type1=1
!!$    elseif(hit_first) then
!!$       type1=2
!!$    elseif(hit) then
!!$       type1=3
!!$    else
!!$       type1=0
!!$    endif      
!!$
!!$  end subroutine check_type
!!$end subroutine setup_decomp_type

module pol_decomp_module
  use grid_class  
      
  !!PRIVATE MEMBER FUNCTIONS:
  private ceil2    
  private pair     
  private partition_opt, partition_opt_f0
  private partition_eval, partition_eval2
  private local_weights_eval
  private pol_decomp_wt_calc
  
!  real (8) ,allocatable :: gvid0_pid(:) ! starting global vertex(node) ID of pid_th processor
  integer ,allocatable :: gvid0_pid(:) ! starting global vertex(node) ID of pid_th processor
  integer ,allocatable :: gvid0_pid_old(:)

  logical :: pol_decomp_new
  
! variables needed for load balancing; f0 variables declared here to avoid
! circular dependence with f0_module  
  integer  :: f0_inode1, f0_inode2
  real (8) :: f0_ratio
  real (8), allocatable :: f0_node_cost(:)
  real (8) :: f0_grid_ptl_imbal
  logical  :: f0_grid_load_balance
  logical  :: f0_node_cost_updated

contains
  !
  ! set gvid0_pid
  subroutine set_weights(grid,spall)
    use grid_class
    use ptl_module
    use sml_module
    use col_module
    use omp_module, only: split_indices
    use perf_monitor
    implicit none
    include 'mpif.h'
    type(grid_type) :: grid
    type(species_type) :: spall(0:ptl_nsp_max)
    !     locals
    integer :: isp
    integer :: i, ierror
    integer :: pid, new_min_max_num_i, new_maxnum
    real (8) :: base_particle_number, ptl_constraint
!
    real (8) :: wt_gvid(grid%nnode), wt_gvid2(grid%nnode) ! weight to determine decomposition
    real (8) :: f0_wt_gvid(grid%nnode), f0_wt_gvid2(grid%nnode) ! weight to determine decomposition
    real (8) :: sum_wt_gvid, sum_f0_wt_gvid, new_min_max_num_r, new_min_max_cost_r

    ! particle number to determine a computational load of an empty node
    ! 10 is ideal for collisionless simulations with 1,000 ~ 10,000 particle per cell simulation case
    ! when nonlinear collision is on this number should be higher
    ! Some heuristic formula will determine this number later.
    base_particle_number=10D0



    if(sml_pol_decomp_simple .or. .not. sml_pol_decomp) then ! no decomposition, but using shift_ie
       return 
    endif

    ! set weight of cells
    ! use number of particles

    if(sml_electron_on) then
       isp=0
    else
       isp=1
    endif

    ! set weights for each node
    call t_startf("LOCAL_WTS_EVAL")
    call local_weights_eval(grid,spall(isp),wt_gvid2,grid%nnode)
    call t_stopf("LOCAL_WTS_EVAL")

    ! calculate weights, with final sum on process 0
    ! to guarantee that all processes use the same partition 
    ! (required since reduction is over a real vector)
    call t_startf("POL_DECOMP_WT_CALC_SUM")
    call pol_decomp_wt_calc(wt_gvid2, wt_gvid, grid%nnode, MPI_SUM)
    call t_stopf("POL_DECOMP_WT_CALC_SUM")

    if (sml_f0_grid .and. f0_grid_load_balance .and. (f0_ratio>0.D0)) then
       ! using max for f0 node cost since follows partition exactly
       f0_wt_gvid2(:) = 0.0D0
       f0_wt_gvid2(f0_inode1:f0_inode2) = f0_node_cost(f0_inode1:f0_inode2)
       call t_startf("POL_DECOMP_WT_CALC_MAX")
       call pol_decomp_wt_calc(f0_wt_gvid2, f0_wt_gvid, grid%nnode, MPI_MAX)
       call t_stopf("POL_DECOMP_WT_CALC_MAX")
    endif

    ! calculate partition (on process 0)
    if (sml_mype == 0) then
       ! to avoid too large amount of boundary grid on a proc, add base_particle_number particles per cell
       wt_gvid = wt_gvid + base_particle_number

       sum_wt_gvid = sum(wt_gvid)
       wt_gvid=wt_gvid/sum_wt_gvid

       if (sml_f0_grid .and. f0_grid_load_balance .and. (f0_ratio>0.D0)) then
          call partition_opt(wt_gvid, grid%nnode, gvid0_pid, sml_pe_per_plane, new_min_max_num_r, output=.false.)
          ptl_constraint = f0_grid_ptl_imbal*new_min_max_num_r

          sum_f0_wt_gvid = sum(f0_wt_gvid)
          f0_wt_gvid=f0_wt_gvid/sum_f0_wt_gvid

          call partition_opt_f0(f0_wt_gvid, wt_gvid, grid%nnode, ptl_constraint, &
                                gvid0_pid, sml_pe_per_plane, new_min_max_cost_r)
       else
          call partition_opt(wt_gvid, grid%nnode, gvid0_pid, sml_pe_per_plane, new_min_max_num_r)
       endif
    endif

    ! broadcast new partition
    call MPI_Bcast(gvid0_pid, sml_pe_per_plane+1, MPI_INTEGER, 0, sml_comm, ierror)
     call assert( ierror .eq. MPI_SUCCESS, &
                  'set_weights: MPI_Bcast problem in line ',__LINE__)

    ! for each particle type, determine new maximum particle count and determine whether
    ! need to increase local particle storage

    ! first, do current isp. if not load balancing sml_f0_grid work, 
    ! then use upper bound generated by partition_opt. otherwise, compute directly
    if (sml_f0_grid .and. f0_grid_load_balance .and. (f0_ratio>0.D0)) then
       if (sml_mype == 0) then
          wt_gvid2(1)=wt_gvid(1)
          do i=1+1, grid%nnode
             wt_gvid2(i)=wt_gvid2(i-1)+wt_gvid(i)
          enddo
          !
          new_min_max_num_r = wt_gvid2(gvid0_pid(1)-1)
          do pid=1, sml_pe_per_plane-2
             if (new_min_max_num_r < (wt_gvid2(gvid0_pid(pid+1)-1) - wt_gvid2(gvid0_pid(pid)-1))) then
                new_min_max_num_r = wt_gvid2(gvid0_pid(pid+1)-1) - wt_gvid2(gvid0_pid(pid)-1)
             endif
          enddo
          pid = sml_pe_per_plane-1
          if (new_min_max_num_r < (1.0D0 - wt_gvid2(gvid0_pid(pid)-1))) then
             new_min_max_num_r = (1.0D0 - wt_gvid2(gvid0_pid(pid)-1))
          endif
          !
          new_min_max_num_i = (new_min_max_num_r*sum_wt_gvid)/(sml_totalpe/sml_pe_per_plane)
       endif
    else
       if (sml_mype == 0) then
          !  min_max_num_i == approximate maximum number of particles per process after load balancing
          new_min_max_num_i = (new_min_max_num_r*sum_wt_gvid)/(sml_totalpe/sml_pe_per_plane)
       endif
    endif

    call mpi_bcast(new_min_max_num_i,1,mpi_integer,0,sml_comm,ierror)
     call assert( ierror .eq. MPI_SUCCESS, &
                  'set_weights: MPI_Bcast problem in line ',__LINE__)
    spall(isp)%min_max_num = new_min_max_num_i

    ! use sml_max_imblance*sml_max_imbalance as expansion factor for upper bound on particle memory
    new_maxnum = (sml_max_imbalance*sml_max_imbalance)*spall(isp)%min_max_num
    if (new_maxnum > spall(isp)%maxnum) then
       if (sml_mype == 0) then
          if (isp == 0) then
             write(6,*) "set_weights: increased elec. maxnum from,to:", spall(isp)%maxnum, new_maxnum
          else
             write(6,*) "set_weights: increased ion   maxnum from,to:", spall(isp)%maxnum, new_maxnum
             if (col_mode == 2) then
                call t_startf("SET_WTS_COL2_MEM_REALLOC")
                call col2_mem_reallocation(spall(isp)%maxnum, new_maxnum)
                call t_stopf("SET_WTS_COL2_MEM_REALLOC")
             endif
          endif
!pw       call flush(6)
       endif
       call t_startf("SET_WTS_PTL_MEM_REALLOC")
       call ptl_mem_reallocation(spall(isp), new_maxnum, spall(isp)%lost_nummax)
       call t_stopf("SET_WTS_PTL_MEM_REALLOC")
    endif
     
    ! if new partition chosen for electrons, next determine ion distribution for
    ! new partition and adjust local particle storage, as needed
    if (isp == 0) then
       isp = 1
       call local_weights_eval(grid,spall(isp),wt_gvid2,grid%nnode)

       call t_startf("POL_DECOMP_WT_CALC_SUM")
       call pol_decomp_wt_calc(wt_gvid2, wt_gvid, grid%nnode, MPI_SUM)
       call t_stopf("POL_DECOMP_WT_CALC_SUM")
       if (sml_mype == 0) then
          sum_wt_gvid = sum(wt_gvid)
          wt_gvid=wt_gvid/sum_wt_gvid
          !
          wt_gvid2(1)=wt_gvid(1)
          do i=1+1, grid%nnode
             wt_gvid2(i)=wt_gvid2(i-1)+wt_gvid(i)
          enddo
          !
          new_min_max_num_r = wt_gvid2(gvid0_pid(1)-1)
          do pid=1, sml_pe_per_plane-2
             if (new_min_max_num_r < (wt_gvid2(gvid0_pid(pid+1)-1) - wt_gvid2(gvid0_pid(pid)-1))) then
                new_min_max_num_r = wt_gvid2(gvid0_pid(pid+1)-1) - wt_gvid2(gvid0_pid(pid)-1)
             endif
          enddo
          pid = sml_pe_per_plane-1
          if (new_min_max_num_r < (1.0D0 - wt_gvid2(gvid0_pid(pid)-1))) then
             new_min_max_num_r = (1.0D0 - wt_gvid2(gvid0_pid(pid)-1))
          endif
          !
          new_min_max_num_i = (new_min_max_num_r*sum_wt_gvid)/(sml_totalpe/sml_pe_per_plane)
       endif
       call mpi_bcast(new_min_max_num_i,1,mpi_integer,0,sml_comm,ierror)
        call assert( ierror .eq. MPI_SUCCESS, &
                     'set_weights: MPI_Bcast problem in line ',__LINE__)
       spall(isp)%min_max_num = new_min_max_num_i

       new_maxnum = (sml_max_imbalance*sml_max_imbalance)*spall(isp)%min_max_num
       if (new_maxnum > spall(isp)%maxnum) then
          if (sml_mype == 0) then
             write(6,*) "set_weights: increased ion   maxnum from,to:", spall(isp)%maxnum, new_maxnum
             if (col_mode == 2) then
                call t_startf("SET_WTS_COL2_MEM_REALLOC")
                call col2_mem_reallocation(spall(isp)%maxnum, new_maxnum)
                call t_stopf("SET_WTS_COL2_MEM_REALLOC")
             endif
!pw          call flush(6)
          endif
          call t_startf("SET_WTS_PTL_MEM_REALLOC")
          call ptl_mem_reallocation(spall(isp), new_maxnum, spall(isp)%lost_nummax)
          call t_stopf("SET_WTS_PTL_MEM_REALLOC")
       endif
    endif

  end subroutine set_weights
!
  ! determine distribution of particles in each node after particles assigned to this
  ! process are shifted one step
  subroutine local_weights_eval(grid, sp, weights, wts_size)
    use grid_class
    use ptl_module
    use sml_module
    use omp_module, only: split_indices
    implicit none
    type(grid_type) :: grid
    type(species_type):: sp    
    integer , intent(in)  :: wts_size
    real (8), intent(out) :: weights(wts_size)

    ! local variables
    integer, allocatable :: cell_part(:)
    integer :: i_beg(sml_nthreads), i_end(sml_nthreads)
    integer :: ith, i
    integer :: iphi_g, max_node, itr, ml(1)

    real (8) :: phi_mid, x(3), xff(2), p(3)

    allocate(cell_part(sp%num))

    call split_indices(sp%num, sml_nthreads, i_beg, i_end)

    ! search index

!$OMP PARALLEL DO &
!$OMP PRIVATE( ITH, I, IPHI_G, PHI_MID, &
!$OMP          X, XFF, P, ITR, ML, MAX_NODE )
    do ith=1,sml_nthreads
!    do i=1, sp%num
       do i=i_beg(ith),i_end(ith)
          if(sp%ptl(i)%gid<=0) then
             cell_part(i)=-1
             cycle
          endif
          ! Modulo operation
          if(sp%ptl(i)%ph(3)>= sml_2pi_wedge_n .or. sp%ptl(i)%ph(3)< 0D0 ) then
             sp%ptl(i)%ph(3)=modulo(sp%ptl(i)%ph(3),sml_2pi_wedge_n)
          endif
          
          ! global iphi index
          iphi_g=floor(sp%ptl(i)%ph(3)/grid%delta_phi)

          !get field following position
          phi_mid=(real(iphi_g,8)+0.5D0)*grid%delta_phi
          x=sp%ptl(i)%ph(1:3)

          call field_following_pos2(x(1:2),x(3),phi_mid,xff)
          
          ! search triangle
          call search_tr2(grid,xff,itr,p)
          
          if(itr>0) then
             ml=maxloc(p)
             max_node= grid%nd(ml(1),itr)
             cell_part(i) = max_node
          else
             if(sml_sheath_mode==0 .or. sml_gstep <=0 ) then
                cell_part(i)=-1
             else
                ! do again with phase0
                x=sp%phase0(1:3,i)
                call field_following_pos2(x(1:2),x(3),phi_mid,xff)             
                call search_tr2(grid,xff,itr,p)
                if(itr>0) then
                   ml=maxloc(p)
                   max_node= grid%nd(ml(1),itr)
                   cell_part(i) = max_node
                else
                   ! give up for load balanching
                   cell_part(i) = -1
                endif
             endif
          endif
       enddo
    enddo

    !set weights for each node

    weights = 0.0d0
    do i=1, sp%num       
       max_node = cell_part(i)
       if(max_node<=0) then
          cycle
       endif
       weights(max_node) = weights(max_node) + 1.d0
    enddo
    deallocate(cell_part)

  end subroutine local_weights_eval
!
  ! determine particle weights, using sum, and f0 weights, using max,
  ! with results on process 0, when working with a poloidal decomposition
  subroutine pol_decomp_wt_calc(l_wt_gvid, g_wt_gvid, nnode, operator)
    use sml_module
    implicit none
    include 'mpif.h'
    
    integer, intent(in) :: nnode
    integer, intent(in) :: operator
    real (8), intent(in) :: l_wt_gvid(nnode) ! locally determined weights
    real (8), intent(out):: g_wt_gvid(nnode) ! weights

    ! local variables
    integer :: i, ierror
    integer :: min_node, max_node
    integer :: l_node_range(2), intpl_node_range(2)
    integer :: g_node_ranges(2,0:sml_plane_totalpe-1)
    integer :: g_node_counts(0:sml_plane_totalpe-1)
    integer :: g_node_displs(0:sml_plane_totalpe-1)
    integer :: intpl_node_count
    logical :: use_gather
!
    real (8) :: intpl_wt_gvid(nnode) ! work space for sum over interplanes

    ! initialize output array
    g_wt_gvid = 0.D0

    ! determine first node with particles or f0 cost
    min_node = nnode
    do i=1,nnode
       if (l_wt_gvid(i) > 0.0d0) then
          min_node = i
          exit
       endif
    enddo

    ! determine last node with particles or f0 cost
    max_node = 1
    do i=nnode,min_node,-1
       if (l_wt_gvid(i) > 0.0d0) then
          max_node = i
          exit
       endif
    enddo

    ! determine min and max across the interplane
    l_node_range(1) = -min_node
    l_node_range(2) = max_node
    intpl_node_range(:) = 0
    call MPI_Allreduce(l_node_range, intpl_node_range, &
                       2, MPI_INTEGER, MPI_MAX, &
                       sml_intpl_comm, ierror)
     call assert( ierror .eq. MPI_SUCCESS, &
                  'pol_decomp_wt_calc: MPI_Allreduce problem in line ',__LINE__)
    intpl_node_range(1) = -intpl_node_range(1)
    min_node = intpl_node_range(1)
    max_node = intpl_node_range(2)

    ! determine max of particles (or f0 cost) for each node in the interplane, 
    ! in intpl root
    if (max_node >= min_node) then
       intpl_node_count = (max_node - min_node + 1)
       call MPI_Reduce(l_wt_gvid(min_node:max_node), g_wt_gvid(min_node:max_node), &
                       intpl_node_count, MPI_REAL8, operator, 0, &
                       sml_intpl_comm, ierror)
        call assert( ierror .eq. MPI_SUCCESS, &
                     'pol_decomp_wt_calc: MPI_Reduce problem in line ',__LINE__)
    endif

    ! finish determination of weights within plane 0
    if (sml_intpl_mype == 0) then

       intpl_wt_gvid(:) = 0.0d0
       if (max_node >= min_node) then
          intpl_wt_gvid(min_node:max_node) = g_wt_gvid(min_node:max_node)
       endif

       ! determine min and max within plane 0
       g_node_ranges(:,:) = 0
       call MPI_Allgather(intpl_node_range, 2, MPI_INTEGER, &
                          g_node_ranges, 2, MPI_INTEGER, &
                          sml_plane_comm, ierror)
        call assert( ierror .eq. MPI_SUCCESS, &
                     'pol_decomp_wt_calc: MPI_Allgather problem in line ',__LINE__)

       ! if ranges are nonempty and nonoverlapping, then use gather,
       ! else use reduce. 
       use_gather = .true.

       ! all nonempty?
       do i=0,sml_plane_totalpe-1
          if (g_node_ranges(2,i) < g_node_ranges(1,i)) then
             use_gather = .false.
             exit
          endif
       enddo

       ! nonoverlapping?
       ! remember that displs(i)+1 is the first node for process i
       if (use_gather) then
          g_node_counts(0) = (g_node_ranges(2,0) - g_node_ranges(1,0) + 1)
          g_node_displs(0) = g_node_ranges(1,0) - 1
          do i=0,sml_plane_totalpe-2
             if (g_node_ranges(2,i) >= g_node_ranges(1,i+1)) then
                use_gather = .false.
                exit
             else
                g_node_counts(i+1) = (g_node_ranges(2,i+1) - g_node_ranges(1,i+1) + 1)
                g_node_displs(i+1) = g_node_ranges(1,i+1) - 1
             endif
          enddo
       endif

       if (use_gather) then
          call MPI_Gatherv(intpl_wt_gvid(min_node:max_node), intpl_node_count, MPI_REAL8, &
                           g_wt_gvid, g_node_counts, g_node_displs, &
                           MPI_REAL8, 0, sml_plane_comm, ierror)
           call assert( ierror .eq. MPI_SUCCESS, &
                        'pol_decomp_wt_calc: MPI_Gatherv problem in line ',__LINE__)
       else
          call MPI_Reduce(intpl_wt_gvid, g_wt_gvid, nnode, &
                          MPI_REAL8, operator, 0, sml_plane_comm, ierror)
           call assert( ierror .eq. MPI_SUCCESS, &
                        'pol_decomp_wt_calc: MPI_Reduce problem in line ',__LINE__)
       endif

    endif

  end subroutine pol_decomp_wt_calc
!
  ! find partition that minimizes the maximum f0_grid weight for a single partition element 
  ! when subject to a constraint on the maximum particle load imbalance
  subroutine partition_opt_f0(weights_f0, weights_ptl, wts_size, ptl_constraint, &
                              partition, part_size, max_part_wt)
    implicit none
    integer , intent(in)  :: wts_size
    real (8), intent(in)  :: weights_f0(wts_size)
    real (8), intent(in)  :: weights_ptl(wts_size)
    real (8), intent(in)  :: ptl_constraint
    integer , intent(in)  :: part_size
    integer , intent(out) :: partition(part_size+1)
    real (8), intent(out) :: max_part_wt

    ! local variables
    real (8) :: a, a_out, ptl_out
    real (8) :: b, c, b_out, c_out

    ! use bisection to determine max f0_grid weight subject to the particle weight constraint

    ! Lower bound on maximum weight is maxval(weights_f0))
    a = maxval(weights_f0)

    ! Evaluate partition for maxval(weights_f0) and check whether it generates the solution
    call partition_eval2(weights_f0, weights_ptl, wts_size, a, ptl_constraint, partition, &
                         part_size, a_out, ptl_out)

    if (a_out > 0.0d0) then
       ! a is a strict lower bound

       ! find an upper bound
       if (a < 1.0d0/real(part_size,8)) then
          c = 1.0d0/real(part_size,8)
          call partition_eval2(weights_f0, weights_ptl, wts_size, c, ptl_constraint, &
                               partition, part_size, c_out, ptl_out)
       else
          c = a
          c_out = a_out
       endif
       do while ((c_out > 0.0d0) .or. (ptl_out > 0.0d0))
          a = c
          c = 2*c
          call partition_eval2(weights_f0, weights_ptl, wts_size, c, ptl_constraint, &
                               partition, part_size, c_out, ptl_out)
       enddo
    
       if (c_out /= 0.0d0) then
          ! c is a strict upper bound, so use bisection-like algorithm to find a "good"
          ! feasible solution. Termination criterion is 1% of lower bound
          do while ((c-a) > 0.01*a)
             b = 0.5*(a+c)
             call partition_eval2(weights_f0, weights_ptl, wts_size, b, ptl_constraint, &
                                  partition, part_size, b_out, ptl_out)
             if      ((b_out > 0.0d0) .or. (ptl_out > 0.0d0)) then
                ! new lower bound
                a = b
             else if ((b_out < 0.0d0) .and. (ptl_out < 0.0d0)) then
                ! new upper bound
                c = b
             else
                ! exact solution
                exit
             endif
          enddo

          if (b == a) then
             b = c
          endif

       endif  

    else

       b = a
       b_out = a_out

    endif  

    ! regenerate partition
    call partition_eval2(weights_f0, weights_ptl, wts_size, b, ptl_constraint, &
                            partition, part_size, b_out, ptl_out)

    !debug
    !if(sml_mype==0) 
    write(*,216) b,b_out+b,ptl_constraint,ptl_out+ptl_constraint
216 format('partition_opt_f0 (f0):',e11.5,1x,e11.5,' (ptl:)',e11.5,1x,e11.5)
    !endif

    max_part_wt = b

    ! end condition
    partition(part_size+1)=wts_size+1   ! next (imaginary) gid on (imaginary) pid 

  end subroutine partition_opt_f0
!
  ! find partition that minimizes the maximum weight for a single partition element
  subroutine partition_opt(weights, wts_size, partition, part_size, max_part_wt, output)
    implicit none
    integer , intent(in)  :: wts_size
    real (8), intent(in)  :: weights(wts_size)
    integer , intent(in)  :: part_size
    integer , intent(out) :: partition(part_size+1)
    real (8), intent(out) :: max_part_wt
    logical , optional    :: output

    ! local variables
    logical :: alg2_worked
    logical :: l_output

    integer  :: i, pid
    integer  :: func_eval1, func_eval2
    integer  :: alg2_failure_mode

    real (8) :: weights2(wts_size)
    real (8) :: a, b, c, a_out, b_out, c_out
    real (8) :: pw, part_wt

    if (present(output)) then
       l_output = output
    else
       l_output = .true.
    endif

    ! Algorithm 1:
    ! Determine part_size element partition that minimizes the maximum weight
    ! while still assigning at least 1 node per partition.
    ! Goal is to find smallest maximum weight goal which is satisfied by the last partition
    ! after forcing it to be satisfied by all other partition elements. This is approx. a
    ! zero-finding problem (for a discontinuous but monotonic function) 
    ! where the approximation is required to be nonnegative (feasible).
    ! Using a bisection-like search method.

    ! Lower bound on maximum weight is maxval(wt_gvid)
    a = maxval(weights)

    weights2(1)=weights(1)
    do i=1+1,wts_size
       weights2(i)=weights2(i-1)+weights(i)
    enddo

    ! Evaluate partition for maxval(weights) and check whether it generates the solution
    call partition_eval(weights, wts_size, a, partition, part_size, a_out)
    func_eval1 = 1

    if (a_out > 0.0d0) then
       ! a is a strict lower bound

          ! find an upper bound
          if (a < 1.0d0/real(part_size,8)) then
             c = 1.0d0/real(part_size,8)
             call partition_eval(weights, wts_size, c, partition, part_size, c_out)
             func_eval1 = func_eval1 + 1
          else
             c = a
             c_out = a_out
          endif
          do while (c_out > 0.0d0)
             a = c
             c = 2*c
             call partition_eval(weights, wts_size, c, partition, part_size, c_out)
             func_eval1 = func_eval1 + 1
          enddo
    
          if (c_out /= 0.0d0) then
             ! c is a strict upper bound, so use bisection-like algorithm to find a "good"
             ! feasible solution. Termination criterion is 1% of lower bound
             func_eval2 = 0
             do while ((c-a) > 0.01*a)
                b = 0.5*(a+c)
                call partition_eval(weights, wts_size, b, partition, part_size, b_out)
                func_eval2 = func_eval2 + 1
                if      (b_out > 0.0d0) then
                   ! new lower bound
                   a = b
                else if (b_out < 0.0d0) then
                   ! new upper bound
                   c = b
                else
                   ! exact solution
                   exit
                endif
             enddo

             if (b == a) then
                b = c
                ! partition for approximate solution regenerated later, if necessary
                !  call partition_eval(weights, wts_size, b, partition, part_size, b_out)
                !  func_eval2 = func_eval2 + 1
             endif

          endif  

       else

          b = a
          b_out = a_out
          func_eval2 = 0

       endif  

       ! Algorithm 2:
       ! distribute weights as equally as possible over all processes
       alg2_worked = .true.
       alg2_failure_mode = 0
       pw  = 1.0d0/real(part_size,8)
       pid = 1
       partition(pid) = 1
       i = 1
       do while (i <= wts_size)
          if (weights2(i) > pw*real(pid)) then ! find weight
             if (pid >= part_size+1) then
                alg2_worked = .false.
                alg2_failure_mode = 1
                exit
             endif
             pid = pid+1
             partition(pid) = i
             ! do the same checking with increased pid
             i = i-1
          endif
          i=i+1
       enddo

       if (pid/=part_size) then
          if (pid/=part_size+1) then
             alg2_worked = .false.
             alg2_failure_mode = pid - part_size
          endif
       endif

       if (alg2_worked) then
          ! determine max cost
          max_part_wt = weights2(partition(2)-1)
          do pid=2,part_size-1
             part_wt = weights2(partition(pid+1)-1) - weights2(partition(pid)-1)
             if (max_part_wt < part_wt) then
                max_part_wt = part_wt
             endif
          enddo
          part_wt = (1.0D0 - weights2(partition(part_size)-1))
          if (max_part_wt < part_wt) then
             max_part_wt = part_wt
          endif
       else
          max_part_wt = 2.0d0
       endif

       if (max_part_wt > 1.01*b) then
          ! using Algorithm 1 result, so regenerate partition
          call partition_eval(weights, wts_size, b, partition, part_size, b_out)
          func_eval2 = func_eval2 + 1

          !debug
          if (l_output) then
             !if(sml_mype==0) 
             write(*,116) b,b_out+b,func_eval1+func_eval2
116          format('partition_opt: (using Alg. 1)',e13.5,1x,e13.5,1x,i6)
             if (alg2_worked) then
                write(*,117) max_part_wt
117             format('partition_opt: (alt. Alg. 2)',1x,e13.5,1x,i1)
             else
               write(*,118) alg2_failure_mode
118            format('partition_opt: (alt. Alg. 2 failed)',1x,i1)
             endif
             !endif
          endif

          max_part_wt = b
       else
          !debug
          if (l_output) then
             !if(sml_mype==0) 
             write(*,119) max_part_wt
119          format('partition_opt: (using Alg. 2)',e13.5)
             write(*,120) b,b_out+b,func_eval1+func_eval2
120          format('partition_opt: (alt. Alg. 1)',1x,e13.5,1x,e13.5,1x,i6)
             !endif
          endif
       endif

       ! end condition
       !partition(part_size+1)=real(wts_size+1)   ! next (imaginary) gid on (imaginary) pid 
       partition(part_size+1)=wts_size+1   ! next (imaginary) gid on (imaginary) pid 

       ! debug
       !if(sml_mype==0) then 
       !write(999,*) 1, partition(1), partition(2) - partition(1), weights2(partition(2)-1)
       !do pid=2, part_size-1
       !   write(999,*) pid, partition(pid), partition(pid+1) - partition(pid), &
       !                     weights2(partition(pid+1)-1) - weights2(partition(pid)-1)
       !enddo
       !pid = part_size
       !write(999,*) pid, partition(pid), partition(pid+1) - partition(pid), &
       !                  (1.0D0 - weights2(partition(pid)-1))
       !close(999)
       !endif

  end subroutine partition_opt
!
  ! evaluate the 'goodness' of the partition of the nodes across processes
  ! when given two weight functions and two goals
  subroutine partition_eval2(weights1, weights2, wts_size, goal1, goal2, partition, part_size, &
                             remainder1, remainder2)
    implicit none
    integer , intent(in)  :: wts_size
    real (8), intent(in)  :: weights1(wts_size)
    real (8), intent(in)  :: weights2(wts_size)
    real (8), intent(in)  :: goal1
    real (8), intent(in)  :: goal2
    integer , intent(in)  :: part_size
    integer , intent(out) :: partition(part_size)
    real (8), intent(out) :: remainder1
    real (8), intent(out) :: remainder2

    ! local variables
    integer  :: pid, i, partid
    real (8) :: cur_sum1, cur_sum2

    pid=1
    partition(pid)=1
    cur_sum1 = weights1(1)
    cur_sum2 = weights2(1)
    i=2
    ! assign weights until run out of weights or run out 
    ! of partition elements, leaving enough weights to assign 
    ! a single weight to any remaining elements
    do while ((i <= (wts_size-(part_size-pid))) .and. (pid < part_size))
       if ((cur_sum1 + weights1(i) > goal1) .or. (cur_sum2 + weights2(i) > goal2)) then
          pid=pid+1
          partition(pid)=i
          cur_sum1 = weights1(i)
          cur_sum2 = weights2(i)
       else
          cur_sum1 = cur_sum1 + weights1(i)
          cur_sum2 = cur_sum2 + weights2(i)
       endif
       i=i+1
    enddo
    do partid = pid+1, part_size
       partition(partid)=i
       i=i+1
    enddo
    cur_sum1 = 0.0d0
    cur_sum2 = 0.0d0
    do i=partition(part_size),wts_size
       cur_sum1 = cur_sum1 + weights1(i)
       cur_sum2 = cur_sum2 + weights2(i)
    enddo
    remainder1 = cur_sum1 - goal1
    remainder2 = cur_sum2 - goal2
  end subroutine partition_eval2
!
  ! evaluate the 'goodness' of the partition of the nodes across processes
  subroutine partition_eval(weights, wts_size, goal, partition, part_size, remainder)
    implicit none
    integer , intent(in)  :: wts_size
    real (8), intent(in)  :: weights(wts_size)
    real (8), intent(in)  :: goal
    integer , intent(in)  :: part_size
    integer , intent(out) :: partition(part_size)
    real (8), intent(out) :: remainder

    ! local variables
    integer  :: pid, i, partid
    real (8) :: cur_sum

    pid=1
    partition(pid)=1
    cur_sum = weights(1)
    i=2
    ! assign weights until run out of weights or run out 
    ! of partition elements, leaving enough weights to assign 
    ! a single weight to any remaining elements
    do while ((i <= (wts_size-(part_size-pid))) .and. (pid < part_size))
       if (cur_sum + weights(i) > goal) then
          pid=pid+1
          partition(pid)=i
          cur_sum = weights(i)
       else
          cur_sum = cur_sum + weights(i)
       endif
       i=i+1
    enddo
    do partid = pid+1, part_size
       partition(partid)=i
       i=i+1
    enddo
    cur_sum = 0.0d0
    do i=partition(part_size),wts_size
       cur_sum = cur_sum + weights(i)
    enddo
    remainder = cur_sum - goal
  end subroutine partition_eval
!
  subroutine setup_pol_decomp(grid,spall)
    use grid_class
    use ptl_module
    use sml_module
    implicit none
    include 'mpif.h'
    type(grid_type) :: grid
    type(species_type) :: spall(0:ptl_nsp_max)

    ! allocation -- side effect
    allocate(gvid0_pid(0:sml_pe_per_plane))
    allocate(gvid0_pid_old(0:sml_pe_per_plane))

    !call set_weights(grid,spall)

   end subroutine setup_pol_decomp
!
  subroutine finalize_pol_decomp
    implicit none
    if(allocated(gvid0_pid)) deallocate(gvid0_pid)
    if(allocated(gvid0_pid_old)) deallocate(gvid0_pid_old)
  end subroutine finalize_pol_decomp
!
  subroutine shift_ie(grid,psn,sp,phase0,ptl,shift_opt)
    use grid_class
    use psn_class
    use ptl_module
    use sml_module
    use omp_module , only : split_indices
    use perf_monitor
    implicit none
    include 'mpif.h'
    type(grid_type) :: grid
    type(psn_type) :: psn
    type(species_type):: sp    
    real (kind=8) :: phase0(ptl_nphase,sp%maxnum)  ! sp%phase0
    type(ptl_type) :: ptl(sp%maxnum)   ! sp%ptl
    integer, optional :: shift_opt(num_shift_ie_opts)

    ! shift_ie communication options. Defined in module.F90
    ! integer, parameter :: num_shift_ie_opts              = 10
    ! integer, parameter :: index_shift_ie_max_nthreads    =  1
    ! integer, parameter :: index_shift_ie_use_alltoall    =  2
    ! integer, parameter :: index_shift_ie_use_hs_barrier0 =  3
    ! integer, parameter :: index_shift_ie_large_limit     =  4
    ! integer, parameter :: index_shift_ie_handshake       =  5
    ! integer, parameter :: index_shift_ie_use_hs_barrier1 =  6
    ! integer, parameter :: index_shift_ie_use_sendrecv    =  7
    ! integer, parameter :: index_shift_ie_all_sendrecv    =  8
    ! integer, parameter :: index_shift_ie_use_isend       =  9
    ! integer, parameter :: index_shift_ie_use_rsend       = 10
    ! integer, parameter :: use_def_shift_ie_opt           = -2
    ! integer, parameter :: true_shift_ie_opt              =  1
    ! integer, parameter :: false_shift_ie_opt             =  0

    !
    integer :: alloc_stat
    integer :: send_tot, cnt_g, rearr_steps
    integer :: jindex, nrecv2, i, s, pid, itag, isize, ierror

    integer :: iphs, iphe, icts, icte, igid, iph0s, iph0e! index order of buffer
#ifdef PURE_RK4
    integer :: idphs, idphe, idphts, idphte
#endif
    integer :: m, mtop, last, arr_start, arr_end, pc_index
    integer :: ipe, max_steps, steps, step, recv_tot, r_offset
    integer :: signal, maxreq, maxreqh, rstep
    integer :: ith, jth, l_index, rb_index

    integer :: l_cnt_g(sml_nthreads)
    integer :: i_beg(sml_nthreads), i_end(sml_nthreads)
    integer :: s_beg(sml_nthreads), s_end(sml_nthreads)
    integer :: r_beg(sml_nthreads), r_end(sml_nthreads)
    integer, dimension(MPI_STATUS_SIZE) :: istatus

    integer, allocatable :: pid_part(:), l_hole(:,:), l_nsend(:,:), nsend(:)
    integer, allocatable :: hole(:), rearr_src(:), startid(:), l_startid(:,:)
    integer, allocatable :: nrecv(:), rrequest(:), hs_rrequest(:), srequest(:)
    integer, allocatable :: nspid(:), nrxxx(:), pids(:)

    real (8), allocatable :: sbuffer(:,:), rbuffer(:,:)

    !
    !for sheath -- temp variable
    integer :: itr
    real (8) :: p(3)

    !
    ! support for both simple toroidal decomposition and toroidal/poloidal decomposition
    integer :: comm, my_pid, total_pid, pid_type

    ! The following parameters control the parallel algorithm options.

    ! (a) max_nthreads limits the amount of OpenMP parallelism to
    !     exploit. There are no computational reasons to limit the
    !     number of threads, but memory requirements increase with
    !     the number of threads. If memory is limited, max_nthreads
    !     provides a mechanism to control the memory usage.
    !     If max_nthreads <=0, then there is no limit.
    integer, parameter :: def_max_nthreads = -1
    integer :: use_nthreads, max_nthreads

    ! (b) The first step of the communication algorithm determines
    !     the communication pattern. Two options are provided: using
    !     mpi_alltoall to communicate this information, or using 
    !     mpi_allreduce followed by targeted point-to-point messages
    !     For this second option, a barrier can also be used to
    !     guarantee that all receive requests are posted before the
    !     sends commence. This is meant to improve robustness when 
    !     there are concerns.
    logical, parameter :: def_use_alltoall = .false.
    logical, parameter :: def_use_hs_barrier0 = .false.
    logical :: use_alltoall, use_hs_barrier0

    ! (c) The second step of the communication algorithm moves
    !     the particle information. There are a number of different
    !     communication protocol options:

    !     (i) The communication algorithm preposts up to 'large_limit'
    !         mpi_irecv requests. Performance is generally improved by
    !         making large_limit as large as possible (<= 0 implies no
    !         limit), but too many outstanding requests can also hurt
    !         performance in some circumstances. Note that if the 
    !         mpi_allreduce option in (b) indicates more than 
    !         large_limit number of messages being received in the
    !         point-to-point phase, then mpi_alltoall is used
    !         instead
    integer, parameter :: def_large_limit = -1
    integer :: large_limit

    !     (ii) Due to the irregular communication pattern, some processes
    !         can receive many more messages than others. The handshake
    !         option forces the sending process to wait to send the
    !         information until the receiving process has indicated that
    !         it is ready to receive it. This can be important for 
    !         robustness. By preposting receive requests, the overhead 
    !         of these extra communications can be controlled.
    !         When large_limit <= 0, an alternative is to simply separate
    !         the receive and send requests with a barrier request. This can
    !         be used in other situations as well, though the motivation 
    !         is missing
    logical, parameter :: def_handshake = .false.
    logical, parameter :: def_use_hs_barrier1 = .false.
    logical :: handshake, use_hs_barrier1

    !     (iii) A complementary way to control the number of messages
    !         a single process receives unexpectedly is to use the XOR
    !         ordering and a 'logical' sendrecv, i.e. two processes
    !         swap messages even if only one them has any data to send.
    logical, parameter :: def_use_sendrecv = .false.
    logical :: use_sendrecv

    !     (iv) An even more extreme option is to force all processes
    !         to swap messages even if neither has any data to send to
    !         the other. This requires that both use_sendrecv and
    !         all_sendrecv be true. Without use_sendrecv, all_sendrecv
    !         does nothing of importance
    logical, parameter :: def_all_sendrecv = .false.
    logical :: all_sendrecv

    !     (v) The communication algorithm uses either mpi_irecv/mpi_send
    !         or mpi_irecv/mpi_isend. When using the handshake option,
    !         the ready send (mpi_rsend, mpi_irsend) is also an option
    logical, parameter :: def_use_isend = .false.
    logical, parameter :: def_use_rsend = .true.
    logical :: use_isend, use_rsend

    call t_startf("SHIFT_IE_INIT")

    ! set parallel algorithm options
    max_nthreads    = def_max_nthreads
    use_alltoall    = def_use_alltoall
    use_hs_barrier0 = def_use_hs_barrier0
    large_limit     = def_large_limit
    handshake       = def_handshake
    use_hs_barrier1 = def_use_hs_barrier1
    use_sendrecv    = def_use_sendrecv
    all_sendrecv    = def_all_sendrecv
    use_isend       = def_use_isend
    use_rsend       = def_use_rsend

    if (present(shift_opt)) then

       if (shift_opt(index_shift_ie_max_nthreads) .ne. use_def_shift_ie_opt) then
          max_nthreads = shift_opt(index_shift_ie_max_nthreads)
       endif

       if (shift_opt(index_shift_ie_use_alltoall) .ne. use_def_shift_ie_opt) then
          if (shift_opt(index_shift_ie_use_alltoall) == true_shift_ie_opt) then
             use_alltoall = .true.
          else
             use_alltoall = .false.
          endif
       endif

       if (shift_opt(index_shift_ie_use_hs_barrier0) .ne. use_def_shift_ie_opt) then
          if (shift_opt(index_shift_ie_use_hs_barrier0) == true_shift_ie_opt) then
             use_hs_barrier0 = .true.
          else
             use_hs_barrier0 = .false.
          endif
       endif

       if (shift_opt(index_shift_ie_large_limit) .ne. use_def_shift_ie_opt) then
          large_limit = shift_opt(index_shift_ie_large_limit)
          if (large_limit == 1) large_limit = 2
       endif

       if (shift_opt(index_shift_ie_handshake) .ne. use_def_shift_ie_opt) then
          if (shift_opt(index_shift_ie_handshake) == true_shift_ie_opt) then
             handshake = .true.
          else
             handshake = .false.
          endif
       endif

       if (shift_opt(index_shift_ie_use_hs_barrier1) .ne. use_def_shift_ie_opt) then
          if (shift_opt(index_shift_ie_use_hs_barrier1) == true_shift_ie_opt) then
             use_hs_barrier1 = .true.
          else
             use_hs_barrier1 = .false.
          endif
       endif

       if (shift_opt(index_shift_ie_use_sendrecv) .ne. use_def_shift_ie_opt) then
          if (shift_opt(index_shift_ie_use_sendrecv) == true_shift_ie_opt) then
             use_sendrecv = .true.
          else
             use_sendrecv = .false.
          endif
       endif

       if (shift_opt(index_shift_ie_all_sendrecv) .ne. use_def_shift_ie_opt) then
          if (shift_opt(index_shift_ie_all_sendrecv) == true_shift_ie_opt) then
             all_sendrecv = .true.
          else
             all_sendrecv = .false.
          endif
       endif

       if (shift_opt(index_shift_ie_use_isend) .ne. use_def_shift_ie_opt) then
          if (shift_opt(index_shift_ie_use_isend) == true_shift_ie_opt) then
             use_isend = .true.
          else
             use_isend = .false.
          endif
       endif

       if (shift_opt(index_shift_ie_use_rsend) .ne. use_def_shift_ie_opt) then
          if (shift_opt(index_shift_ie_use_rsend) == true_shift_ie_opt) then
             use_rsend = .true.
          else
             use_rsend = .false.
          endif
       endif

    endif

    if ( sml_pol_decomp_simple .or. .not. sml_pol_decomp ) then 
       comm = sml_intpl_comm
       my_pid = sml_intpl_mype
       total_pid = sml_intpl_totalpe
       pid_type = 1
    else
       comm = sml_comm
       my_pid = sml_mype
       total_pid = sml_totalpe
       pid_type = 0
    endif

    if (max_nthreads > 0) then
       use_nthreads = min(sml_nthreads,max_nthreads)
    else
       use_nthreads = sml_nthreads
    endif

    allocate(pid_part(sp%num), stat=alloc_stat)
     call assert( alloc_stat .eq. 0, &
                  'shift_ie: alloc(pid_part) return istat=',alloc_stat)
     pid_part(:) = 0

    allocate(l_hole((sp%num/use_nthreads)+1,use_nthreads), stat=alloc_stat)
     call assert( alloc_stat .eq. 0, &
                  'shift_ie: alloc(l_hole) return istat=',alloc_stat)
     l_hole(:,:) = 0

    allocate(l_nsend(0:total_pid-1,use_nthreads), stat=alloc_stat)
     call assert( alloc_stat .eq. 0, &
                  'shift_ie: alloc(l_nsend) return istat=',alloc_stat)
     l_nsend(:,:) = 0

    allocate(nsend(0:total_pid-1), stat=alloc_stat)
     call assert( alloc_stat .eq. 0, &
                  'shift_ie: alloc(nsend) return istat=',alloc_stat)
     nsend(:) = 0

    ! calculate data size and indices for a single particle
    iphs=1
    iphe=ptl_nphase
    
    icts=iphe+1
    icte=iphe+ptl_nconst  
    
    iph0s=icte+1
    iph0e=icte+ptl_nphase
    
    isize=iph0e

#ifdef PURE_RK4  
    !Extra data to transfer
    idphs=iph0e+1
    idphe=iph0e+ptl_nphase

    idphts=idphe+1
    idphte=idphe+ptl_nphase

    isize=idphte
#endif
    isize=isize+1 ! gid as real value

    ! find particles to be sent off process
    call split_indices(sp%num, use_nthreads, i_beg, i_end)

!$OMP PARALLEL DO &
!$OMP PRIVATE( ITH, I, PID, ITR, P )
    do ith=1,use_nthreads
       call t_startf("SHIFT_IE_LOOP1") 

       l_nsend(:,ith) = 0
       l_cnt_g(ith) = 0

       do i=i_beg(ith), i_end(ith)

          ! Modulo operation
          if(ptl(i)%ph(3)>= sml_2pi_wedge_n .or. ptl(i)%ph(3)< 0D0 ) then
             ptl(i)%ph(3)=modulo(ptl(i)%ph(3),sml_2pi_wedge_n)
          endif

          !search destination proc ID
          if ( ptl(i)%gid > 0 ) then
             call search_pid(grid,ptl(i)%ph(1:3),pid,pid_type)
             if ( pid == -2 ) then
                call sheath_calculation(grid,psn,sp,i,sp%type,itr,p,ith)
                if ( itr > 0 ) then
                   ! for sheath reflected particle
                   call search_pid(grid,ptl(i)%ph(1:3),pid,pid_type) 
                else
                   pid=-1
                endif
             endif
          else
             pid=-1
          endif

          if ( pid >= 0 ) then
             if ( pid /= my_pid ) then
                l_nsend(pid,ith)=l_nsend(pid,ith)+1
                l_cnt_g(ith) = l_cnt_g(ith) + 1
                l_hole(l_cnt_g(ith),ith) = i
             endif
             pid_part(i)=pid
          else
             call remove_particle(sp,i,-1,ith)
             pid_part(i)=my_pid
          endif
       enddo
       
       call t_stopf("SHIFT_IE_LOOP1") 
    enddo

    ! determine number of particles moving off process
    nsend(:) = l_nsend(:,1)
    do ith=2,use_nthreads
       do pid=0, total_pid-1
          nsend(pid) = nsend(pid) + l_nsend(pid,ith)
       enddo
    enddo
    send_tot=sum(nsend)

    ! Save location of particles moving off process (creating 'holes').
    !  Per thread lists (l_hole) are consecutive, since the list of particles 
    !  was traversed in a sequentual fashion, so just concatenate them.
    ! Also determine which threads will be handling which particles
    !  during buffer packing.
    call split_indices(send_tot, use_nthreads, s_beg, s_end)
    l_nsend = 0

    allocate(hole(send_tot), stat=alloc_stat)
     call assert( alloc_stat .eq. 0, &
                  'shift_ie: alloc(hole) return istat=',alloc_stat)
     hole(:) = 0

    allocate(rearr_src(send_tot), stat=alloc_stat)
     call assert( alloc_stat .eq. 0, &
                  'shift_ie: alloc(rearr_src) return istat=',alloc_stat)
     rearr_src(:) = 0

    allocate(startid(0:total_pid-1), stat=alloc_stat)
     call assert( alloc_stat .eq. 0, &
                  'shift_ie: alloc(startid) return istat=',alloc_stat)
     startid(:) = 0

    allocate(l_startid(0:total_pid-1,use_nthreads), stat=alloc_stat)
     call assert( alloc_stat .eq. 0, &
                  'shift_ie: alloc(l_startid) return istat=',alloc_stat)
     l_startid(:,:) = 0

    cnt_g = 0
    jth = 1
    do ith=1,use_nthreads
       do s=1,l_cnt_g(ith)
          ! save hole
          i = l_hole(s,ith)
          cnt_g = cnt_g + 1
          hole(cnt_g) = i

          ! calculate per thread information
          pid = pid_part(i)
          l_nsend(pid,jth) = l_nsend(pid,jth)  + 1
          if (s_end(jth) .eq. cnt_g) jth = jth+1
       enddo
    enddo
    deallocate(l_hole, stat=alloc_stat)
     call assert( alloc_stat .eq. 0, &
                  'shift_ie: dealloc(l_hole) return istat=',alloc_stat)

    ! SEND
    ! prepare send buffer indexing scheme
    jindex=1 
    startid(0)=1

    l_index=jindex
    l_startid(0,1)=startid(0)
    do ith=2,use_nthreads
       l_index=l_index+l_nsend(0,ith-1)
       l_startid(0,ith)=l_index
    enddo

    do pid=1, total_pid-1
       jindex=jindex+nsend(pid-1)
       startid(pid)=jindex

       l_index=jindex
       l_startid(pid,1)=startid(pid)
       do ith=2,use_nthreads
          l_index=l_index+l_nsend(pid,ith-1)
          l_startid(pid,ith)=l_index
       enddo

    enddo
    
    !allocate send buffer array
    allocate(sbuffer(isize,send_tot), stat=alloc_stat)
     call assert( alloc_stat .eq. 0, &
                  'shift_ie: alloc(sbuffer) return istat=',alloc_stat)
     sbuffer(:,:) = 0.D0

    ! pack to send buffer
!$OMP PARALLEL DO &
!$OMP PRIVATE( ITH, S, I, PID, JINDEX )
    do ith=1,use_nthreads
       call t_startf("SHIFT_IE_LOOP2") 

       l_nsend(:,ith) = 0

       do s=s_beg(ith), s_end(ith)
          i = hole(s)
          pid = pid_part(i)
          jindex = l_startid(pid,ith) + l_nsend(pid,ith)

          ! write particle info into buffer
          sbuffer(iphs:iphe,jindex)   =ptl(i)%ph
          sbuffer(icts:icte,jindex)   =ptl(i)%ct
          sbuffer(isize,jindex)       =ptl(i)%gid
          sbuffer(iph0s:iph0e,jindex) =phase0(:,i)
#ifdef PURE_RK4
          sbuffer(idphs:idphe,jindex)  =ptl(i)%dph
          sbuffer(idphts:idphte,jindex)=ptl(i)%dpht
#endif
          ! increase counter
          l_nsend(pid,ith)=l_nsend(pid,ith)+1
       enddo

       call t_stopf("SHIFT_IE_LOOP2") 
    enddo

    deallocate(pid_part,l_nsend,l_startid, stat=alloc_stat)
     call assert( alloc_stat .eq. 0, &
                  'shift_ie: dealloc(pid_part,l_nsend,l_startid) return istat=',alloc_stat)

    allocate(nrecv(0:total_pid-1), stat=alloc_stat)
     call assert( alloc_stat .eq. 0, &
                  'shift_ie: alloc(nrecv) return istat=',alloc_stat)
     nrecv(:) = 0

    allocate(rrequest(total_pid), stat=alloc_stat)
     call assert( alloc_stat .eq. 0, &
                  'shift_ie: alloc(rrequest) return istat=',alloc_stat)
     rrequest(:) = MPI_REQUEST_NULL

    if (handshake) then
       allocate(hs_rrequest(total_pid), stat=alloc_stat)
        call assert( alloc_stat .eq. 0, &
                     'shift_ie: alloc(hs_rrequest) return istat=',alloc_stat)
        hs_rrequest(:) = MPI_REQUEST_NULL
    endif

    if (use_alltoall) then

       ! notify processes how much data sending to each of them
       call t_startf("SHIFT_IE_A2A")
       ! nrecv(:) = 0 (set above)
       call mpi_alltoall(nsend, 1, MPI_INTEGER, &
                         nrecv, 1, MPI_INTEGER, &
                         comm, ierror)
        call assert( ierror .eq. MPI_SUCCESS, &
                     'shift_ie: MPI_Alltoall problem in line ',__LINE__)
       call t_stopf("SHIFT_IE_A2A")

    else

       allocate(nspid(0:total_pid-1), stat=alloc_stat)
        call assert( alloc_stat .eq. 0, &
                     'shift_ie: alloc(nspid) return istat=',alloc_stat)
        nspid(:) = 0

       allocate(nrxxx(0:total_pid-1), stat=alloc_stat)
        call assert( alloc_stat .eq. 0, &
                     'shift_ie: alloc(nrxxx) return istat=',alloc_stat)
        nrxxx(:) = 0

       ! find number of processes to send particles to
       ! nspid = 0 (set above)
       do pid=0, total_pid-1
          if (nsend(pid) > 0) nspid(pid) = 1
       enddo
    
       ! find number of processes sending particles to me
       call t_startf("SHIFT_IE_RED")
       call mpi_allreduce(nspid, nrxxx, total_pid, MPI_INTEGER, &
                          MPI_SUM, comm, ierror)
        call assert( ierror .eq. MPI_SUCCESS, &
                     'shift_ie: MPI_Allreduce problem in line ',__LINE__)
       call t_stopf("SHIFT_IE_RED")

       deallocate(nspid, stat=alloc_stat)
        call assert( alloc_stat .eq. 0, &
                     'shift_ie: dealloc(nspid) return istat=',alloc_stat)

       max_steps = maxval(nrxxx)

       if ((large_limit > 1) .and. (max_steps > large_limit)) then

          ! notify processes how much data sending to each of them
          ! using mpi_alltoall
          call t_startf("SHIFT_IE_A2A")
          ! nrecv(:) = 0 (set above)
          call mpi_alltoall(nsend, 1, MPI_INTEGER, &
                            nrecv, 1, MPI_INTEGER, &
                            comm, ierror)
           call assert( ierror .eq. MPI_SUCCESS, &
                        'shift_ie: MPI_Alltoall problem in line ',__LINE__)
          call t_stopf("SHIFT_IE_A2A")

       else

          ! notify processes how much data sending to each of them
          ! using point-=to-point messages
          steps = nrxxx(my_pid)

          call t_startf("SHIFT_IE_SR0")
          ! post receive requests
          itag = total_pid + my_pid
          do step=1, steps
             call mpi_irecv(nrxxx(step), 1, MPI_INTEGER, MPI_ANY_SOURCE, &
                            itag, comm, rrequest(step), ierror)
              call assert( ierror .eq. MPI_SUCCESS, &
                           'shift_ie: MPI_Irecv problem in line ',__LINE__)
          enddo

          if (use_hs_barrier0) then
             call t_startf("SHIFT_IE_BAR0")
             call mpi_barrier(comm, ierror)
              call assert( ierror .eq. MPI_SUCCESS, &
                           'shift_ie: MPI_Barrier problem in line ',__LINE__)
             call t_stopf("SHIFT_IE_BAR0")
          endif

          ! send message lengths to destinations using XOR ordering
          do ipe=1, ceil2(total_pid)-1
             pid = pair(total_pid,ipe,my_pid)
             if (pid >= 0) then
                if (nsend(pid) > 0) then
                   itag = total_pid + pid
                   call mpi_send(nsend(pid), 1, MPI_INTEGER, pid, itag, &
                                 comm, ierror)
                    call assert( ierror .eq. MPI_SUCCESS, &
                                 'shift_ie: MPI_Send problem in line ',__LINE__)
                endif
             endif
          enddo

          ! get size of message from each source
          nrecv = 0
          do step=1, steps
             call mpi_wait(rrequest(step),istatus,ierror)
              call assert( ierror .eq. MPI_SUCCESS, &
                           'shift_ie: MPI_Wait problem in line ',__LINE__)
             nrecv(istatus(MPI_SOURCE)) = nrxxx(step)
          enddo

          call t_stopf("SHIFT_IE_SR0")
       endif

       deallocate(nrxxx, stat=alloc_stat)
        call assert( alloc_stat .eq. 0, &
                     'shift_ie: dealloc(nrxxx) return istat=',alloc_stat)
    endif

    allocate(pids(total_pid), stat=alloc_stat)
     call assert( alloc_stat .eq. 0, &
                  'shift_ie: alloc(pids) return istat=',alloc_stat)
     pids(:) = -1

    ! precalculate XOR communication ordering and find amount of data 
    ! to be received
    recv_tot = 0
    ! pids = -1 (set above)
    steps = 0
    do ipe=1, ceil2(total_pid)-1
       pid = pair(total_pid,ipe,my_pid)
       if (pid >= 0) then
          if ((nrecv(pid) > 0) .or. (nsend(pid) > 0) .or. (all_sendrecv)) then
             steps = steps + 1
             pids(steps) = pid
             recv_tot = recv_tot + nrecv(pid)
          end if
       end if
    end do

    ! determine maximum number of outstanding send/receive requests
    if ((large_limit > 1) .and. (large_limit .le. steps)) then
       maxreq = large_limit
       maxreqh = maxreq/2
    else
       maxreq  = steps
       maxreqh = steps
    endif

    ! allocate receive buffer array
    allocate(rbuffer(isize,recv_tot), stat=alloc_stat)
     call assert( alloc_stat .eq. 0, &
                  'shift_ie: alloc(rbuffer) return istat=',alloc_stat)
     rbuffer(:,:) = 0.D0

    allocate(srequest(total_pid), stat=alloc_stat)
     call assert( alloc_stat .eq. 0, &
                  'shift_ie: alloc(sreqeust) return istat=',alloc_stat)
     srequest(:) = MPI_REQUEST_NULL

    call t_stopf("SHIFT_IE_INIT")

    call t_startf("SHIFT_IE_SR1")
    ! Initialize hs variable
    signal = 1

    ! Post initial handshake receive requests
    if (handshake) then
       do step=1, maxreq
          pid = pids(step)
          if ((nsend(pid) > 0) .or. (use_sendrecv)) then
             call mpi_irecv(signal, 1, MPI_INTEGER, pid, my_pid, &
                            comm, hs_rrequest(step), ierror)
              call assert( ierror .eq. MPI_SUCCESS, &
                           'shift_ie: MPI_Irecv problem in line ',__LINE__)
          endif
       enddo
    endif

    ! Post initial receive requests
    r_offset = 0
    do step=1, maxreq
       pid = pids(step)
       if ((nrecv(pid) > 0) .or. (use_sendrecv)) then
          if (handshake) then
             call mpi_irecv(rbuffer(1,r_offset+1), nrecv(pid)*isize, &
                            MPI_REAL8, pid, pid, &
                            comm, rrequest(step), ierror)
              call assert( ierror .eq. MPI_SUCCESS, &
                           'shift_ie: MPI_Irecv problem in line ',__LINE__)
             call mpi_send(signal, 1, MPI_INTEGER, pid, pid, &
                           comm, ierror)
              call assert( ierror .eq. MPI_SUCCESS, &
                           'shift_ie: MPI_Send problem in line ',__LINE__)
          else
             call mpi_irecv(rbuffer(1,r_offset+1), nrecv(pid)*isize, &
                            MPI_REAL8, pid, pid, &
                            comm, rrequest(step), ierror)
              call assert( ierror .eq. MPI_SUCCESS, &
                           'shift_ie: MPI_Irecv problem in line ',__LINE__)
          endif
          r_offset = r_offset + nrecv(pid)
       endif
    enddo
    rstep = maxreq
    call t_stopf("SHIFT_IE_SR1")

    call t_startf("SHIFT_IE_REARR")
    ! re-arrange particle array
    mtop=sp%num
    ! # of particles remain on local PE
    sp%num=sp%num-send_tot
    ! fill the hole
    last=send_tot
    !        print *, '[',my_pid,']','pnum,mtop,last',pnum,mtop,last

    ! determine what is to be copied where in rearrangement (before copying)
    rearr_steps = 0
    FILLHOLE : do i=1, send_tot
       m=hole(i)
       if( m > sp%num )  exit FILLHOLE
       !when empty space in the end
       do while( mtop == hole(last) )
          mtop=mtop-1
          last=last-1
       enddo
       rearr_steps = rearr_steps + 1
       rearr_src(rearr_steps) = mtop
       mtop=mtop-1
       if( mtop == sp%num) exit FILLHOLE
    enddo FILLHOLE

    ! do the copy
    if (rearr_steps .ge. use_nthreads) then
       call split_indices(rearr_steps, use_nthreads, s_beg, s_end)
!$OMP PARALLEL DO &
!$OMP PRIVATE( ITH, I, M, MTOP )
       do ith=1, use_nthreads

          do i=s_beg(ith),s_end(ith)
             m = hole(i)
             mtop = rearr_src(i)
             ptl(m)  = ptl(mtop)
             phase0(1:ptl_nphase,m) = phase0(1:ptl_nphase,mtop)
          enddo

       enddo
    else
       do i=1,rearr_steps
          m = hole(i)
          mtop = rearr_src(i)
          ptl(m)  = ptl(mtop)
          phase0(1:ptl_nphase,m) = phase0(1:ptl_nphase,mtop)
       enddo
    endif
    call t_stopf("SHIFT_IE_REARR")

    if (use_hs_barrier1) then
       call t_startf("SHIFT_IE_BAR1")
       call mpi_barrier(comm, ierror)
        call assert( ierror .eq. MPI_SUCCESS, &
                     'shift_ie: MPI_Barrier problem in line ',__LINE__)
       call t_stopf("SHIFT_IE_BAR1")
    endif

    call t_startf("SHIFT_IE_SR2")
    ! SEND
    do step=1, steps
       pid = pids(step)
       if ((nsend(pid) > 0) .or. (use_sendrecv)) then
          if (handshake) then
             call mpi_wait (hs_rrequest(step), istatus, ierror )
              call assert( ierror .eq. MPI_SUCCESS, &
                           'shift_ie: MPI_Wait problem in line ',__LINE__)
          endif
          if ((handshake) .and. (use_rsend)) then
             if (use_isend) then
                call mpi_irsend(sbuffer(:,startid(pid)), nsend(pid)*isize, &
                                MPI_REAL8, pid, my_pid, &
                                comm, srequest(step), ierror)
                 call assert( ierror .eq. MPI_SUCCESS, &
                              'shift_ie: MPI_Irsend problem in line ',__LINE__)
             else
                call mpi_rsend(sbuffer(:,startid(pid)), nsend(pid)*isize, &
                               MPI_REAL8, pid, my_pid, &
                               comm, ierror)
                 call assert( ierror .eq. MPI_SUCCESS, &
                              'shift_ie: MPI_Rsend problem in line ',__LINE__)
             endif
          else
             if (use_isend) then
                call mpi_isend(sbuffer(:,startid(pid)), nsend(pid)*isize, &
                               MPI_REAL8, pid, my_pid, &
                               comm, srequest(step), ierror)
                 call assert( ierror .eq. MPI_SUCCESS, &
                              'shift_ie: MPI_Isend problem in line ',__LINE__)
             else
                call mpi_send(sbuffer(:,startid(pid)), nsend(pid)*isize, &
                              MPI_REAL8, pid, my_pid, &
                              comm, ierror)
                 call assert( ierror .eq. MPI_SUCCESS, &
                              'shift_ie: MPI_Send problem in line ',__LINE__)
             endif
          endif
       endif

       if (step > maxreqh) then
          ! Wait for oldest irecv request to complete
          pid = pids(step-maxreqh)
          if ((nrecv(pid) > 0) .or. (use_sendrecv)) then
             call mpi_wait( rrequest(step-maxreqh), istatus, ierror )
              call assert( ierror .eq. MPI_SUCCESS, &
                           'shift_ie: MPI_Wait problem in line ',__LINE__)
          endif

          if (rstep < steps) then

             rstep = rstep + 1
             pid = pids(rstep)

             if (handshake) then
                ! Submit a new handshake irecv request
                if ((nsend(pid) > 0) .or. (use_sendrecv)) then
                   call mpi_irecv(signal, 1, MPI_INTEGER, pid, my_pid, &
                                  comm, hs_rrequest(rstep), ierror)
                    call assert( ierror .eq. MPI_SUCCESS, &
                                 'shift_ie: MPI_Irecv problem in line ',__LINE__)
                endif
             endif

             ! Submit a new irecv request
             if ((nrecv(pid) > 0) .or. (use_sendrecv)) then
                if (handshake) then
                   call mpi_irecv(rbuffer(1,r_offset+1), nrecv(pid)*isize, &
                                  MPI_REAL8, pid, pid, &
                                  comm, rrequest(rstep), ierror)
                    call assert( ierror .eq. MPI_SUCCESS, &
                                 'shift_ie: MPI_Irecv problem in line ',__LINE__)
                   call mpi_send(signal, 1, MPI_INTEGER, pid, pid, &
                                 comm, ierror)
                    call assert( ierror .eq. MPI_SUCCESS, &
                                 'shift_ie: MPI_Send problem in line ',__LINE__)
                else
                   call mpi_irecv(rbuffer(1,r_offset+1), nrecv(pid)*isize, &
                                  MPI_REAL8, pid, pid, &
                                  comm, rrequest(rstep), ierror)
                    call assert( ierror .eq. MPI_SUCCESS, &
                                 'shift_ie: MPI_Irecv problem in line ',__LINE__)
                endif
                r_offset = r_offset + nrecv(pid)
             endif

          endif

          if (use_isend) then
             ! Wait for outstanding i(r)send request to complete
             pid = pids(step-maxreqh)
             if ((nsend(pid) > 0) .or. (use_sendrecv)) then
                call mpi_wait( srequest(step-maxreqh), istatus, ierror )
                 call assert( ierror .eq. MPI_SUCCESS, &
                              'shift_ie: MPI_Wait problem in line ',__LINE__)
             endif
          endif

       endif
       
    enddo

    ! RECEIVE
    r_offset = 0
    do step=1, steps
       pid = pids(step)
       if ((nrecv(pid) > 0) .or. (use_sendrecv)) then
          if (step > steps-maxreqh) then
             call mpi_wait(rrequest(step),istatus,ierror)
              call assert( ierror .eq. MPI_SUCCESS, &
                           'shift_ie: MPI_Wait problem in line ',__LINE__)
          endif
          !unpack received particle
          if( sp%num + nrecv(pid) > sp%maxnum ) then 
             print * , 'Error :sp%num + nrecv(pid) > sp%maxnum: [', &
                       sml_mype, my_pid,'] step=',step,', nums:', &
                       sp%num+nrecv(pid),sp%maxnum
             stop
          endif
          if (nrecv(pid) .ge. use_nthreads) then
             call split_indices(nrecv(pid), use_nthreads, r_beg, r_end)
!$OMP PARALLEL DO &
!$OMP PRIVATE( ITH, M, L_INDEX, RB_INDEX )
             do ith=1,use_nthreads
                
                do m=r_beg(ith),r_end(ith)
                   l_index = sp%num+m
                   rb_index = r_offset+m
                   
                   ptl(l_index)%ph  =rbuffer(iphs:iphe  ,rb_index)
                   ptl(l_index)%ct  =rbuffer(icts:icte  ,rb_index)
                   ptl(l_index)%gid =rbuffer(isize      ,rb_index)
                   phase0(:,l_index)=rbuffer(iph0s:iph0e,rb_index)
#ifdef PURE_RK4
                   ptl(l_index)%dph =rbuffer(idphs:idphe  ,rb_index)
                   ptl(l_index)%dpht=rbuffer(idphts:idphte,rb_index)
#endif
                enddo

             enddo
          else
             do m=1,nrecv(pid)
                ptl(sp%num+m)%ph    =rbuffer(iphs:iphe  ,r_offset+m)
                ptl(sp%num+m)%ct    =rbuffer(icts:icte  ,r_offset+m)
                ptl(sp%num+m)%gid   =rbuffer(isize      ,r_offset+m)
                phase0(:,sp%num+m)  =rbuffer(iph0s:iph0e,r_offset+m)
#ifdef PURE_RK4
                ptl(sp%num+m)%dph   =rbuffer(idphs:idphe  ,r_offset+m)
                ptl(sp%num+m)%dpht  =rbuffer(idphts:idphte,r_offset+m)
#endif
             enddo
          endif
          r_offset = r_offset + nrecv(pid)
          sp%num = sp%num + nrecv(pid)
       endif
    enddo

    ! wait for sends to complete
    if (use_isend) then
       do step=steps-maxreqh+1,steps
          pid = pids(step)
          if ((nsend(pid) > 0) .or. (use_sendrecv)) then
             call mpi_wait(srequest(step),istatus,ierror)
              call assert( ierror .eq. MPI_SUCCESS, &
                           'shift_ie: MPI_Wait problem in line ',__LINE__)
          endif
       enddo
    endif
    call t_stopf("SHIFT_IE_SR2")

    if (handshake) then
       deallocate(hs_rrequest, stat=alloc_stat)
        call assert( alloc_stat .eq. 0, &
                     'shift_ie: dealloc(hs_rrequest) return istat=',alloc_stat)
    endif
    deallocate(nsend, hole, rearr_src, startid, sbuffer, &
               nrecv, rbuffer, rrequest, srequest, pids, &
               stat=alloc_stat)
     call assert( alloc_stat .eq. 0, &
                  'shift_ie: dealloc(shift_ie buffers) return istat=',alloc_stat)

  end subroutine shift_ie
!
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
!
  
end module pol_decomp_module

subroutine search_pid(grid,x,pid,pid_type)
  use grid_class
  use pol_decomp_module
  use sml_module
  implicit none
  type(grid_type) :: grid
  real (8) :: x(3)
  integer :: pid
  integer :: pid_type
  !
  integer :: itr, iphi_g, ml(1)
  real (8) :: xff(2), p(3), pw, x0(3), phi_mid
  integer, external :: pid_gvid
!  real (8) :: max_node
  integer :: max_node, intpl_pid, plane_pid
  logical :: intpl_only

  if (pid_type == 1) then
     intpl_only = .true.
  else
     intpl_only = .false.
  endif

  ! follow field line to poloidal plane
  iphi_g = floor(x(3)/grid%delta_phi)
  intpl_pid = iphi_g

  ! simple decomposition -- (no poloidal decomposition)
  if ( sml_pol_decomp_simple .or. .not. sml_pol_decomp ) then 
     ! returning interplane pid
     if (intpl_only) then
        pid = intpl_pid
     else
        pid = sml_plane_mype + (sml_plane_totalpe*intpl_pid)
     endif
     return
  endif

  ! find field following position
  phi_mid=(real(iphi_g,8)+0.5D0)*grid%delta_phi
  call field_following_pos2(x(1:2),x(3),phi_mid,xff)

  ! search 
  call search_tr2(grid,xff,itr,p)
  if(itr>0) then
     ml=maxloc(p)
     max_node= grid%nd(ml(1),itr) ! + (3D0*p(ml(1))-1D0)/2D0  

     ! pseudo binary search
     plane_pid=pid_gvid(max_node,grid%nnode,sml_pe_per_plane)
  
     ! calculate global pid from plane and intpl coordinates
     pid = plane_pid + (sml_plane_totalpe*intpl_pid)

  else
     if(sml_sheath_mode==0 .or. sml_gstep <=0 ) then
        pid=-1
     else
        pid=-2
     endif
  endif
     
end subroutine search_pid

integer function pid_gvid(gvid,nnodes,nproc)
  use pol_decomp_module
  implicit none
!  real (8) :: gvid
  integer :: gvid
  integer ::  nnodes,nproc
  !
   integer :: pid
   
!  integer :: mid, left, right

!  right=nproc
!  left=0

!  do while ( right > left+1)
!     mid=(right+left)/2
!     if(gvid0_pid(mid)<=gvid ) then
!        left=mid
!     else
!        right=mid
!     endif
!  enddo

  pid=floor(real(gvid)*real(nproc)/real(nnodes)) ! initial guess using linear assumption
  !print *, 'pid gvid gvid0_pid(pid)', pid, gvid, gvid0_pid(pid)

  do while(gvid0_pid(pid) > gvid)
     pid=pid-1
  end do
  
  do while(gvid0_pid(pid+1) <= gvid)
     pid=pid+1
  enddo

  pid_gvid=min(nproc-1,max(0,pid))
!  pid_gvid=left
  
end function  pid_gvid


  
