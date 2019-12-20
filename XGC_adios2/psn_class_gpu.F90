module psn_class_gpu
use dimensions_mod_gpu
use precision_mod_gpu
use cudafor
!use boundary_class_gpu, only : boundary2_type

!sa : enabling all these texture pointers won't do anything
!sa : right now, and there aren't any checks that they
!sa : are within the max texture size

#if 0
  real (kind=work_p),target, allocatable :: psn_pot_rho_ff(:,:,:) 
  real (kind=work_p),target, allocatable :: psn_E_rho_ff(:,:,:,:)
  real (kind=work_p),target, allocatable :: psn_ddpotdt_phi(:,:,:)
  real (kind=work_p),target, allocatable :: psn_pot0(:)

  real (kind=work_p),texture, pointer :: t_psn_pot_rho_ff(:,:,:) 
  real (kind=work_p),texture, pointer :: t_psn_E_rho_ff(:,:,:,:)
  real (kind=work_p),texture, pointer :: t_psn_ddpotdt_phi(:,:,:)
  real (kind=work_p),texture, pointer :: t_psn_pot0(:)

#else
  real (kind=work_p),allocatable :: psn_pot_rho_ff(:,:,:) 
  real (kind=work_p),allocatable :: psn_E_rho_ff(:,:,:,:)
  real (kind=work_p),allocatable :: psn_ddpotdt_phi(:,:,:)
  real (kind=work_p),allocatable :: psn_pot0(:)
#endif
  real (kind=work_p),allocatable :: psn_sheath_lost(:)
  real (kind=work_p),allocatable :: psn_sheath_pot(:)
  integer, allocatable :: psn_node_to_wall(:)
  integer, allocatable :: psn_wall_nodes(:)

#ifdef USE_CALC_GRADIENT
  real (kind=work_p), allocatable :: psn_pot_phi_real(:,:)
  integer, allocatable :: psn_ff_hdp_tr(:,:)
  integer, allocatable :: psn_ff_1dp_tr(:,:)
  real (kind=work_p), allocatable ::    psn_ff_hdp_p(:,:,:)
  real (kind=work_p), allocatable ::    psn_ff_1dp_p(:,:,:)
  real (kind=work_p), allocatable ::    psn_bfollow_1_dx(:)

  attributes(device) :: psn_pot_phi_real,            &
                        psn_ff_hdp_tr, psn_ff_hdp_p, &
                        psn_ff_1dp_tr, psn_ff_1dp_p, &
                        psn_bfollow_1_dx
#else

!sa: These textures actually work, and there is verification that
!sa: there are within the allowable texture size before binding. 
#ifdef USE_PSN_TEXTURE
  real (kind=work_p),texture,pointer :: t0_psn_E_phi_ff(:,:,:)
  real (kind=work_p),texture,pointer :: t1_psn_E_phi_ff(:,:,:)  
!sa: Split this into 2 arrays for better performance, and texture binding
  real (kind=work_p),target,allocatable :: psn_E_phi_ff0(:,:,:)
  real (kind=work_p),target,allocatable :: psn_E_phi_ff1(:,:,:)  
  private :: t0_psn_E_phi_ff, t1_psn_E_phi_ff
#else
!sa: Split this into 2 arrays for better performance
  real (kind=work_p),allocatable :: psn_E_phi_ff0(:,:,:)
  real (kind=work_p),allocatable :: psn_E_phi_ff1(:,:,:)  
#endif

  attributes(device) :: psn_E_phi_ff0, psn_E_phi_ff1
#endif
 
  attributes(device) :: psn_pot_rho_ff, psn_E_rho_ff,       &
                        psn_ddpotdt_phi, psn_pot0, psn_sheath_lost, psn_sheath_pot, &
                        psn_node_to_wall, psn_wall_nodes !, psn_E00_ff, psn_pot_phi_ff,psn_ddpotdt

  attributes(constant) :: psn_nwall

  ! Aliases for grid_nd texture access
#ifdef USE_GRID_TEXTURE
  integer, pointer, texture :: t_alias_grid_nd(:,:)
#else
  integer, pointer :: t_alias_grid_nd(:,:)
#endif
  contains

  attributes(host) &
  subroutine update_device_psn_type(psn)
  use sml_module, only : sml_mype
  use psn_class, only :  psn_type
#ifdef USE_GRID_TEXTURE
  use grid_class_gpu, only : grid_nd_
#else
  use grid_class_gpu, only : grid_nd
#endif
  use cudafor
  implicit none

  type(psn_type) :: psn
!  real (kind=work_p), parameter :: undefined = huge(0.0_work_p)
!  integer :: i3,i4
!  integer :: lb1,lb2,lb3,lb4, ub1,ub2,ub3,ub4
  integer :: icount_pr,icount_er, icount_e0, icount_pp, icount_ep, icount_dd
  integer :: ierr
#ifdef USE_GPU_EMU
  logical :: use_cudaMemcpy = .false.
#else
  logical :: use_cudaMemcpy = .true.
#endif
  integer, parameter :: idebug = 0 

  integer :: lb1,ub1, lb2,ub2, lb3,ub3, lb4,ub4, r2
  logical :: isvalid
  integer :: icount
  !! Temporary to reshape psn_E_phi_ff
  real (kind=work_p),allocatable :: reshape_E_phi_ff0(:,:,:)
  real (kind=work_p),allocatable :: reshape_E_phi_ff1(:,:,:)  

  psn_nwall = psn%nwall

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!! copy psn data from CPU to GPU !!!
  
!!  allocate(psn_pot_rho_ff(0:1,0:nrho_dim,nnode_dim), 
!!           psn_E_rho_ff(3,0:1,0:nrho_dim,nnode_dim), 
!!           psn_ddpotdt(nnode_dim,0:1), &
!!           psn_E00_ff(2,0:1,nnode_dim), 
!!           psn_pot_phi_ff(0:1,nnode_dim,0:(nphi_dim-1)),                                      &
!!           psn_E_phi_ff(3,0:1,nnode_dim,0:(nphi_dim-1)),
!!           psn_ddpotdt_phi(nnode_dim,0:1,0:(nphi_dim-1)))


! ========================
! allocate psn_pot_rho_ff(:,:,:)
! ========================
  if (.not.allocated(psn_pot_rho_ff)) then

    lb1 = lbound(psn%pot_rho_ff,1)
    lb2 = lbound(psn%pot_rho_ff,2)
    lb3 = lbound(psn%pot_rho_ff,3)

    ub1 = ubound(psn%pot_rho_ff,1)
    ub2 = ubound(psn%pot_rho_ff,2)
    ub3 = ubound(psn%pot_rho_ff,3)

    allocate( psn_pot_rho_ff(lb1:ub1,lb2:ub2,lb3:ub3), stat=ierr)
    call assert(ierr.eq.0,'alloc(psn_pot_rho_ff)',ierr)
    psn_pot_rho_ff = 0
  endif

  isvalid = size(psn%pot_rho_ff).eq.size(psn_pot_rho_ff)
  call assert(isvalid,'invalid size(psn_pot_rho_ff)',size(psn_pot_rho_ff))

! ======================
! allocate psn_E_rho_ff(:,:,:,:)
! ======================
  if (.not.allocated(psn_E_rho_ff)) then

    lb1 = lbound(psn%E_rho_ff,1)
    lb2 = lbound(psn%E_rho_ff,2)
    lb3 = lbound(psn%E_rho_ff,3)
    lb4 = lbound(psn%E_rho_ff,4)

    ub1 = ubound(psn%E_rho_ff,1)
    ub2 = ubound(psn%E_rho_ff,2)
    ub3 = ubound(psn%E_rho_ff,3)
    ub4 = ubound(psn%E_rho_ff,4)

    allocate( psn_E_rho_ff(lb1:ub1, lb2:ub2, lb3:ub3, lb4:ub4),stat=ierr)
    call assert(ierr.eq.0,'alloc(psn_E_rho_ff)',ierr)

    psn_E_rho_ff = 0
  endif

  isvalid = size(psn_E_rho_ff).eq.size(psn%E_rho_ff)
  call assert(isvalid,'invalid size(psn_E_rho_ff)',size(psn_E_rho_ff))

! ======================
! allocate psn_pot0(:)
! ======================
  if (.not.allocated(psn_pot0)) then

    lb1 = lbound(psn%pot0,1)

    ub1 = ubound(psn%pot0,1)

    allocate( psn_pot0(lb1:ub1),stat=ierr)
    call assert(ierr.eq.0,'alloc(psn_pot0)',ierr)

    psn_pot0 = 0
  endif

  isvalid = size(psn_pot0).eq.size(psn%pot0)
  call assert(isvalid,'invalid size(psn_pot0)',size(psn_pot0))

! ======================
! allocate psn_sheath_lost(:)
! ======================
  if (.not.allocated(psn_sheath_lost)) then

    lb1 = lbound(psn%sheath_lost(:,1),1)

    ub1 = ubound(psn%sheath_lost(:,1),1)
!print *, 'sheath_lost', lb1, ub1, psn%nwall, size(psn%sheath_lost)
    if(ub1>=lb1) allocate( psn_sheath_lost(lb1:ub1),stat=ierr)
    call assert(ierr.eq.0,'alloc(psn_sheath_lost)',ierr)

    psn_sheath_lost = 0
  endif

  isvalid = size(psn_sheath_lost).eq.size(psn%sheath_lost(:,1))
  call assert(isvalid,'invalid size(psn_sheath_lost)',size(psn_sheath_lost))

! ======================
! allocate psn_sheath_pot(:)
! ======================
  if (.not.allocated(psn_sheath_pot)) then

    lb1 = lbound(psn%sheath_pot,1)

    ub1 = ubound(psn%sheath_pot,1)
!print *, 'sheath_lost', lb1, ub1, psn%nwall, size(psn%sheath_lost)
    if(ub1>=lb1) allocate( psn_sheath_pot(lb1:ub1),stat=ierr)
    call assert(ierr.eq.0,'alloc(psn_sheath_pot)',ierr)

    psn_sheath_pot = 0
  endif

  isvalid = size(psn_sheath_pot).eq.size(psn%sheath_pot)
  call assert(isvalid,'invalid size(psn_sheath_pot)',size(psn_sheath_pot))

! ======================
! allocate psn_wall_nodes(:)
! ======================
  if (.not.allocated(psn_wall_nodes)) then

    lb1 = lbound(psn%wall_nodes,1)

    ub1 = ubound(psn%wall_nodes,1)

    if(ub1>=lb1) allocate( psn_wall_nodes(lb1:ub1),stat=ierr)
    call assert(ierr.eq.0,'alloc(psn_wall_nodes)',ierr)

    psn_wall_nodes = 0
  endif

  isvalid = size(psn_wall_nodes).eq.size(psn%wall_nodes)
  call assert(isvalid,'invalid size(psn_wall_nodes)',size(psn_wall_nodes))

! ======================
! allocate psn_node_to_wall(:)
! ======================
  if (.not.allocated(psn_node_to_wall)) then

    lb1 = lbound(psn%node_to_wall,1)

    ub1 = ubound(psn%node_to_wall,1)

    if(ub1>=lb1) allocate( psn_node_to_wall(lb1:ub1),stat=ierr)
    call assert(ierr.eq.0,'alloc(psn_node_to_wall)',ierr)

    psn_node_to_wall = 0
  endif

  isvalid = size(psn_node_to_wall).eq.size(psn%node_to_wall)
  call assert(isvalid,'invalid size(psn_node_to_wall)',size(psn_node_to_wall))


! =====================
! allocate  psn_ddpotdt(:,:)
! =====================

!  if (.not.allocated(psn_ddpotdt)) then
  
!    lb1 = lbound(psn%ddpotdt,1)
!    lb2 = lbound(psn%ddpotdt,2)

!    ub1 = ubound(psn%ddpotdt,1)
!    ub2 = ubound(psn%ddpotdt,2)

!    allocate( psn_ddpotdt(lb1:ub1, lb2:ub2),stat=ierr)
!    call assert(ierr.eq.0,'alloc(psn_ddpotdt)',ierr)

!    psn_ddpotdt = 0
!  endif

!  isvalid = size(psn%ddpotdt).eq.size(psn_ddpotdt)
!  call assert(isvalid,'invalid size(psn_ddpotdt)',size(psn_ddpotdt))


! =====================
! allocate  psn_E00_ff(:,:,:)
! =====================

!  if (.not.allocated(psn_E00_ff)) then

!    lb1 = lbound(psn%E00_ff,1)
!    lb2 = lbound(psn%E00_ff,2)
!    lb3 = lbound(psn%E00_ff,3)

!    ub1 = ubound(psn%E00_ff,1)
!    ub2 = ubound(psn%E00_ff,2)
!    ub3 = ubound(psn%E00_ff,3)

!    allocate( psn_E00_ff(lb1:ub1, lb2:ub2, lb3:ub3), stat=ierr)
!    call assert(ierr.eq.0,'alloc(psn_E00_ff)',ierr)

!    psn_E00_ff = 0
!  endif

!  isvalid = size(psn%E00_ff).eq.size(psn_E00_ff)
!  call assert(isvalid,'invalid size(psn_E00_ff)',size(psn_E00_ff))


! =====================
! allocate  psn_pot_phi_ff(:,:,:)
! =====================

!  if (.not.allocated(psn_pot_phi_ff)) then

!    lb1 = lbound(psn%pot_phi_ff,1)
!    lb2 = lbound(psn%pot_phi_ff,2)
!    lb3 = lbound(psn%pot_phi_ff,3)

!    ub1 = ubound(psn%pot_phi_ff,1)
!    ub2 = ubound(psn%pot_phi_ff,2)
!    ub3 = ubound(psn%pot_phi_ff,3)

!    allocate( psn_pot_phi_ff(lb1:ub1, lb2:ub2, lb3:ub3), stat=ierr)
!    call assert(ierr.eq.0,'alloc(psn_pot_phi_ff)',ierr)

!    psn_pot_phi_ff = 0
!  endif

!  isvalid = size(psn%pot_phi_ff).eq.size(psn_pot_phi_ff)
!  call assert(isvalid,'invalid size(psn_pot_phi_ff)',size(psn_pot_phi_ff))

#ifdef USE_CALC_GRADIENT

! =====================
! allocate  psn_pot_phi_real(:,:)
! =====================

  if (.not.allocated(psn_pot_phi_real)) then

    lb1 = lbound(psn%pot_phi_real,1)
    lb2 = lbound(psn%pot_phi_real,2)

    ub1 = ubound(psn%pot_phi_real,1)
    ub2 = ubound(psn%pot_phi_real,2)

    allocate( psn_pot_phi_real(lb1:ub1, lb2:ub2), stat=ierr)
    call assert(ierr.eq.0,'alloc(psn_pot_phi_real)',ierr)

    psn_pot_phi_real = 0
  endif

  isvalid = size(psn%pot_phi_real).eq.size(psn_pot_phi_real)
  call assert(isvalid,'invalid size(psn_pot_phi_real)',size(psn_pot_phi_real))


! =====================
! allocate  psn_ff_hdp_tr(:,:)
! =====================

  if (.not.allocated(psn_ff_hdp_tr)) then

    lb1 = lbound(psn%ff_hdp_tr,1)
    lb2 = lbound(psn%ff_hdp_tr,2)

    ub1 = ubound(psn%ff_hdp_tr,1)
    ub2 = ubound(psn%ff_hdp_tr,2)

    allocate( psn_ff_hdp_tr(lb1:ub1, lb2:ub2), stat=ierr)
    call assert(ierr.eq.0,'alloc(psn_ff_hdp_tr)',ierr)

    psn_ff_hdp_tr = 0
  endif

  isvalid = size(psn%ff_hdp_tr).eq.size(psn_ff_hdp_tr)
  call assert(isvalid,'invalid size(psn_ff_hdp_tr)',size(psn_ff_hdp_tr))



! =====================
! allocate  psn_ff_hdp_p(:,:,:)
! =====================

  if (.not.allocated(psn_ff_hdp_p)) then

    lb1 = lbound(psn%ff_hdp_p,1)
    lb2 = lbound(psn%ff_hdp_p,2)
    lb3 = lbound(psn%ff_hdp_p,3)

    ub1 = ubound(psn%ff_hdp_p,1)
    ub2 = ubound(psn%ff_hdp_p,2)
    ub3 = ubound(psn%ff_hdp_p,3)

    allocate( psn_ff_hdp_p(lb1:ub1, lb2:ub2, lb3:ub3), stat=ierr)
    call assert(ierr.eq.0,'alloc(psn_ff_hdp_p)',ierr)

    psn_ff_hdp_p = 0
  endif

  isvalid = size(psn%ff_hdp_p).eq.size(psn_ff_hdp_p)
  call assert(isvalid,'invalid size(psn_ff_hdp_p)',size(psn_ff_hdp_p))



! =====================
! allocate  psn_ff_1dp_tr(:,:)
! =====================

  if (.not.allocated(psn_ff_1dp_tr)) then

    lb1 = lbound(psn%ff_1dp_tr,1)
    lb2 = lbound(psn%ff_1dp_tr,2)

    ub1 = ubound(psn%ff_1dp_tr,1)
    ub2 = ubound(psn%ff_1dp_tr,2)

    allocate( psn_ff_1dp_tr(lb1:ub1, lb2:ub2), stat=ierr)
    call assert(ierr.eq.0,'alloc(psn_ff_1dp_tr)',ierr)

    psn_ff_1dp_tr = 0
  endif

  isvalid = size(psn%ff_1dp_tr).eq.size(psn_ff_1dp_tr)
  call assert(isvalid,'invalid size(psn_ff_1dp_tr)',size(psn_ff_1dp_tr))


! =====================
! allocate  psn_ff_1dp_p(:,:,:)
! =====================

  if (.not.allocated(psn_ff_1dp_p)) then

    lb1 = lbound(psn%ff_1dp_p,1)
    lb2 = lbound(psn%ff_1dp_p,2)
    lb3 = lbound(psn%ff_1dp_p,3)

    ub1 = ubound(psn%ff_1dp_p,1)
    ub2 = ubound(psn%ff_1dp_p,2)
    ub3 = ubound(psn%ff_1dp_p,3)

    allocate( psn_ff_1dp_p(lb1:ub1, lb2:ub2, lb3:ub3), stat=ierr)
    call assert(ierr.eq.0,'alloc(psn_ff_1dp_p)',ierr)

    psn_ff_1dp_p = 0
  endif

  isvalid = size(psn%ff_1dp_p).eq.size(psn_ff_1dp_p)
  call assert(isvalid,'invalid size(psn_ff_1dp_p)',size(psn_ff_1dp_p))


! =====================
! allocate  psn_bfollow_1_dx(:,:)
! =====================

  if (.not.allocated(psn_bfollow_1_dx)) then

    lb1 = lbound(psn%bfollow_1_dx,1)
    ub1 = ubound(psn%bfollow_1_dx,1)

    allocate( psn_bfollow_1_dx(lb1:ub1), stat=ierr)
    call assert(ierr.eq.0,'alloc(psn_bfollow_1_dx)',ierr)

    psn_bfollow_1_dx = 0
  endif

  isvalid = size(psn%bfollow_1_dx).eq.size(psn_bfollow_1_dx)
  call assert(isvalid,'invalid size(psn_bfollow_1_dx)',size(psn_bfollow_1_dx))

#else

! =====================
! allocate  psn_E_phi_ff(:,:,:,:)
! =====================

  if (.not.allocated(psn_E_phi_ff0)) then

    lb1 = lbound(psn%E_phi_ff,1)
    lb2 = lbound(psn%E_phi_ff,2)
    lb3 = lbound(psn%E_phi_ff,3)
    lb4 = lbound(psn%E_phi_ff,4)

    ub1 = ubound(psn%E_phi_ff,1)
    ub2 = ubound(psn%E_phi_ff,2)
    ub3 = ubound(psn%E_phi_ff,3)
    ub4 = ubound(psn%E_phi_ff,4)

    !! gpu access is much faster if we make this psn_E_phi_ff(3,node,iphi,0:1)
    !! And spliting it into 0/1 arrays is needed for alignment reasons
    allocate( reshape_E_phi_ff0(lb1:ub1, lb3:ub3, lb4:ub4), stat=ierr)
    allocate( reshape_E_phi_ff1(lb1:ub1, lb3:ub3, lb4:ub4), stat=ierr)    
    allocate( psn_E_phi_ff0(lb1:ub1, lb3:ub3, lb4:ub4), stat=ierr)
    call assert(ierr.eq.0,'alloc(psn_E_phi_ff0)',ierr)
    allocate( psn_E_phi_ff1(lb1:ub1, lb3:ub3, lb4:ub4), stat=ierr)    
    call assert(ierr.eq.0,'alloc(psn_E_phi_ff1)',ierr)

#ifdef USE_PSN_TEXTURE
    icount = ((ub1-lb1+1)*(ub3-lb3+1)*(ub4-lb4+1))
    if (icount >= texture_max) then
      print *, "Texture to big to bind!"
      print *, "Compile without USE_PSN_TEXTURE to run this grid"
    endif
    call assert(icount.lt.2**27,'E_phi_ff too big to bind! size:',icount)
#endif
    psn_E_phi_ff0 = 0
    psn_E_phi_ff1 = 0    

    reshape_E_phi_ff0(:,:,:) = psn%E_phi_ff(:,0,:,:)
    reshape_E_phi_ff1(:,:,:) = psn%E_phi_ff(:,1,:,:)    

  endif

  isvalid = size(psn%E_phi_ff).eq.(size(psn_E_phi_ff0) + size(psn_E_phi_ff1))
  call assert(isvalid,'invalid size(psn_E_phi_ff0 + psn_E_phi_ff1)',size(psn_E_phi_ff0) + size(psn_E_phi_ff1))

#endif



    
! =====================
! allocate  psn_ddpotdt_phi(:,:,:)
! =====================

  if (.not.allocated(psn_ddpotdt_phi)) then

    lb1 = lbound(psn%ddpotdt_phi,1)
    lb2 = lbound(psn%ddpotdt_phi,2)
    lb3 = lbound(psn%ddpotdt_phi,3)

    ub1 = ubound(psn%ddpotdt_phi,1)
    ub2 = ubound(psn%ddpotdt_phi,2)
    ub3 = ubound(psn%ddpotdt_phi,3)

    allocate( psn_ddpotdt_phi(lb1:ub1, lb2:ub2, lb3:ub3), stat=ierr)
    call assert(ierr.eq.0,'alloc(psn_ddpotdt_phi)',ierr)

    psn_ddpotdt_phi = 0
  endif

  isvalid = size(psn%ddpotdt_phi).eq.size(psn_ddpotdt_phi)
  call assert(isvalid,'invalid size(psn_ddpotdt_phi)',size(psn_ddpotdt_phi))







!!  icount_pr = 2*(nrho_dim+1)*nnode_dim
!!  icount_er = 3*2*(nrho_dim+1)*nnode_dim
!!  icount_e0 = nnode_dim*2
!!  icount_pp = 2*nnode_dim*nphi_dim
!!  icount_ep = 3*2*nnode_dim*nphi_dim
!!  icount_dd = nnode_dim*2*nphi_dim
  
  if(use_cudaMemcpy) then
     icount = size(psn_pot_rho_ff)
     ierr = cudaMemcpy( psn_pot_rho_ff, psn%pot_rho_ff,    &
               icount, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0, 'cudaMemcpy(psn_pot_rho_ff)',ierr)

     icount = size(psn_E_rho_ff)
     ierr = cudaMemcpy( psn_E_rho_ff, psn%E_rho_ff,    &
               icount, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0, 'cudaMemcpy(psn_E_rho_ff)',ierr)

     icount = size(psn_pot0)
     ierr = cudaMemcpy( psn_pot0, psn%pot0,    &
               icount, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0, 'cudaMemcpy(psn_pot0)',ierr)

     icount = size(psn_wall_nodes)
     ierr = cudaMemcpy( psn_wall_nodes, psn%wall_nodes,    &
               icount, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0, 'cudaMemcpy(psn_wall_nodes)',ierr)

!     icount = size(psn_sheath_lost)
!     ierr = cudaMemcpy( psn_sheath_lost, psn%sheath_lost,    &
!               icount, cudaMemcpyHostToDevice)
!     call assert(ierr.eq.0, 'cudaMemcpy(psn_sheath_lost)',ierr)

     icount = size(psn_sheath_pot)
     ierr = cudaMemcpy( psn_sheath_pot, psn%sheath_pot,    &
               icount, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0, 'cudaMemcpy(psn_sheath_pot)',ierr)

     icount = size(psn_node_to_wall)
     ierr = cudaMemcpy( psn_node_to_wall, psn%node_to_wall,    &
               icount, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0, 'cudaMemcpy(psn_node_wall)',ierr)


!     icount = size(psn_ddpotdt)
!     ierr = cudaMemcpy( psn_ddpotdt, psn%ddpotdt,    &
!               icount, cudaMemcpyHostToDevice)
!     call assert(ierr.eq.0, 'cudaMemcpy(psn_ddpotdt)',ierr)

!     icount = size(psn_E00_ff)
!     ierr = cudaMemcpy( psn_E00_ff, psn%E00_ff,    &
!               icount, cudaMemcpyHostToDevice)
!     call assert(ierr.eq.0, 'cudaMemcpy(psn_E00_ff)',ierr)



!     icount = size(psn_pot_phi_ff)
!     ierr = cudaMemcpy( psn_pot_phi_ff, psn%pot_phi_ff,    &
!               icount, cudaMemcpyHostToDevice)
!     call assert(ierr.eq.0,'cudaMemcpy(psn_pot_phi_ff)',ierr)


#ifdef USE_CALC_GRADIENT
     icount = size(psn_pot_phi_real)
     ierr = cudaMemcpy( psn_pot_phi_real, psn%pot_phi_real,    &
               icount, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(psn_pot_phi_real)',ierr)

     icount = size(psn_ff_hdp_tr)
     ierr = cudaMemcpy( psn_ff_hdp_tr, psn%ff_hdp_tr,    &
               icount, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(psn_ff_hdp_tr)',ierr)

     icount = size(psn_ff_hdp_p)
     ierr = cudaMemcpy( psn_ff_hdp_p, psn%ff_hdp_p,    &
               icount, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(psn_ff_hdp_p)',ierr)


     icount = size(psn_ff_1dp_tr)
     ierr = cudaMemcpy( psn_ff_1dp_tr, psn%ff_1dp_tr,    &
               icount, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(psn_ff_1dp_tr)',ierr)

     icount = size(psn_ff_1dp_p)
     ierr = cudaMemcpy( psn_ff_1dp_p, psn%ff_1dp_p,    &
               icount, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(psn_ff_1dp_p)',ierr)

     icount = size(psn_bfollow_1_dx)
     ierr = cudaMemcpy( psn_bfollow_1_dx, psn%bfollow_1_dx,    &
               icount, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(psn_bfollow_1_dx)',ierr)


#else
     icount = size(psn_E_phi_ff0)
     ierr = cudaMemcpy( psn_E_phi_ff0, reshape_E_phi_ff0,    &
               icount, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(psn_E_phi_ff0)',ierr)

     icount = size(psn_E_phi_ff1)
     ierr = cudaMemcpy( psn_E_phi_ff1, reshape_E_phi_ff1,    &
               icount, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(psn_E_phi_ff1)',ierr)


#endif


     icount = size(psn_ddpotdt_phi)
     ierr = cudaMemcpy( psn_ddpotdt_phi, psn%ddpotdt_phi,    &
               icount, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(psn_ddpotdt_phi)',ierr)
  else
     psn_pot_rho_ff = psn%pot_rho_ff
     psn_E_rho_ff = psn%E_rho_ff
!     psn_ddpotdt = psn%ddpotdt
!     psn_E00_ff = psn%E00_ff

!     psn_pot_phi_ff = psn%pot_phi_ff

#ifdef USE_CALC_GRADIENT
     psn_pot_phi_real = psn%pot_phi_real
     psn_ff_hdp_tr = psn%ff_hdp_tr
     psn_ff_hdp_p = psn%ff_hdp_p
     psn_ff_1dp_tr = psn%ff_1dp_tr
     psn_ff_1dp_p = psn%ff_1dp_p
     psn_bfollow_1_dx = psn%bfollow_1_dx
#else
     psn_E_phi_ff0 = reshape_E_phi_ff0
     psn_E_phi_ff1 = reshape_E_phi_ff1     
#endif
     psn_ddpotdt_phi = psn%ddpotdt_phi
  endif
 
#if 0
  
  t_psn_pot_rho_ff => psn_pot_rho_ff
  t_psn_E_rho_ff => psn_E_rho_ff
  t_psn_ddpotdt_phi => psn_ddpotdt_phi
  t_psn_pot0 => psn_pot0
#endif

#ifndef USE_CALC_GRADIENT
#ifdef USE_PSN_TEXTURE
  t0_psn_E_phi_ff => psn_E_phi_ff0
  t1_psn_E_phi_ff => psn_E_phi_ff1
#endif
  deallocate(reshape_E_phi_ff0)
  deallocate(reshape_E_phi_ff1)  
#endif


#ifdef USE_GRID_TEXTURE
  t_alias_grid_nd => grid_nd_
#else
  t_alias_grid_nd => grid_nd
#endif
  return
  end subroutine update_device_psn_type



  attributes(host) &
  subroutine update_host_psn_type(psn)
  use psn_class, only : psn_type
  use cudafor
  implicit none
  integer, parameter :: idebug = 0
  integer :: ierr, icount 
  type(psn_type):: psn
  real (8) :: tmp(psn%nwall)
#ifdef USE_GPU_EMU
  logical, parameter :: use_cudaMemcpy = .false.
#else
  logical, parameter :: use_cudaMemcpy = .true.
#endif

  if(use_cudaMemcpy) then

     !icount = size(psn%sheath_lost)
     icount = size(psn_sheath_lost)
     ierr = cudaMemcpy( tmp, psn_sheath_lost,  &
               icount, cudaMemcpyDeviceToHost)
     call assert(ierr.eq.0,'cudaMemcpy(psn%sheath_lost)',ierr)

  else

     tmp = psn_sheath_lost

  endif
  psn%sheath_lost(:,1) = psn%sheath_lost(:,1) + tmp

  if(allocated(psn_pot_rho_ff)) then
    deallocate( psn_pot_rho_ff,stat=ierr)
    call assert(ierr.eq.0,'dealloc(psn_pot_rho_ff)',ierr)
#if 0
    nullify( t_psn_pot_rho_ff )
#endif
  endif

  if(allocated(psn_E_rho_ff)) then
    deallocate( psn_E_rho_ff,stat=ierr)
    call assert(ierr.eq.0,'dealloc(psn_E_rho_ff)',ierr)
#if 0
    nullify( t_psn_E_rho_ff )
#endif
  endif

  if(allocated(psn_ddpotdt_phi)) then
    deallocate( psn_ddpotdt_phi,stat=ierr)
    call assert(ierr.eq.0,'dealloc(psn_ddpotdt_phi)',ierr)
#if 0
    nullify( t_psn_ddpotdt_phi )
#endif

  endif

  if(allocated(psn_node_to_wall)) then
    deallocate( psn_node_to_wall,stat=ierr)
    call assert(ierr.eq.0,'dealloc(psn_node_to_wall)',ierr)
  endif

  if(allocated(psn_wall_nodes)) then
    deallocate( psn_wall_nodes,stat=ierr)
    call assert(ierr.eq.0,'dealloc(psn_wall_nodes)',ierr)
  endif

  if(allocated(psn_sheath_lost)) then
    deallocate( psn_sheath_lost,stat=ierr)
    call assert(ierr.eq.0,'dealloc(psn_sheath_lost)',ierr)
  endif

  if(allocated(psn_sheath_pot)) then
    deallocate( psn_sheath_pot,stat=ierr)
    call assert(ierr.eq.0,'dealloc(psn_sheath_pot)',ierr)
  endif

  if(allocated(psn_pot0)) then
    deallocate( psn_pot0,stat=ierr)
    call assert(ierr.eq.0,'dealloc(psn_pot0)',ierr)
#if 0
    nullify( t_psn_pot0 )
#endif
  endif

#ifdef USE_CALC_GRADIENT
  if(allocated(psn_pot_phi_real)) then
    deallocate( psn_pot_phi_real,stat=ierr)
    call assert(ierr.eq.0,'dealloc(psn_pot_phi_real)',ierr)
  endif

  if(allocated(psn_ff_hdp_tr)) then
    deallocate( psn_ff_hdp_tr,stat=ierr)
    call assert(ierr.eq.0,'dealloc(psn_ff_hdp_tr)',ierr)
  endif

  if(allocated(psn_ff_hdp_p)) then
    deallocate( psn_ff_hdp_p,stat=ierr)
    call assert(ierr.eq.0,'dealloc(psn_ff_hdp_p)',ierr)
  endif

  if(allocated(psn_ff_1dp_tr)) then
    deallocate( psn_ff_1dp_tr,stat=ierr)
    call assert(ierr.eq.0,'dealloc(psn_ff_1dp_tr)',ierr)
  endif

  if(allocated(psn_ff_1dp_p)) then
    deallocate( psn_ff_1dp_p,stat=ierr)
    call assert(ierr.eq.0,'dealloc(psn_ff_1dp_p)',ierr)
  endif

  if(allocated(psn_bfollow_1_dx)) then
    deallocate( psn_bfollow_1_dx,stat=ierr)
    call assert(ierr.eq.0,'dealloc(psn_bfollow_1_dx)',ierr)
  endif

#else

  if(allocated(psn_E_phi_ff0)) then
#ifdef USE_PSN_TEXTURE
    nullify( t0_psn_E_phi_ff )
    nullify( t1_psn_E_phi_ff )    
#endif
    deallocate( psn_E_phi_ff0,stat=ierr)
    call assert(ierr.eq.0,'dealloc(psn_E_phi_ff0)',ierr)
    deallocate( psn_E_phi_ff1,stat=ierr)
    call assert(ierr.eq.0,'dealloc(psn_E_phi_ff1)',ierr)
  endif

#endif

  nullify(t_alias_grid_nd)

  return
  end subroutine update_host_psn_type

attributes(device) &
subroutine efield_gk_elec_gpu(i,fld,itr,p) !require fld%(r,z,phi,br,bz,bphi)
  use sml_module_gpu, only : sml_turb_efield, sml_mype, sml_turb_efield
  use grid_class_gpu, only : grid_delta_phi, grid_rhomax, grid_drho, grid_nrho                             

  use fld_module, only : fld_type
!  use ptl_module_gpu, only :  type_gpu, rhoi_gpu
  use precision_mod_gpu
  implicit none
  integer, intent(in) :: i,itr
  type(fld_type), intent(inout) :: fld
  real (kind=work_p), intent(in) :: p(3)

  integer :: ip,node,irho, iphi
  real (kind=work_p) :: wphi(0:1), wp, rho, rhon, wrho(2)
  real (kind=work_p) :: E(3),  B  !, D(3), pot,  E00(2),  ddpotdt

#ifdef USE_CALC_GRADIENT
  real (kind=work_p) :: E_phi_ff(1:3,0:1)
#endif

  iphi=floor(fld%phi/grid_delta_phi)
!  wphi(1)=(fld%phi/grid%delta_phi)  - grid%iphi_offset
  wphi(1)=(fld%phi/grid_delta_phi) - iphi
  wphi(0)=1D0 - wphi(1)

!  pot=0D0
  E=0D0
!  E00=0D0
  !D=0D0 !debug
!  ddpotdt=0D0

  !get E
 
  if(itr>0) then        
     do ip=1, 3
        node=t_alias_grid_nd(ip,itr)
        wp=p(ip)
        
           ! for electron -- rho=0 case, optimized.           
!           pot = pot + wp*wphi(0)*psn_pot_phi_ff(0,node,iphi) 
!           pot = pot + wp*wphi(1)*psn_pot_phi_ff(1,node,iphi)

#ifdef USE_CALC_GRADIENT
          call calc_E_phi_ff_gpu(node,iphi,E_phi_ff )
          E = E + wp*wphi(0)*E_phi_ff(:,0)
          E = E + wp*wphi(1)*E_phi_ff(:,1)
#else


#ifdef USE_PSN_TEXTURE
           E   = E   + wp*wphi(0)*t0_psn_E_phi_ff(:,node,iphi)
           E   = E   + wp*wphi(1)*t1_psn_E_phi_ff(:,node,iphi)
#else
           E   = E   + wp*wphi(0)*psn_E_phi_ff0(:,node,iphi)
           E   = E   + wp*wphi(1)*psn_E_phi_ff1(:,node,iphi)
#endif
#endif


!           E00 = E00 + wp*wphi(0)*psn_E00_ff(:,0,node)
!           E00 = E00 + wp*wphi(1)*psn_E00_ff(:,1,node)

!           ddpotdt = ddpotdt + wp*wphi(0)*psn_ddpotdt_phi(node,0,iphi)
!           ddpotdt = ddpotdt + wp*wphi(1)*psn_ddpotdt_phi(node,1,iphi)



     enddo
  else
!      print *, 'E-field ion invalid tr : (i,itr,mype,gid)=', i,itr,sml_mype,sp%ptl(i)%gid
  endif
    
  !E(3) was para E-field and becomes Ephi
  if(sml_turb_efield) then
     B=sqrt(fld%br**2 + fld%bz**2 + fld%bphi**2)
     E(3)=(E(3)*B- E(1)*fld%Br - E(2)*fld%Bz)/fld%Bphi

     !debug
     !D(3)=(D(3)*B- D(1)*fld%Br - D(2)*fld%Bz)/fld%Bphi
  else
     E(3)=0D0
  endif
  
  fld%Er=E(1)
  fld%Ez=E(2)
  fld%Ephi=E(3)
!  fld%Epot=pot
!  fld%Er00=E00(1)
!  fld%Ez00=E00(2) 
!  fld%ddpotdt=ddpotdt


end subroutine efield_gk_elec_gpu

end module psn_class_gpu
