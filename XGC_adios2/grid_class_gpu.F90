module grid_class_gpu
  use dimensions_mod_gpu
  use precision_mod_gpu
  use grid_class, only :grid_type

  use cudafor
  implicit none

  integer :: grid_nnode      !! number of nodes in a grid system
  integer :: grid_ntriangle  !! number of trianble in a grid system
  integer :: grid_iphi_offset  !! toroidal angle indexing offset
  integer :: grid_nrho

#ifdef USE_GRID_TEXTURE


  integer,texture,pointer :: grid_nd(:,:)  !! 3 node numbers that each triangle has
  integer,texture,pointer :: grid_adj(:,:) !! 3 adjacent triangle
  integer,texture,pointer :: grid_guess_table(:,:)
  integer,texture,pointer :: grid_guess_list(:)
  integer,texture,pointer :: grid_guess_xtable(:,:)
  integer,texture,pointer :: grid_guess_count(:,:)

  integer, texture,pointer :: grid_rgn(:)   !! region value for each node point
  real (kind=work_p),texture,pointer :: grid_x(:,:)
 
  real (kind=work_p),texture,pointer :: grid_mapping(:,:,:)  !! shape function coefficient



  integer,target,allocatable :: grid_nd_(:,:)  !! 3 node numbers that each triangle has
  integer,target,allocatable :: grid_adj_(:,:) !! 3 adjacent triangle
  integer,target,allocatable :: grid_guess_table_(:,:)
  integer,target,allocatable :: grid_guess_list_(:)
  integer,target,allocatable :: grid_guess_xtable_(:,:)
  integer,target,allocatable :: grid_guess_count_(:,:)

  integer, target,allocatable :: grid_rgn_(:)   !! region value for each node point
  real (kind=work_p),target,allocatable :: grid_x_(:,:)
 
  real (kind=work_p),target,allocatable :: grid_mapping_(:,:,:)  !! shape function coefficient



#define GRID_ND grid_nd_
#define GRID_ADJ grid_adj_
#define GRID_GUESS_TABLE grid_guess_table_
#define GRID_GUESS_LIST grid_guess_list_
#define GRID_GUESS_XTABLE grid_guess_xtable_
#define GRID_GUESS_COUNT grid_guess_count_
#define GRID_RGN grid_rgn_
#define GRID_X grid_x_
#define GRID_MAPPING grid_mapping_

#else
  ! grid_nd is always a target, so it can be aliased in psn_class_gpu
  integer,allocatable,target :: grid_nd(:,:)  !! 3 node numbers that each triangle has
  integer,allocatable :: grid_adj(:,:) !! 3 adjacent triangle
  integer,allocatable :: grid_guess_table(:,:)
  integer,allocatable :: grid_guess_list(:)
  integer,allocatable :: grid_guess_xtable(:,:)
  integer,allocatable :: grid_guess_count(:,:)

  integer, allocatable :: grid_rgn(:)   !! region value for each node point
  real (kind=work_p),allocatable :: grid_x(:,:)
 
  real (kind=work_p),allocatable :: grid_mapping(:,:,:)  !! shape function coefficient

#define GRID_ND grid_nd
#define GRID_ADJ grid_adj
#define GRID_GUESS_TABLE grid_guess_table
#define GRID_GUESS_LIST grid_guess_list
#define GRID_GUESS_XTABLE grid_guess_xtable
#define GRID_GUESS_COUNT grid_guess_count
#define GRID_RGN grid_rgn
#define GRID_X grid_x
#define GRID_MAPPING grid_mapping

#endif

  integer,dimension(2) :: grid_guess_n
  real (kind=work_p),dimension(2) :: grid_guess_min 
  real (kind=work_p),dimension(2) :: grid_guess_max
  real (kind=work_p),dimension(2) :: grid_guess_d 
  real (kind=work_p),dimension(2) :: grid_inv_guess_d


#ifdef USE_CALC_GRADIENT

  integer :: grid_max_v
#ifdef USE_GRID_TEXTURE
  integer, texture,pointer :: grid_v_node(:,:)
  real (kind=work_p), texture,pointer :: grid_gradx(:,:)
  real (kind=work_p), texture,pointer :: grid_grady(:,:)
  integer, texture,pointer :: grid_num_v_node(:)

  integer, target,allocatable :: grid_v_node_(:,:)
  real (kind=work_p), target,allocatable :: grid_gradx_(:,:)
  real (kind=work_p), target,allocatable :: grid_grady_(:,:)
  integer, target,allocatable :: grid_num_v_node_(:)

#define GRID_V_NODE grid_v_node_
#define GRID_GRADX grid_gradx_
#define GRID_GRADY grid_grady_
#define GRID_NUM_V_NODE grid_num_v_node_

#else
  integer, allocatable :: grid_v_node(:,:)
  real (kind=work_p), allocatable :: grid_gradx(:,:)
  real (kind=work_p), allocatable :: grid_grady(:,:)
  integer, allocatable :: grid_num_v_node(:)

#define GRID_V_NODE grid_v_node
#define GRID_GRADX grid_gradx
#define GRID_GRADY grid_grady
#define GRID_NUM_V_NODE grid_num_v_node


#endif
#endif


  real (kind=work_p) :: grid_delta_phi
  real (kind=work_p) :: grid_rhomax
  real (kind=work_p) :: grid_drho
  integer, parameter :: idebug = 0


  attributes(device) :: GRID_ND, GRID_ADJ 
  attributes(device) :: GRID_GUESS_TABLE, GRID_GUESS_LIST
  attributes(device) :: GRID_GUESS_XTABLE,  GRID_GUESS_COUNT
  attributes(device) :: GRID_RGN, GRID_X, GRID_MAPPING

  attributes(constant) :: grid_guess_n
  attributes(constant) :: grid_iphi_offset, grid_nrho  
  attributes(constant) :: grid_nnode, grid_ntriangle, grid_guess_min, &
                        grid_guess_max, grid_guess_d, grid_inv_guess_d,  &
                        grid_delta_phi, grid_rhomax, grid_drho 
#ifdef USE_CALC_GRADIENT
  attributes(constant) :: grid_max_v


  attributes(device) :: GRID_V_NODE, GRID_GRADX,                       &
                        GRID_GRADY, GRID_NUM_V_NODE 
#endif

  contains

  attributes(host) &
  subroutine check_grid_tex_size(arsize)
    use dimensions_mod_gpu
    implicit none
    integer, intent(in) :: arsize
    if (arsize >= texture_max) then
      print *, "Please recompile without USE_GRID_TEXTURE!"
    endif
    call assert(arsize < texture_max,"grid array too large to bind! size:",arsize)
    return
  end subroutine check_grid_tex_size

  attributes(host) &
  subroutine update_device_grid_type(grid )
  use assert_mod
  use sml_module, only : sml_mype
  use grid_class, only : grid_type
  type(grid_type), intent(in) :: grid
#ifdef USE_GPU_EMU
  logical, parameter :: use_cudaMemcpy = .false.
#else
  logical, parameter :: use_cudaMemcpy = .true.
#endif
  integer :: icount_nd, icount_adj, icount_gn, icount_map, icount_gt, icount_gx, icount_gc, &
             icount_list, icount_gmin, icount_gmax, icount_gd, icount_igd, icount_rg, icount_x, ierr
  integer, parameter :: iundefined = huge(0)

  integer :: icount 
  integer :: lb1,ub1,  lb2,ub2, lb3,ub3
  logical :: isvalid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!! copy grid data from CPU to GPU !!!

!  icount_nd = 3*ntriangle_dim
!  icount_adj = 3*ntriangle_dim
!  icount_map = 2*3*ntriangle_dim
!  icount_gt = guess_table_n1_dim*guess_table_n2_dim
!  icount_list = guess_list_size
!  icount_gx = guess_table_n1_dim*guess_table_n2_dim 
!  icount_gc = guess_table_n1_dim*guess_table_n2_dim
   
!!   allocate(grid_nd(3,ntriangle_dim), &
!!            grid_adj(3,ntriangle_dim),
!!            grid_guess_table(guess_table_n1_dim,guess_table_n2_dim), &
!!            grid_guess_list(guess_list_size), &
!!            grid_guess_xtable(guess_table_n1_dim,guess_table_n2_dim), &
!!            grid_guess_count(guess_table_n1_dim,guess_table_n2_dim), &
!!            grid_mapping(2,3,ntriangle_dim), &

  
! =====================
! allocate grid_nd(:,:)
! =====================

  if (.not.allocated(GRID_ND)) then
    lb1 = lbound(grid%nd,1)
    lb2 = lbound(grid%nd,2)

    ub1 = ubound(grid%nd,1)
    ub2 = ubound(grid%nd,2)

    allocate( GRID_ND( lb1:ub1, lb2:ub2 ), stat=ierr)
    call assert(ierr.eq.0,'alloc(GRID_ND)',ierr)

  endif

  isvalid = size(GRID_ND).eq.size(grid%nd)
  call assert(isvalid,'invalid size(GRID_ND)',size(GRID_ND))



! =====================
! allocate GRID_ADJ(:,:)
! =====================

  if (.not.allocated(GRID_ADJ)) then
    lb1 = lbound(grid%adj,1)
    lb2 = lbound(grid%adj,2)

    ub1 = ubound(grid%adj,1)
    ub2 = ubound(grid%adj,2)

    allocate( GRID_ADJ( lb1:ub1, lb2:ub2 ), stat=ierr)
    call assert(ierr.eq.0,'alloc(GRID_ADJ)',ierr)

  endif

  isvalid = size(GRID_ADJ).eq.size(grid%adj)
  call assert(isvalid,'invalid size(GRID_ADJ)',size(GRID_ADJ))



! =====================
! allocate GRID_GUESS_TABLE(:,:)
! =====================

  if (.not.allocated(GRID_GUESS_TABLE)) then
    lb1 = lbound(grid%guess_table,1)
    lb2 = lbound(grid%guess_table,2)

    ub1 = ubound(grid%guess_table,1)
    ub2 = ubound(grid%guess_table,2)

    allocate( GRID_GUESS_TABLE( lb1:ub1, lb2:ub2 ), stat=ierr)
    call assert(ierr.eq.0,'alloc(GRID_GUESS_TABLE)',ierr)

  endif

  isvalid = size(GRID_GUESS_TABLE).eq.size(grid%guess_table)
  call assert(isvalid,'invalid size(GRID_GUESS_TABLE)',size(GRID_GUESS_TABLE))



! =====================
! allocate GRID_GUESS_XTABLE(:,:)
! =====================

  if (.not.allocated(GRID_GUESS_XTABLE)) then
    lb1 = lbound(grid%guess_xtable,1)
    lb2 = lbound(grid%guess_xtable,2)

    ub1 = ubound(grid%guess_xtable,1)
    ub2 = ubound(grid%guess_xtable,2)

    allocate( GRID_GUESS_XTABLE( lb1:ub1, lb2:ub2 ), stat=ierr)
    call assert(ierr.eq.0,'alloc(GRID_GUESS_XTABLE)',ierr)

  endif

  isvalid = size(GRID_GUESS_XTABLE).eq.size(grid%guess_xtable)
  call assert(isvalid,'invalid size(GRID_GUESS_XTABLE)',size(GRID_GUESS_XTABLE))




! =====================
! allocate GRID_GUESS_COUNT(:,:)
! =====================

  if (.not.allocated(GRID_GUESS_COUNT)) then
    lb1 = lbound(grid%guess_count,1)
    lb2 = lbound(grid%guess_count,2)

    ub1 = ubound(grid%guess_count,1)
    ub2 = ubound(grid%guess_count,2)

    allocate( GRID_GUESS_COUNT( lb1:ub1, lb2:ub2 ), stat=ierr)
    call assert(ierr.eq.0,'alloc(GRID_GUESS_COUNT)',ierr)

  endif

  isvalid = size(GRID_GUESS_COUNT).eq.size(grid%guess_count)
  call assert(isvalid,'invalid size(GRID_GUESS_COUNT)',size(GRID_GUESS_COUNT))

! =====================
! allocate GRID_GUESS_LIST(:)
! =====================

  if (.not.allocated(GRID_GUESS_LIST)) then
    lb1 = lbound(grid%guess_list,1)

    ub1 = ubound(grid%guess_list,1)

    allocate( GRID_GUESS_LIST( lb1:ub1 ), stat=ierr)
    call assert(ierr.eq.0,'alloc(GRID_GUESS_LIST)',ierr)

  endif

  isvalid = size(GRID_GUESS_LIST).eq.size(grid%guess_list)
  call assert(isvalid,'invalid size(GRID_GUESS_LIST)',size(GRID_GUESS_LIST))


! =====================
! allocate GRID_MAPPING(:,:,:)
! =====================

  if (.not.allocated(GRID_MAPPING)) then
    lb1 = lbound(grid%mapping,1)
    lb2 = lbound(grid%mapping,2)
    lb3 = lbound(grid%mapping,3)

    ub1 = ubound(grid%mapping,1)
    ub2 = ubound(grid%mapping,2)
    ub3 = ubound(grid%mapping,3)

    allocate( GRID_MAPPING( lb1:ub1, lb2:ub2, lb3:ub3 ), stat=ierr)
    call assert(ierr.eq.0,'alloc(GRID_MAPPING)',ierr)

  endif

  isvalid = size(GRID_MAPPING).eq.size(grid%mapping)
  call assert(isvalid,'invalid size(GRID_MAPPING)',size(GRID_MAPPING))







! =====================
! allocate GRID_RGN(:)
! =====================

  if (.not.allocated(GRID_RGN)) then
    lb1 = lbound(grid%rgn,1)

    ub1 = ubound(grid%rgn,1)

    allocate( GRID_RGN( lb1:ub1 ), stat=ierr)
    call assert(ierr.eq.0,'alloc(GRID_RGN)',ierr)

  endif

  isvalid = size(GRID_RGN).eq.size(grid%rgn)
  call assert(isvalid,'invalid size(GRID_RGN)',size(GRID_RGN))

! =====================
! allocate GRID_X(:)
! =====================

  if (.not.allocated(GRID_X)) then
    lb1 = lbound(grid%x,1)
    lb2 = lbound(grid%x,2)

    ub1 = ubound(grid%x,1)
    ub2 = ubound(grid%x,2)

    allocate( GRID_X( lb1:ub1, lb2:ub2 ), stat=ierr)
    call assert(ierr.eq.0,'alloc(GRID_X)',ierr)

  endif

  isvalid = size(GRID_X).eq.size(grid%x)
  call assert(isvalid,'invalid size(GRID_X)',size(GRID_X))

#ifdef USE_CALC_GRADIENT

! =====================
! allocate GRID_V_NODE(:,:)
! =====================

  if (.not.allocated(GRID_V_NODE)) then
    lb1 = lbound(grid%v_node,1)
    lb2 = lbound(grid%v_node,2)

    ub1 = ubound(grid%v_node,1)
    ub2 = ubound(grid%v_node,2)

    allocate( GRID_V_NODE( lb1:ub1, lb2:ub2 ), stat=ierr)
    call assert(ierr.eq.0,'alloc(GRID_V_NODE)',ierr)

  endif

  isvalid = size(GRID_V_NODE).eq.size(grid%v_node)
  call assert(isvalid,'invalid size(GRID_V_NODE)',size(GRID_V_NODE))


! =====================
! allocate GRID_GRADX(:,:)
! =====================

  if (.not.allocated(GRID_GRADX)) then
    lb1 = lbound(grid%gradx,1)
    lb2 = lbound(grid%gradx,2)

    ub1 = ubound(grid%gradx,1)
    ub2 = ubound(grid%gradx,2)

    allocate( GRID_GRADX( lb1:ub1, lb2:ub2 ), stat=ierr)
    call assert(ierr.eq.0,'alloc(GRID_GRADX)',ierr)

  endif

  isvalid = size(GRID_GRADX).eq.size(grid%gradx)
  call assert(isvalid,'invalid size(GRID_GRADX)',size(GRID_GRADX))


! =====================
! allocate GRID_GRADY(:,:)
! =====================

  if (.not.allocated(GRID_GRADY)) then
    lb1 = lbound(grid%grady,1)
    lb2 = lbound(grid%grady,2)

    ub1 = ubound(grid%grady,1)
    ub2 = ubound(grid%grady,2)

    allocate( GRID_GRADY( lb1:ub1, lb2:ub2 ), stat=ierr)
    call assert(ierr.eq.0,'alloc(GRID_GRADY)',ierr)

  endif

  isvalid = size(GRID_GRADY).eq.size(grid%grady)
  call assert(isvalid,'invalid size(GRID_GRADY)',size(GRID_GRADY))


! =====================
! allocate GRID_NUM_V_NODE(:)
! =====================

  if (.not.allocated(GRID_NUM_V_NODE)) then
    lb1 = lbound(grid%num_v_node,1)
    ub1 = ubound(grid%num_v_node,1)

    allocate( GRID_NUM_V_NODE( lb1:ub1 ), stat=ierr)
    call assert(ierr.eq.0,'alloc(GRID_NUM_V_NODE)',ierr)

  endif

  isvalid = size(GRID_NUM_V_NODE).eq.size(grid%num_v_node)
  call assert(isvalid,'invalid size(GRID_NUM_V_NODE)',size(GRID_NUM_V_NODE))


#endif



  grid_nnode = grid%nnode
  grid_ntriangle = grid%ntriangle
  grid_iphi_offset = grid%iphi_offset
  grid_nrho = grid%nrho
  grid_rhomax = grid%rhomax
  grid_drho = grid%drho
  grid_delta_phi = grid%delta_phi

#ifdef USE_CALC_GRADIENT
  grid_max_v = size(grid%v_node,dim=1)
#endif

  if(use_cudaMemcpy) then


     icount_nd = size(GRID_ND)
     ierr = cudaMemcpy( GRID_ND, grid%nd,    &
               icount_nd, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(GRID_ND)',ierr)


     icount_adj = size(GRID_ADJ)
     ierr = cudaMemcpy( GRID_ADJ, grid%adj,  &
               icount_adj, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(GRID_ADJ)',ierr)

     icount_map = size(GRID_MAPPING)
     ierr = cudaMemcpy( GRID_MAPPING, grid%mapping,    &
              icount_map, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(GRID_MAPPING)',ierr)


     icount_rg = size(GRID_RGN)
     ierr = cudaMemcpy( GRID_RGN, grid%rgn,    &
               icount_rg, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(GRID_RGN)',ierr)

     icount_x = size(GRID_X)
     ierr = cudaMemcpy( GRID_X, grid%x,    &
               icount_x, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(GRID_X)',ierr)


     icount_gt = size(GRID_GUESS_TABLE)
     ierr = cudaMemcpy( GRID_GUESS_TABLE, grid%guess_table,    &
               icount_gt, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(GRID_GUESS_TABLE)',ierr)


     icount_gx = size(GRID_GUESS_XTABLE)
     ierr = cudaMemcpy( GRID_GUESS_XTABLE, grid%guess_xtable,    &
               icount_gx, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(GRID_GUESS_XTABLE)',ierr)


     icount_gc = size(GRID_GUESS_COUNT)
     ierr = cudaMemcpy( GRID_GUESS_COUNT, grid%guess_count,    &
               icount_gc, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(GRID_GUESS_COUNT)',ierr)


     icount_list = size(GRID_GUESS_LIST)
     ierr = cudaMemcpy( GRID_GUESS_LIST, grid%guess_list,  &
               icount_list, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(GRID_GUESS_LIST)',ierr)



     icount = size(grid_guess_n)
     ierr = cudaMemcpy( grid_guess_n, grid%guess_n,    &
               icount, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(grid_guess_n)',ierr)

     icount = size(grid_guess_min)
     ierr = cudaMemcpy( grid_guess_min, grid%guess_min,    &
               icount, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(grid_guess_min)',ierr)


     icount = size(grid_guess_max)
     ierr = cudaMemcpy( grid_guess_max, grid%guess_max,    &
               icount, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(grid_guess_max)',ierr)


     icount = size(grid_guess_d)
     ierr = cudaMemcpy( grid_guess_d, grid%guess_d,    &
               icount, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(grid_guess_d)',ierr)


     icount = size(grid_inv_guess_d)
     ierr = cudaMemcpy( grid_inv_guess_d, grid%inv_guess_d,    &
               icount, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(grid_inv_guess_d)',ierr)

#ifdef USE_CALC_GRADIENT
     icount = size(GRID_V_NODE)
     ierr = cudaMemcpy( GRID_V_NODE, grid%v_node,    &
               icount, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(GRID_V_NODE)',ierr)

     icount = size(GRID_GRADX)
     ierr = cudaMemcpy( GRID_GRADX, grid%gradx,    &
               icount, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(GRID_GRADX)',ierr)

     icount = size(GRID_GRADY)
     ierr = cudaMemcpy( GRID_GRADY, grid%grady,    &
               icount, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(GRID_GRADY)',ierr)

     icount = size(GRID_NUM_V_NODE)
     ierr = cudaMemcpy( GRID_NUM_V_NODE, grid%num_v_node,    &
               icount, cudaMemcpyHostToDevice)
     call assert(ierr.eq.0,'cudaMemcpy(GRID_NUM_V_NODE)',ierr)


#endif

  else
     GRID_ND = grid%nd
     GRID_ADJ = grid%adj

     GRID_MAPPING = grid%mapping
     GRID_X = grid%x
     GRID_RGN = grid%rgn

     GRID_GUESS_TABLE = grid%guess_table
     GRID_GUESS_XTABLE = grid%guess_xtable
     GRID_GUESS_COUNT = grid%guess_count
     GRID_GUESS_LIST = grid%guess_list


     grid_guess_n = grid%guess_n
     grid_guess_min = grid%guess_min
     grid_guess_max = grid%guess_max
     grid_guess_d = grid%guess_d
     grid_inv_guess_d = grid%inv_guess_d

#ifdef USE_CALC_GRADIENT

     GRID_V_NODE = grid%v_node
     GRID_GRADX = grid%gradx
     GRID_GRADY = grid%grady
     GRID_NUM_V_NODE = grid%num_v_node

#endif

  endif

#ifdef USE_GRID_TEXTURE
   call check_grid_tex_size(size(grid_nd_))
   call check_grid_tex_size(size(grid_adj_))
   call check_grid_tex_size(size(grid_guess_table_))
   call check_grid_tex_size(size(grid_guess_list_))
   call check_grid_tex_size(size(grid_guess_xtable_))   
   call check_grid_tex_size(size(grid_guess_count_))
   call check_grid_tex_size(size(grid_rgn_))
   call check_grid_tex_size(size(grid_x_))   
   call check_grid_tex_size(size(grid_mapping_))
   grid_nd => grid_nd_
   grid_adj => grid_adj_
   grid_guess_table => grid_guess_table_
   grid_guess_list => grid_guess_list_
   grid_guess_xtable => grid_guess_xtable_
   grid_guess_count => grid_guess_count_
   grid_rgn => grid_rgn_
   grid_x => grid_x_
   grid_mapping => grid_mapping_

#ifdef USE_CALC_GRADIENT
   grid_v_node => grid_v_node_
   grid_gradx => grid_gradx_
   grid_grady => grid_grady_
   grid_num_v_node => grid_num_v_node_
#endif


#endif

  return
  end subroutine update_device_grid_type


  attributes(host) &
  subroutine update_host_grid_type()

  integer :: ierr

  if (allocated(GRID_ND)) then
    deallocate(GRID_ND,stat=ierr)
    call assert(ierr.eq.0,'dealloc(GRID_ND)',ierr)
  endif


  if (allocated(GRID_ADJ)) then
    deallocate(GRID_ADJ,stat=ierr)
    call assert(ierr.eq.0,'dealloc(GRID_ADJ)',ierr)
  endif



  if (allocated(GRID_GUESS_TABLE)) then
    deallocate(GRID_GUESS_TABLE,stat=ierr)
    call assert(ierr.eq.0,'dealloc(GRID_GUESS_TABLE)',ierr)
  endif


  if (allocated(GRID_GUESS_XTABLE)) then
    deallocate(GRID_GUESS_XTABLE,stat=ierr)
    call assert(ierr.eq.0,'dealloc(GRID_GUESS_XTABLE)',ierr)
  endif


  if (allocated(GRID_GUESS_COUNT)) then
    deallocate(GRID_GUESS_COUNT,stat=ierr)
    call assert(ierr.eq.0,'dealloc(GRID_GUESS_COUNT)',ierr)
  endif



  if (allocated(GRID_GUESS_LIST)) then
    deallocate(GRID_GUESS_LIST,stat=ierr)
    call assert(ierr.eq.0,'dealloc(GRID_GUESS_LIST)',ierr)
  endif

  if (allocated(GRID_MAPPING)) then
    deallocate(GRID_MAPPING,stat=ierr)
    call assert(ierr.eq.0,'dealloc(GRID_MAPPING)',ierr)
  endif


  if (allocated(GRID_RGN)) then
    deallocate(GRID_RGN,stat=ierr)
    call assert(ierr.eq.0,'dealloc(GRID_RGN)',ierr)
  endif


  if (allocated(GRID_X)) then
    deallocate(GRID_X,stat=ierr)
    call assert(ierr.eq.0,'dealloc(GRID_X)',ierr)
  endif

#ifdef USE_GRID_TEXTURE
  nullify(grid_nd)
  nullify(grid_adj)

  nullify(grid_x)
  nullify(grid_rgn)
  nullify(grid_mapping)

  nullify(grid_guess_table)
  nullify(grid_guess_xtable)
  nullify(grid_guess_count)
  nullify(grid_guess_list)
#endif

  return
  end subroutine update_host_grid_type


  !! These search routines rely on texture pointers for
  !! performance, and thus need to be placed in the
  !! same module where the pointers are declared
  attributes(device) &
  subroutine t_coeff_gpu(itr,x,p)
    use precision_mod_gpu
    implicit none
    integer, intent(in) :: itr
    real (kind=work_p), intent(in) :: x(2)
    real (kind=work_p), intent(out) :: p(3)

    integer :: nd
    real (kind=work_p) :: dx(2)

    ! dx=x - grid%x(:,grid%nd(3,itr))
    dx(1:2) = x(1:2) - grid_mapping(1:2,3,itr)
    p(1:2)= grid_mapping(1:2,1,itr)*dx(1) + grid_mapping(1:2,2,itr)*dx(2)
    p(3)=1.0d0 - p(1) - p(2)


  end subroutine t_coeff_gpu

  attributes(device) &
  subroutine search_tr_check_guess_gpu(x,init,itr,p)
    use precision_mod_gpu
    use cudafor
    implicit none
    real (kind=work_p), intent(in) :: x(2)
    integer, intent(in) :: init
    integer, intent(out) :: itr
    real (kind=work_p), intent(out) :: p(3)
    !debug
    real (kind=work_p), parameter :: zero = 0.0d0
    real (kind=work_p) :: mnvp, mxvp

    itr=init
    if(init<=0 .or. init>grid_ntriangle) then
!       print *, 'invalid guess. init=1 is used instead',init
       itr=1
    endif


    call t_coeff_gpu(itr,x,p)

    mnvp = min(p(1),min(p(2),p(3)))
    mxvp = max(p(1),max(p(2),p(3)))
    if ( mnvp < zero .or.  mxvp > 1D0) then
      itr = -1
    endif
    return
  end subroutine search_tr_check_guess_gpu

attributes(device) &
subroutine search_tr2_gpu(  xy, itr, p )
  use precision_mod_gpu
  implicit none
  real(kind=work_p) :: xy(2)
  integer :: itr
  real(kind=work_p) :: p(3)

  real(kind=work_p), parameter :: zero = 0.0d0
  real(kind=work_p), parameter :: eps = 10.0d0*epsilon(zero)
  integer :: ij(2), istart,iend, k, itrig
  integer :: i,j,  ilo,ihi,  jlo,jhi
  real(kind=work_p) ::  dx(2), pmin, pmax, dp,  minval_p
  logical :: is_found


  ilo = lbound( grid_guess_table, 1 )
  ihi = ubound( grid_guess_table, 1 )

  jlo = lbound( grid_guess_table, 2 )
  jhi = ubound( grid_guess_table, 2 )

  ij = (xy - grid_guess_min)*grid_inv_guess_d + 1
  i = max(ilo, min(ihi, ij(1)) )
  j = max(jlo, min(jhi, ij(2)) )


  istart = grid_guess_xtable(i,j)
  iend = istart + grid_guess_count(i,j) - 1


  itr = -1
  do k=istart,iend
     itrig = grid_guess_list(k)
     ! call t_coeff( grid, itrig, xy, p )

    dx(1:2) = xy(1:2) - grid_mapping(1:2,3,itrig)
    p(1:2)= grid_mapping(1:2,1,itrig)*dx(1) +                          &
            grid_mapping(1:2,2,itrig)*dx(2)
    p(3)=1.0d0 - p(1) - p(2)

    minval_p = min(p(1),min(p(2),p(3)))
    if (minval_p .ge. -eps) then
        itr = itrig
        exit
     endif
  enddo

  return
end subroutine search_tr2_gpu


end module grid_class_gpu

#undef GRID_V_NODE 
#undef GRID_GRADX 
#undef GRID_GRADY 
#undef GRID_NUM_V_NODE 


#undef GRID_ND 
#undef GRID_ADJ 
#undef GRID_GUESS_TABLE 
#undef GRID_GUESS_LIST 
#undef GRID_GUESS_XTABLE 
#undef GRID_GUESS_COUNT 
#undef GRID_RGN 
#undef GRID_X 
#undef GRID_MAPPING 

