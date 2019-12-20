subroutine tri_pol_decomp(grid,ptl,part,use_ion_in)  
  use grid_class, only : grid_type
  use ptl_module, only : ptl_type, species_type
  use sml_module, only : sml_pe_per_plane  
  use assert_mod, only : assert
  implicit none
  
  include 'mpif.h'
  
  type(grid_type),intent(in) :: grid
  type(ptl_type), target :: ptl
  integer, intent(inout) :: part(:)
  logical, optional, intent(in) :: use_ion_in
  
  
  integer, allocatable, dimension(:) :: tri_weight
  integer, allocatable, dimension(:) :: recvbuf
  
  
  type(species_type), pointer :: sp
  
  real(8) :: x2(2)
  logical :: is_valid
  integer :: icount,nparts 
  integer :: i, itr 
  real(8) :: p(3)
  logical :: use_ion
  
  
  integer :: ntriangle,lsize
  integer, pointer, dimension(:,:) :: nd
  
  integer :: comm_reduceall
  integer :: comm,ierror
  logical :: is_ok
  

  ntriangle = grid%ntriangle
  nd => grid%nd


  !       ------------
  !	extra checks
  !       ------------
  is_ok = (size(part) >= ntriangle)
  call assert(is_ok,'size(part) < ntriangle ',size(part))



  !       -------------------------------
  !       compute the weight of triangle
  !
  !       default is to use number of ion particles
  !       -------------------------------
  use_ion = .true.
  if (present(use_ion_in)) use_ion = use_ion_in

  if (use_ion) then
     sp => ptl%ion
  else
     sp => ptl%elec
  endif

  !       ---------------------------------
  !       use number of particles as weight
  !       ---------------------------------
  lsize = max(1, ntriangle)
  allocate( tri_weight(lsize), stat=ierror )
  is_ok = (ierror.eq.0)
  call assert( is_ok,'allocate( tri_weight ) ',ierror)

  tri_weight(:) = 0
  do i=1,sp%num
     if (sp%gid(i) <= 0) cycle

     x2 = sp%phase(1:2,i)
     call search_tr2(grid,x2,itr,p)
     is_valid = (1 <= itr) .and. (itr <= ntriangle)
     if (is_valid) then
        tri_weight(itr) = tri_weight(itr) + 1
     endif
  enddo

  !       -------------------------------
  !       global sum to get global weight
  !       -------------------------------

  !       -----------------------------------------------------------
  !       interface: mpi_allreduce( sendbuf, recvbuf, count, 
  !                                 datatype, op, comm)
  !       -----------------------------------------------------------


  !       ============================================
  !       !!! NOTE: !!!
  !       global sum over ALL processors on ALL planes
  !       ============================================

  comm_reduceall = sml_comm

  icount = ntriangle
  allocate( recvbuf(max(1,icount)), stat=ierror )
  is_ok = (ierror.eq.0)
  call assert(is_ok,'allocate(recvbuf) ',ierror)

  call MPI_ALLREDUCE( tri_weight, recvbuf, icount,                    &
       &           MPI_INTEGER, MPI_SUM, comm_reduceall, ierror)
  is_ok = (ierror.eq.MPI_SUCCESS)
  call assert(is_ok,'mpi_allreduce(tri_weight)',ierror)

  tri_weight(1:icount) = recvbuf(1:icount)

  deallocate( recvbuf, stat=ierror )
  is_ok = (ierror.eq.0)
  call assert(is_ok,'deallocate(recvbuf) ',ierror)




  !       ---------------------------
  !       perform graph partitioning
  !       ---------------------------

  !       =============================================== 
  !       !!! NOTE: !!!!
  !       NOTE: nparts is number of processors per plane
  !       =============================================== 
  nparts = sml_pe_per_plane

  !       -------------------------------------------------------------
  !	each processor independently run a "serial" copy of parmetis
  !	this is the simplest and seems to give the best quality
  !       -------------------------------------------------------------
  comm  = MPI_COMM_SELF

  call tri_pol_decomp_part(ntriangle,nd,nparts,tri_weight,          &
       &            comm,part)

  !       -----------------------------------------------
  !       optional: 
  !       perform a broad cast to make sure all
  !       processors agree on the initial partition
  !       -----------------------------------------------

  comm = sml_comm
  call MPI_BCAST( part, ntriangle, MPI_INTEGER, 0, comm,ierror)
  is_ok = (ierror.eq.MPI_SUCCESS)
  call assert(is_ok,'mpi_bcast(part) ',ierror)


  deallocate( tri_weight, stat=ierror)
  is_ok = (ierror.eq.0)
  call assert(is_ok,'deallocate(tri_weight) ',ierror)

  return
end subroutine tri_pol_decomp
