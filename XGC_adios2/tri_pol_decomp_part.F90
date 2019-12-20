subroutine tri_pol_decomp_part(ntriangle,nd,nparts,tri_weight,comm,part)
#if defined(F2003)
  use iso_c_binding, only : c_float
#endif
  use assert_mod, only : assert
  implicit none
  include 'mpif.h'

  integer, parameter :: idebug = 0

  integer, intent(in) :: ntriangle, nd(3,ntriangle)
  integer, intent(in) :: tri_weight(ntriangle)
  integer, intent(in) :: comm,nparts
  integer, intent(inout) :: part(ntriangle)

  integer :: ip, ncommonnodes
  integer, allocatable, dimension(:) :: elmdist,eptr,eind
  integer, allocatable, dimension(:,:) :: elmwgt
  integer :: edgecut,numflag,wgtflag

  integer, parameter :: dp = selected_real_kind(15,30)
  real(kind=dp) :: avg_weight, dscale
#if !defined(F2003)
  integer, parameter :: c_float = selected_real_kind(6,30)
#endif
  integer, allocatable, dimension(:) :: local_part, recvbuf

  !       ---------------------
  !       setup element weights
  !       set ncon = 1 if use just the number of particles as weight
  !       set ncon = 2 if also consider number of triangles as constraint
  !       ---------------------
  integer, parameter :: ncon = 2

  real(c_float), dimension(ncon,nparts) :: tpwgts
  real(c_float), dimension(ncon) :: ubvec
  integer, dimension(0:2) :: options

  integer :: i,i1,i2,ival,icount, ierror
  integer :: lsize, ielm, iremain, ielm_start
  integer, allocatable, dimension(:) ::  isize

  integer :: nproc,irank
  integer, parameter :: root  = 0
  logical ::  is_ok


#if defined(F2003)
  !       ----------------------
  !       interface to ParMETIS
  !       see parmetis.h
  !       ----------------------
  interface

     subroutine ParMETIS_V3_PartMeshKway( elmdist, eptr, eind, elmwgt,     &
          &           wgtflag, numflag, ncon, ncommonnodes, nparts,                &
          &           tpwgts, ubvec, options, edgecut, part, comm )                &
          &           bind(c,name='ParMETIS_V3_PartMeshKway')
       use iso_c_binding, only : c_float, c_int
       implicit none

       !       -------------------------------------------------
       !       see parmetis.h, in particular, whether idxtype is
       !       defined as "short" or "int"
       !
       !       integer, parameter :: indxtype = c_short
       !       -------------------------------------------------
       integer, parameter :: indxtype = c_int

       integer(kind=indxtype) :: elmdist(*),eptr(*),eind(*),elmwgt(*)
       integer(kind=c_int) :: wgtflag,numflag,ncon,ncommonnodes,nparts
       real(c_float) :: tpwgts(*), ubvec(*)
       integer(kind=c_int) :: options(*), edgecut, part(*)
       integer(kind=c_int) :: comm
     end subroutine ParMETIS_V3_PartMeshKway

  end interface
#endif


  call MPI_COMM_SIZE( comm, nproc, ierror)
  is_ok = (ierror .eq. MPI_SUCCESS)
  call assert(is_ok,'mpi_comm_size ', ierror)

  call MPI_COMM_RANK( comm, irank, ierror)
  is_ok  = (ierror .eq. MPI_SUCCESS)
  call assert(is_ok,'mpi_comm_rank ',ierror)


  !	--------------------------------------------------
  !	roughly equal number of triangle/elements on each
  !	local processor
  !	--------------------------------------------------
  allocate(isize(nproc),elmdist(nproc+1), stat=ierror)
  is_ok = (ierror.eq.0)
  call assert(is_ok,'allocate(elmdist) ',ierror)


  lsize = int(ntriangle/nproc)
  isize(:) = lsize

  iremain = ntriangle - lsize*nproc
  is_ok = (0 .le. iremain).and.(iremain .le. nproc)
  call assert(is_ok,'invalid iremain,iremain=',iremain)


  do i=1,iremain
     isize(i) = isize(i) + 1
  enddo
  is_ok = (sum(isize(1:nproc)) .eq. ntriangle) 
  call assert(is_ok,'invalid sum(isize) ',                         &
       &              ntriangle-sum(isize(1:nproc)))


  !	---------------------------
  !	setup element distribution 
  !	---------------------------
  elmdist(1) = 1
  do i=1,nproc
     elmdist(i+1) = elmdist(i) + isize(i)
  enddo

  !       --------------------------------------------
  !       create element list eptr(:),eind(:) compatible with ParMETIS
  !       --------------------------------------------
  lsize = isize(irank+1)
  allocate( eptr(lsize+1), eind(3*lsize+1), stat=ierror)
  is_ok = (ierror.eq.0)
  call assert(is_ok,'allocate(eptr) ',ierror)

  ielm_start = elmdist(irank+1)
  do i=1,lsize
     ip = (i-1)*3 + 1
     eptr(i) = ip

     ielm = ielm_start + (i-1)
     eind(ip)   = nd(1,ielm)
     eind(ip+1) = nd(2,ielm)
     eind(ip+2) = nd(3,ielm)
  enddo
  eptr(lsize+1) = 3*lsize+1


  lsize = isize(irank+1)
  allocate( elmwgt(ncon,lsize), stat=ierror)
  is_ok = (ierror.eq.0)
  call assert(is_ok,'allocate(elmwgt) ', ierror)

  !       -------------------------------------------------------
  !       in case the weights are too large, this may lead to
  !       partitions with no triangles
  !       scale constraints by average weight to get some balance
  !       -------------------------------------------------------
  avg_weight =  0
  do i=1,ntriangle
     avg_weight = avg_weight + dble( tri_weight(i) )
  enddo
  avg_weight = avg_weight/dble(ntriangle)

  dscale = 0.1d0
  elmwgt = max(1,int(dscale*avg_weight))

  lsize = isize(irank+1)
  ielm_start = elmdist(irank+1)
  do i=1,lsize
     elmwgt(1,i) = tri_weight(ielm_start + i-1)
  enddo


  !       ---------------------------
  !       set wgtflag
  !
  !       0: no weights (vwgt and adjwgt are both NULL)
  !       1: weights on the edges only (vwgt is NULL)
  !       2: weights on the vertices only (adjwgt is NULL)
  !       3: weights on both the vertices and edges
  !       ---------------------------------

  !       ----------------------------
  !       use weights on the vertices only
  !       ----------------------------
  wgtflag = 2

  !       ----------------------------------------------------
  !       degree of connectivity, 
  !       ncommonnodes=2 triangles are connected if they share an edge
  !       ncommonnodes=1 triangles are connected if they share a vertex
  !       ----------------------------------------------------
  ncommonnodes = 2

  !       -------------------------------------------
  !       Fortran style numbering that starts from 1
  !       -------------------------------------------
  numflag = 1

  !       ------------------------------------------------------
  !       tpwgts array size of (ncon by nparts)
  !       set (1/nparts) if all sub domains are to be same size
  !       ------------------------------------------------------

  tpwgts = real(dble(1)/dble(nparts),kind=kind(tpwgts))

  !       -------------------------------------------------------------
  !       ubvec array size of ncon to specify the imbalance tolerance
  !       1 being perfect, a value of 1.05 is recommended
  !       -------------------------------------------------------------

  ubvec = 1.05

  !       -------------------------------------------------------
  !       options array of integers to pass additional parameters
  !       options(0) mean default values are used
  !       options(1) level of information to be returned, default 0
  !       options(2) random number seed, default 15
  !       -------------------------------------------------------

  options = 0
  if (idebug >= 1) then
     options(0) = 1
     options(1) = 1
  endif


  !       ---------------------------------------------------
  !       local_part is an array of size equal to the number of locally stored
  !       vertices. Upon successful completion, the partiion vector of the
  !       locally stored vertices is written to this array
  !       ---------------------------------------------------

  lsize = isize(irank+1)
  allocate( local_part( lsize ), stat=ierror)
  is_ok = (ierror .eq. 0)
  call assert( is_ok, 'allocate( local_part )', ierror)

  edgecut = 0
  call ParMETIS_V3_PartMeshKway(elmdist, eptr, eind, elmwgt,           &
       &      wgtflag, numflag, ncon, ncommonnodes, nparts,                    &
       &      tpwgts, ubvec, options, edgecut, local_part,                     &
       &      comm)


  !       -------------------------------------
  !       collect results into a global vector
  !       -------------------------------------
  part(:) = 0

  i1 = elmdist(irank+1)
  i2 = i1 + lsize-1
  part( i1:i2 ) = local_part(1:lsize)

  deallocate( local_part, stat=ierror)
  is_ok = (ierror.eq.0)
  call assert( is_ok,'deallocate(local_part) ', ierror)

  if (nproc >= 1) then
     !       -----------------------------------------------------------
     !       interface: mpi_allreduce( sendbuf, recvbuf, count,
     !                                 datatype, op, comm)
     !       -----------------------------------------------------------
     icount = ntriangle
     allocate(  recvbuf(icount), stat=ierror)
     is_ok = (ierror.eq.0)
     call assert(is_ok,'allocate(recvbuf) ', ierror)

     call MPI_ALLREDUCE( part, recvbuf, icount,                       &
          &            MPI_INTEGER, MPI_SUM, comm, ierror)
     is_ok = (ierror.eq.MPI_SUCCESS)
     call assert(is_ok,'mpi_allreduce(part) ',ierror)

     part(1:icount) = recvbuf(1:icount)

     deallocate( recvbuf, stat=ierror)
     is_ok = (ierror.eq.0)
     call assert( is_ok,'deallocate(recvbuf) ', ierror)

  endif


  !       ---------------------------------------------------------
  !       make sure numbering of part(1:ntriangle) is 0..(nparts-1)
  !       ---------------------------------------------------------
  ival = minval( part(1:ntriangle) )
  if (ival .ne. 0) then
     part(1:ntriangle) = part(1:ntriangle) - ival
  endif

  do i=1,ntriangle
     is_ok = (0 <= part(i)).and.(part(i) <= (nparts-1))
     call assert(is_ok,'tri_pol_decomp_part: invalid part(i) ',     &
          &               part(i) )
  enddo

  !       -------
  !       cleanup
  !       -------
  deallocate( elmwgt, stat=ierror)
  is_ok = (ierror.eq.0)
  call assert(is_ok,'deallocate(elmwgt) ',ierror)

  deallocate( eptr, eind, stat=ierror) 
  is_ok = (ierror.eq.0)
  call assert(is_ok,'deallocate(eptr) ',ierror)

  deallocate( isize, elmdist, stat=ierror) 
  is_ok = (ierror.eq.0)
  call assert(is_ok,'deallocate(isize) ',ierror)

  return

end subroutine tri_pol_decomp_part
