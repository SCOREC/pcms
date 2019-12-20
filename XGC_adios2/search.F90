module grid_class
  use boundary_class
  use mat_class
  integer, parameter :: grid_rgn_wall = 100
  !****** poloidal decomposition **********
  ! defined here for dependancies
  type decomp_type
     integer :: ntr               ! local triangle number
     integer :: nn1, nn2, nn3 ! local node number (1: inner, 2: +owned, 3: +ghost)
     integer, allocatable :: l2g_tr(:), l2g_node(:)
     !
     integer, allocatable :: gtr_2_pe_ltr(:,:)
     integer, allocatable :: gnd_2_pe_lnd(:,:)
     integer, allocatable :: lnd_2_pe_lnd(:,:)
     ! list of PEs that contains a node as a ghost node (except self) -- need?
     integer, allocatable :: npe_ghost(:) ! (l_nn1+1:l_nn2) -->
     integer :: npe_ghost_max
     integer, allocatable :: pes_ghost(:,:)   !  (l_npe_max, l_nn1+1:l_nn2)
  end type decomp_type


  type grid_type
     integer :: nnode      !! number of nodes in a grid system
     integer :: ntriangle  !! number of trianble in a grid system

     !for nodes
     integer, allocatable :: rgn(:)   !! region value for each node point
     real (kind=8), allocatable :: psi(:)  !! psi value of node point (mabye need to change to a integer value)
     real (kind=8), allocatable :: x(:,:)  !! R-Z Positon of each node point

     !for triangle
     integer, allocatable :: nd(:,:)  !! 3 node numbers that each triangle has
     integer, allocatable :: adj(:,:) !! 3 adjacent triangle

     real (kind=8), allocatable :: mapping(:,:,:)  !! shape function coefficient

     !for node->triangle structure
     integer, allocatable :: tr_node(:,:)  !! triangle index that has a given node point
     integer, allocatable :: num_t_node(:) !! number of trangle for each node, triangle index for each node

     ! --------------------------------------------------------------------------
     ! gradient computed as  weighted sum of potential values from neighbor nodes
     !   potv(1:k) = pot( v_node(1:k,i) )
     !  Ex = sum( gradx( 1:k, i) * potv(1:k) )
     !  Ey = sum( grady( 1:k, i) * potv(1:k) )
     !  where k = num_v_node(i)
     ! -------------------------------------------------------------------------
     !for node->node structure and computing gradient at vertex
     integer, allocatable :: v_node(:,:) !!  v_node(1:max_t, 1:nnodes)
     integer, allocatable :: num_v_node(:) !! num_v_node(1:nnodes)
     real(kind=8), allocatable :: gradx(:,:), grady(:,:)  !!  gradx(1:max_t,1:nnodes)


     !rh For FD gradient in (R,Z) plane --->
     ! For electric field calculation
     real (kind=8), allocatable :: unit_vecs(:,:,:)  !! for the Psi and theta unit vectors in (R,Z) coords
     integer, allocatable :: basis(:)  !! 0 for Psi-theta, 1 for R-Z
     type (mat_type) :: gradientx, gradienty  !! For the gradient operation

     !rh For poloidal smoothing
     type (mat_type) :: smooth_pol


     !guess table
     integer :: guess_n(2)             !! number of array size that
     real (kind=8) :: guess_min(2), guess_max(2), guess_d(2), inv_guess_d(2)  !! min, max, delta for guess table in RZ space
     integer, allocatable :: guess_table(:,:)        !! guess value for triangle search

     !phi direction
     integer :: iphi_offset  !! toroidal angle indexing offset
     integer :: nphi      !! number of toroidal cross section of one PE
     real (kind=8) :: delta_phi,phimin,phimax       !! Delta toroidal angle


     !rho
     integer :: nrho
     real (8) :: rhomax, drho


     integer :: npsi00
     real (8) :: psi00min, psi00max, dpsi00

     ! 2d <--> 1d00 conversion
     type(mat_type) :: cnv_2d_00
     real (8), allocatable :: cnv_norm_2d(:),cnv_norm_1d00(:)

     !node volume
     real (kind=8),allocatable ::  surf_vol(:), inv_node_vol(:),node_area(:),node_vol(:)
     !triangle volume
     real (kind=8),allocatable ::  tr_vol(:),tr_area(:)
     !temporory variables for calculations - number of node point
     real (kind=8), allocatable  :: rtmp1(:),rtmp2(:)

     ! field following
     real (8), allocatable :: node_vol_ff(:,:)
     real (8), allocatable :: node_vol_nearest(:)

     !
     integer, pointer, dimension(:) :: guess_list
     integer, pointer, dimension(:,:) :: guess_xtable, guess_count

     real (kind=8), allocatable :: bfield(:,:)
     real (kind=8), allocatable :: v_curv(:,:), v_gradb(:,:)
     real (kind=8), allocatable :: nb_curl_nb(:)

#ifdef XGC1_EM
     real (kind=8), allocatable :: gradpsi(:,:), absgradpsi(:), dene(:), tempe(:)
     real (kind=8), allocatable :: B_cross_gradB(:,:), j0_par(:), curl_nb(:,:)
     real (kind=8), allocatable :: tearing_drive(:,:), tearing_drive2(:)
     type(mat_type) :: kink_mat  !! (b x grad(j0/eB)).grad for kink drive
#endif

     type(decomp_type) :: p

     ! More geometric quantities:
     real (kind=8), allocatable :: epspar(:,:)  ! Neoclassical minor/major radius and inv. aspect ratio
     real (kind=8), allocatable :: qsafety(:), qsafety2(:)   ! Safety factor
     real (kind=8), allocatable :: trapped(:)   ! Trapped particle fraction
     integer :: nsurf
     integer, allocatable :: itheta0(:), ntheta0(:)

  end type grid_type

contains

  !! initialize adjacent triangle
  subroutine init_triangle(grid)
    implicit none
    type(grid_type) :: grid
    integer :: i,j,k,l,kp,lp
    real (kind=8) :: dx1(2),dx2(2),det
    integer :: maxnum,nd
!
    integer :: ndsize, jj, njj, iv, kk, jtrig
    integer, allocatable, dimension(:) :: jlist
    logical, allocatable, dimension(:) :: is_seen
    integer, dimension(4) :: ilist4, jlist4
    logical :: isfound
    integer :: max_v, max_t, nnode, ntri

    ! mapping matrix init
    do i=1, grid%ntriangle

       dx1=grid%x(:,grid%nd(1,i)) - grid%x(:,grid%nd(3,i))
       dx2=grid%x(:,grid%nd(2,i)) - grid%x(:,grid%nd(3,i))

       det=1./( dx1(1)*dx2(2) - dx2(1)*dx1(2))

       grid%mapping(1,1,i)=dx2(2)*det
       grid%mapping(1,2,i)=-dx2(1)*det
       grid%mapping(2,1,i)=-dx1(2)*det
       grid%mapping(2,2,i)=dx1(1)*det


       grid%mapping(:,3,i)=grid%x(:,grid%nd(3,i))
    enddo

    ! set node->triangle array
    allocate(grid%num_t_node(grid%nnode))
    do i=1,2
       grid%num_t_node(:)=0
       do j=1, grid%ntriangle
          do k=1,3
             nd=grid%nd(k,j)
             grid%num_t_node( nd )= grid%num_t_node( nd ) + 1
             if(i==2) grid%tr_node(grid%num_t_node(nd),nd)=j
          enddo
       enddo
       if(i==1) then
          maxnum=maxval(grid%num_t_node(:))
          allocate( grid%tr_node(maxnum,grid%nnode) )
       endif
    enddo

!   ----------------------------------------------
!   use temporary storage jlist(:) and is_seen(:)
!   to speed up search for adjacent triangle
!   ----------------------------------------------

    allocate( is_seen(grid%ntriangle), jlist( 3*maxnum + 1 ) )
    is_seen(:) = .false.
    jlist(:) = 0

    grid%adj=0
    ! find adjacent triangle  ->To do :  make it efficient using tr_node later
    do i=1, grid%ntriangle
       ! tr1=>grid%triangle(i)


!      -------------------------------------------
!      store potential list of neighbor triangles
!      in "jlist". Use is_seen(:) to avoid duplicates
!      -------------------------------------------
       is_seen(i) = .true.

       njj = 0
       do k=1,3
         iv = grid%nd(k,i)
         ndsize = grid%num_t_node(iv)

         do kk=1,ndsize
          jtrig = grid%tr_node(kk,iv)
          if (.not.is_seen( jtrig )) then
            njj = njj + 1
            jlist(njj) = jtrig
            is_seen(jtrig) = .true.
          endif
         enddo
       enddo

!      ---------------------------------------
!      reset boolean vector for next iteration
!      ---------------------------------------
       is_seen(i) = .false.
       is_seen( jlist(1:njj) ) = .false.

       ilist4(1:3) = grid%nd(1:3,i)
       ilist4(4) = ilist4(1)

       do jj=1,njj
          j = jlist(jj)

          jlist4(1:3) = grid%nd(1:3,j)
          jlist4(4) = jlist4(1)

          if(i/=j) then
             do k=1,3
                kp= k+1
                do l=1,3
                   lp=l+1
                   if ((ilist4(k)+ilist4(kp)).ne.                              &
     &                 (jlist4(l)+jlist4(lp))) then
                      cycle
                   endif

                   isfound = ((ilist4(k)  .eq.jlist4(l)).and.                 &
     &                        (ilist4(kp).eq.jlist4(lp))) .or.                &
     &                       ((ilist4(k)  .eq.jlist4(lp)).and.                &
     &                        (ilist4(kp).eq.jlist4(l)))
                   if (isfound) then
                      grid%adj(mod(kp,3)+1,i)=j
                      exit
                   endif
                enddo
             enddo
          end if
       enddo
    enddo

  deallocate( is_seen, jlist )

! ----------------
! setup node->node
! ----------------
  max_t = size(grid%tr_node,dim=1)
  max_v = max_t + 2
  nnode = grid%nnode
  ntri = grid%ntriangle

  allocate( grid%num_v_node( nnode ), &
            grid%v_node(max_v,nnode) )
  call gen_v_node(nnode, ntri, max_t, max_v,       &
        grid%nd, grid%num_t_node, grid%tr_node,     &
        grid%num_v_node, grid%v_node )



  end subroutine init_triangle

  !! Coefficient calculation
  subroutine t_coeff(grid,itr,x,p)
    implicit none
    type(grid_type) , intent(in) :: grid
    integer, intent(in) :: itr
    real (kind=8), intent(in) :: x(2)
    real (kind=8), intent(out) :: p(3)

    integer :: nd
    real (kind=8) :: dx(2)

    ! dx=x - grid%x(:,grid%nd(3,itr))
    dx(1:2) = x(1:2) - grid%mapping(1:2,3,itr)
    p(1:2)= grid%mapping(1:2,1,itr)*dx(1) + grid%mapping(1:2,2,itr)*dx(2)
    p(3)=1.0d0 - p(1) - p(2)


  end subroutine t_coeff

  !! search triangle with guess table
  subroutine search_tr(grid,x,itr,p)
    implicit none
    type (grid_type) , intent(in) :: grid  ! grid system
    real (kind=8), intent(in) :: x(2)    ! input position
    integer, intent(out) :: itr !triangle number which contains x
    real (kind=8), intent(out) :: p(3)   !weiting value
    integer :: init,count

    integer :: itr2
    real(kind=8), dimension(3) :: p2

    ! initial guess
    call guess(grid,x,init)

    ! search
    call search_tr_with_guess(grid,x,init,itr,p,count)
    ! error message for debug only
    !if(itr<0) then
    !   print *,"Search failed", x,init,p
    !endif


  end subroutine search_tr

  !!return initial guess value for triangle search
  subroutine guess(grid,x,init)
    implicit none
    type(grid_type), intent(in) :: grid
    real (kind=8), intent(in) :: x(2)
    integer, intent(out) :: init
    integer :: i(2)

    i= (x-grid%guess_min)*grid%inv_guess_d +1
    !error message for debug only
!    if(i(1)<=0 .or. i(1)>grid%guess_n(1)) then
!       print *, 'Invaild number for guess table- R',i(1),x(1),grid%guess_min(1),grid%guess_max(1)
!    endif
!    if(i(2)<=0 .or. i(2)>grid%guess_n(2)) then
!       print *, 'Invaild number for guess table- Z',i(2),x(2),grid%guess_min(2),grid%guess_max(2)
!    endif

    i=min(max(i,1),grid%guess_n)
    init=grid%guess_table(i(1),i(2))

  end subroutine guess

  !! initialize guess table
  subroutine init_guess_table(grid)
    implicit none
    type(grid_type) :: grid
    integer :: n(2),i,j,k(2),l,init,itr,count
    real (kind=8) :: x(2),p(3)
    real (kind=8) :: dist, min_dist
    integer ::  min_dist_node,nd,itr2,find
    ! grid%guess%n, min, max  should be initialized

    logical, parameter :: use_original = .false.
    integer, parameter :: itr_init = 1
    logical :: is_found
    integer, dimension(2) :: ij,ijmin,ijmax
    real (kind=8), dimension(2) :: xy,xy1,xy2,xy3,xyc,xymin,xymax
    real (kind=8) :: xmin,xmax,ymin,ymax
    integer :: ilo,ihi,jlo,jhi,i1,j1,i2,j2
    logical, parameter :: do_refine = .false.

    logical :: isok
    integer :: it,ierr
    real(kind=8), allocatable, dimension(:) :: parea
    integer, allocatable, dimension(:) :: iperm

    real(kind=8), parameter :: zero = 0.0d0



    n = grid%guess_n
    grid%guess_d = (grid%guess_max - grid%guess_min)/n
    grid%inv_guess_d = 1D0/grid%guess_d

    allocate( grid%guess_table(n(1),n(2)) )
    grid%guess_table=0 ! added for safety 2006/10/11

    call init_guess_list(grid)

    if (use_original) then

       ! initialize
       init=1

       do i=1, n(1)
          do j=1, n(2)
             ! find most probable triangle index
             k(1)=i
             k(2)= mod(i,2)*j + mod(i+1,2)*(n(2)-j+1)
             x=grid%guess_d*(real(k)-0.5)+grid%guess_min


             !find nearest node point
             min_dist=1D50
             do nd=1, grid%nnode
                dist=(grid%x(1,nd)-x(1))**2 + (grid%x(2,nd)-x(2))**2
                if(min_dist > dist) then
                   min_dist=dist
                   min_dist_node=nd
                endif
             enddo

             !search the triangles that has this node point
             find=0
             do l=1, grid%num_t_node(min_dist_node)
                itr=grid%tr_node(l,min_dist_node)
                call t_coeff(grid,itr,x,p)
                ! if( minval(p) >= 0D0 .and. maxval(p) <= 1D0 ) then
                if (minval(p) >= zero) then
                   find=1
                   exit
                endif
             enddo

             if(find/=1) then
                itr2=itr
                call search_tr_with_guess(grid,x,init,itr2,p,count)
             endif

             if(itr2<0) then
                grid%guess_table(k(1),k(2))=itr
             else
                grid%guess_table(k(1),k(2))=itr2
             endif


             !search all triangle
             if(itr<0) then
                do l=1, grid%ntriangle
                   !get mapping value
                   call t_coeff(grid,l,x,p)
                   !check inside
                   !if( minval(p) >= 0. .and. maxval(p) <=1 ) then
                   if (minval(p) >= zero) then
                      itr=l
                   endif
                enddo
             endif

             grid%guess_table(k(1),k(2))=itr

             !          if(itr>0) init=itr ! set new initial triangle value to old finding one
          enddo
       enddo

    else


       !   -----------------------------------
       !   The guess_table only need a close enough
       !   starting triangle. This should affect
       !   efficiency but not correctness
       !   --------------------------------------

       !   -----------------------------------
       !   default value for starting triangle
       !   -----------------------------------
       ilo = lbound(grid%guess_table,1)
       ihi = ubound(grid%guess_table,1)
       jlo = lbound(grid%guess_table,2)
       jhi = ubound(grid%guess_table,2)

       itr = itr_init
       grid%guess_table(:,:) = itr




       !   --------------------------------
       !   simple but sloppy initialization
       !   --------------------------------

       !   -------------------------------------
       !   sort the triangles by the patch size
       !   handle the smaller patches last
       !   -------------------------------------
       allocate( parea(grid%ntriangle),iperm(grid%ntriangle),stat=ierr )
       call assert(ierr.eq.0,'allocate(parea),ntriangle ',grid%ntriangle)

       do itr=1,grid%ntriangle
          xy1(1:2)= grid%x(1:2,grid%nd(1,itr))
          xy2(1:2)= grid%x(1:2,grid%nd(2,itr))
          xy3(1:2)= grid%x(1:2,grid%nd(3,itr))

          !       -----------------------------------------
          !       determine the bottom left and upper right
          !       corners of enclosing rectangle
          !       -----------------------------------------
          xmin = min( xy1(1), min(xy2(1), xy3(1)) )
          xmax = max( xy1(1), max(xy2(1), xy3(1)) )
          ymin = min( xy1(2), min(xy2(2), xy3(2)) )
          ymax = max( xy1(2), max(xy2(2), xy3(2)) )

          parea(itr) = (ymax-ymin)*(xmax-xmin)
       enddo

       call dshelldec(grid%ntriangle,parea,iperm)
       !   --------------------------------------
       !   check parea(iperm(:)) decreasing order
       !   --------------------------------------
       do itr=1,grid%ntriangle-1
          isok = parea(iperm(itr)).ge.parea(iperm(itr+1))
          call assert( isok,'parea not in sorted order',itr)
       enddo


       do it=1,grid%ntriangle
          itr = iperm(it)

          xy1(1:2)= grid%x(1:2,grid%nd(1,itr))
          xy2(1:2)= grid%x(1:2,grid%nd(2,itr))
          xy3(1:2)= grid%x(1:2,grid%nd(3,itr))

          !       -----------------------------------------
          !       determine the bottom left and upper right
          !       corners of enclosing rectangle
          !       -----------------------------------------
          xmin = min( xy1(1), min(xy2(1), xy3(1)) )
          xmax = max( xy1(1), max(xy2(1), xy3(1)) )
          ymin = min( xy1(2), min(xy2(2), xy3(2)) )
          ymax = max( xy1(2), max(xy2(2), xy3(2)) )

          xymin(1) = xmin
          xymin(2) = ymin
          ijmin(1:2)= (xymin(1:2)-grid%guess_min(1:2))*grid%inv_guess_d +1

          xymax(1) = xmax
          xymax(2) = ymax
          ijmax(1:2)= (xymax(1:2)-grid%guess_min(1:2))*grid%inv_guess_d +1


          i1 = max(ilo,min(ihi, ijmin(1)))
          i2 = max(ilo,min(ihi, ijmax(1)))
          j1 = max(jlo,min(jhi, ijmin(2)))
          j2 = max(jlo,min(jhi, ijmax(2)))
          grid%guess_table( i1:i2, j1:j2 ) = itr

       enddo

       deallocate( parea, iperm, stat=ierr)
       call assert(ierr.eq.0,'deallocate(parea)',ierr)


       !   -----------------------------------------
       !   refine triangle assignment  in guess_table
       !   -----------------------------------------

       do itr=1,grid%ntriangle
          xy1(1:2)= grid%x(1:2,grid%nd(1,itr))
          xy2(1:2)= grid%x(1:2,grid%nd(2,itr))
          xy3(1:2)= grid%x(1:2,grid%nd(3,itr))


          xy(1:2) = xy1(1:2)
          ij(1:2) = (xy(1:2) - grid%guess_min(1:2))*grid%inv_guess_d + 1
          i1 = max(ilo,min(ihi, ij(1)) )
          j1 = max(jlo,min(jhi, ij(2)) )
          grid%guess_table( i1,j1 ) = itr

          xy(1:2) = xy2(1:2)
          ij(1:2) = (xy(1:2) - grid%guess_min(1:2))*grid%inv_guess_d + 1
          i1 = max(ilo,min(ihi, ij(1)) )
          j1 = max(jlo,min(jhi, ij(2)) )
          grid%guess_table( i1,j1 ) = itr

          xy(1:2) = xy3(1:2)
          ij(1:2) = (xy(1:2) - grid%guess_min(1:2))*grid%inv_guess_d + 1
          i1 = max(ilo,min(ihi, ij(1)) )
          j1 = max(jlo,min(jhi, ij(2)) )
          grid%guess_table( i1,j1 ) = itr

          !       --------------------
          !       centroid of triangle
          !       --------------------
          xy(1:2) = (xy1(1:2) + xy2(1:2) + xy3(1:2) )/3.0d0
          ij(1:2) = (xy(1:2) - grid%guess_min(1:2))*grid%inv_guess_d + 1
          i1 = max(ilo,min(ihi, ij(1)) )
          j1 = max(jlo,min(jhi, ij(2)) )
          grid%guess_table( i1,j1 ) = itr

       enddo


       !   ------------------------------------------
       !   use triangle determined by init_guess_list
       !   ------------------------------------------
       do j=1,ubound(grid%guess_count,2)
          do i=1,ubound(grid%guess_count,1)
             if (grid%guess_count(i,j) .ge. 1) then
                itr = grid%guess_list( grid%guess_xtable(i,j) )
                grid%guess_table(i,j) = itr
             endif
          enddo
       enddo





    endif

    !debug-begin
    !   print*,'grid%ntriangle ', grid%ntriangle
    !   print*,'end of init_guess_table'
    !debug-end

    return
  end subroutine init_guess_table



 !! search triangle with initial guess
  subroutine search_tr_with_guess(grid,x,init,itr,p,count)
    implicit none
    type(grid_type),target :: grid
    real (kind=8), intent(in) :: x(2)
    integer, intent(in) :: init
    integer, intent(out) :: itr
    real (kind=8), intent(out) :: p(3)
    integer :: find,min_index(1),next_itr,current_itr
    !debug
    integer :: tmp,i,count,mloc
    real (kind=8) :: mval
    real (kind=8), parameter :: zero = 0.0d0
    real (kind=8), parameter :: eps = 1.0d-10


    find=0
    count=0

    itr=init
    if(init<=0 .or. init>grid%ntriangle) then
       print *, 'invalid guess. init=1 is used instead',init
       itr=1
    endif


    do while(find/=1 .and. count<grid%nnode)
       count=count+1
       call t_coeff(grid,itr,x,p)
       ! if( minval(p) >= 0D0 .and. maxval(p) <=1D0 ) then
       if (minval(p) >= zero) then
          find=1
       else
          min_index=minloc(p)
          next_itr=grid%adj(min_index(1),itr)
          if(next_itr>0) then
             itr=next_itr
          else
             itr=-1
             return
          endif
       endif
    enddo
    if(count>=grid%nnode) then
       call t_coeff(grid,itr,x,p)
       ! if( minval(p) > -1D-10 .and. maxval(p) <1+1D-10) then
       if (minval(p) > -eps) then
          find=1
          mval=-1D0
          do i=1, 3
             if(mval < p(i)) then
                mval=p(i)
                mloc=i
             endif
          enddo
          p(:)=0D0
          p(mloc)=1D0
       else
          itr=-1
          print *, 'search error : p=',p
!       print *, 'Too large loop number in search tr'
       endif
    endif

  end subroutine search_tr_with_guess


end module grid_class

include 'init_guess_list.F90'
include 'check_guess_table.F90'
include 'dshelldec.F90'
include 'checkoverlap.F90'
include 'search_vtr.F90'


!! Return MPI processor information
subroutine get_pe_info(mype,totalpe)
  use sml_module
  implicit none
  integer :: mype, totalpe

  mype=sml_mype
  totalpe=sml_totalpe

end subroutine get_pe_info

!! return number toroidal cross section
subroutine get_nphi(nphi_each)
  use sml_module
  implicit none
  integer :: nphi_each

  nphi_each=sml_plane_per_pe
end subroutine get_nphi

subroutine get_node_vol(grid)
  use grid_class
  use sml_module
  use eq_module, only: eq_axis_r
  implicit none
  type(grid_type) :: grid
!  real (kind=8) , allocatable :: tr_area(:)  !,tr_vol(:)
  real (kind=8) :: vol_tmp,dx1(2),dx2(2),com_r,area,area_tmp
  integer :: i,j,tr,nd(3)
  integer :: max_v, max_t, nnode

  do i=1, grid%ntriangle
     ! find triangle area - 1/2 * vec1 x vec2
     nd(:)=grid%nd(:,i)
     dx1(:)=grid%x(:,nd(1))- grid%x(:,nd(3))
     dx2(:)=grid%x(:,nd(2))- grid%x(:,nd(3))
     area= 0.5D0 * abs( dx1(1)*dx2(2) - dx1(2)*dx2(1))
     ! find triangle center of mass - Radius
     if (sml_cylindrical) then
       ! Cylindrical limit with periodic boundary conditions
       com_r = eq_axis_r
     else
       ! Toroidal geometry
       com_r = 1D0/3D0 * ( grid%x(1,nd(1)) + grid%x(1,nd(2)) +grid% x(1,nd(3)) )
     endif
     !get volume - area * 2pi * Radius / nphi
     grid%tr_vol(i)= area * sml_2pi * com_r / real(sml_nphi_2pi)
     grid%tr_area(i)=area
  enddo
  do i=1, grid%nnode
!     vol_tmp =0D0
     area_tmp=0D0
     do j=1, grid%num_t_node(i)
!        vol_tmp =vol_tmp +tr_vol( grid%tr_node(j,i))/3D0
        area_tmp=area_tmp+grid%tr_area(grid%tr_node(j,i))/3D0
     enddo
     grid%node_area(i)=area_tmp
     grid%node_vol(i)=area_tmp*sml_2pi*grid%x(1,i)/real(sml_nphi_2pi)
     grid%inv_node_vol(i)=1D0/grid%node_vol(i)
  enddo

  !debug
  if(sml_mype==0) print *, 'sum(tr_vol)',sum(grid%tr_vol)
  if(sml_mype==0) print *, 'sum(node_area)',sum(grid%node_area)
  if(sml_mype==0) print *, 'sum(node_vol)',sum(grid%node_vol)



  ! --------------------------------------------------------------------------
  ! compute weights for computing 2D gradient at vertices
  !
  ! gradient computed as  weighted sum of potential values from neighbor nodes
  !   potv(1:k) = pot( v_node(1:k,i) )
  !  Ex = sum( gradx( 1:k, i) * potv(1:k) )
  !  Ey = sum( grady( 1:k, i) * potv(1:k) )
  !  where k = num_v_node(i)
  ! -------------------------------------------------------------------------
  max_t = size( grid%tr_node,dim=1)
  max_v = 2 + max_t
  nnode = grid%nnode
  allocate( grid%gradx(max_v,nnode),  &
            grid%grady(max_v,nnode)  )

  call gen_gradxy(grid%nnode,grid%ntriangle, max_t,max_v, &
        grid%nd, grid%mapping, grid%tr_area, &
        grid%num_t_node, grid%tr_node, &
        grid%num_v_node, grid%v_node, &
        grid%gradx, grid%grady,   sml_mype.eq.0 )


end subroutine get_node_vol


subroutine search_tr2( grid, xy, itr, p )
  use grid_class
  implicit none
  type(grid_type) :: grid
  real(kind=8) :: xy(2)
  integer :: itr
  real(kind=8) :: p(3)

  real(kind=8), parameter :: zero = 0.0d0
  real(kind=8), parameter :: eps = 10.0d0*epsilon(zero)
  integer :: ij(2), istart,iend, k, itrig
  integer :: i,j,  ilo,ihi,  jlo,jhi
  real(kind=8) ::  dx(2), pmin, pmax, dp
  logical :: is_found


  ilo = lbound( grid%guess_table, 1 )
  ihi = ubound( grid%guess_table, 1 )

  jlo = lbound( grid%guess_table, 2 )
  jhi = ubound( grid%guess_table, 2 )

  ij = (xy - grid%guess_min)*grid%inv_guess_d + 1
  i = max(ilo, min(ihi, ij(1)) )
  j = max(jlo, min(jhi, ij(2)) )


  istart = grid%guess_xtable(i,j)
  iend = istart + grid%guess_count(i,j) - 1


  itr = -1
  do k=istart,iend
     itrig = grid%guess_list(k)
     ! call t_coeff( grid, itrig, xy, p )

    dx(1:2) = xy(1:2) - grid%mapping(1:2,3,itrig)
    p(1:2)= grid%mapping(1:2,1,itrig)*dx(1) +                          &
            grid%mapping(1:2,2,itrig)*dx(2)
    p(3)=1.0d0 - p(1) - p(2)

     if (minval(p) .ge. -eps) then
        itr = itrig
        exit
     endif
  enddo

  return
end subroutine search_tr2

!#ifdef OLD_CONVERT_GRID
subroutine convert_001d_2_grid(grid,v1d,v2d)
  use eq_module
  use grid_class
  implicit none
  type(grid_type), intent(in) :: grid
  real (8), intent(in)  :: v1d(grid%npsi00)
  real (8), intent(out) :: v2d(grid%nnode)
  !
  integer :: i, ip
  real (8) :: pn, wp

  do i=1, grid%nnode
     pn=(grid%psi(i)-grid%psi00min)/grid%dpsi00
     ip=floor(pn)+1
     if(0 < ip .and. ip < grid%npsi00 .and. is_rgn12(grid%x(1,i),grid%x(2,i),grid%psi(i)) ) then
        wp=1D0 - ( pn - real(ip-1,8) )
     elseif( ip <= 0 ) then
        ip=1
        wp=1D0
     else
        ip=grid%npsi00-1
        wp=0D0
     endif

     v2d(i)=v1d(ip)*wp  + v1d(ip+1)*(1D0-wp)
  end do

end subroutine convert_001d_2_grid

#ifdef OLD_CONVERT_GRID
! this routine may cause NAN when v1d is too fine -- does not matter if we use convert_001d_2_grid
subroutine convert_grid_2_001d(grid,v2d,v1d)
  use eq_module
  use grid_class
  implicit none
  type(grid_type), intent(in) :: grid
  real (8), intent(in) :: v2d(grid%nnode)
  real (8), intent(out)  :: v1d(grid%npsi00)

  !
  integer :: i, ip
  real (8) :: pn, wp
  real (8) :: dum(grid%nnode)

  v1d=0D0
  dum=0D0

  do i=1, grid%nnode
     pn=(grid%psi(i)-grid%psi00min)/grid%dpsi00
     ip=floor(pn)+1
     if(0 < ip .and. ip < grid%npsi00 .and. is_rgn12(grid%x(1,i),grid%x(2,i),grid%psi(i) ) then
        wp=1D0 - ( pn - real(ip-1,8) )
     elseif( ip <= 0 ) then
        ip=1
        wp=1D0
     else
        ip=grid%npsi00-1
        wp=0D0
     endif


     v1d(ip  )=v1d(ip  ) + v2d(i)*grid%node_vol(i)* wp
     v1d(ip+1)=v1d(ip+1) + v2d(i)*grid%node_vol(i)*(1D0-wp)
     dum(ip  )=dum(ip  ) + grid%node_vol(i)* wp
     dum(ip+1)=dum(ip+1) + grid%node_vol(i)*(1D0-wp)

  end do
  v1d=v1d/dum


end subroutine convert_grid_2_001d
#else
subroutine convert_grid_init(grid,nsub)
  use eq_module
  use grid_class
  use sml_module
  implicit none
  type(grid_type), intent(inout) :: grid
  integer,intent(in) :: nsub
  !
  integer :: itr, nd(3)
  real (8) :: dx1(2),dx2(2),area
  integer :: j,k,ip
  real (8) :: c1,c2,c3,xc(2),volume,psi,pn,wp
  integer :: itmp
  real (8), external :: psi_interpol

  ! matrix initialization
#ifdef ITER_GRID
  call new_mat(grid%cnv_2d_00,grid%nnode,60) ! 30 may be too large --> make it dynamic later
#else
  call new_mat(grid%cnv_2d_00,grid%nnode,30) ! 30 may be too large --> make it dynamic later
#endif
  call set_non_square_mat(grid%cnv_2d_00,grid%npsi00)


  ! assign matrix elements
  do itr=1, grid%ntriangle ! for all triangle

     nd=grid%nd(:,itr)
     dx1= ( grid%x(:,nd(1)) - grid%x(:,nd(3)) )/real(nsub,8)
     dx2= ( grid%x(:,nd(2)) - grid%x(:,nd(3)) )/real(nsub,8)
     area= 0.5D0*abs( dx1(1)*dx2(2) - dx2(1)*dx1(2) )

     do j=1, nsub       ! for all subtriangle
        do k=1, 2*j-1 ! for all subtriangle
           itmp=  (k-1)/2   ! make integer if not integer
           c1=(j-1)- itmp + real(mod(k,2)*2-1,8)/3D0
           itmp= k/2
           c2=itmp + real(mod(k,2)*2-1,8)/3D0
           c3=real(nsub,8) - c1 - c2

           xc=grid%x(:,nd(3)) + c1*dx1 + c2*dx2 ! center of subtriangle
           volume=area*xc(1)                             ! volume of subtriangle (2pi is missing)

           psi=psi_interpol(xc(1),xc(2),0,0)  ! find psi of center


           pn=(psi-grid%psi00min)/grid%dpsi00  ! parameter for linear interpolation
           ip=floor(pn)+1

           if(0 < ip .and. ip < grid%npsi00 .and. is_rgn12(xc(1),xc(2),psi) &
                  .and. (.not. sml_00_xz_up .or. xc(2) > eq_x_z) ) then  ! only when psi is inside and exclude private region
              wp=1D0 - ( pn - real(ip-1,8) )

              !add matrix comp
              ! factor of nsub is missing but ok (normalization later)
              call set_value(grid%cnv_2d_00,nd(1),ip  ,wp*c1*volume,1)
              call set_value(grid%cnv_2d_00,nd(2),ip  ,wp*c2*volume,1)
              call set_value(grid%cnv_2d_00,nd(3),ip  ,wp*c3*volume,1)

              call set_value(grid%cnv_2d_00,nd(1),ip+1 ,(1D0-wp)*c1*volume,1)
              call set_value(grid%cnv_2d_00,nd(2),ip+1 ,(1D0-wp)*c2*volume,1)
              call set_value(grid%cnv_2d_00,nd(3),ip+1 ,(1D0-wp)*c3*volume,1)

           endif
        enddo
     enddo
  enddo

  !normalize -- average opertion
  call normalize_const

contains
  subroutine normalize_const
    implicit none
    real (8), allocatable :: tmp(:)


    allocate(grid%cnv_norm_2d(grid%nnode),grid%cnv_norm_1d00(grid%npsi00))
    allocate(tmp(grid%npsi00))


    ! obtain cnv_norm_1d00
    grid%cnv_norm_2d=1D0  ! as a tmp variable here
    call mat_transpose_mult(grid%cnv_2d_00,grid%cnv_norm_2d,grid%cnv_norm_1d00)

    ! obtain cnv_norm_2d
    tmp=1D0
    call mat_mult(grid%cnv_2d_00,tmp,grid%cnv_norm_2d)
    grid%cnv_norm_2d = grid%cnv_norm_2d + 1D-99 ! prevent NaN

    deallocate(tmp)

  end subroutine normalize_const

!!$  subroutine normalize_row(mat)
!!$    implicit none
!!$    type(mat_type) :: mat
!!$    real (8), allocatable :: tmp1(:), tmp2(:)
!!$    integer :: i,j,k
!!$
!!$    allocate(tmp1(mat%n),tmp2(mat%m))
!!$    tmp1=1D0
!!$    call mat_transpose_mult(mat,tmp1,tmp2)
!!$
!!$    do i=1, mat%n
!!$       do j=1, mat%nelement(i)
!!$          k = mat%eindex(j,i)
!!$          mat%value(j,i)=mat%value(j,i)/tmp2(k)
!!$       enddo
!!$    enddo
!!$    deallocate(tmp1,tmp2)
!!$  end subroutine normalize_row
end subroutine convert_grid_init



!subroutine convert_001d_2_grid(grid,v1d,v2d)
!  use eq_module
!  use grid_class
!  implicit none
!  type(grid_type), intent(in) :: grid
!  real (8), intent(in)  :: v1d(grid%npsi00)
!  real (8), intent(out) :: v2d(grid%nnode)
!  !
!  v2d=0D0
!  call mat_mult(grid%cnv_2d_00,v1d,v2d)
!  v2d=v2d/grid%cnv_norm_2d
!
!end subroutine convert_001d_2_grid


subroutine convert_grid_2_001d(grid,v2d,v1d)
  use grid_class
  implicit none
  type(grid_type), intent(in) :: grid
  real (8), intent(in) :: v2d(grid%nnode)
  real (8), intent(out)  :: v1d(grid%npsi00)


  call mat_transpose_mult(grid%cnv_2d_00,v2d,v1d)
  v1d=v1d/grid%cnv_norm_1d00

end subroutine convert_grid_2_001d
#endif

#include "gen_v_node.F90"
#include "calc_gradxy.F90"
#include "gen_gradxy.F90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!rh 2nd order FD (R,Z) gradient operator
!rh Things that can be improved:
!rh ---> 4th order gradient
!rh ---> generalize for Laplacian
!rh ---> Use proper metric for distance calculations (not done yet due to small step size)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Initialization of matrices for gradient-operator
subroutine init_gradient_mat(grid)
  use grid_class
  use mat_class
  use eq_module, only: eq_x_r, eq_x_z, eq_axis_r, eq_axis_z, is_rgn12
  use sml_module, only: sml_gradpsi_limit, sml_inpsi, sml_outpsi
  implicit none
  type (grid_type), intent(inout) :: grid
  real (kind=8) :: d1, d2, tmp, x(2), dum(2), area, area_sum
  real (kind=8) :: det, dx1(2), dx2(2), r, z, psi
  integer :: i, j, nnode, tnum, nd(3), self, a, b
  real (kind=8), parameter :: eps=1.D-8
  real (kind=8), external :: psi_interpol

  nnode=grid%nnode

  ! Set up mapping for (Psi,theta) derivative
  ! matrix initialization
  j=size(grid%tr_node,1) ! maximum number of triangles per node
  ! The number of entries in each line of the matrix is 2*j+1;
  !"+1" for the diagonal element, 2*j for the off-diagonal elements
  call new_mat(grid%gradientx,nnode,j*2+1)
  call new_mat(grid%gradienty,nnode,j*2+1)

  allocate(grid%unit_vecs(2,2,nnode))
  allocate(grid%basis(nnode))
  grid%basis(:)=0
  grid%unit_vecs=0.D0

  do i=1,nnode
    x=grid%x(:,i)
    r=x(1)
    z=x(2)
    psi=grid%psi(i)
    d1=sqrt((grid%x(1,i)-eq_axis_r)**2+(grid%x(2,i)-eq_axis_z)**2)
    d2=sqrt((grid%x(1,i)-eq_x_r)**2+(grid%x(2,i)-eq_x_z)**2)
    dum(1)=psi_interpol(x(1),x(2),1,0)
    dum(2)=psi_interpol(x(1),x(2),0,1)
    tmp=sqrt(sum(dum**2))
    if (d1 .gt. eps .and. d2 .gt. eps .and. tmp .gt. sml_gradpsi_limit .and. &
        sml_inpsi .le. psi .and. psi .le. sml_outpsi) then
      ! Condition is distance from X point and axis and lower limit for |grad(Psi)|
      ! Unit vector in Psi direction (in (R,Z) space)
      grid%unit_vecs(1,1,i)=dum(1)
      grid%unit_vecs(1,2,i)=dum(2)
      grid%unit_vecs(1,:,i)=grid%unit_vecs(1,:,i)/tmp
      ! Unit vector in poloidal direction (in (R,Z) space) (counter-clockwise!)
      ! This should always be counter-clockwise
      grid%unit_vecs(2,1,i)=-grid%unit_vecs(1,2,i)
      grid%unit_vecs(2,2,i)=grid%unit_vecs(1,1,i)

      grid%basis(i)=0

      !if (is_rgn12(r,z,psi) .and. sml_inpsi .le. grid%psi(i) .and. grid%psi(i) .le. sml_outpsi) then
        call set_gradient_mat_flux_aligned(grid,i)
      !endif
    elseif (d1 .le. eps .or. (psi .lt. sml_inpsi .or. psi .gt. sml_outpsi) ) then
      ! grad(X) can be only in grad(Psi) direction, which cannot be
      ! described in cylindrical coordinates.
      ! Use R-Z derivatives and make the derivative symmetric
      grid%unit_vecs(1,1,i)=1.D0
      grid%unit_vecs(1,2,i)=0.D0
      grid%unit_vecs(2,1,i)=0.D0
      grid%unit_vecs(2,2,i)=1.D0
      grid%basis(i)=1
      call set_value(grid%gradientx,i,i,0D0,1)
      call set_value(grid%gradienty,i,i,0D0,1)
      !call set_gradient_mat_triangle(grid,i)
    else
      ! In all other cases, use (R,Z) unit vectors
      grid%unit_vecs(1,1,i)=1.D0
      grid%unit_vecs(1,2,i)=0.D0
      grid%unit_vecs(2,1,i)=0.D0
      grid%unit_vecs(2,2,i)=1.D0
      grid%basis(i)=1

      !if (is_rgn12(r,z,psi) .and. sml_inpsi .le. grid%psi(i) .and. grid%psi(i) .le. sml_outpsi) then
        call set_gradient_mat_triangle(grid,i)
      !endif
    endif
  enddo

end subroutine init_gradient_mat

subroutine set_gradient_mat_flux_aligned(grid,ip)
  use grid_class
  use mat_class
  use eq_module, only: eq_x_psi
  implicit none
  type (grid_type), intent(inout) :: grid
  integer, intent(in) :: ip
  real (kind=8) :: x(2), psi, dum(3), tmp
  real (kind=8) :: p1(3), p2(3), p3(3), p4(3), dl1, dl2, dl3, dl4, dist_psi, dist_theta
  integer :: i, j, nnode, tnum, nd(3), self, a, b, itr1, itr2, itr3, itr4, dir

  nnode=grid%nnode
  tnum=grid%num_t_node(ip)
  dist_psi=0D0
  dist_theta=0D0
  x=grid%x(:,ip)
  psi=grid%psi(ip)

  ! Step 1: Determine characteristic distance in grad(psi) and grad(theta) direction
  call get_char_length(grid,ip,dum)
  ! Use 80% of the characteristic distance for the generation of the finite difference points
  ! because, in the calculation of dist_psi and dist_theta, distance to some points
  ! may be zero (e.g. dist_psi for 2 points on the same flux-surface).
  dist_psi = dum(2)*0.8D0  !! Use dum(2) for characteristic psi-difference
  !dist_psi = dum(1)*0.8D0 !! Use dum(1) for characteristic length in grad(psi)-direction
  dist_theta = dum(3)*0.8D0

  ! Step 2: Follow grad(Psi) to find the points and weights for finite difference
  !         Psi-derivative
  dir=1  !! positive grad(psi)-direction <--> outwards
  call get_node_psi(grid,x,psi,dist_psi,dir,dl1,itr1,p1)
  dir=-1  !! negative grad(psi)-direction <--> inwards
  call get_node_psi(grid,x,psi,dist_psi,dir,dl2,itr2,p2)

  ! Step 3: Follow grad(theta) to find the points and weights for finite difference
  !         theta-derivative
  dir=1  !! positive grad(psi)-direction <--> outwards
  call get_node_theta(grid,x,dist_theta,dir,dl3,itr3,p3)
  dir=-1  !! negative grad(psi)-direction <--> inwards
  call get_node_theta(grid,x,dist_theta,dir,dl4,itr4,p4)

  if (itr1 < 1 .or. itr2 < 1 .or. itr3 < 1 .or. itr4 < 1) then
    ! If any of the above operations did not return a valid point
    ! activate fallback --> use triangle gradient operation
    grid%unit_vecs(1,1,ip)=1.D0
    grid%unit_vecs(1,2,ip)=0.D0
    grid%unit_vecs(2,1,ip)=0.D0
    grid%unit_vecs(2,2,ip)=1.D0
    grid%basis(ip)=1
    call set_gradient_mat_triangle(grid,ip)
    return
  endif

  ! If we have found 4 valid points, set up the finite difference matrix
  ! using the flux-surface aligned

  ! I) Psi-derivative
  ! a) positive grad(Psi) direction
  nd=grid%nd(:,itr1)
  do i=1,3
    tmp=p1(i)*dl2/(dl1*(dl1+dl2))
    call set_value(grid%gradientx,ip,nd(i),tmp,1)
  enddo
  ! b) negative grad(Psi) direction
  nd=grid%nd(:,itr2)
  do i=1,3
    tmp=-p2(i)*dl1/(dl2*(dl1+dl2))
    call set_value(grid%gradientx,ip,nd(i),tmp,1)
  enddo
  ! c) mid-point
  tmp=(dl1**2-dl2**2)/(dl1*dl2*(dl1+dl2))
  call set_value(grid%gradientx,ip,ip,tmp,1)


  ! II) theta-derivative
  ! a) positive grad(Psi) direction
  nd=grid%nd(:,itr3)
  do i=1,3
    tmp=p3(i)*dl4/(dl3*(dl3+dl4))
    call set_value(grid%gradienty,ip,nd(i),tmp,1)
  enddo
  ! b) negative grad(Psi) direction
  nd=grid%nd(:,itr4)
  do i=1,3
    tmp=-p4(i)*dl3/(dl4*(dl3+dl4))
    call set_value(grid%gradienty,ip,nd(i),tmp,1)
  enddo
  ! c) mid-point
  tmp=(dl3**2-dl4**2)/(dl3*dl4*(dl3+dl4))
  call set_value(grid%gradienty,ip,ip,tmp,1)

  ! All finished!

end subroutine set_gradient_mat_flux_aligned

subroutine set_gradient_mat_triangle(grid,ip)
  use grid_class
  use mat_class
  implicit none
  type (grid_type), intent(inout) :: grid
  integer, intent(in) :: ip
  real (kind=8) :: tmp, dum(2), area, area_sum
  real (kind=8) :: det, dx1(2), dx2(2)
  integer :: i, j, nnode, tnum, nd(3), self, a, b

  nnode=grid%nnode
  tnum=grid%num_t_node(ip)
  area_sum=0D0

  do j=1, tnum
    area_sum = area_sum + grid%tr_area(grid%tr_node(j,ip))
  enddo

  do j=1,tnum

    area=grid%tr_area(grid%tr_node(j,ip))
    nd=grid%nd(:,grid%tr_node(j,ip))

    if (ip==nd(1)) then
      self=ip
      a=nd(2)
      b=nd(3)
    elseif(ip==nd(2)) then
      self=ip
      a=nd(1)
      b=nd(3)
    else
      self=ip
      a=nd(1)
      b=nd(2)
    endif

    dum=grid%x(:,a) - grid%x(:,self)
    dx1(1)=sum(dum*grid%unit_vecs(1,:,self))
    dx1(2)=sum(dum*grid%unit_vecs(2,:,self))

    dum=grid%x(:,b) - grid%x(:,self)
    dx2(1)=sum(dum*grid%unit_vecs(1,:,self))
    dx2(2)=sum(dum*grid%unit_vecs(2,:,self))

    det=1./( dx1(1)*dx2(2) - dx2(1)*dx1(2))

    tmp=(dx1(2)-dx2(2))*det*area/area_sum
    call set_value(grid%gradientx,self,self,tmp,1)
    tmp=dx2(2)*det*area/area_sum
    call set_value(grid%gradientx,self,a,tmp,1)
    tmp=-dx1(2)*det*area/area_sum
    call set_value(grid%gradientx,self,b,tmp,1)
    tmp=(dx2(1)-dx1(1))*det*area/area_sum
    call set_value(grid%gradienty,self,self,tmp,1)
    tmp=-dx2(1)*det*area/area_sum
    call set_value(grid%gradienty,self,a,tmp,1)
    tmp=dx1(1)*det*area/area_sum
    call set_value(grid%gradienty,self,b,tmp,1)

  enddo

end subroutine set_gradient_mat_triangle

! This routine uses get_next_point_gl to find a valid point in grad(psi)
! direction using a psi-difference. It reduces the input distance until it finds a valid point
! or the target distance is less than 10% of the input distance
subroutine get_node_psi(grid,x,psi_in,dpsi_in,direction,dist_out,itr,p)
  use grid_class
  implicit none
  type (grid_type), intent(in) :: grid
  real (kind=8),intent(in) :: x(2), dpsi_in, psi_in
  integer, intent(in) :: direction
  real (kind=8), intent(inout) :: p(3), dist_out
  integer, intent(inout) :: itr
  real (kind=8) :: x2(2), dpsi, psi
  real (kind=8), external :: psi_interpol

  dpsi=dpsi_in
  itr=-1

  do while (itr<1 .and. dpsi .ge. 0.1D0*dpsi_in)
    call get_next_point_gl(x,psi_in,dpsi,direction,x2,dist_out)
    psi=psi_interpol(x2(1),x2(2),0,0)
    !rh call search_tr2(grid,x2,psi,itr,p)
    call search_tr2(grid,x2,itr,p)
    if (itr<1) then
      dpsi=0.4D0*dpsi
    else
      exit
    endif
  enddo
  if (itr<1) then
    ! This is not an error, although it should not happen.
    ! This is likely a node close to the wall.
    dist_out=0.D0
!    if (sml_mype==0) then
!      print *,'init_grid_deriv: No base point found in Psi-direction'
!      print '(a11,i3)', 'Direction: ', direction
!      print '(a13,2f16.6)', 'Coordinates: ',x(1),x(2)
!    endif
  endif

end subroutine get_node_psi

! This routine uses get_next_point_perp to find a valid point in grad(psi)
! direction using a distance in m. It reduces the input distance until it finds a valid point
! or the target distance is less than 10% of the input distance
subroutine get_node_perp(grid,x,dlperp_in,direction,dist_out,itr,p)
  use grid_class
  implicit none
  type (grid_type), intent(in) :: grid
  real (kind=8),intent(in) :: x(2), dlperp_in
  integer, intent(in) :: direction
  real (kind=8), intent(inout) :: p(3), dist_out
  integer, intent(inout) :: itr
  real (kind=8) :: x2(2), dlperp, psi
  real (kind=8), external :: psi_interpol

  dlperp=dlperp_in
  itr=-1

  do while (itr<1 .and. dlperp .ge. 0.1D0*dlperp_in)
    call get_next_point_perp(x,dlperp,direction,x2,dist_out)
    psi=psi_interpol(x2(1),x2(2),0,0)
    !rh call search_tr2(grid,x2,psi,itr,p)
    call search_tr2(grid,x2,itr,p)
    if (itr<1) then
      dlperp=0.4D0*dlperp
    else
      exit
    endif
  enddo
  if (itr<1) then
    ! This is not an error, although it should not happen.
    ! This is likely a node close to the wall.
    dist_out=0.D0
!    if (sml_mype==0) then
!      print *,'init_grid_deriv: No base point found in Psi-direction'
!      print '(a11,i3)', 'Direction: ', direction
!      print '(a13,2f16.6)', 'Coordinates: ',x(1),x(2)
!    endif
  endif

end subroutine get_node_perp

! This routine uses get_next_point_tang to find a valid point in grad(theta)
! direction using a distance in m. It reduces the input distance until it finds a valid point
! or the target distance is less than 10% of the input distance
subroutine get_node_theta(grid,x,dltheta_in,direction,dist_out,itr,p)
  use grid_class
  implicit none
  type (grid_type), intent(in) :: grid
  real (kind=8),intent(in) :: x(2), dltheta_in
  integer, intent(in) :: direction
  real (kind=8), intent(inout) :: p(3), dist_out
  integer, intent(inout) :: itr
  real (kind=8) :: x2(2), dltheta, psi
  real (kind=8), external :: psi_interpol

  dltheta=dltheta_in
  itr=-1

  do while (itr<1 .and. dltheta .ge. 0.1D0*dltheta_in)
    call get_next_point_tang(x,dltheta,direction,x2,dist_out)
    psi=psi_interpol(x2(1),x2(2),0,0)
    !rh call search_tr2(grid,x2,psi,itr,p)
    call search_tr2(grid,x2,itr,p)
    if (itr<1) then
      dltheta=0.4D0*dltheta
    else
      exit
    endif
  enddo
  if (itr<1) then
    ! This is not an error, although it should not happen.
    ! This is likely a node close to the wall.
    dist_out=0.D0
!    if (sml_mype==0) then
!      print *,'init_grid_deriv: No base point found in theta-direction'
!      print '(a11,i3)', 'Direction: ', direction
!      print '(a13,2f16.6)', 'Coordinates: ',x(1),x(2)
!    endif
  endif

end subroutine get_node_theta

! RK2 routine to follow the Psi gradient a given distance dPsi in Psi space
subroutine get_next_point_gl(x_in,psi_in,dpsi,direction,x_out,dist)
  implicit none
  real (kind=8), intent(in) :: x_in(2), psi_in, dpsi
  real (kind=8), intent(inout) :: x_out(2), dist
  integer, intent(in) :: direction
  real (kind=8), dimension(2) :: gradpsi
  real (kind=8) :: x1(2),x2(2), absgradpsi2, dp, psi
  real (kind=8), external :: psi_interpol
  integer :: i
  integer, parameter :: maxsteps=20

  x2=x_in
  ! Plan for maxsteps/2 iterations -->
  dp=dpsi/real(maxsteps/2,8)
  dist=0.D0

  do i=1,maxsteps
    x1=x2
    gradpsi(1)=psi_interpol(x1(1),x1(2),1,0)
    gradpsi(2)=psi_interpol(x1(1),x1(2),0,1)
    ! absgradpsi is actually |grad(Psi)|^2
    absgradpsi2=(gradpsi(1)**2+gradpsi(2)**2)
    ! I need to implement a fallback in case we are on the separatrix:
    if (absgradpsi2 .lt. 1d-13) then
       print *, 'get_next_point_gl: Warning, gradient smaller than threshold!'
    endif
    ! Predictor step: dl_perp=dpsi/|grad(Psi)|
    ! ==> dl_perp uvector_perp = dpsi*grad(Psi)/|grad(Psi)|^2
    x2=x1+direction*(0.5D0*dp)*gradpsi/absgradpsi2
    gradpsi(1)=psi_interpol(x2(1),x2(2),1,0)
    gradpsi(2)=psi_interpol(x2(1),x2(2),0,1)
    absgradpsi2=(gradpsi(1)**2+gradpsi(2)**2)
    ! Corrector step
    x2=x1+direction*dp*gradpsi/absgradpsi2
    psi=psi_interpol(x2(1),x2(2),0,0)
    if (dabs(psi-psi_in) .ge. dpsi) then
      if (dabs(psi-psi_in) .gt. 1.5D0*dpsi) then
        ! We have gone too far, go one step back
        x2=x1
      else
        dist=dist+dsqrt(sum((x2-x1)**2))
      endif
      ! exit the loop
      x_out=x2
      exit
    endif
    dist=dist+dsqrt(sum((x2-x1)**2))
  enddo
end subroutine get_next_point_gl

! RK2 routine to follow the Psi gradient a given distance dlperp
subroutine get_next_point_perp(x_in,dlperp,direction,x_out,dist)
  implicit none
  real (kind=8), intent(in) :: x_in(2), dlperp
  real (kind=8), intent(inout) :: x_out(2), dist
  integer, intent(in) :: direction
  real (kind=8), dimension(2) :: gradpsi
  real (kind=8) :: x1(2),x2(2), absgradpsi, dlp, dlc
  real (kind=8), external :: psi_interpol
  integer :: i
  integer, parameter :: maxsteps=20

  x2=x_in
  ! Plan for maxsteps/2 iterations -->
  dlp=dlperp/real(maxsteps/2,8)
  dist=0.D0
  dlc=0

  do i=1,maxsteps
    x1=x2
    gradpsi(1)=psi_interpol(x1(1),x1(2),1,0)
    gradpsi(2)=psi_interpol(x1(1),x1(2),0,1)
    absgradpsi=dsqrt(gradpsi(1)**2+gradpsi(2)**2)
    ! I need to implement a fallback in case we are on the separatrix:
    if (absgradpsi .lt. 1d-13) then
       print *, 'get_next_point_gl: Warning, gradient smaller than threshold!'
    endif
    ! Predictor step
    x2=x1+direction*(0.5D0*dlp)*gradpsi/absgradpsi
    gradpsi(1)=psi_interpol(x2(1),x2(2),1,0)
    gradpsi(2)=psi_interpol(x2(1),x2(2),0,1)
    absgradpsi=dsqrt(gradpsi(1)**2+gradpsi(2)**2)
    ! Corrector step
    x2=x1+direction*dlp*gradpsi/absgradpsi
    dlc=dlc+dsqrt(sum((x2-x1)**2))
    if (dlc .ge. dlperp) then
      if (dlc .gt. 1.5D0*dlperp) then
        ! We have gone too far, go one step back
        dist=dlc-dsqrt(sum((x2-x1)**2))
        x2=x1
      else
        dist=dlc
      endif
      ! exit the loop
      x_out=x2
      exit
    endif
    dist=dlc
  enddo
end subroutine get_next_point_perp

! RK4 routine to follow the flux surface
subroutine get_next_point_tang(x_in,dltheta,direction,x_out,dist)
  implicit none
  real (kind=8), intent(in) :: x_in(2), dltheta
  real (kind=8), intent(inout) :: x_out(2), dist
  integer, intent(in) :: direction
  real (kind=8), dimension(2) :: polv
  real (kind=8) :: abspolv, dlt, dlc
  real (kind=8) :: xp(2,4), x(2,5)
  real (kind=8), external :: psi_interpol
  integer :: i
  integer, parameter :: maxsteps=50

  x(:,5)=x_in
  ! Plan for maxsteps/2 iterations -->
  dlt=dltheta/real(maxsteps/2,8)
  dist=0.D0
  dlc=0

  do i=1,maxsteps
    x(:,1)=x(:,5)
    ! Direction is counter-clockwise!
    polv(1)=-psi_interpol(x(1,1),x(2,1),0,1)
    polv(2)=psi_interpol(x(1,1),x(2,1),1,0)
    abspolv=dsqrt(polv(1)**2+polv(2)**2)
    xp(:,1)=polv/abspolv
    ! I need to implement a fallback in case we are on the separatrix:
    !if (abspolv .lt. 1d-13) then
    !   print *, 'get_next_point_gl: Warning, gradient smaller than threshold!'
    !endif
    ! Predictor steps
    x(:,2)=x(:,1)+direction*(0.5D0*dlt)*xp(:,1)
    polv(1)=-psi_interpol(x(1,2),x(2,2),0,1)
    polv(2)=psi_interpol(x(1,2),x(2,2),1,0)
    abspolv=dsqrt(polv(1)**2+polv(2)**2)
    xp(:,2)=polv/abspolv

    x(:,3)=x(:,1)+direction*(0.5D0*dlt)*xp(:,2)
    polv(1)=-psi_interpol(x(1,3),x(2,3),0,1)
    polv(2)=psi_interpol(x(1,3),x(2,3),1,0)
    abspolv=dsqrt(polv(1)**2+polv(2)**2)
    xp(:,3)=polv/abspolv

    x(:,4)=x(:,1)+direction*dlt*xp(:,3)
    polv(1)=-psi_interpol(x(1,4),x(2,4),0,1)
    polv(2)=psi_interpol(x(1,4),x(2,4),1,0)
    abspolv=dsqrt(polv(1)**2+polv(2)**2)
    xp(:,4)=polv/abspolv

    ! Corrector step
    x(:,5)=x(:,1)+direction*(dlt/6D0)*(xp(:,1)+2D0*(xp(:,2)+xp(:,3))+xp(:,4))
    dlc=dlc+dsqrt(sum((x(:,5)-x(:,1))**2))
    if (dlc .ge. dltheta) then
      if (dlc .gt. 1.5D0*dltheta) then
        ! We have gone too far, go one step back
        dist=dlc-dsqrt(sum((x(:,5)-x(:,1))**2))
        x(:,5)=x(:,1)
      else
        dist=dlc
      endif
      ! exit the loop
      x_out=x(:,5)
      exit
    endif
    dist=dlc
  enddo
end subroutine get_next_point_tang

subroutine get_char_length(grid,ip,dl)
  use grid_class
  use mat_class
  use eq_module, only: eq_x_psi
  implicit none
  type (grid_type), intent(inout) :: grid
  integer, intent(in) :: ip
  real (kind=8), intent(out) :: dl(3)
  real (kind=8) :: x(2), psi, norm, dum(2), tmp
  real (kind=8) :: dist_psipsi, dist_psi, dist_theta
  integer :: i, j, nnode, tnum, nd(3), self, a, b
  integer :: count1, count2
  !! Assume 100 flux-surfaces in region 1 and 1% error
  !! --> criterion for two grid points being on the same flux-surface:
  real (kind=8) :: eps1

  eps1=eq_x_psi*1D-4 ! maybe even 1D-5?
  count1=0
  count2=0
  nnode=grid%nnode
  tnum=grid%num_t_node(ip)
  dist_psi=0D0
  dist_psipsi=0D0
  dist_theta=0D0
  x=grid%x(:,ip)
  psi=grid%psi(ip)

  ! Step 1: Determine characteristic distance in grad(psi) and grad(theta) direction
  do j=1,tnum

    nd=grid%nd(:,grid%tr_node(j,ip))

    if (ip==nd(1)) then
      self=ip
      a=nd(2)
      b=nd(3)
    elseif(ip==nd(2)) then
      self=ip
      a=nd(1)
      b=nd(3)
    else
      self=ip
      a=nd(1)
      b=nd(2)
    endif

    dum=grid%x(:,a) - grid%x(:,self)
    if (abs(grid%psi(a)-grid%psi(self)) .gt. eps1) then
      ! a and ip are on different surfaces
      dist_psi = dist_psi + abs(sum(dum*grid%unit_vecs(1,:,self)))  !! when using distance in grad(Psi) direction
      dist_psipsi = dist_psipsi + abs(grid%psi(a) - grid%psi(self))  !! when using Psi difference in grad(Psi) direction
      count1 = count1+1
    else
      ! a and ip are on the same surface
      dist_theta = dist_theta + abs(sum(dum*grid%unit_vecs(2,:,self)))
      count2 = count2+1
    endif

    dum=grid%x(:,b) - grid%x(:,self)
    if (abs(grid%psi(b)-grid%psi(self)) .gt. eps1) then
      ! a and ip are on different surfaces
      dist_psi = dist_psi + abs(sum(dum*grid%unit_vecs(1,:,self)))  !! when using distance in grad(Psi) direction
      dist_psipsi = dist_psipsi + abs(grid%psi(b) - grid%psi(self))  !! when using Psi difference in grad(Psi) direction
      count1 = count1+1
    else
      ! a and ip are on the same surface
      dist_theta = dist_theta + abs(sum(dum*grid%unit_vecs(2,:,self)))
      count2 = count2+1
    endif
  enddo

  dist_psi = dist_psi/real(count1,8)
  dist_psipsi = dist_psipsi/real(count1,8)
  if (count2 .gt. 0) then
    dist_theta = dist_theta/real(count2,8)
  else
    dist_theta=dist_psi
  endif

  dl(1)=dist_psi
  dl(2)=dist_psipsi
  dl(3)=dist_theta

end subroutine get_char_length

! Generates matrix for v.grad(X) from psi-theta derivatives
subroutine make_v_dot_grad_mat(grid,mat,v)
  use grid_class
  use mat_class
  implicit none
  type(grid_type), intent(inout) :: grid
  real (kind=8), intent(in) :: v(grid%nnode,3)
  type(mat_type), intent(out) :: mat
  type(mat_type) :: ddr, ddz
  real (kind=8) :: v2(2), val
  integer :: i, j, k, width

  j=size(grid%tr_node,1) ! maximum number of triangles per node

  ! 1) Create matrices for R and Z derivatives
  !call new_mat(ddr,grid%nnode,2*(2*j+1))
  !call new_mat(ddz,grid%nnode,2*(2*j+1))

  ! 2) Calculate elements of ddr and ddz
  !do i=1,grid%nnode
  !  width=grid%gradientx%nelement(i)
  !  do j=1,width
  !    val = grid%unit_vecs(1,1,i)*grid%gradientx%value(j,i)
  !    k   = grid%gradientx%eindex(j,i)
  !    call set_value(ddr,i,k,val,1)
  !    val = grid%unit_vecs(1,2,i)*grid%gradientx%value(j,i)
  !    k   = grid%gradientx%eindex(j,i)
  !    call set_value(ddz,i,k,val,1)
  !  enddo
  !  width=grid%gradienty%nelement(i)
  !  do j=1,width
  !    val = grid%unit_vecs(2,1,i)*grid%gradienty%value(j,i)
  !    k   = grid%gradienty%eindex(j,i)
  !    call set_value(ddr,i,k,val,1)
  !    val = grid%unit_vecs(2,2,i)*grid%gradienty%value(j,i)
  !    k   = grid%gradienty%eindex(j,i)
  !    call set_value(ddz,i,k,val,1)
  !  enddo
  !enddo

  call new_mat(mat,grid%nnode,2*(2*j+1))

  ! Calculate
  !do i=1,grid%nnode
  !  width=ddr%nelement(i)
  !  do j=1,width
  !    val = v(i,1)*ddr%value(j,i)
  !    k   = ddr%eindex(j,i)
  !    call set_value(mat,i,k,val,1)
  !  enddo
  !  width=ddz%nelement(i)
  !  do j=1,width
  !    val = v(i,2)*ddz%value(j,i)
  !    k   = ddz%eindex(j,i)
  !    call set_value(mat,i,k,val,1)
  !  enddo
  !enddo
#ifndef OPTIM_GYRO_AVG_MAT
  do i=1,grid%nnode
    if (grid%basis(i)==0) then
      ! psi and theta components of (b x grad(j0/eB))
      v2(1) = v(i,1)*grid%unit_vecs(1,1,i) + v(i,2)*grid%unit_vecs(1,2,i)
      v2(2) = v(i,1)*grid%unit_vecs(2,1,i) + v(i,2)*grid%unit_vecs(2,2,i)
    else
      v2(1) = v(i,1)
      v2(2) = v(i,2)
    endif
    width=grid%gradientx%nelement(i)
    do j=1,width
      val = v2(1)*grid%gradientx%value(j,i)
      k   = grid%gradientx%eindex(j,i)
      call set_value(mat,i,k,val,1)
    enddo
    width=grid%gradienty%nelement(i)
    do j=1,width
      val = v2(2)*grid%gradienty%value(j,i)
      k   = grid%gradienty%eindex(j,i)
      call set_value(mat,i,k,val,1)
    enddo
  enddo
#else
  do i=1,grid%nnode
    if (grid%basis(i)==0) then
      ! psi and theta components of (b x grad(j0/eB))
      v2(1) = v(i,1)*grid%unit_vecs(1,1,i) + v(i,2)*grid%unit_vecs(1,2,i)
      v2(2) = v(i,1)*grid%unit_vecs(2,1,i) + v(i,2)*grid%unit_vecs(2,2,i)
    else
      v2(1) = v(i,1)
      v2(2) = v(i,2)
    endif
    width=grid%gradientx%nelement(i)
    do j=1,width
      val = v2(1)*grid%gradientx%value(grid%gradientx%a(i)+j)
      k   = grid%gradientx%eindex(grid%gradientx%a(i)+j)
      call set_value(mat,i,k,val,1)
    enddo
    width=grid%gradienty%nelement(i)
    do j=1,width
      val = v2(2)*grid%gradienty%value(grid%gradienty%a(i)+j)
      k   = grid%gradienty%eindex(grid%gradienty%a(i)+j)
      call set_value(mat,i,k,val,1)
    enddo
  enddo
#endif
  !call del_mat(ddr)
  !call del_mat(ddz)

end subroutine make_v_dot_grad_mat


! Computes the gradient of grid based quantities
subroutine grid_deriv(grid,qty,qty_deriv_x,qty_deriv_y)
  use mat_class
  use grid_class
  use omp_module
  use sml_module, only: sml_nthreads
  implicit none
  type(grid_type), intent(in) :: grid
  real (kind=8), intent(in) :: qty(grid%nnode)
  real (kind=8), intent(inout) :: qty_deriv_x(grid%nnode), qty_deriv_y(grid%nnode)
  real (kind=8) :: dumx(grid%nnode), dumy(grid%nnode)
  integer :: i, ith, i_beg(sml_nthreads), i_end(sml_nthreads)

  qty_deriv_x=0.D0
  qty_deriv_y=0.D0

  call mat_mult(grid%gradientx,qty,dumx(:))
  call mat_mult(grid%gradienty,qty,dumy(:))

  qty_deriv_x=dumx
  qty_deriv_y=dumy

  call split_indices(grid%nnode, sml_nthreads, i_beg, i_end)

!rh For conversion to R-Z coordinates
!$OMP PARALLEL DO &
!$OMP PRIVATE( I, ITH)
  do ith=1,sml_nthreads
    do i=i_beg(ith),i_end(ith)
      qty_deriv_x(i)=dumx(i)*grid%unit_vecs(1,1,i)+dumy(i)*grid%unit_vecs(2,1,i)
      qty_deriv_y(i)=dumx(i)*grid%unit_vecs(1,2,i)+dumy(i)*grid%unit_vecs(2,2,i)
    enddo
  enddo

end subroutine grid_deriv

subroutine v_dot_grad_perp(grid,v,x,output)
  use grid_class
  use omp_module
  use sml_module, only: sml_nthreads
  implicit none

  type(grid_type), intent(in) :: grid
  real (kind=8), intent(in) :: v(grid%nnode,2), x(grid%nnode)
  real (kind=8), intent(out) :: output(grid%nnode)
  real (kind=8), dimension(grid%nnode) :: dxdr, dxdz
  integer :: i, ith, i_beg(sml_nthreads), i_end(sml_nthreads)

  output(:)=0D0
  dxdr(:)=0D0
  dxdz(:)=0D0

  call grid_deriv(grid,x,dxdr,dxdz)

  call split_indices(grid%nnode, sml_nthreads, i_beg, i_end)

!rh For conversion to R-Z coordinates
!$OMP PARALLEL DO &
!$OMP PRIVATE( I, ITH)
  do ith=1,sml_nthreads
    do i=i_beg(ith),i_end(ith)
      output(i)=v(i,1)*dxdr(i) + v(i,2)*dxdz(i)
    enddo
  enddo

end subroutine v_dot_grad_perp

! Computes the gradient of grid based quantities
subroutine grid_theta_deriv(grid,qty,qty_deriv)
  use mat_class
  use grid_class
  implicit none
  type(grid_type), intent(in) :: grid
  real (kind=8), intent(in) :: qty(grid%nnode)
  real (kind=8), intent(inout) :: qty_deriv(grid%nnode)

  qty_deriv=0.D0

  call mat_mult(grid%gradienty,qty,qty_deriv)

end subroutine grid_theta_deriv

! Initialization of poloidal smoothing matrix
subroutine init_pol_smooth_mat(grid)
  use grid_class
  use mat_class
  use smooth_module
  use eq_module, only: eq_x_r, eq_x_z, eq_axis_r, eq_axis_z, is_rgn12
  use sml_module
  implicit none
  type (grid_type), intent(inout) :: grid
  real (kind=8) :: d1, d2, tmp, x(2), dum(3), dltheta_sum, dltheta1, fac
  real (kind=8) :: x_in(2), x_out(2), r, z, psi, p(3), norm, smooth_len
  integer :: width1, width2, dir
  integer :: i, j, k, nnode, itr, nd(3), dum1(1), maxp, dltheta_count
  real (kind=8), allocatable :: weights(:,:), dl(:), dist(:), vol(:), dltheta(:)
  integer, allocatable :: tr_indx(:)
  real (kind=8), parameter :: eps=1.D-8
  real (kind=8), external :: psi_interpol

  nnode=grid%nnode
  width1=smooth_pol_width
  width2=2*width1+1

  allocate(weights(3,width2),tr_indx(width2),dltheta(width2))
  allocate(dl(width2), dist(width2), vol(width2))

  ! Set up matrix for poloidal smoothing
  ! matrix initialization
  call new_mat(grid%smooth_pol,nnode,3*width2)

  do i=1,nnode

    weights(:,:)=0D0
    tr_indx(:)=-1
    dl(:)=0D0
    dist(:)=0D0
    vol(:)=0D0
    dltheta(:)=0D0
    dltheta_sum=0D0
    dltheta_count=0

    x=grid%x(:,i)
    r=x(1)
    z=x(2)
    psi=grid%psi(i)
    d1=sqrt((grid%x(1,i)-eq_axis_r)**2+(grid%x(2,i)-eq_axis_z)**2)
    d2=sqrt((grid%x(1,i)-eq_x_r)**2+(grid%x(2,i)-eq_x_z)**2)
    dum(1)=psi_interpol(x(1),x(2),1,0)
    dum(2)=psi_interpol(x(1),x(2),0,1)
    tmp=sqrt(sum(dum**2))
    if (d1 .gt. eps .and. d2 .gt. eps .and. tmp .gt. sml_gradpsi_limit .and. is_rgn12(r,z,psi)   &
        .and. sml_inpsi .le. psi .and. psi .le. sml_outpsi) then
      ! Condition is distance from X point and axis, lower limit for |grad(Psi)|
      ! point in regions 1/2, and point inside the simulation boundaries
      ! Outside of [sml_inpsi,sml_outpsi], we have no charge density --> potential
      ! should at least be smooth

      call get_char_length(grid,i,dum)
      dltheta(width1+1)=dum(3)
      dltheta_sum=dltheta_sum+dum(3)
      dltheta_count=dltheta_count+1
      !dltheta=min(dltheta,0.07D0)

      ! The point itself
      tr_indx(width1+1)=grid%tr_node(1,i)
      nd=grid%nd(:,tr_indx(width1+1))
      do k=1,3
        if (nd(k)==i) then
          weights(k,width1+1)=1D0
        else
          weights(k,width1+1)=0D0
        endif
      enddo


      ! Look for points in the positive theta-direction
      dir=1
      x_in=x
      do j=1,width1
        fac=1D0
        itr=-1
        do while(fac>0.1D0 .and. itr==-1)
          call get_next_point_tang(x_in,fac*dltheta(width1+j),dir,x_out,dltheta1)
          !rh call search_tr2(grid,x_out,-1D10,itr,p)
          call search_tr2(grid,x_out,itr,p)
          fac=fac*0.5D0
        enddo
        if (itr<1) then
          exit
        else
          weights(:,width1+1+j)=p

          dum1=maxloc(p)
          maxp=grid%nd(dum1(1),itr)
          call get_char_length(grid,maxp,dum)
          dltheta(width1+1+j)=dum(3)
          dltheta_sum=dltheta_sum+dum(3)
          dltheta_count=dltheta_count+1

          tr_indx(width1+1+j)=itr
          dist(width1+1+j)=dist(width1+j)+real(dir,8)*dltheta1
          dl(width1+1+j)=dltheta1
          x_in=x_out
        endif
      enddo

      if (j .lt. width1+1) then
        ! We walked out of the grid
        ! Fill up weight array
        dl(width1+1+j:)=0D0
        dist(width1+1+j:)=dist(width1+j)
        tr_indx(width1+1+j:)=tr_indx(width1+j)
        do k=j,width1
          weights(:,width1+1+k)=weights(:,width1+j)
        enddo
      endif

      ! Look for points in the negative theta-direction
      dir=-1
      x_in=x
      do j=1,width1
        fac=1D0
        itr=-1
        do while(fac>0.1D0 .and. itr==-1)
          call get_next_point_tang(x_in,fac*dltheta(width1+2-j),dir,x_out,dltheta1)
          !rh call search_tr2(grid,x_out,-1D10,itr,p)
          call search_tr2(grid,x_out,itr,p)
          fac=fac*0.5D0
        enddo
        if (itr<1) then
          exit
        else
          weights(:,width1+1-j)=p

          dum1=maxloc(p)
          maxp=grid%nd(dum1(1),itr)
          call get_char_length(grid,maxp,dum)
          dltheta(width1+1-j)=dum(3)
          dltheta_sum=dltheta_sum+dum(3)
          dltheta_count=dltheta_count+1

          tr_indx(width1+1-j)=itr
          dist(width1+1-j)=dist(width1+2-j)+real(dir,8)*dltheta1
          dl(width1+1-j)=dltheta1
          x_in=x_out
        endif
      enddo

      if (j .lt. width1+1) then
        ! We walked out of the grid
        ! Fill up weight array
        dl(1:width1+1-j)=0D0
        dist(1:width1+1-j)=dist(width1+2-j)
        tr_indx(1:width1+1-j)=tr_indx(width1+2-j)
        do k=j,width1
          weights(:,width1+1-k)=weights(:,width1+2-j)
        enddo
      endif

      dum1=minval(tr_indx)

      smooth_len=smooth_pol_d0*dltheta_sum/real(dltheta_count,8)

      ! Now set-up the matrix
      ! a) Determine the normalization factor
      norm=0D0
      vol=0D0
      do j=1,width2
        nd=grid%nd(:,tr_indx(j))
        do k=1,3
          !norm = norm + weights(k,j)*exp(-(dist(j)/smooth_len)**2)
          vol(j)=vol(j)+weights(k,j)*grid%node_vol(nd(k))
        enddo
        norm = norm+exp(-(dist(j)/smooth_len)**2)*vol(j)
      enddo

      ! b) Set the matrix elements
      do j=1,width2
        nd=grid%nd(:,tr_indx(j))
        do k=1,3
          tmp=weights(k,j)*exp(-(dist(j)/smooth_len)**2)*grid%node_vol(nd(k))/norm
          call set_value(grid%smooth_pol,i,nd(k),tmp,1)
        enddo
      enddo

    else
      ! Identity
      call set_value(grid%smooth_pol,i,i,1D0,1)
    endif
  enddo

  deallocate(weights,tr_indx,vol,dl,dist,dltheta)

end subroutine init_pol_smooth_mat


#ifndef OLD_SMOOTH_POL
! Computes the gradient of grid based quantities
subroutine smooth_pol(grid,qty,qty_smooth)
  use mat_class
  use grid_class
  implicit none
  type(grid_type), intent(in) :: grid
  real (kind=8), intent(in) :: qty(grid%nnode)
  real (kind=8), intent(inout) :: qty_smooth(grid%nnode)

  qty_smooth=0.D0

  call mat_mult(grid%smooth_pol,qty,qty_smooth)

end subroutine smooth_pol
#endif




