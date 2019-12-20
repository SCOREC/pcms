subroutine sort_particles(grid,sp)
  use grid_class
  use ptl_module
  use perf_monitor
  use sml_module
  use omp_module , only : split_indices
  implicit none

  type(grid_type) :: grid
  type(species_type) :: sp

  integer :: i, ith, i_beg(sml_nthreads), i_end(sml_nthreads)
  integer :: ilo,ihi,jlo,jhi
  integer :: nx, ny, n, xydim
  real(kind=8) :: xmin,ymin, inv_dx,inv_dy
  integer*8, allocatable, dimension(:) :: gid
  real(kind=8), allocatable, dimension(:) :: xcoord,ycoord 
  real(kind=8), allocatable, dimension(:,:) :: tphase0
  type(ptl_type), allocatable, dimension(:) :: tptl
  integer, allocatable, dimension(:) :: iperm

  if(sp%num==0) return  ! nothing to permute

  call split_indices(sp%num, sml_nthreads, i_beg, i_end)

! --------------
! sort particles
! --------------
  ilo = lbound( grid%guess_table,1)
  ihi = ubound( grid%guess_table,1)
  jlo = lbound( grid%guess_table,2)
  jhi = ubound( grid%guess_table,2)

  allocate( iperm(sp%num))
  allocate( gid(sp%num), xcoord(sp%num), ycoord(sp%num) )

!$omp parallel do  private(ith,i)
  do ith=1,sml_nthreads
     do i=i_beg(ith),i_end(ith)
        gid(i) = sp%ptl(i)%gid
        xcoord(i) = sp%ptl(i)%ph(1)
        ycoord(i) = sp%ptl(i)%ph(2)

        iperm(i) = i
     enddo
  enddo

  call gen_perm( ilo,ihi,jlo,jhi,   &
                 grid%guess_min, grid%inv_guess_d, &
                 sp%num, gid, xcoord, ycoord, iperm )

  deallocate( gid, xcoord, ycoord )
! ------------------------
! rearrange particle ptl and phase0 data in sp
! ------------------------

! --------------
! permute sp%ptl
! --------------
  allocate(tptl(sp%num))
!$omp parallel do private(ith,i)
  do ith=1,sml_nthreads
     do i=i_beg(ith),i_end(ith)
        tptl(i) = sp%ptl( iperm(i) )
     enddo
 enddo

!$omp parallel do private(ith,i)
  do ith=1,sml_nthreads
     do i=i_beg(ith),i_end(ith)
        sp%ptl(i) = tptl(i)
     enddo
  enddo
  deallocate( tptl )

! -----------------
! permute sp%phase0
! -----------------
  allocate(tphase0(lbound(sp%phase0,1):ubound(sp%phase0,1),sp%num))

!$omp parallel do private(ith,i)
  do ith=1,sml_nthreads
     do i=i_beg(ith),i_end(ith)
        tphase0(:,i) = sp%phase0(:, iperm(i) )
     enddo
  enddo

!$omp parallel do private(ith,i)
  do ith=1,sml_nthreads
     do i=i_beg(ith),i_end(ith)
        sp%phase0(:,i) = tphase0(:,i)
     enddo
  enddo
  deallocate( tphase0 )

  deallocate( iperm )

  return
end subroutine sort_particles
