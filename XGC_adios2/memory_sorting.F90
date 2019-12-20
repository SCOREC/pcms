! need to be updated for handling ptl(i) data structure

subroutine memory_sorting(spall,use_permute_inplace)
  use ptl_module
  use eq_module
#if defined(CAM_TIMERS)
  use perf_mod
#endif
  use sml_module, only : sml_electron_on
  use perm_module, permute_copy => permscopy
  implicit none
  type(species_type):: spall(0:ptl_nsp_max)
  logical  :: use_permute_inplace 
  integer :: isp


!!$  do isp=ptl_isp, ptl_nsp
!!$     call mem_sort_one_sp(spall(isp))
!!$  enddo
!!$#if defined(CAM_TIMERS)
!!$  if (use_permute_inplace) then
!!$     call t_stopf('memory_sorting:inplace')
!!$  else
!!$     call t_stopf('memory_sorting:copy')
!!$  endif
!!$#endif
!!$  return
!!$contains
!!$  subroutine mem_sort_one_sp(sp)
!!$    implicit none
!!$    type(species_type),target :: sp
!!$
!!$    real(kind=8), dimension(:,:), pointer :: rz
!!$    integer, dimension( (eq_mr-1)*(eq_mz-1)+1) :: xstart
!!$    integer :: lb,ub, lb2,ub2, ierr, i,n
!!$    integer, allocatable, dimension(:) :: iperm
!!$    integer :: ncycle, ldim
!!$    integer, allocatable, dimension(:) :: leader
!!$
!!$#if defined(CAM_TIMERS)
!!$    if (use_permute_inplace) then
!!$       call t_startf('memory_sorting:inplace')
!!$    else
!!$       call t_startf('memory_sorting:copy')
!!$    endif
!!$#endif
!!$
!!$    n = sp%num
!!$
!!$    allocate( iperm(1:n), leader(1:n), stat=ierr)
!!$    call assert( ierr.eq.0,                                          &
!!$         &         'memory_sorting: allocate( iperm ) ',ierr )
!!$
!!$    rz => sp%phase(1:2, 1:sp%num )
!!$    call geomhash_rz( n, rz, xstart, iperm )
!!$
!!$
!!$    if (use_permute_inplace) then
!!$       ldim = n
!!$       allocate( leader(ldim), stat=ierr )
!!$       call assert( ierr.eq.0,                                          &
!!$            &         'memory_sorting: allocate( leader ) ',ierr )
!!$
!!$#if defined(CAM_TIMERS)
!!$       call t_startf('gencycle')
!!$       call gencycle( n, iperm, ncycle, leader, ldim )
!!$       call t_stopf('gencycle')
!!$    endif
!!$#endif 
!!$
!!$    !         ---------------
!!$    !         copy phase(:,:)
!!$    !         ---------------
!!$    if (use_permute_inplace) then
!!$       call permute_inplace(n,iperm,ncycle,leader, sp%phase )
!!$    else
!!$       call permute_copy(n,iperm,    sp%phase)
!!$    endif
!!$
!!$
!!$    !         ----------------
!!$    !         copy phase0(:,:)
!!$    !         ----------------
!!$    if (use_permute_inplace) then
!!$       call permute_inplace(n,iperm,ncycle,leader, sp%phase0 )
!!$    else
!!$       call permute_copy(n,iperm,    sp%phase0)
!!$    endif
!!$
!!$
!!$
!!$    !         -----------
!!$    !         copy gid(:)
!!$    !         -----------
!!$    if (use_permute_inplace) then
!!$       call permute_inplace(n,iperm,ncycle,leader, sp%gid )
!!$    else
!!$       call permute_copy(n,iperm,    sp%gid )
!!$    endif
!!$
!!$    !         -----------
!!$    !         copy derivs(:,:,1:n)
!!$    !         -----------
!!$    if (use_permute_inplace) then
!!$       call permute_inplace(n,iperm,ncycle,leader, sp%derivs )
!!$    else
!!$       call permute_copy(n,iperm,    sp%derivs )
!!$    endif
!!$
!!$    if (use_permute_inplace) then
!!$       call permute_inplace(n,iperm,ncycle,leader, sp%ct )
!!$    else
!!$       call permute_copy(n,iperm,    sp%derivs )
!!$    endif
!!$
!!$    if (use_permute_inplace) then
!!$       call permute_inplace(n,iperm,ncycle,leader, sp%dph )
!!$    else
!!$       call permute_copy(n,iperm,    sp%derivs )
!!$    endif
!!$
!!$    if (use_permute_inplace) then
!!$       call permute_inplace(n,iperm,ncycle,leader, sp%dpht )
!!$    else
!!$       call permute_copy(n,iperm,    sp%derivs )
!!$    endif
!!$
!!$
!!$
!!$
!!$
!!$
!!$#ifdef NEED_RESET_DERIVS
!!$    !         --------------------
!!$    !         copy reset_derivs(:)
!!$    !         --------------------
!!$    if (use_permute_inplace) then
!!$       call permute_inplace(n,iperm,ncycle,leader, sp%reset_derivs)
!!$    else
!!$       call permute_copy(n,iperm,   sp%reset_derivs)
!!$    endif
!!$#endif
!!$    !         -------------
!!$    !         copy angle(:)
!!$    !         -------------
!!$
!!$    deallocate( iperm, leader, stat=ierr)
!!$    call assert( ierr.eq.0,                                            &
!!$         &         'memory_sorting: deallocate( iperm ) ',ierr )
!!$  end subroutine mem_sort_one_sp


end subroutine memory_sorting
                              
