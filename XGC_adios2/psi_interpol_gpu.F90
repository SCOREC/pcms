!cdir$r unroll 
#ifdef USE_GPU
  attributes(device) &
#endif
real (kind=8) function psi_interpol_gpu(r_in,z_in,r_der,z_der)
  use itp_module_gpu
  use bicub_mod_gpu, only : bicub_interpol1
  use precision_mod_gpu

  implicit none
  real (kind=work_p), intent(in) :: r_in, z_in
  integer , intent(in) :: r_der, z_der
!  real (kind=8) :: r,z
!  real (kind=8) ,external:: psi_interpol_pspline
  real (kind=work_p) :: dpsi_dr, dpsi_dz, psi


  call bicub_interpol1(r_in,z_in,psi,dpsi_dr,dpsi_dz)
  if(r_der==0 .and. z_der==0) then
     psi_interpol_gpu=psi
  elseif(r_der==1 .and. z_der==0) then
     psi_interpol_gpu=dpsi_dr
  elseif(r_der==0 .and. z_der==1) then
     psi_interpol_gpu=dpsi_dz
  else
     !print *, 'Do not use second derivative of psi_interpol with USE_BICUB_MOD.'
     !print *, 'Replace it with psi_interpol_bicub(r,z,psi,dpsi_dr,dpsi_dz,d2psi_drdz,d2psi_dr2,d2psi_dz2).'
     !stop 
  endif

  if (r_der==0 .and. z_der==0) then
     psi_interpol_gpu = max(1D-99, psi_interpol_gpu)
  endif
  
end function psi_interpol_gpu

!sa : Unrolled evals in this file are only kept for archival purposes.
!sa : They don't actually seem to be used on the GPU.

!cdir$r unroll 
#ifdef USE_GPU
  attributes(device) &
#endif
  subroutine eval_bicub_2_loop(x,y,xc,yc,acoef,f00,f10,f01,f11,f20,f02)
    use bicub_mod_gpu, only : ndeg
    use precision_mod_gpu
    implicit none
    ! ----------------------------------
    ! evaluate bicubic polynomial f(x,y)
    ! and high order derivatives
    !
    ! note (xc,yc) is offset or center of box
    !
    ! f00 = f(x,y)
    ! f10 = df/dx
    ! f01 = df/dy
    ! f11 = df^2/dx/dy
    ! f20 = df^2/dx/dx
    ! f02 = df^2/dy/dy
    ! ----------------------------------
    !integer, parameter :: ndeg = 3
    real (kind=work_p), intent(in) :: x,y,xc,yc
    real (kind=work_p), intent(in) :: acoef(0:ndeg,0:*)
    real (kind=work_p), intent(inout) :: f00, f10,f01,f11,f20,f02
    integer, parameter :: idebug = 0
    real (kind=work_p), parameter :: tol = 1.0d-5
    real (kind=work_p) :: ftmp, abserr, dx,dy
    logical :: isok
    integer :: i,j
    real (kind=work_p), dimension(0:ndeg) :: xv,yv
    real (kind=work_p), dimension(0:ndeg) :: fx, fy
    real (kind=work_p), dimension(0:ndeg) :: dfx, dfy
    real (kind=work_p), dimension(0:ndeg) :: dfx2, dfy2
    logical, parameter :: use_fy = .false.
    f00 = 0.0d0
    f01 = 0.0d0
    f10 = 0.0d0
    f11 = 0.0d0
    f20 = 0.0d0
    f02 = 0.0d0
    dx = (x-xc)
    dy = (y-yc)
    xv(0) = 1.0d0
    xv(1) = dx
    xv(2) = dx*dx
    xv(3) = (dx*dx)*dx
    yv(0) = 1.0d0
    yv(1) = dy
    yv(2) = dy*dy
    yv(3) = (dy*dy)*dy
    if (use_fy) then
       fy(0:ndeg) = 0.0d0
       dfy(0:ndeg) = 0.0d0
       dfy2(0:ndeg) = 0.0d0
       do i=0,ndeg
          do j=0,ndeg
             fy(i) = fy(i) + acoef(i,j)*yv(j)
          enddo
          do j=1,ndeg
             dfy(i) = dfy(i) + acoef(i,j)*dble(j)*yv(j-1)
          enddo
          do j=2,ndeg
             dfy2(i) = dfy2(i) + acoef(i,j)*dble(j*(j-1))*yv(j-2)
          enddo
       enddo
       ! ----------------
       ! f00 = f(x,y)
       ! f10 = df/dx
       ! f01 = df/dy
       ! f11 = df^2/dx/dy
       ! f20 = df^2/dx/dx
       ! f02 = df^2/dy/dy
       ! ----------------
       do i=0,ndeg
          f00 = f00 + xv(i)*fy(i)
          f01 = f01 + xv(i)*dfy(i)
          f02 = f02 + xv(i)*dfy2(i)
       enddo
       do i=1,ndeg
          dfx(i) = dble(i)*xv(i-1)

          f10 = f10 + dfx(i)*fy(i)
          f11 = f11 + dfx(i)*dfy(i)
       enddo

       do i=2,ndeg
          dfx2(i) = dble(i*(i-1))*xv(i-2)

          f20 = f20 + dfx2(i)*fy(i)
       enddo


    else

       fx(0:ndeg) = 0.0d0
       dfx(0:ndeg) = 0.0d0
       dfx2(0:ndeg) = 0.0d0


       do j=0,ndeg
          do i=0,ndeg
             fx(j) = fx(j) + xv(i)*acoef(i,j)
          enddo
          do i=1,ndeg
             dfx(j) = dfx(j) + dble(i)*xv(i-1)*acoef(i,j)
          enddo
          do i=2,ndeg
             dfx2(j) = dfx2(j) + dble(i*(i-1))*xv(i-2)*acoef(i,j)
          enddo
       enddo

       do j=0,ndeg
          f00 = f00 + fx(j)*yv(j)
          f10 = f10 + dfx(j)*yv(j)
          f20 = f20 + dfx2(j)*yv(j)
       enddo



       do j=1,ndeg
          dfy(j) = dble(j)*yv(j-1)

          f01 = f01 + fx(j)*dfy(j)
          f11 = f11 + dfx(j)*dfy(j)
       enddo

       do j=2,ndeg
          dfy2(j) = dble(j*(j-1))*yv(j-2)

          f02 = f02 + fx(j)*dfy2(j)
       enddo




    endif

#if (0)
    if (idebug.ge.1) then

       ! ---------
       ! check f00
       ! ---------
       ftmp = 0.0d0
       do j=0,ndeg
          do i=0,ndeg
             ftmp = ftmp + acoef(i,j)*(dx**i)*(dy**j)
          enddo
       enddo
       abserr = abs(ftmp - f00)
       isok = (abserr.le.tol) .or. &
            & (abserr .le. tol*max(abs(f00),abs(ftmp)))
       if (.not.isok) then
!           write(*,*) 'ftmp, f00 ', ftmp, f00
       endif

       ! ---------
       ! check f01
       ! ---------

       ftmp = 0.0d0
       do j=1,ndeg
          do i=0,ndeg
             ftmp = ftmp + acoef(i,j)*( dx**i)*(j*dy**(j-1))
          enddo
       enddo
       abserr = abs(ftmp - f01)
       isok = (abserr.le.tol) .or. &
            & (abserr .le. tol*max(abs(f01),abs(ftmp)))
       if (.not.isok) then
!           write(*,*) 'ftmp,f01 ', ftmp,f01
       endif

       ! ---------
       ! check f10
       ! ---------
       ftmp = 0.0d0
       do j=0,ndeg
          do i=1,ndeg
             ftmp = ftmp + acoef(i,j)*(i*dx**(i-1))*dy**j
          enddo
       enddo
       abserr = abs(ftmp - f10)
       isok = (abserr.le.tol).or. &
            & (abserr .le. tol*max(abs(f10),abs(ftmp)))
       if (.not.isok) then
!           write(*,*) 'ftmp,f10,abserr ',ftmp,f10,abserr
       endif


       ! ---------
       ! check f11
       ! ---------
       ftmp = 0.0d0
       do j=1,ndeg
          do i=1,ndeg
             ftmp = ftmp + acoef(i,j)*(i*dx**(i-1))*(j*dy**(j-1))
          enddo
       enddo
       abserr = abs(ftmp - f11)
       isok = (abserr.le.tol).or. &
            & (abserr .le. tol*max(abs(f11),abs(ftmp)))
       if (.not.isok) then
!           write(*,*) 'ftmp,f11 ',ftmp,f11
       endif


       ! ---------
       ! check f20
       ! ---------
       ftmp = 0.0d0
       do j=0,ndeg
          do i=2,ndeg
             ftmp = ftmp + acoef(i,j)*(i*(i-1))*dx**(i-2)*dy**j
          enddo
       enddo
       abserr = abs(ftmp - f20)
       isok = (abserr.le.tol).or. &
            & (abserr .le. tol*max(abs(f20),abs(ftmp)))
       if (.not.isok) then
!          write(*,*) 'ftmp,f20 ',ftmp,f20
       endif



       ! ---------
       ! check f02
       ! ---------
       ftmp = 0.0d0
       do j=2,ndeg
          do i=0,ndeg
             ftmp = ftmp + acoef(i,j)*dx**i*(j*(j-1)*dy**(j-2))
          enddo
       enddo
       abserr = abs(ftmp - f02)
       isok = (abserr.le.tol).or. &
            & (abserr .le. tol*max(abs(f02),abs(ftmp)))
       if (.not.isok) then
!           write(*,*) 'ftmp,f02 ',ftmp,f02
       endif

    endif
#endif

    return
  end subroutine eval_bicub_2_loop

!cdir$r unroll 
#ifdef USE_GPU
  attributes(device) &
#endif
  subroutine eval_bicub_1_loop(x,y,xc,yc,acoef,f00,f10,f01)
    use bicub_mod_gpu, only : ndeg
    use precision_mod_gpu
    implicit none
    ! ----------------------------------
    ! evaluate bicubic polynomial f(x,y)
    ! and high order derivatives
    !
    ! note (xc,yc) is offset or center of box
    !
    ! f00 = f(x,y)
    ! f10 = df/dx
    ! f01 = df/dy
    ! f11 = df^2/dx/dy
    ! f20 = df^2/dx/dx
    ! f02 = df^2/dy/dy
    ! ----------------------------------
    !integer, parameter :: ndeg = 3
    real (kind=work_p), intent(in) :: x,y,xc,yc
    real (kind=work_p), intent(in) :: acoef(0:ndeg,0:*)
    real (kind=work_p), intent(inout) :: f00, f10,f01
    integer, parameter :: idebug = 0
    real (kind=work_p), parameter :: tol = 1.0d-5
    real (kind=work_p) :: ftmp, abserr, dx,dy
    logical :: isok
    integer :: i,j
    real (kind=work_p), dimension(0:ndeg) :: xv,yv
    real (kind=work_p), dimension(0:ndeg) :: fx, fy
    real (kind=work_p), dimension(0:ndeg) :: dfx, dfy
    logical, parameter :: use_fy = .false.
    f00 = 0.0d0
    f01 = 0.0d0
    f10 = 0.0d0
    dx = (x-xc)
    dy = (y-yc)
    xv(0) = 1.0d0
    xv(1) = dx
    xv(2) = dx*dx
    xv(3) = (dx*dx)*dx
    yv(0) = 1.0d0
    yv(1) = dy
    yv(2) = dy*dy
    yv(3) = (dy*dy)*dy
    if (use_fy) then
       fy(0:ndeg) = 0.0d0
       dfy(0:ndeg) = 0.0d0
       do i=0,ndeg
          do j=0,ndeg
             fy(i) = fy(i) + acoef(i,j)*yv(j)
          enddo
          do j=1,ndeg
             dfy(i) = dfy(i) + acoef(i,j)*dble(j)*yv(j-1)
          enddo
       enddo
       ! ----------------
       ! f00 = f(x,y)
       ! f10 = df/dx
       ! f01 = df/dy
       ! f11 = df^2/dx/dy
       ! f20 = df^2/dx/dx
       ! f02 = df^2/dy/dy
       ! ----------------
       do i=0,ndeg
          f00 = f00 + xv(i)*fy(i)
          f01 = f01 + xv(i)*dfy(i)
       enddo
       do i=1,ndeg
          dfx(i) = dble(i)*xv(i-1)
          f10 = f10 + dfx(i)*fy(i)
       enddo
    else
       fx(0:ndeg) = 0.0d0
       dfx(0:ndeg) = 0.0d0
       do j=0,ndeg
          do i=0,ndeg
             fx(j) = fx(j) + xv(i)*acoef(i,j)
          enddo
          do i=1,ndeg
             dfx(j) = dfx(j) + dble(i)*xv(i-1)*acoef(i,j)
          enddo
       enddo
       do j=0,ndeg
          f00 = f00 + fx(j)*yv(j)
          f10 = f10 + dfx(j)*yv(j)
       enddo
       do j=1,ndeg
          dfy(j) = dble(j)*yv(j-1)
          f01 = f01 + fx(j)*dfy(j)
       enddo
    endif

#if (0)
    if (idebug.ge.1) then
       ! ---------
       ! check f00
       ! ---------
       ftmp = 0.0d0
       do j=0,ndeg
          do i=0,ndeg
             ftmp = ftmp + acoef(i,j)*(dx**i)*(dy**j)
          enddo
       enddo
       abserr = abs(ftmp - f00)
       isok = (abserr.le.tol) .or. &
            & (abserr .le. tol*max(abs(f00),abs(ftmp)))
       if (.not.isok) then
!           write(*,*) 'ftmp, f00 ', ftmp, f00
       endif
       ! ---------
       ! check f01
       ! ---------
       ftmp = 0.0d0
       do j=1,ndeg
          do i=0,ndeg
             ftmp = ftmp + acoef(i,j)*( dx**i)*(j*dy**(j-1))
          enddo
       enddo
       abserr = abs(ftmp - f01)
       isok = (abserr.le.tol) .or. &
            & (abserr .le. tol*max(abs(f01),abs(ftmp)))
       if (.not.isok) then
!           write(*,*) 'ftmp,f01 ', ftmp,f01
       endif
       ! ---------
       ! check f10
       ! ---------
       ftmp = 0.0d0
       do j=0,ndeg
          do i=1,ndeg
             ftmp = ftmp + acoef(i,j)*(i*dx**(i-1))*dy**j
          enddo
       enddo
       abserr = abs(ftmp - f10)
       isok = (abserr.le.tol).or. &
            & (abserr .le. tol*max(abs(f10),abs(ftmp)))
       if (.not.isok) then
!           write(*,*) 'ftmp,f10,abserr ',ftmp,f10,abserr
       endif
    endif
#endif
    return
  end subroutine eval_bicub_1_loop

!cdir$r unroll 
#ifdef USE_GPU
  attributes(device) &
#endif
  subroutine eval_bicub_00_loop(x,y,xc,yc,acoef,f00 )
    use bicub_mod_gpu, only : ndeg
    use precision_mod_gpu
    implicit none
    ! ----------------------------------
    ! evaluate bicubic polynomial f(x,y)
    ! and high order derivatives
    !
    ! note (xc,yc) is offset or center of box
    !
    ! f00 = f(x,y)
    ! f10 = df/dx
    ! f01 = df/dy
    ! f11 = df^2/dx/dy
    ! f20 = df^2/dx/dx
    ! f02 = df^2/dy/dy
    ! ----------------------------------
    !integer, parameter :: ndeg = 3
    real (kind=work_p), intent(in) :: x,y,xc,yc
    real (kind=work_p), intent(in) :: acoef(0:ndeg,0:*)
    real (kind=work_p), intent(inout) :: f00
    integer, parameter :: idebug = 0
    real (kind=work_p), parameter :: tol = 1.0d-5
    real (kind=work_p) :: ftmp, abserr
    logical :: isok
    integer :: i,j
    real (kind=work_p), dimension(0:ndeg) :: xv,yv, fy
    real (kind=work_p) :: dx, dy
    dx = (x-xc)
    xv(0) = 1.0d0
    xv(1) = dx
    xv(2) = dx*dx
    xv(3) = (dx*dx)*dx
    dy = (y-yc)
    yv(0) = 1.0d0
    yv(1) = dy
    yv(2) = dy*dy
    yv(3) = (dy*dy)*dy

       fy(0:ndeg) = 0.0d0
       do i=0,ndeg
          do j=0,ndeg
             fy(i) = fy(i) + acoef(i,j) * yv(j)
          enddo
       enddo
       f00 = 0.0d0
       do i=0,ndeg
          f00 = f00 + xv(i) * fy(i)
       enddo
    return
  end subroutine eval_bicub_00_loop


!cdir$r unroll 
#ifdef USE_GPU
   attributes(device) &
#endif
  subroutine eval_bicub_2_unroll(x,y,xc,yc,A,f00,f10,f01,f11,f20,f02)
    use bicub_mod_gpu
    use precision_mod_gpu
    implicit none
    ! ----------------------------------
    ! evaluate bicubic polynomial f(x,y)
    ! and high order derivatives
    !
    ! note (xc,yc) is offset or center of box
    !
    ! f00 = f(x,y)
    ! f10 = df/dx
    ! f01 = df/dy
    ! f11 = df^2/dx/dy
    ! f20 = df^2/dx/dx
    ! f02 = df^2/dy/dy
    ! ----------------------------------
    !integer, parameter :: ndeg = 3
    real (kind=work_p), intent(in) :: x,y,xc,yc
    real (kind=work_p), intent(in) :: A(0:ndeg,0:*)
    real (kind=work_p), intent(inout) :: f00, f10,f01,f11,f20,f02
    integer, parameter :: idebug = 0
    real (kind=work_p), parameter :: tol = 1.0d-5
    real (kind=work_p) :: ftmp, abserr, dx,dy
    logical :: isok
    integer :: i,j
!     real (kind=work_p), dimension(0:ndeg) :: xv,yv
!     real (kind=work_p), dimension(0:ndeg) :: fx 
!     real (kind=work_p), dimension(0:ndeg) :: dfx
!     real (kind=work_p), dimension(0:ndeg) :: dfx2

   real(kind=work_p), dimension(0:ndeg) :: fx,dfx,dfx2



!    f00 = 0.0d0
!    f01 = 0.0d0
!    f10 = 0.0d0
!    f11 = 0.0d0
!    f20 = 0.0d0
!    f02 = 0.0d0
    dx = (x-xc)
    dy = (y-yc)

!     xv(0) = 1.0d0
!     xv(1) = dx
!     xv(2) = dx*dx
!     xv(3) = (dx*dx)*dx
!     yv(0) = 1.0d0
!     yv(1) = dy
!     yv(2) = dy*dy
!     yv(3) = (dy*dy)*dy


!       fx(0:ndeg) = 0.0d0
!       dfx(0:ndeg) = 0.0d0
!       dfx2(0:ndeg) = 0.0d0


!       do j=0,ndeg
!          do i=0,ndeg
!             fx(j) = fx(j) + xv(i)*A(i,j)
!          enddo
!       enddo

        j = 0
        fx(j) = A(0,j) + dx*A(1,j) + (dx*dx)*A(2,j) + ((dx*dx)*dx)*A(3,j)
        j = 1
        fx(j) = A(0,j) + dx*A(1,j) + (dx*dx)*A(2,j) + ((dx*dx)*dx)*A(3,j)
        j = 2
        fx(j) = A(0,j) + dx*A(1,j) + (dx*dx)*A(2,j) + ((dx*dx)*dx)*A(3,j)
        j = 3
        fx(j) = A(0,j) + dx*A(1,j) + (dx*dx)*A(2,j) + ((dx*dx)*dx)*A(3,j)

!!          dfy(j) = dble(j)*yv(j-1)
!!          dfy2(j) = dble(j*(j-1))*yv(j-2)
!          f01 = f01 + fx(j)*dfy(j), j=1:ndeg
!          f02 = f02 + fx(j)*dfy2(j), j=2:ndeg

       f00 = fx(0) + dy*fx(1) + (dy*dy)*fx(2) + ((dy*dy)*dy)*fx(3)
       f01 =            fx(1) + 2*dy*fx(2) +   3*(dy*dy)*fx(3)
       f02 =                    2*fx(2)    +  3*2*dy*fx(3)


!       do j=0,ndeg
!          do i=1,ndeg
!             dfx(j) = dfx(j) + dble(i)*xv(i-1)*A(i,j)
!          enddo
!       enddo

        j = 0
        dfx(j) = A(1,j) + 2*dx*A(2,j) + 3*(dx*dx)*A(3,j)
        j = 1
        dfx(j) = A(1,j) + 2*dx*A(2,j) + 3*(dx*dx)*A(3,j)
        j = 2
        dfx(j) = A(1,j) + 2*dx*A(2,j) + 3*(dx*dx)*A(3,j)
        j = 3
        dfx(j) = A(1,j) + 2*dx*A(2,j) + 3*(dx*dx)*A(3,j)



!       do j=1,ndeg
!!          dfy(j) = dble(j)*yv(j-1)
!
!          f11 = f11 + dfx(j)*dfy(j)
!        enddo


!          f11 = f11 + dfx(j)*dfy(j), j=1:ndeg
!          f10 = f10 + dfx(j)*yv(j), j=0:ndeg
       f10 = dfx(0) + dy*dfx(1) + (dy*dy)*dfx(2) + ((dy*dy)*dy)*dfx(3)
       f11 =             dfx(1) + 2*dy*dfx(2)    + 3*dy*dy*dfx(3)


!       do j=0,ndeg
!          do i=2,ndeg
!             dfx2(j) = dfx2(j) + dble(i*(i-1))*xv(i-2)*A(i,j)
!          enddo
!       enddo

       j = 0
       dfx2(j) = (3*2*A(3,j)*dx + 2*1*A(2,j))
       j = 1
       dfx2(j) = (3*2*A(3,j)*dx + 2*1*A(2,j))
       j = 2
       dfx2(j) = (3*2*A(3,j)*dx + 2*1*A(2,j))
       j = 3
       dfx2(j) = (3*2*A(3,j)*dx + 2*1*A(2,j))

!        do j=0,ndeg
!           f20 = f20 + dfx2(j)*yv(j)
!        enddo

       f20 = dfx2(0) + dy*dfx2(1) + (dy*dy)*dfx2(2) + ((dy*dy)*dy)*dfx2(3)


!       do j=1,ndeg
!!          dfy(j) = dble(j)*yv(j-1)
!
!          f01 = f01 + fx(j)*dfy(j)
!          f11 = f11 + dfx(j)*dfy(j)
!       enddo

!       do j=2,ndeg
!!          dfy2(j) = dble(j*(j-1))*yv(j-2)
!
!          f02 = f02 + fx(j)*dfy2(j)
!       enddo



    return
  end subroutine eval_bicub_2_unroll






!cdir$r unroll 
#ifdef USE_GPU
   attributes(device) &
#endif
  subroutine eval_bicub_1_unroll(x,y,xc,yc,A,f00,f10,f01)
    use bicub_mod_gpu
    use precision_mod_gpu
    implicit none
    ! ----------------------------------
    ! evaluate bicubic polynomial f(x,y)
    ! and high order derivatives
    !
    ! note (xc,yc) is offset or center of box
    !
    ! f00 = f(x,y)
    ! f10 = df/dx
    ! f01 = df/dy
    ! f11 = df^2/dx/dy
    ! f20 = df^2/dx/dx
    ! f02 = df^2/dy/dy
    ! ----------------------------------
    !integer, parameter :: ndeg = 3
    real (kind=work_p), intent(in) :: x,y,xc,yc
    real (kind=work_p), intent(in) :: A(0:ndeg,0:*)
    real (kind=work_p), intent(inout) :: f00, f10,f01
    integer, parameter :: idebug = 0
    real (kind=work_p), parameter :: tol = 1.0d-5
    real (kind=work_p) :: ftmp, abserr, dx,dy
    logical :: isok
    integer :: i,j
!     real (kind=work_p), dimension(0:ndeg) :: xv,yv
!     real (kind=work_p), dimension(0:ndeg) :: fx, fy
!     real (kind=work_p), dimension(0:ndeg) :: dfx, dfy

    real (kind=work_p), dimension(0:ndeg) :: fx,dfx
 
!    f00 = 0.0d0
!    f01 = 0.0d0
!    f10 = 0.0d0
    dx = (x-xc)
    dy = (y-yc)

!     xv(0) = 1.0d0
!     xv(1) = dx
!     xv(2) = dx*dx
!     xv(3) = (dx*dx)*dx
!     yv(0) = 1.0d0
!     yv(1) = dy
!     yv(2) = dy*dy
!     yv(3) = (dy*dy)*dy

!       fx(0:ndeg) = 0.0d0
!       dfx(0:ndeg) = 0.0d0
!       do j=0,ndeg
!          do i=0,ndeg
!             fx(j) = fx(j) + xv(i)*A(i,j)
!          enddo
!          do i=1,ndeg
!             dfx(j) = dfx(j) + dble(i)*xv(i-1)*A(i,j)
!          enddo
!       enddo


       j = 0
       fx(j) = A(0,j) + dx*A(1,j) + (dx*dx)*A(2,j) + ((dx*dx)*dx)*A(3,j)
       j = 1
       fx(j) = A(0,j) + dx*A(1,j) + (dx*dx)*A(2,j) + ((dx*dx)*dx)*A(3,j)
       j = 2
       fx(j) = A(0,j) + dx*A(1,j) + (dx*dx)*A(2,j) + ((dx*dx)*dx)*A(3,j)
       j = 3
       fx(j) = A(0,j) + dx*A(1,j) + (dx*dx)*A(2,j) + ((dx*dx)*dx)*A(3,j)

       f00 = fx(0) + dy*fx(1) + (dy*dy)*fx(2) + ((dy*dy)*dy)*fx(3)
       f01 =            fx(1) + 2*dy*fx(2)    + 3*(dy*dy)*fx(3)

       j = 0
       dfx(j) =         A(1,j) + 2*dx*A(2,j)   + 3*(dx*dx)*A(3,j)
       j = 1
       dfx(j) =         A(1,j) + 2*dx*A(2,j)   + 3*(dx*dx)*A(3,j)
       j = 2
       dfx(j) =         A(1,j) + 2*dx*A(2,j)   + 3*(dx*dx)*A(3,j)
       j = 3
       dfx(j) =         A(1,j) + 2*dx*A(2,j)   + 3*(dx*dx)*A(3,j)

       f10 =  dfx(0) + dy*dfx(1) + (dy*dy)*dfx(2) + ((dy*dy)*dy)*dfx(3)

!       do j=0,ndeg
!          f00 = f00 + fx(j)*yv(j)
!          f10 = f10 + dfx(j)*yv(j)
!       enddo
!
!       do j=1,ndeg
!          dfy(j) = dble(j)*yv(j-1)
!          f01 = f01 + fx(j)*dfy(j)
!       enddo

    return
  end subroutine eval_bicub_1_unroll




!cdir$r unroll 
#ifdef USE_GPU
   attributes(device) &
#endif
  subroutine eval_bicub_00_unroll(x,y,xc,yc,A,f00 )
    use bicub_mod_gpu, only : ndeg
    use precision_mod_gpu
    implicit none
    ! ----------------------------------
    ! evaluate bicubic polynomial f(x,y)
    ! and high order derivatives
    !
    ! note (xc,yc) is offset or center of box
    !
    ! f00 = f(x,y)
    ! f10 = df/dx
    ! f01 = df/dy
    ! f11 = df^2/dx/dy
    ! f20 = df^2/dx/dx
    ! f02 = df^2/dy/dy
    ! ----------------------------------
    !integer, parameter :: ndeg = 3
    real (kind=work_p), intent(in) :: x,y,xc,yc
    real (kind=work_p), intent(in) :: A(0:ndeg,0:*)
    real (kind=work_p), intent(inout) :: f00
    integer :: i
!     real (kind=work_p), dimension(0:ndeg) :: xv,yv
    real (kind=work_p), dimension(0:ndeg) :: fy
    real (kind=work_p) :: dx, dy
    dx = (x-xc)
    dy = (y-yc)

!     xv(0) = 1.0d0
!     xv(1) = dx
!     xv(2) = dx*dx
!     xv(3) = (dx*dx)*dx
!     yv(0) = 1.0d0
!     yv(1) = dy
!     yv(2) = dy*dy
!     yv(3) = (dy*dy)*dy

!       fy(0:ndeg) = 0.0d0
!       do i=0,ndeg
!          do j=0,ndeg
!             fy(i) = fy(i) + A(i,j) * yv(j)
!          enddo
!       enddo

        i = 0
        fy(i) = A(i,0) + A(i,1)*dy + A(i,2)*(dy*dy) + A(i,3)*((dy*dy)*dy)
        i = 1
        fy(i) = A(i,0) + A(i,1)*dy + A(i,2)*(dy*dy) + A(i,3)*((dy*dy)*dy)
        i = 2
        fy(i) = A(i,0) + A(i,1)*dy + A(i,2)*(dy*dy) + A(i,3)*((dy*dy)*dy)
        i = 3
        fy(i) = A(i,0) + A(i,1)*dy + A(i,2)*(dy*dy) + A(i,3)*((dy*dy)*dy)

!       f00 = 0.0d0
!       do i=0,ndeg
!          f00 = f00 + xv(i) * fy(i)
!       enddo

    f00 = fy(0) + dx*fy(1) + (dx*dx)*fy(2) + ((dx*dx)*dx)*fy(3)
    return
  end subroutine eval_bicub_00_unroll

