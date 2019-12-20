! This routine is written by E. D'Azevedo originally.
! generalization (making class) by S. Ku and J. Lang (3/8/2012)
!
! double precision is replaced by real (8), to be consistent with other XGC1 subroutines.


module bicub_mod
  implicit none

  type bicub_type
     integer :: nr, nz
     real (8), allocatable, dimension(:) :: rc_cub, zc_cub
     real (8), allocatable, dimension(:,:,:,:) :: acoef_cub
     real (8) :: dr_inv, dz_inv, rmin, zmin
  end type bicub_type


  type(bicub_type) :: psi_bicub, bmod_bicub
  

  interface eval_bicub
     module procedure eval_bicub_2, eval_bicub_1, eval_bicub_0
     ! 0 : function value only
     ! 1 : + 1st derivatives
     ! 2 : + 2nd derivatives
  end interface

  interface bicub_interpol
     module procedure bicub_interpol1,bicub_interpol2, bicub_interpol0
  end interface

contains
  subroutine eval_bicub_2(x,y,xc,yc,acoef,f00,f10,f01,f11,f20,f02)
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
    integer, parameter :: ndeg = 3
    real (8), intent(in) :: x,y,xc,yc
    real (8), intent(in) :: acoef(0:ndeg,0:ndeg)
    real (8), intent(inout) :: f00, f10,f01,f11,f20,f02
    integer, parameter :: idebug = 0
    real (8), parameter :: tol = 1.0d-5
    real (8) :: ftmp, abserr, dx,dy
    logical :: isok
    integer :: i,j
    real (8), dimension(0:ndeg) :: xv,yv
    real (8), dimension(0:ndeg) :: fx, fy
    real (8), dimension(0:ndeg) :: dfx, dfy
    real (8), dimension(0:ndeg) :: dfx2, dfy2
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
       fy(:) = 0.0d0
       dfy(:) = 0.0d0
       dfy2(:) = 0.0d0
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

       fx(:) = 0.0d0
       dfx(:) = 0.0d0
       dfx2(:) = 0.0d0


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
          write(*,*) 'ftmp, f00 ', ftmp, f00
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
          write(*,*) 'ftmp,f01 ', ftmp,f01
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
          write(*,*) 'ftmp,f10,abserr ',ftmp,f10,abserr
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
          write(*,*) 'ftmp,f11 ',ftmp,f11
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
          write(*,*) 'ftmp,f20 ',ftmp,f20
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
          write(*,*) 'ftmp,f02 ',ftmp,f02
       endif

    endif

    return
  end subroutine eval_bicub_2
  subroutine eval_bicub_1(x,y,xc,yc,acoef,f00,f10,f01)
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
    integer, parameter :: ndeg = 3
    real (8), intent(in) :: x,y,xc,yc
    real (8), intent(in) :: acoef(0:ndeg,0:ndeg)
    real (8), intent(inout) :: f00, f10,f01
    integer, parameter :: idebug = 0
    real (8), parameter :: tol = 1.0d-5
    real (8) :: ftmp, abserr, dx,dy
    logical :: isok
    integer :: i,j
    real (8), dimension(0:ndeg) :: xv,yv
    real (8), dimension(0:ndeg) :: fx, fy
    real (8), dimension(0:ndeg) :: dfx, dfy
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
       fy(:) = 0.0d0
       dfy(:) = 0.0d0
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
       fx(:) = 0.0d0
       dfx(:) = 0.0d0
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
          write(*,*) 'ftmp, f00 ', ftmp, f00
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
          write(*,*) 'ftmp,f01 ', ftmp,f01
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
          write(*,*) 'ftmp,f10,abserr ',ftmp,f10,abserr
       endif
    endif
    return
  end subroutine eval_bicub_1
  subroutine eval_bicub_0(x,y,xc,yc,acoef,f00 )
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
    integer, parameter :: ndeg = 3
    real (8), intent(in) :: x,y,xc,yc
    real (8), intent(in) :: acoef(0:ndeg,0:ndeg)
    real (8), intent(inout) :: f00
    integer, parameter :: idebug = 0
    real (8), parameter :: tol = 1.0d-5
    real (8) :: ftmp, abserr
    logical :: isok
    integer :: i,j
    real (8), dimension(0:ndeg) :: xv,yv, fy
    real (8) :: dx, dy
    logical, parameter :: use_intrinsic = .false.
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
    if (use_intrinsic) then
       fy = matmul( acoef(:,:), yv(:) )
       f00 = dot_product( xv(:), fy(:) )
    else
       fy(:) = 0.0d0
       do i=0,ndeg
          do j=0,ndeg
             fy(i) = fy(i) + acoef(i,j) * yv(j)
          enddo
       enddo
       f00 = 0.0d0
       do i=0,ndeg
          f00 = f00 + xv(i) * fy(i)
       enddo
    endif
    if (idebug.ge.1) then
       ! ---------
       ! check f00
       ! ---------
       dx = (x-xc)
       dy = (y-yc)
       ftmp = 0.0d0
       do j=0,ndeg
          do i=0,ndeg
             ftmp = ftmp + acoef(i,j)*(dx**i)*(dy**j)
          enddo
       enddo
       abserr = abs(ftmp - f00)
       isok = (abserr.le.tol).or. &
            & (abserr .le. tol*max(abs(f00),abs(ftmp)))
       if (.not.isok) then
          write(*,*) 'ftmp, f00,abserr ', ftmp, f00,abserr
       endif
    endif
    return
  end subroutine eval_bicub_0
  subroutine gencoef_bicub( xv, yv, xc,yc, fval, acoef )
    implicit none
    ! generate coefficients of bicubic polynomial
    ! using interpolation on a tensor product grid
    !
    ! fval(i,j) = f( xv(i), yv(j) )
    !
    ! note (xc,yc) is offset to improve conditioning
    !
    ! bicubic(x,y) = sum( acoef(i,j)*x^i * y^j, i=0:ndeg,j=0:ndeg)
    !
    integer, parameter :: ndeg = 3
    logical, parameter :: use_inverse = .false.
    real (8),dimension(0:ndeg),intent(in) :: xv
    real (8),dimension(0:ndeg),intent(in) :: yv
    real (8), intent(in) :: xc,yc
    real (8),dimension(0:ndeg,0:ndeg),intent(in)::fval
    real (8),dimension(0:ndeg,0:ndeg),intent(inout)::acoef
    real (8), dimension(0:ndeg,0:ndeg) :: Tx,Ty
    real (8), dimension(0:ndeg,0:ndeg) :: Txinv,Tyinv
    real (8), dimension(0:ndeg,0:ndeg) :: tmpmat, tmpmat2
    integer, dimension(0:ndeg) :: ipivx,ipivy
    integer :: m,n,nrhs, info, ld1, ld2
    integer :: i
    interface
       subroutine dgetrs(trans,n,nrhs,A,lda,ipiv,B,ldb,info)
         character trans
         integer n,nrhs,lda,ldb,info
         integer ipiv(*)
         real (8) A(lda,*),B(ldb,*)
       end subroutine dgetrs
       subroutine dgetrf(m,n,A,lda,ipiv,info)
         integer m,n,lda,info
         integer ipiv(*)
         real (8) A(lda,*)
       end subroutine dgetrf
    end interface
    ipivx(:) = 0
    ipivy(:) = 0
    info = 0
    ! -------------------------------------------------
    ! Tx(:,:), Ty(:,:) initialized with
    ! basis functions (x^j) evaluated at xv(:) or yv(:)
    ! -------------------------------------------------
    Tx(:,0) = 1.0d0
    Tx(:,1) = (xv(:)-xc)
    Tx(:,2) = Tx(:,1)*Tx(:,1)
    Tx(:,3) = Tx(:,1)*Tx(:,2)
    Ty(:,0) = 1.0d0
    Ty(:,1) = (yv(:)-yc)
    Ty(:,2) = Ty(:,1)*Ty(:,1)
    Ty(:,3) = Ty(:,1)*Ty(:,2)
    ! ---------------------------
    ! setup Txinv(:,:), Tyinv(:,:) as identity
    ! ---------------------------
    if (use_inverse) then
       Txinv(:,:) = 0.0d0
       Tyinv(:,:) = 0.0d0
       do i=0,ndeg
          Txinv(i,i) = 1.0d0
          Tyinv(i,i) = 1.0d0
       enddo
    endif
    ! ----------------
    ! generate inverse for Tx, Ty
    ! ----------------
    m = ndeg + 1
    n = ndeg + 1
    ld1 = ndeg+1
    info = 0
    call dgetrf(m,n,Tx,ld1,ipivx,info)
    if (info.ne.0) then
       write(*,*) 'dgetrf Tx return info=',info
       stop '** error in gencoef_bicub '
    endif
    if (use_inverse) then
       n = ndeg + 1
       nrhs = ndeg + 1
       ld1 = size(Tx,1)
       ld2 = size(Txinv,1)
       info = 0
       call dgetrs('N',n,nrhs,Tx,ld1, ipivx,Txinv,ld2,info)
       if (info.ne.0) then
          write(*,*) 'dgetrs Tx return info=',info
          stop '** error in gencoef_bicub'
       endif
    endif
    info = 0
    call dgetrf(m,n,Ty,ld1,ipivy,info)
    if (info.ne.0) then
       write(*,*) 'dgetrf Ty return info=',info
       stop '** error in gencoef_bicub '
    endif
    if (use_inverse) then
       n = ndeg + 1
       nrhs = ndeg + 1
       ld1 = size(Ty,1)
       ld2 = size(Tyinv,1)
       info = 0
       call dgetrs('N',n,nrhs,Ty,ld1, ipivy,Tyinv,ld2,info)
       if (info.ne.0) then
          write(*,*) 'dgetrs Ty return info=',info
          stop '** error in gencoef_bicub'
       endif
    endif
    ! ---------------------------------------
    ! Y = kron(A,B)*X = B*X*transpose(A)
    ! inv(kron(A,B)) = kron( inv(A), inv(B) )
    !
    ! kron( Ty, Tx)*acoef = fval
    !
    ! acoef = inv( kron(Ty,Tx) )
    ! acoef = kron( inv(Ty), inv(Tx) )
    ! acoef = kron( Tyinv, Txinv)*fval
    ! acoef = Txinv * fval * transpose(Tyinv)
    ! ---------------------------------------
    if (use_inverse) then
       tmpmat(:,:) = matmul( Txinv(:,:), fval(:,:) )
       acoef(:,:) = matmul( tmpmat(:,:), transpose( Tyinv(:,:) ) )
    else
       ! ----------------------------------------------
       ! acoef = (Tx \ fval)/Ty
       ! tmpmat = (Tx \ fval)
       !
       ! acoef = tmpmat * inv( transpose(Ty) )
       ! transpose(acoef) = inv(Ty) * transpose(tmpmat)
       !
       ! tmpmat2 = Ty \ (transpose(tmpmat)
       ! acoef = transpose(tmpmat2)
       ! ----------------------------------------------
       n = ndeg + 1
       nrhs = ndeg + 1
       ld1 = size(Tx,1)
       ld2 = size(tmpmat,1)
       tmpmat(:,:) = fval(:,:)
       info = 0
       call dgetrs('N',n,nrhs,Tx,ld1,ipivx,tmpmat,ld2,info)
       if (info.ne.0) then
          write(*,*) 'dgetrs Tx tmpmat return info=',info
          stop '** error in gencoef_bicub'
       endif
       tmpmat2(:,:) = transpose(tmpmat(:,:))
       ld1 = size(Ty,1)
       ld2 = size(tmpmat2,1)
       info = 0
       call dgetrs( 'N', n,nrhs, Ty,ld1,ipivy,tmpmat2,ld2,info)
       if (info.ne.0) then
          write(*,*) 'dgetrs Ty tmpmat2 info=',info
          stop '** error in gencoef_bicub'
       endif
       acoef(:,:) = transpose( tmpmat2(:,:) )
    endif
    return
  end subroutine gencoef_bicub
  subroutine check_bicub()
    implicit none
    integer, parameter :: ndeg = 3
    integer, parameter :: idebug = 0
    real (8), dimension(0:ndeg,0:ndeg) :: acoef
    real (8), dimension(0:ndeg,0:ndeg) :: acoef_org
    real (8), dimension(0:ndeg,0:ndeg) :: fval
    real (8), dimension(0:ndeg) :: xv,yv
    real (8), parameter :: tol = 1.0d-6
    real (8) :: abserr, maxerr,x,y,xc,yc
    real (8) :: f00,f10,f01,f11,f20,f02
    real (8) :: g00,g10,g01
    integer :: i,j, ii,jj
    logical :: isok
    call random_number(acoef)
    acoef = 2.0d0*acoef - 1.0d0
    acoef_org = acoef
    ! ------------------
    ! fit to [0,1]x[0,1]
    ! ------------------
    xv(0) = 0.0d0
    xv(1) = dble(1)/dble(3)
    xv(2) = dble(2)/dble(3)
    xv(3) = 1.0d0
    yv(:) = xv(:)
    xc = dble(1)/dble(2)
    yc = dble(1)/dble(2)
    do j=0,ndeg
       do i=0,ndeg
          call eval_bicub( xv(i),yv(j),xc,yc, acoef, &
               & f00,f10,f01,f11,f20,f02)
          fval(i,j) = f00
       enddo
    enddo
    ! ------------------------------
    ! check 2 versions of eval_bicub
    ! ------------------------------
    do j=0,ndeg
       do i=0,ndeg
          call eval_bicub( xv(i),yv(j),xc,yc, acoef, &
               & f00,f10,f01,f11,f20,f02)
          call eval_bicub( xv(i),yv(j),xc,yc, acoef, &
               & g00,g10,g01 )
          ! ---------
          ! check g00
          ! ---------
          abserr = abs(g00-f00)
          isok = (abserr.le.tol).or. &
               & (abserr.le.tol*max(abs(g00),abs(f00)))
          if (.not.isok) then
             write(*,*) 'g00,f00,abserr ',g00,f00,abserr
             stop '** error in check_bicub'
          endif
          ! ---------
          ! check g01
          ! ---------
          abserr = abs(g01-f01)
          isok = (abserr.le.tol).or. &
               & (abserr.le.tol*max(abs(g01),abs(f01)))
          if (.not.isok) then
             write(*,*) 'g01,f01,abserr ',g01,f01,abserr
             stop '** error in check_bicub'
          endif
          ! ---------
          ! check g10
          ! ---------
          abserr = abs(g10-f10)
          isok = (abserr.le.tol).or. &
               & (abserr.le.tol*max(abs(g10),abs(f10)))
          if (.not.isok) then
             write(*,*) 'g10,f10,abserr ',g10,f10,abserr
             stop '** error in check_bicub'
          endif
       enddo
    enddo
    ! ---------------------------------------
    ! generate coefficients for interpolation
    ! ---------------------------------------
    call gencoef_bicub( xv,yv,xc,yc, fval, acoef )
    maxerr = 0.0d0
    do j=0,ndeg
       do i=0,ndeg
          abserr = abs(acoef(i,j) - acoef_org(i,j))
          maxerr = max( maxerr, abserr )
          if (idebug.ge.2) then
             write(*,*) 'i,j,acoef,acoef_org ', &
                  & i,j,acoef(i,j),acoef_org(i,j)
          endif
       enddo
    enddo
    write(*,*) 'maxerr = ', maxerr
    call random_number(xv)
    call random_number(yv)
    do j=0,ndeg
       do i=0,ndeg
          call eval_bicub( xv(i),yv(j),xc,yc, acoef, &
               & f00,f10,f01,f11,f20,f02)
          fval(i,j) = f00
       enddo
    enddo
    do jj=0,ndeg
       do ii=0,ndeg
          x = xv(ii)
          y = yv(jj)
          f00 = 0.0d0
          do j=0,ndeg
             do i=0,ndeg
                f00 = f00 + acoef(i,j)*((x-xc)**i)*((y-yc)**j)
             enddo
          enddo
          abserr = abs(f00 - fval(ii,jj))
          maxerr = max( maxerr, abserr )
          isok = abserr .le. tol*max(abs(f00),abs(fval(ii,jj)))
          if ((.not.isok).or.(idebug.ge.2)) then
             write(*,*) 'x,y ', x,y
             write(*,*) 'ii,jj,f00,fval(ii,jj) ',ii,jj,f00,fval(ii,jj)
          endif
          call eval_bicub_0( x,y,xc,yc,acoef,f00)
          abserr = abs(f00 - fval(ii,jj))
          maxerr = max( maxerr, abserr )
          isok = abserr.le.max(abs(f00),abs(fval(ii,jj)))
          if ((.not.isok).or.(idebug.ge.2)) then
             write(*,*) 'eval_bicub_0: x,y ', x,y
             write(*,*) 'ii,jj,f00,fval(ii,jj) ',ii,jj,f00,fval(ii,jj)
          endif
       enddo
    enddo
    write(*,*) 'random evaluation maxerr ', maxerr
    return
  end subroutine check_bicub
  subroutine setup_bicub(var,eq_mr,eq_rgrid, eq_mz, eq_zgrid,interpol_pspline)
    implicit none
    integer, intent(in) :: eq_mr, eq_mz
    real (8), dimension(eq_mr), intent(in) :: eq_rgrid
    real (8), dimension(eq_mz), intent(in) :: eq_zgrid
    integer, parameter :: idebug = 0
    integer, parameter :: ndeg = 3
    real (8), parameter :: tol = 1.0d-5
    logical, parameter :: use_psi_interpol_bicub = .true.
    real (8), dimension(0:ndeg) :: rv, zv
    real (8), dimension(0:ndeg,0:ndeg) :: acoef
    real (8), dimension(0:ndeg,0:ndeg) :: fval
    real (8) :: dr, dz
    real (8) :: r,z, f_00, f00,f10,f01,f11,f20,f02
    real (8) :: abserr, rc, zc
    real (8) :: maxabserr
    integer :: nr,nz, i,j, ierr, i0,j0
    logical :: isok
    real (8) :: psi, dpsi_dr,dpsi_dz
    real (8) :: d2psi_d2r,d2psi_drdz,d2psi_d2z
    integer, parameter :: n = ndeg + 1
    real (8), dimension(n) :: dweight
    real (8) :: pi
    type(bicub_type) :: var
    intrinsic :: atan, cos
    interface
       real (8) function interpol_pspline(r,z,ider,idez)
         real (8) r, z
         integer ider, idez
       end function interpol_pspline
    end interface
    nr = eq_mr - 1
    nz = eq_mz - 1
    var%nr=nr
    var%nz=nz
    ! ---------------------------------------------------
    ! assume equal spacing in eq_rgrid(:) and eq_zgrid(:)
    ! ---------------------------------------------------
    dr = eq_rgrid(2)-eq_rgrid(1)
    dz = eq_zgrid(2)-eq_zgrid(1)
    var%dr_inv = dble(1)/dr
    var%dz_inv = dble(1)/dz
    var%rmin = eq_rgrid(1)
    var%zmin = eq_zgrid(1)
    allocate( var%rc_cub(nr), var%zc_cub(nz), &
         & var%acoef_cub( 0:ndeg,0:ndeg, nr, nz) , stat=ierr)
    isok = (ierr.eq.0)
    if (.not.isok) then
       write(*,*) 'allocate psi_acoef return ierr ',ierr
       stop '** error in setup_psi_bicub '
    endif
    ! -----------------------------
    ! setup interpolation positions
    ! over [0,1]
    ! use chebyshev points
    ! -----------------------------
    pi = atan( dble(1) ) * dble(4)
    do i=1,n
       dweight(i) = cos( dble(i)/dble(n+1) * pi )
    enddo
    dweight(:) = (dweight(:) + dble(1))/dble(2)
    do i=1,nr
       var%rc_cub(i) = (eq_rgrid(i) + eq_rgrid(i+1))/2.0d0
    enddo
    do j=1,nz
       var%zc_cub(j) = (eq_zgrid(j) + eq_zgrid(j+1))/2.0d0
    enddo
    do j=1,nz
       do i=1,nr
          rv(0) = eq_rgrid(i) + dweight(1)*dr
          rv(1) = eq_rgrid(i) + dweight(2)*dr
          rv(2) = eq_rgrid(i) + dweight(3)*dr
          rv(3) = eq_rgrid(i) + dweight(4)*dr
          zv(0) = eq_zgrid(j) + dweight(1)*dz
          zv(1) = eq_zgrid(j) + dweight(2)*dz
          zv(2) = eq_zgrid(j) + dweight(3)*dz
          zv(3) = eq_zgrid(j) + dweight(4)*dz
          rc = var%rc_cub(i)
          zc = var%zc_cub(j)
          do j0=0,ndeg
             do i0=0,ndeg
                fval(i0,j0) = interpol_pspline( rv(i0), zv(j0), 0, 0)
             enddo
          enddo
          call gencoef_bicub(rv,zv,rc,zc, fval, acoef )
          var%acoef_cub(:,:, i,j) = acoef(:,:)
       enddo
    enddo
    ! ------------------------------------
    ! extra check for psi_interpol_bicub_2
    ! ------------------------------------
    if (idebug.ge.1) then
       maxabserr = dble(0)
       ierr = 0
       do j=1,nz
          do i=1,nr
             rc = var%rc_cub(i)
             zc = var%zc_cub(j)
             r = var%rc_cub(i)
             z = var%zc_cub(j)
             psi = interpol_pspline(r,z,0,0)
             dpsi_dr = interpol_pspline(r,z,1,0)
             dpsi_dz = interpol_pspline(r,z,0,1)
             d2psi_d2r = interpol_pspline(r,z,2,0)
             d2psi_drdz = interpol_pspline(r,z,1,1)
             d2psi_d2z = interpol_pspline(r,z,0,2)
             if (use_psi_interpol_bicub) then
                call bicub_interpol(var,r,z, f00,f10,f01,f11,f20,f02)
             else
                acoef(:,:) = var%acoef_cub(:,:,i,j)
                call eval_bicub( r,z, rc,zc, acoef, &
                     & f00,f10,f01,f11,f20,f02)
             endif
             ! ---------
             ! check psi
             ! ---------
             abserr = abs(f00 - psi)
             maxabserr = max( abserr, maxabserr )
             isok = (abserr.le.tol).or. &
                  & (abserr .le. tol * max( abs(f00),abs(psi) ))
             if (.not.isok) then
                write(*,*) 'i,j, f00,psi,abserr ', &
                     & i,j, f00,psi,abserr
                ierr = ierr + 1
             endif
             !f_00 = bicub_interpol0(var,r,z)
             call bicub_interpol0(var,r,z,f_00)  ! from  function to subroutine --pathscale
             abserr = abs(f_00 - psi)
             maxabserr = max( abserr, maxabserr )
             isok = (abserr.le.tol).or. &
                  & (abserr .le. tol * max( abs(f_00),abs(psi) ))
             if (.not.isok) then
                write(*,*) 'i,j, f_00,psi,abserr ', &
                     & i,j, f_00,psi,abserr
                ierr = ierr + 1
             endif
             ! -------------
             ! check dpsi_dr
             ! -------------
             abserr = abs(f10 - dpsi_dr)
             maxabserr = max( abserr, maxabserr )
             isok = (abserr.le.tol).or. &
                  & (abserr .le. tol*max(abs(f10),abs(dpsi_dr)))
             if (.not.isok) then
                write(*,*) 'i,j, f10,dpsi_dr,abserr ', &
                     & i,j, f10,dpsi_dr,abserr
                ierr = ierr + 1
             endif
             ! -------------
             ! check dpsi_dz
             ! -------------
             abserr = abs(f01 - dpsi_dz)
             maxabserr = max( abserr, maxabserr )
             isok = (abserr.le.tol).or. &
                  & (abserr .le. tol*max( abs(f01),abs(dpsi_dz)))
             if (.not.isok) then
                write(*,*) 'i,j, f01,dpsi_dz,abserr', &
                     & i,j, f01,dpsi_dz,abserr
                ierr = ierr + 1
             endif
             ! ----------------
             ! check d2psi_drdz
             ! ----------------
             abserr = abs(f11 - d2psi_drdz)
             maxabserr = max( abserr, maxabserr )
             isok = (abserr.le.tol).or. &
                  & (abserr .le. tol*max(abs(f11),abs(d2psi_drdz)))
             if (.not.isok) then
                write(*,*) 'i,j,f11,d2psi_drdz,abserr', &
                     & i,j,f11,d2psi_drdz,abserr
                ierr = ierr + 1
             endif
             ! ---------------
             ! check d2psi_d2r
             ! ---------------
             abserr = abs(f20 - d2psi_d2r)
             maxabserr = max( abserr, maxabserr )
             isok = (abserr.le.tol).or. &
                  & (abserr .le. tol*max(abs(f20),abs(d2psi_d2r)))
             if (.not.isok) then
                write(*,*) 'i,j, f20,d2psi_d2r,abserr ', &
                     & i,j, f20,d2psi_d2r,abserr
                ierr = ierr + 1
             endif
             ! ---------------
             ! check d2psi_d2z
             ! ---------------
             abserr = abs(f02 - d2psi_d2z)
             maxabserr = max( abserr, maxabserr )
             isok = (abserr.le.tol).or. &
                  & (abserr .le. tol*max(abs(f02),abs(d2psi_d2z)))
             if (.not.isok) then
                write(*,*) 'i,j, f02,d2psi_d2z,abserr ', &
                     & i,j, f02,d2psi_d2z,abserr
                ierr = ierr + 1
             endif
          enddo
       enddo
       isok = (ierr.eq.0)
       if (.not.isok) then
          stop '** error in setup_psi_bicub'
       endif
       write(*,*) 'psi_interpol_bicub_2 maxabserr = ', maxabserr
    endif
    ! ------------------------------------
    ! extra check for psi_interpol_bicub_1
    ! ------------------------------------
    if (idebug.ge.1) then
       maxabserr = dble(0)
       ierr = 0
       do j=1,nz
          do i=1,nr
             rc = var%rc_cub(i)
             zc = var%zc_cub(j)
             r = var%rc_cub(i)
             z = var%zc_cub(j)
             psi = interpol_pspline(r,z,0,0)
             dpsi_dr = interpol_pspline(r,z,1,0)
             dpsi_dz = interpol_pspline(r,z,0,1)
             if (use_psi_interpol_bicub) then
                call bicub_interpol(var,r,z, f00,f10,f01)
             else
                acoef(:,:) = var%acoef_cub(:,:,i,j)
                call eval_bicub( r,z, rc,zc, acoef, &
                     & f00,f10,f01)
             endif
             ! ---------
             ! check psi
             ! ---------
             abserr = abs(f00 - psi)
             maxabserr = max( abserr, maxabserr )
             isok = (abserr.le.tol).or. &
                  & (abserr .le. tol * max( abs(f00),abs(psi) ))
             if (.not.isok) then
                write(*,*) 'i,j, f00,psi,abserr ', &
                     & i,j, f00,psi,abserr
                ierr = ierr + 1
             endif
             !f_00 = bicub_interpol0(var,r,z)
             call bicub_interpol0(var,r,z,f_00)
             abserr = abs(f_00 - psi)
             maxabserr = max( abserr, maxabserr )
             isok = (abserr.le.tol).or. &
                  & (abserr .le. tol * max( abs(f_00),abs(psi) ))
             if (.not.isok) then
                write(*,*) 'i,j, f_00,psi,abserr ', &
                     & i,j, f_00,psi,abserr
                ierr = ierr + 1
             endif
             ! -------------
             ! check dpsi_dr
             ! -------------
             abserr = abs(f10 - dpsi_dr)
             maxabserr = max( abserr, maxabserr )
             isok = (abserr.le.tol).or. &
                  & (abserr .le. tol*max(abs(f10),abs(dpsi_dr)))
             if (.not.isok) then
                write(*,*) 'i,j, f10,dpsi_dr,abserr ', &
                     & i,j, f10,dpsi_dr,abserr
                ierr = ierr + 1
             endif
             ! -------------
             ! check dpsi_dz
             ! -------------
             abserr = abs(f01 - dpsi_dz)
             maxabserr = max( abserr, maxabserr )
             isok = (abserr.le.tol).or. &
                  & (abserr .le. tol*max( abs(f01),abs(dpsi_dz)))
             if (.not.isok) then
                write(*,*) 'i,j, f01,dpsi_dz,abserr', &
                     & i,j, f01,dpsi_dz,abserr
                ierr = ierr + 1
             endif
          enddo
       enddo
       isok = (ierr.eq.0)
       if (.not.isok) then
          stop '** error in setup_bicub'
       endif
       write(*,*) 'psi_interpol_bicub_1 maxabserr = ', maxabserr
    endif
    return
  end subroutine setup_bicub

  subroutine bicub_interpol2(var,r,z, f00,f10,f01,f11,f20,f02)
    implicit none
    type(bicub_type) :: var
    real (8), intent(in) :: r,z
    real (8), intent(inout) :: f00,f10,f01,f11,f20,f02
    integer, parameter :: idebug = 0
    real (8), parameter :: tol = 1.0d-5
    real (8) :: abserr, g00,g10,g01,g11,g20,g02
    integer :: ierr
    logical :: isok
    interface
       ! only for psi debugging
       real (8) function psi_interpol_pspline(r,z,ider,idez)
         real (8) r,z
         integer ider,idez
       end function psi_interpol_pspline
    end interface
    integer :: i,j
    ! ---------------
    ! perform hashing
    ! ---------------
    !nr = ubound(psi_rc,1)
    !nz = ubound(psi_zc,1)
    i = max(1,min(var%nr, 1 + int( (r-var%rmin)*var%dr_inv ) ))
    j = max(1,min(var%nz, 1 + int( (z-var%zmin)*var%dz_inv ) ))
    call eval_bicub(r,z, var%rc_cub(i), var%zc_cub(j), var%acoef_cub(:,:,i,j), &
         & f00,f10,f01,f11,f20,f02 )
    if (idebug.ge.1) then
       g00 = psi_interpol_pspline(r,z,0,0)
       g10 = psi_interpol_pspline(r,z,1,0)
       g01 = psi_interpol_pspline(r,z,0,1)
       g11 = psi_interpol_pspline(r,z,1,1)
       g20 = psi_interpol_pspline(r,z,2,0)
       g02 = psi_interpol_pspline(r,z,0,2)
       ierr = 0
       abserr = abs(g00-f00)
       isok = (abserr.le.tol).or. &
            & (abserr .le. tol*max( abs(g00), abs(f00) ))
       if (.not.isok) then
          write(*,*) 'r,z,g00,f00,abserr', r,z,g00,f00,abserr
          ierr = ierr + 1
       endif
       abserr = abs(g10-f10)
       isok = (abserr.le.tol).or. &
            & (abserr .le. tol*max( abs(g10), abs(f10) ))
       if (.not.isok) then
          write(*,*) 'r,z,g10,f10,abserr', r,z,g10,f10,abserr
          ierr = ierr + 1
       endif
       abserr = abs(g01-f01)
       isok = (abserr.le.tol).or. &
            & (abserr .le. tol*max( abs(g01), abs(f01) ))
       if (.not.isok) then
          write(*,*) 'r,z,g01,f01,abserr', r,z,g01,f01,abserr
          ierr = ierr + 1
       endif
       abserr = abs(g11-f11)
       isok = (abserr.le.tol).or. &
            & (abserr .le. tol*max( abs(g11), abs(f11) ))
       if (.not.isok) then
          write(*,*) 'r,z,g11,f11,abserr', r,z,g11,f11,abserr
          ierr = ierr + 1
       endif
       abserr = abs(g20-f20)
       isok = (abserr.le.tol).or. &
            & (abserr .le. tol*max( abs(g20), abs(f20) ))
       if (.not.isok) then
          write(*,*) 'r,z,g20,f20,abserr', r,z,g20,f20,abserr
          ierr = ierr + 1
       endif
       abserr = abs(g02-f02)
       isok = (abserr.le.tol).or. &
            & (abserr .le. tol*max( abs(g02), abs(f02) ))
       if (.not.isok) then
          write(*,*) 'r,z,g02,f02,abserr', r,z,g02,f02,abserr
          ierr = ierr + 1
       endif
       isok = (ierr.eq.0)
       if (.not.isok) then
          stop '** error in psi_interpol_bicub_2'
       endif
    endif
    return
  end subroutine bicub_interpol2
  subroutine bicub_interpol1(var,r,z, f00,f10,f01)
    implicit none
    type(bicub_type) :: var
    real (8), intent(in) :: r,z
    real (8), intent(inout) :: f00,f10,f01
    integer, parameter :: idebug = 0
    real (8), parameter :: tol = 1.0d-5
    real (8) :: abserr, g00,g10,g01
    integer :: ierr
    logical :: isok
    interface
       ! only for psi debugging
       real (8) function psi_interpol_pspline(r,z,ider,idez)
         real (8) r,z
         integer ider,idez
       end function psi_interpol_pspline
    end interface
    integer ::  i,j
    ! ---------------
    ! perform hashing
    ! ---------------
    !nr = ubound(psi_rc,1)
    !nz = ubound(psi_zc,1)
    i = max(1,min(var%nr, 1 + int( (r-var%rmin)*var%dr_inv ) ))
    j = max(1,min(var%nz, 1 + int( (z-var%zmin)*var%dz_inv ) ))
    call eval_bicub(r,z, var%rc_cub(i), var%zc_cub(j), var%acoef_cub(:,:,i,j), &
         & f00,f10,f01)
    if (idebug.ge.1) then
       g00 = psi_interpol_pspline(r,z,0,0)
       g10 = psi_interpol_pspline(r,z,1,0)
       g01 = psi_interpol_pspline(r,z,0,1)
       ierr = 0
       abserr = abs(g00-f00)
       isok = (abserr.le.tol).or. &
            & (abserr .le. tol*max( abs(g00), abs(f00) ))
       if (.not.isok) then
          write(*,*) 'r,z,g00,f00,abserr', r,z,g00,f00,abserr
          ierr = ierr + 1
       endif
       abserr = abs(g10-f10)
       isok = (abserr.le.tol).or. &
            & (abserr .le. tol*max( abs(g10), abs(f10) ))
       if (.not.isok) then
          write(*,*) 'r,z,g10,f10,abserr', r,z,g10,f10,abserr
          ierr = ierr + 1
       endif
       abserr = abs(g01-f01)
       isok = (abserr.le.tol).or. &
            & (abserr .le. tol*max( abs(g01), abs(f01) ))
       if (.not.isok) then
          write(*,*) 'r,z,g01,f01,abserr', r,z,g01,f01,abserr
          ierr = ierr + 1
       endif
       isok = (ierr.eq.0)
       if (.not.isok) then
          stop '** error in psi_interpol_bicub_1'
       endif
    endif
    return
  end subroutine bicub_interpol1
  subroutine bicub_interpol0(var,r,z,f00)
    implicit none
    type(bicub_type):: var
    real (8), intent(in) :: r,z
    real (8), intent(out) :: f00
    interface
       ! only for psi debugging
       real (8) function psi_interpol_pspline(r,z,ider,idez)
         real (8) r,z
         integer ider,idez
       end function psi_interpol_pspline
    end interface
    integer, parameter :: idebug = 0
    real (8), parameter :: tol = 1.0d-5
    real (8) :: g00, abserr
    logical :: isok
    !real (8) :: f00
    integer :: i,j
    ! ---------------
    ! perform hashing
    ! ---------------
    !nr = ubound(psi_rc,1)
    !nz = ubound(psi_zc,1)
    i = max(1,min(var%nr, 1 + int( (r-var%rmin)*var%dr_inv ) ))
    j = max(1,min(var%nz, 1 + int( (z-var%zmin)*var%dz_inv ) ))
    call eval_bicub_0(r,z, var%rc_cub(i), var%zc_cub(j), &
         & var%acoef_cub(:,:,i,j), f00)
    if (idebug.ge.1) then
       g00 = psi_interpol_pspline(r,z,0,0)
       abserr = abs(g00 - f00)
       isok = (abserr .le. tol * max( abs(g00), abs(f00) )).or. &
            & (abserr.le.tol)
       if (.not.isok) then
          write(*,*) 'r,z, f00,g00,abserr ', r,z,f00,g00,abserr
          stop '** error bicub_interpol0 '
       endif
    endif
    !bicub_interpol0 = f00
    return
  end subroutine bicub_interpol0
end module bicub_mod
