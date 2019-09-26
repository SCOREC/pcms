#include "redef.h"
#include "intrinsic_sizes.h"
#include "switches.h"
Module diag_Gyro_LES_real_3d
  Use all_rhs_terms, only: this_nonlinear_term
  Use par_mod
  Use file_io, only: get_unit_nr
  use communications
  use compute_f, only: f_, h_
  use antenna, only: antenna_type
  use collisions, only: equ_collisions
  use diagnostics_energy, only: get_cfgamma_f,get_cfgamma_h,get_cfgamma_ant,&
       &energy_integral
  use dfdzv_terms
  use dgdxy_terms
  use dzv_terms
  use dfdxy_terms, only: add_dfdxy
  use dchidz_term, only: add_dchidz
  use dchidxy_terms, only: add_dchidxy
  use spatial_averages, only: sum_int_3d, sum_int_z
  use prefactors
  use numerical_damping
  USE calc_rhs, only: this_nonlinear_term
  use diag_Gyro_LES_common
  use fourier

  Implicit None

  public:: initialize_diag_spectral_threeD, finalize_diag_spectral_threeD, diag_spectral_threeD

  !****** Quantities for 3D spectral diags
  private
  Character(Len=8) :: filestat='replace', filepos='rewind'
  integer, dimension(4) :: fnumber ! Will store the file id's for all quantity/species couples
  integer, dimension(:,:),allocatable:: fnumber_spec ! Will store the file id's for all quantity/species couples
  Complex, Dimension(:,:,:,:,:,:), Allocatable :: temp_rhs

contains

  !******************************************************************************************
  !************************************ entropy/electrostatic balance ***********************
  !******************************************************************************************

  !****************************************************************************************************
  !**************************************  Diag 3D in real space **************************************
  !****************************************************************************************************

  Subroutine initialize_diag_spectral_threeD

#ifdef with_extended_diags

    Character(len=20), dimension(4) :: label
    integer :: n, m

    ! Arrays to cumulate for time average
    allocate(fnumber_spec(ln1:ln2,1:4))

    call init_kperp

    fnumber=0
    fnumber_spec=0

    ! We initialize the label to be appended to the name of file
    label(1)='collision'
    label(2)='nonlinear'
    label(3)='parallel'
    label(4)='fe'

     do m=1,4
       if (mype.eq.0) then
           call get_unit_nr(fnumber(m))
           OPEN(fnumber(m), file=trim(diagdir)//'/Real_3D_'&
               //trim(adjustl(label(m)))//trim(file_extension),form = 'FORMATTED', &
               status=filestat, position=filepos)
       endif
     enddo

     if (mype.eq.pexyzvwspec(0,0,0,0,0,my_pespec)) then
        do n=ln1,ln2
         do m=1,4
           call get_unit_nr(fnumber_spec(n,m))
           OPEN(fnumber_spec(n,m), file=trim(diagdir)//'/Real_3D_'&
               //trim(adjustl(label(m)))//'_'//trim(spec(n)%name)//trim(file_extension),form = 'FORMATTED', &
               status=filestat, position=filepos)
         enddo
        enddo
     endif

#endif

  End Subroutine initialize_diag_spectral_threeD

  Subroutine finalize_diag_spectral_threeD

#ifdef with_extended_diags
    Implicit None
    Integer:: n,m

    if (mype.eq.0) then
       do m=1,4
          close(fnumber(m))
       enddo
    endif

    If (mype.eq.pexyzvwspec(0,0,0,0,0,my_pespec)) then
        do n=ln1,ln2
            do m=1,4
                close(fnumber_spec(n,m))
            enddo
        enddo
    endif

    deallocate(fnumber_spec)
#endif

  End Subroutine finalize_diag_spectral_threeD


!!!******************************************************************************************
!!! Subroutine that goes from 6d to 3d     ***************************************************
!!!*******************************************************************************************

  Subroutine diag_spectral_threeD

#ifdef with_extended_diags
    implicit none

    integer :: n,lbg1,lbg2
    complex, dimension(:,:,:),pointer :: ptr2
    complex, dimension(:,:,:,:),pointer :: ptr1, ptr3
    complex, dimension(:,:,:,:,:,:), pointer :: ptr4
    real, dimension(0:ly0da-1, 0:li0da/n_procs_y-1 ,lk1:lk2, ln1:ln2) :: v4d_real
    real, dimension(0:ly0da-1, 0:li0da/n_procs_y-1 ,lk1:lk2) :: v3d_real


    Allocate(cfgamma(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
    Allocate(cfgamma_real(0:ly0da-1, 0:li0da/n_procs_y-1, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
    Allocate(temp_rhs(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
    Allocate(temp_rhs_real(0:ly0da-1, 0:li0da/n_procs_y-1, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
    Allocate(v6d_real(0:ly0da-1, 0:li0da/n_procs_y-1, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))

    ptr1 => null()
    ptr2 => null()
    ptr3 => null()
    ptr4 => null()
    lbg1 = 1
    lbg2 = lklmn0

    if (arakawa_zv) then
       call get_cfgamma_h(h_,cfgamma)
    else
       call get_cfgamma_f(f_,emfields,cfgamma) !because we compute f
    endif

    !------------------------------------------------------------------

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!  for the collision term  / 1
    temp_rhs = cmplx(0.0,0.0)
    call equ_collisions(f_,temp_rhs,replace_rhs=.true.)

    cfgamma_real = 0.0
    temp_rhs_real = 0.0
    v6d_real = 0.0
    v4d_real = 0.0
    v3d_real = 0.0
    call nl_to_direct(temp_rhs,temp_rhs_real)
    call nl_to_direct(conjg(cfgamma),cfgamma_real)
    v6d_real = cfgamma_real*temp_rhs_real
    ! total for check
    call threeD_array_real(v6d_real,v3d_real)
    call write3d_real(v3d_real,1)
    ! each species independent
    do n=ln1,ln2
       call threeD_array_real_spec(v6d_real(:,:,:,:,:,n),v4d_real(:,:,:,n))
       call write3d_real_spec(v4d_real(:,:,:,n),n,1)
    enddo


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! for the nonlinear transfer / 2
    temp_rhs = cmplx(0.0,0.0)
    if (.not.nonlin_h) then
        call this_nonlinear_term%add(g_1, ptr1, emfields,&
        &ptr2, ptr3, temp_rhs, lbg1, lbg2, 0)
    else
        call this_nonlinear_term%add(h_, ptr1, emfields,&
        &ptr2, ptr3, temp_rhs, lbg1, lbg2, 0)
    endif

    cfgamma_real = 0.0
    temp_rhs_real = 0.0
    v6d_real = 0.0
    v4d_real = 0.0
    v3d_real = 0.0
    call nl_to_direct(temp_rhs,temp_rhs_real)
    call nl_to_direct(conjg(cfgamma),cfgamma_real)
    v6d_real = cfgamma_real*temp_rhs_real
    ! total for check
    call threeD_array_real(v6d_real,v3d_real)
    call write3d_real(v3d_real,2)
   ! each species independent
    do n=ln1,ln2
       call threeD_array_real_spec(v6d_real(:,:,:,:,:,n),v4d_real(:,:,:,n))
       call write3d_real_spec(v4d_real(:,:,:,n),n,2)
    enddo


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! for the parallel term /3
    temp_rhs = cmplx(0.0,0.0)
    if (arakawa_zv) then
        call equ_dzv(h_,temp_rhs,lbg1,lbg2)
    if (hyp_on_h) then
        if (hypz_compensation) call equ_comp_hypz(emfields,ptr4,temp_rhs,lbg1,lbg2)
        endif
    else
        call equ_dfdzv(f_,temp_rhs,lbg1,lbg2)
        call add_dchidz(emfields, ptr4, temp_rhs, lbg1, lbg2)
    end if

    temp_rhs_real = 0.0
    v6d_real = 0.0
    v4d_real = 0.0
    v3d_real = 0.0
    call nl_to_direct(temp_rhs,temp_rhs_real)
    v6d_real = cfgamma_real*temp_rhs_real
    ! total for check
    call threeD_array_real(v6d_real,v3d_real)
    call write3d_real(v3d_real,3)
    ! each species independent
    do n=ln1,ln2
       call threeD_array_real_spec(v6d_real(:,:,:,:,:,n),v4d_real(:,:,:,n))
       call write3d_real_spec(v4d_real(:,:,:,n),n,3)
    enddo

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Free energy / 4
    temp_rhs = cmplx(0.0,0.0)
    temp_rhs(:,:,:,:,:,:) = g_1/2.

    temp_rhs_real = 0.0
    v6d_real = 0.0
    v4d_real = 0.0
    v3d_real = 0.0
    call nl_to_direct(temp_rhs,temp_rhs_real)
    v6d_real = cfgamma_real*temp_rhs_real
    ! total for check
    call threeD_array_real(v6d_real,v3d_real)
    call write3d_real(v3d_real,4)
    ! each species independent
    do n=ln1,ln2
       call threeD_array_real_spec(v6d_real(:,:,:,:,:,n),v4d_real(:,:,:,n))
       call write3d_real_spec(v4d_real(:,:,:,n),n,4)
    enddo



    deallocate(temp_rhs)
    deallocate(temp_rhs_real)
    deallocate(cfgamma)
    deallocate(cfgamma_real)
    deallocate(v6d_real)

#endif

  End subroutine diag_spectral_threeD

!!!*****************************************************************************************
!!!***********************  tools  *********************************************************
#ifdef with_extended_diags

  subroutine threeD_array_real_spec(v5d,v3d)

    implicit none
    real, dimension(0:ly0da-1,0:li0da/n_procs_y-1,lk1:lk2,ll1:ll2,lm1:lm2), intent(in) :: v5d
    real, dimension(0:ly0da-1,0:li0da/n_procs_y-1,lk1:lk2), intent(out) :: v3d
    integer :: k,l,m
    integer :: nla

    nla = li0da*ly0da/n_procs_y

    v3d = 0.0
    do m=lm1,lm2
        do l=ll1,ll2
            do k=lk1,lk2
               call real_axpy_ij(nla,mat_00(pi1,pj1,k,l,m),v5d(:,:,k,l,m),&
                        &v3d(:,:,k))
            enddo
        enddo
    enddo
    call my_real_sum_vw(v3d,size(v3d))

  end subroutine threeD_array_real_spec

  subroutine threeD_array_real(v6d,v3d)

    implicit none
    real, dimension(0:ly0da-1,0:li0da/n_procs_y-1,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2), intent(in) :: v6d
    real, dimension(0:ly0da-1,0:li0da/n_procs_y-1,lk1:lk2), intent(out) :: v3d
    integer :: k,l,m,n
    integer :: nla

    nla = li0da*ly0da/n_procs_y

    v3d = 0.0
    do n=ln1,ln2
        do m=lm1,lm2
            do l=ll1,ll2
                do k=lk1,lk2
                   call real_axpy_ij(nla,mat_00(pi1,pj1,k,l,m),v6d(:,:,k,l,m,n),&
                        &v3d(:,:,k))
                enddo
            enddo
        enddo
    enddo
    call my_real_sum_vwspec(v3d,size(v3d))

  end subroutine threeD_array_real


!######################################################################################

  subroutine write3d_real(arr3d,m)

    Real, Dimension(0:ly0da-1,0:li0da/n_procs_y-1,lk1:lk2), intent(in):: arr3d
    Real, Dimension(0:ly0da-1,0:li0da/n_procs_y-1,lk1:lk2) :: arrx
    Real, Dimension(0:ly0da-1,0:li0da-1,lk1:lk2) :: arrxy
    Real, Dimension(0:ly0da-1,0:li0da-1, 0:nz0-1)  ::  dest3d
    Integer :: i,j,k,m,ierr,wfile

    wfile = fnumber(m)

    if ((my_pev+my_pew+my_pespec).eq.0) then


          ! Gather energy over x-Distribution on pey==x.
          Do k = lk1, lk2
             Do j = 0, li0da/n_procs_y - 1
                Call mpi_gather(arr3d(0,j,k), ly0da, MPI_REAL_TYPE,&
                     arrx(0,j,k), ly0da, MPI_REAL_TYPE,&
                     0, mpi_comm_x, ierr)
             Enddo
          Enddo

          ! Gather energy over y-Distribution on pey==0.
          Do k = lk1, lk2
             Call mpi_gather(arrx(0,0,k), Size(arrx(:,:,k)), MPI_REAL_TYPE,&
                  arrxy(0,0,k), Size(arrx(:,:,k)), MPI_REAL_TYPE,&
                  0, mpi_comm_y, ierr)
          Enddo

          ! Gather energyf over z-Distribution on pez=0
          CALL mpi_gather(arrxy(0,0,lk1),SIZE(arrxy(:,:,:)), MPI_REAL_TYPE,&
               dest3d(0,0,0),SIZE(arrxy),MPI_REAL_TYPE,&
               0, mpi_comm_z, ierr)


       if (mype.eq.0) then
             do k=0,nz0-1
               do j=0, li0da-1
                  do i = 0, ly0da -1
                     write(wfile,"(3I5,ES12.4)") i,j,k, dest3d(i,j,k)
                  enddo
               enddo
             enddo
           call flush(wfile)
       endif
    endif

    call my_barrier()

  End Subroutine write3d_real

  subroutine write3d_real_spec(arr3d,n,m)

    Real, Dimension(0:ly0da-1,0:li0da/n_procs_y-1,lk1:lk2), intent(in):: arr3d
    Real, Dimension(0:ly0da-1,0:li0da/n_procs_y-1,lk1:lk2) :: arrx
    Real, Dimension(0:ly0da-1,0:li0da-1,lk1:lk2) :: arrxy
    Real, Dimension(0:ly0da-1,0:li0da-1, 0:nz0-1)  ::  dest3d
    Integer :: i,j,k,n,m,ierr,wfile


    wfile = fnumber_spec(n,m)

    if ((my_pev+my_pew).eq.0) then

          ! Gather energy over x-Distribution on pey==x.
          Do k = lk1, lk2
             Do j = 0, li0da/n_procs_y - 1
                Call mpi_gather(arr3d(0,j,k), ly0da, MPI_REAL_TYPE,&
                     arrx(0,j,k), ly0da, MPI_REAL_TYPE,&
                     0, mpi_comm_x, ierr)
             Enddo
          Enddo

          ! Gather energy over y-Distribution on pey==0.
          Do k = lk1, lk2
             Call mpi_gather(arrx(0,0,k), Size(arrx(:,:,k)), MPI_REAL_TYPE,&
                  arrxy(0,0,k), Size(arrx(:,:,k)), MPI_REAL_TYPE,&
                  0, mpi_comm_y, ierr)
          Enddo

          ! Gather energyf over z-Distribution on pez=0
          CALL mpi_gather(arrxy(0,0,lk1),SIZE(arrxy(:,:,:)), MPI_REAL_TYPE,&
               dest3d(0,0,0),SIZE(arrxy),MPI_REAL_TYPE,&
               0, mpi_comm_z, ierr)


       if (mype.eq. pexyzvwspec(0,0,0,0,0,my_pespec) ) then

             do k=0,nz0-1
               do j=0, li0da-1
                  do i = 0, ly0da -1
                     write(wfile,"(3I5,ES12.4)") i,j,k, dest3d(i,j,k)
                  enddo
               enddo
             enddo
           call flush(wfile)
       endif
    endif

    call my_barrier()

  End Subroutine write3d_real_spec

!######################################################################################

  Subroutine nl_to_direct(inarr,outarr)
    Complex,dimension(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2 ), Intent(In):: inarr
    Complex,dimension(li1da:li2da, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2 ) :: inarr_x
    Real, dimension(0:ly0da-1,0:li0da/n_procs_y-1,lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2),   Intent(out):: outarr
    integer :: j,k,l,m,n

    outarr = 0.0

      Do n=ln1,ln2
        Do m=lm1,lm2
          do l=ll1,ll2
            do k=lk1,lk2
              do j=lj1,lj2
                 inarr_x(li1:hkx,j,k,l,m,n) = inarr(li1:hkx,j,k,l,m,n)
                 inarr_x(hkx+1:lkx+kx_offset-1,j,k,l,m,n)= cmplx(0.0,0.0)
                 inarr_x(lkx+kx_offset:li2da,j,k,l,m,n) = inarr(lkx:li2,j,k,l,m,n)
              enddo
            enddo
          enddo
        enddo
      enddo

      Do n=ln1,ln2
        Do m=lm1,lm2
          do l=ll1,ll2
            do k=lk1,lk2
                call nl_to_direct_xy_1d_mod(inarr_x(:,:,k,l,m,n),outarr(:,:,k,l,m,n))
            enddo
          enddo
        enddo
      enddo


  end subroutine nl_to_direct

!######################################################################################
  Subroutine nl_to_direct_xy_1d_mod(inarr,rearr)

    Complex,dimension(li1da:li2da, lj1:lj2), Intent(In):: inarr
    Real, dimension(0:ly0da-1, 0:li0da/n_procs_y-1), Intent(InOut):: rearr

    Complex,dimension(0:ly0da/2, 0:li0da/n_procs_y-1):: temparr
    Complex:: rbuf(0:lj0-1,0:li0da/n_procs_y-1,0:n_procs_y-1)
    Complex:: tempypar(lj1:lj2, li1da:li2da)
    integer:: i, pe, ierr

    !Fourier transform array (first x direction, then transposition, then y)
    if (n_procs_y.gt.1) then
       Call fft_kx_to_x(inarr,tempypar)
       call mpi_alltoall(tempypar(lj1,0),lj0*li0da/n_procs_y,MPI_COMPLEX_TYPE,&
            rbuf(0,0,0),lj0*li0da/n_procs_y,MPI_COMPLEX_TYPE,mpi_comm_y,ierr)
       do i=0,li0da/n_procs_y-1
          do pe=0,n_procs_y-1
             temparr(pe*lj0:(pe+1)*lj0-1,i)=rbuf(:,i,pe)
          enddo
       enddo
    else
       Call fft_kx_to_x(inarr,temparr)
    endif

    temparr(nky0:, :) = 0.

    Call fft_ky_to_y(temparr,rearr)

  End Subroutine nl_to_direct_xy_1d_mod


  subroutine real_axpy_ij(N, A, X, Y)
    INTEGER, INTENT(IN) :: N
    REAL,INTENT(IN) :: A
    REAL,INTENT(IN)    :: X(1:N)
    REAL,INTENT(INOUT) :: Y(1:N)

    Y(:) = Y(:) + A*X(:)

  end subroutine real_axpy_ij


!######################################################################################



#endif


end Module diag_Gyro_LES_real_3d
