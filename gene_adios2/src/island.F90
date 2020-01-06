#include "redef.h"
#include "intrinsic_sizes.h"
!>Module for an external A_par to simulate static magnetic islands
!> For details see: Plama Phys. Control. Fusion 59 034004 (2017)
Module island
  Use par_mod
  Use file_io, only: get_unit_nr
  Use communications
  USE coordinates
  use geometry
  use blockindex
  use vel_space

  implicit none

  public :: omega_island, width_island, damping_island, j_island, xdamp, delta_xdamp,amp_island
  public :: island_zGauss, width_zGauss, island_ky
  public :: Apar_island, island_contrib
  public :: add_Apar_island, initialize_island, finalize_island
  public :: add_Apar_island_nonlocal
  public ::  set_island_defaults
  private
  integer:: init_status=0
  complex, dimension(:,:), allocatable:: Apar_island
  real:: omega_island, width_island, damping_island, xdamp, delta_xdamp,amp_island
  real, dimension(:), allocatable :: Apar_zGauss
  real :: width_zGauss, island_ky
  logical :: island_contrib, island_zGauss
  integer :: j_island, island_kyind

contains

  subroutine set_island_defaults
    integer :: j

    island_contrib  = .false.
    island_zGauss = .false.
    amp_island = 0.0

    ! default settings for standard island
    j_island= 1
    omega_island = 0.0
    width_island  = 40.0
    damping_island = 2.0
    xdamp = 125/3.0
    delta_xdamp = 125/20.0

    ! default settings for Gauss(z) at k_x = 0
    width_zGauss = 1.0 ! FWHM of Gauss in z units
    island_ky = 0.1 ! k_y where Gauss is to be imposed

  end subroutine set_island_defaults


  subroutine initialize_island
    init_status=0

    if (.not. island_zGauss) then
      IF (.not.x_local) then
        call initialize_prefac_island_nonlocal
      else
        call initialize_prefac_island
      endif
    else
      call initialize_island_zGauss
    endif

    init_status=1
    IF (mype.eq.0) write(*,*) "i. init_status", init_status

  end subroutine initialize_island


  subroutine finalize_island
    call finalize_prefac_island

    init_status=0
    IF (mype.eq.0) write(*,*) "f. init_status", init_status
  end subroutine finalize_island


  subroutine add_Apar_island(emfields,a11det_inv_island)
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz,n_fields),intent(inout):: emfields
    real,dimension(li1:li2,lj1:lj2,lk1:lk2),intent(in) :: a11det_inv_island
    integer:: i,j,k

    if (.not. island_zGauss) then ! standard island implementation
      Do j =lj1,lj2
        IF (j_island.eq.j) then
          DO k = lk1, lk2
            DO i = li1, li2
              emfields(i,j,k,2) = emfields(i,j,k,2) + &
                Apar_island(i,k)*EXP(-imag*omega_island*time)&
                *a11det_inv_island(i,j,k)
            END DO
          END DO
        ENDIF
      Enddo
    else ! Gaussian in z implementation
      do j = lj1, lj2
        if (j .eq. island_kyind) then
          do k = lk1, lk2
            emfields(li1,j,k,2) = emfields(li1,j,k,2) + &
              a11det_inv_island(li1,j,k) * Apar_zGauss(k)
          end do
        end if
      end do
    end if

  end subroutine add_Apar_island

  subroutine add_Apar_island_nonlocal(emfields)
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz,1:n_fields),intent(inout):: emfields
    integer:: i,j,k

    Do j =lj1,lj2
      IF (j_island.eq.j) then
      DO k = lk1, lk2
         DO i = li1, li2
            emfields(i,j,k,2) = emfields(i,j,k,2) + &
              Apar_island(i,k)*EXP(-imag*omega_island*time)
         END DO
      END DO
      ENDIF
    Enddo

  end subroutine add_Apar_island_nonlocal

 subroutine initialize_prefac_island
   integer:: i,k, ierr
   real,dimension(0:nx0-1,lk1:lk2)  :: angle
   integer, dimension(0:nx0-1) :: i_kx
   real  :: amplt, minB_loc, minB

   amplt = 0.0
   angle(:,:) = 0.0
   i_kx(:)  = 0

    !Create an index kx array
    do i= 0,hkx+evenx
      i_kx(i) = i
    enddo
    do i=lkx,nx0-1
      i_kx(i) = i-nx0
    enddo

    !Calculate minimum magnetic field
    minB_loc = minval(geom%Bfield)
    call mpi_allreduce(minB_loc,minB,1,MPI_REAL_TYPE,MPI_MIN, mpi_comm_z, ierr)

    if (init_status.eq.0) then
       amplt = minB*shat*width_island**2/(16.0*q0)
       allocate(Apar_island(0:nx0-1,lk1:lk2))
       Apar_island = cmplx(0.E0,0.E0)
          do k=lk1,lk2
             !Loop for the positive kx modes
             do i=0,hkx
                   angle(i,k)= j_island*nexc*zval(k)/(2*pi) + i_kx(i)
                   if (abs(angle(i,k)).lt.1e-10) then
                    Apar_island(i,k) = amplt*cmplx(1.E0,0.E0)
                    if (mype.eq.0) print*, "angle 0.0 condition satisfied for (i, k) =",i,k
                   else
                    Apar_island(i,k) = amplt*(sin(pi*angle(i,k))/(pi*angle(i,k)))&
                             *EXP(-(angle(i,k)/damping_island)**2)*cmplx(1.E0,0.E0)
                   endif
             end do
             ! Loop for the negative kx modes
             do i=lkx,nx0-1
                   angle(i,k)= j_island*nexc*zval(k)/(2*pi) + i_kx(i)
                   if (abs(angle(i,k)).lt.1e-10) then
                    Apar_island(i,k) = amplt*cmplx(1.E0,0.E0)
                    if (mype.eq.0) print*, "angle 0.0 condition satisfied for (i, k) =",i,k
                   else
                    Apar_island(i,k) = amplt*(sin(pi*angle(i,k))/(pi*angle(i,k)))&
                             *EXP(-(angle(i,k)/damping_island)**2)*cmplx(1.E0,0.E0)
                   endif
             enddo
          end do !end k loop
    endif !init status

  end subroutine initialize_prefac_island

 subroutine initialize_prefac_island_nonlocal
   real, dimension(lk1:lk2) :: Bfield_0, theta_0
   integer:: i,k, ierr
   real  :: amplt, minB_loc, minB, C_y0
   real, dimension(0:nx0-1) :: xval_mod
   Character(Len=8) :: filestat='replace', filepos='rewind'

   xval_mod = xval - x0/rhostar
   amplt = 0.0

   !Calculate minimum magnetic field at the center of the domain
   do k=lk1,lk2
    call calc_center_value(geom%Bfield(:,pj1,k),Bfield_0(k))
    call calc_center_value(geom%THETA_pol(:,k),theta_0(k))
   enddo
    call calc_center_value(C_y, C_y0)

   minB_loc = minval(Bfield_0)
   call mpi_allreduce(minB_loc,minB,1,MPI_REAL_TYPE,MPI_MIN, mpi_comm_z, ierr)

   if (init_status.eq.0) then
       amplt = minB*shat*width_island**2/(8.0*q0*major_R)
       allocate(Apar_island(0:nx0-1,lk1:lk2))
       Apar_island = cmplx(0.E0,0.E0)
          do k=lk1,lk2
             do i=li1,li2
                   Apar_island(i,k) = 0.5*(2.0*amp_island*xval_mod(i)/width_island+1.0)&
                     *amplt*(exp(-imag*ky(j_island)*shat*zval(k)*xval_mod(i)))&
                     *0.25*(1+tanh((xval_mod(i)+xdamp)/delta_xdamp))*(1-tanh((xval_mod(i)-xdamp)/delta_xdamp))

                     !if (magn_geometry.eq.'tracer_efit') Apar_island(i,k)=Apar_island(i,k)&
                     !              *exp(imag*(theta_0(k)-zval(k))*q0*n0_global)
                     ! this is not needed since we use definition of straight field angle instead of theta
                     ! for the poloidal mode number of the island
             enddo !end i loop
          end do !end k loop
   endif !init status

  end subroutine initialize_prefac_island_nonlocal


  subroutine initialize_island_zGauss
    integer :: j, k
    real :: offset

    if (init_status .eq. 0) then
      allocate(Apar_zGauss(lk1:lk2))

      island_kyind = -1111
      do j = lj1, lj2
        if (ABS(ky(j)-island_ky) .le. epsilon(island_ky)) island_kyind = j
      enddo

      ! calculate Gauss amplitude at z = +/-pi for offset
      offset = amp_island * EXP(-(pi/width_zGauss)**2.0)

      do k = lk1, lk2
        Apar_zGauss(k) = amp_island * EXP(-(zval(k)/width_zGauss)**2.0) - offset
      enddo
    endif

  end subroutine initialize_island_zGauss


  subroutine finalize_prefac_island

    if (.not. island_zGauss) then
      deallocate(Apar_island)
    else
      deallocate(Apar_zGauss)
    endif

  end subroutine finalize_prefac_island


end module island
