#include "redef.h"
#include "intrinsic_sizes.h"

module df_arakawa_nonlinear_term_mod
  use df_nonlinear_term_mod

  use par_mod, only: imag, pi, ve_max
  use coordinates, only: kjmin
  use communications, only: MPI_COMPLEX_TYPE, MPI_REAL_TYPE, mpi_comm_y,&
       &my_2reals_max_to_all, my_barrier, mpi_comm_xy, reduce_per_thread,&
       &threadlocal_mpi_comm_y
  use discretization
  use blockindex, only: sk,sl,sn
  use fourier
  use prefactors
  USE x_derivatives, only: x_deriv_exc, x_deriv_exc_transp
  use geometry, only: geom, C_xy
  use mpi
  implicit none

  type, public, extends(df_nonlinear_term_t) :: df_arakawa_nonlinear_term_t
   contains
     procedure :: calc => calc_arakawa_nonlinearity_df
     procedure :: getType => getThisType
  end type df_arakawa_nonlinear_term_t


contains

  function getThisType(this)
    class(df_arakawa_nonlinear_term_t) :: this
    character(len=MAX_TYPENAME_LENGTH) :: getThisType

    getThisType = "df_arakawa_nonlinear_term_t"
  end function getThisType

  !>Computes the nonlinearity for the (x) nonlocal version
  !!\todo Get rid of transposition using striding? make some arrays for arakawa allocatable
  !!\todo Check whether the nonlinearity prefactor is important for the CFL criterion
  Subroutine calc_arakawa_nonlinearity_df(this,gy_chi,g_block,vexy,dgdxy,localrhs,first,lb1,lb2)
    class(df_arakawa_nonlinear_term_t) :: this
    Complex, Dimension(li1:li2,lj1:lj2,1:*),Intent(inout) :: gy_chi
    complex, dimension(lbi:ubi,lj1:lj2,1:*),intent(in) :: g_block
    Complex, Dimension(li1:li2,lj1:lj2,2,1:*),Intent(inout) :: vexy, dgdxy
    Complex, Dimension(li1:li2,lj1:lj2,1:*),Intent(inout) :: localrhs
    Logical, Intent(in):: first
    Integer, Intent(in):: lb1,lb2

    ! Local variables
    Complex, Dimension(li1:li2,lj1:lj2,1:this%lbg0) :: nonlin
    Complex, Dimension(0:nj0-1, 0:li0/n_procs_y-1) :: nonlin3,tmp_arr
    complex, dimension(li1:li2, lj1:lj2):: g2d
    Real, Dimension(0:ly0da-1, 0:li0/n_procs_y-1,2) ::  dgdxy_re
    Real, Dimension(0:ly0da-1, 0:li0/n_procs_y-1,2,1) ::  vexy_re !last dimension needed for multiply_max
    Real, Dimension(0:ly0da-1, 0:li0/n_procs_y-1) ::  gy_chi_re, g2d_re
    Real, Dimension(0:ly0da-1, 0:li0/n_procs_y-1) :: nonlin_re
    Real, Dimension(0:ly0da-1,-nib:li0-1+nib):: nonlin2_wb
    Real, Dimension(1:2)  ::  vexy_max_loc
    Integer :: i,j,klmn,ub
#ifdef WITHOMP
    integer :: my_thread, omp_get_thread_num
    my_thread = omp_get_thread_num()
#endif

    ub=li0/n_procs_y-1
    do klmn=1,this%lbg0
       ! anti-aliasing and y fourier transform
       PERFON_I('deal_FFT')
       Call gto_real_nonlinearity_df(this,vexy(:,:,:,klmn),vexy_re(:,:,:,1),2)
       Call gto_real_nonlinearity_df(this,dgdxy(:,:,:,klmn),dgdxy_re,2)

       g2d = g_block(li1:li2,:,klmn)

       !we need the distribution function and the potential in real space
       Call gto_real_nonlinearity_df(this,g2d,g2d_re,1)
       Call gto_real_nonlinearity_df(this,gy_chi(:,:,klmn),gy_chi_re,1)

       PERFOFF_I

       ! get max ExB velocity for timestep estimate
       if (first) then
          call this%prefactor%multiply_max(vexy_max_loc,vexy_re,lb1+klmn-1,lb1+klmn-1)
          ve_max(1)=max(ve_max(1),vexy_max_loc(1)) !ve_x
          ve_max(2)=max(ve_max(2),vexy_max_loc(2)) !ve_y
       end if

       !second Arakawa term
       nonlin2_wb(:,0:ub) = g2d_re*vexy_re(:,:,2,1) - gy_chi_re*dgdxy_re(:,:,2)
       call x_deriv_exc_transp(nonlin2_wb,nonlin_re,li0/n_procs_y)

       !add first Arakawa term (the 'standard' nonlinear term) and transform back to fourier space
       nonlin_re = nonlin_re -vexy_re(:,:,1,1)*dgdxy_re(:,:,2) + vexy_re(:,:,2,1)*dgdxy_re(:,:,1)
       Call to_fourier_y(nonlin_re,tmp_arr)

       !third Arakawa term
       nonlin_re = gy_chi_re(:,:)*dgdxy_re(:,:,1) - g2d_re(:,:)*vexy_re(:,:,1,1)
       ! fourier transform nonlin_re3 back to ky space and compute its y derivative
       Call to_fourier_y(nonlin_re,nonlin3)
       do i=0,ub
          do j=0,nj0-1
             nonlin3(j,i)=imag*kjmin*j*nonlin3(j,i)
          end do
       end do

       ! add all terms and divide by three to calculate the average
       tmp_arr=1/3.*tmp_arr+1/3.*nonlin3

       ! Transpose and remove zeros in y for dealiasing

       Call transpose_cmplx(nj0, li0, tmp_arr, 0, nonlin(:,:,klmn), 0)
    end do


    ! multiply with the prefactor and write to localrhs
    call this%prefactor%multiply_with(nonlin,localrhs,this%lbg0,lb1,lb2)

  End Subroutine calc_arakawa_nonlinearity_df

  !> Go to real space
  !! This routine is equal to the one of df_nonlinear_term, with the
  !! small exception, that here it is clear, that we do not have
  !! explicit x-dealiasing. This saves us one copy operation.
  subroutine gto_real(this,inarr,outarr,howmany)
    class(df_nonlinear_term_t) :: this
    integer, intent(IN) :: howmany
    Complex, Dimension(li1:li2,lj1:lj2,1:howmany),Intent(in) :: inarr
    Real, Dimension(0:ly0da-1, 0:li0/n_procs_y-1,1:howmany), target, Intent(out) :: outarr

    ! Local variables
    !Complex, Dimension(li1:li2,lj1:lj2,1:howmany) :: tmp_arr1
    Complex, Dimension(0:nj0-1, 0:li0/n_procs_y-1) :: tmp_arr2
    real, dimension(:,:), pointer :: p_out
    Integer:: klmn, i_block,number_of_lbg0_blocks
#ifdef WITHOMP
    integer :: my_thread, omp_get_thread_num

    my_thread = omp_get_thread_num()
#endif
    number_of_lbg0_blocks = howmany/this%lbg0
    if (number_of_lbg0_blocks.eq.0) then
       ! Transpose x-y
       Call transpose_cmplx(li0, nj0, inarr(:,:,1), 0,tmp_arr2 , 0)
       ! Fourier transfrom in y (include dealiasing step)
       Call to_real_y(tmp_arr2,outarr(:,:,1))
    else
       do i_block=1,number_of_lbg0_blocks

          ! Transpose x-y
          do klmn=1,this%lbg0
             Call transpose_cmplx(li0, nj0, inarr(:,:,(i_block-1)*this%lbg0+klmn), 0,tmp_arr2 , 0)

             ! Fourier transfrom in y (include dealiasing step)
             p_out => outarr(:,:,(i_block-1)*this%lbg0+klmn)
             Call to_real_y(tmp_arr2,p_out) !outarr(:,:,(i_block-1)*this%lbg0+klmn))
          end do
       end do
    end if

  end subroutine gto_real

end module df_arakawa_nonlinear_term_mod
