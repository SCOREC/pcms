#include "redef.h"
#include "intrinsic_sizes.h"

module poloidal_planes_aux
  use communications
  use par_in, only: file_extension, diagdir, n0_global,write_hac,write_h5,write_std
  use par_mod, only: time, itime
  use discretization, only: li1,li2,li0,lj1,lj2,lj0,lk1,lk2,lk0,&
      & ln1,ln2,lbz,ubz,&
      & pmi1gl,pmi2gl,pmx0,&
      & ky0_ind, nky0,nj0,nz0,&
      & n_procs_x,n_procs_y,n_procs_z,n_procs_s,&
      & my_pex,my_pey,my_pez,my_pev,my_pew,my_pespec,mype
  use par_poloidal_planes, only: nz_out,n_cuts,istep_planes, &
                              &  phi_cut, n_cuts, y_BC, &
                              &  mylj0,mylj1,mylj2, res_fact,&
                              &  sign_phi, isnap_field_in, isnap_mom_in,&
                              &  istep_field_in, istep_mom_in
  use geometry, only: C_y, q_prof
  use par_other, only: pi, equil_par_curr,n_fields
  use fourier, only: to_real_y_arb_nl, to_fourier_y_arb_nl, &
                   & to_real_y_arb_lin, to_fourier_y_arb_lin
  use flr_corr, only: n_flr_corrs,correct_for_flr
#ifdef WITHFUTILS
  use futils
  use par_poloidal_planes, only: fidfield_in_h5
#endif

  implicit none

  public :: define_cuts
  public :: compute_pseudo_inv
  public :: gto_real_n, gto_fourier_n

  private

contains

  subroutine define_cuts
    real :: L_tor
    integer :: o

    if (n_cuts.lt.2*nky0-1) then
		print *, "Increase number of planes"
       stop
    endif

    !poloidal  planes, uniformly.
    allocate(phi_cut(0:n_cuts-1))
    L_tor=sign_phi*2.0*pi/n0_global/(n_cuts)
!    L_tor=2*pi/3/16
    do o=0,n_cuts-1
!       phi_cut(o)=L_tor*(o)
       phi_cut(o)=L_tor*(o+1)!-0.5*L_tor
    end do

  end subroutine define_cuts


  subroutine compute_pseudo_inv(in_mat,out_mat)
    !A is the matrix that for each x,z, gives y. size is [ncuts x 2*nky0]
    !inversion with lapack using LSQR. RHS is unit matrix [ncuts x ncuts]
    !ncuts>=2*ky0
    real, dimension(0:n_cuts-1,0:2*nky0*res_fact-1), intent(IN) :: in_mat !matrix to be inverted
    real, dimension(0:2*nky0*res_fact-1,0:n_cuts-1), intent(OUT) :: out_mat !this is inverse

    real, dimension(0:n_cuts-1,0:n_cuts-1) :: unity_mat
    real, dimension(0:n_cuts-1,0:2*nky0*res_fact-1) :: tmp_mat
    integer :: o

    !internal lapack
    integer :: lwork, info
    real, dimension(:), allocatable :: work
    real :: work_tmp

    call set_unity_matrix(unity_mat,n_cuts)

    call dgels('N',n_cuts,2*nky0*res_fact,n_cuts,in_mat,n_cuts,unity_mat,n_cuts,work_tmp,-1,info)

    lwork=int(work_tmp)
    allocate(work(max(1,lwork)))

    tmp_mat=in_mat
    call dgels('N',n_cuts,2*nky0*res_fact,n_cuts,tmp_mat,n_cuts,unity_mat,n_cuts,work,lwork,info)
    deallocate(work)

    do o=0,n_cuts-1
       out_mat(:,o)=unity_mat(0:2*nky0*res_fact-1,o)
    end do

  end subroutine compute_pseudo_inv


  subroutine set_unity_matrix(unit_mat,n)
    integer, intent(IN) :: n
    real, dimension(0:n-1,0:n-1),intent(INOUT) :: unit_mat
    integer :: i

    unit_mat=0.0

    forall(i=0:n-1) unit_mat(i,i) = 1

  end subroutine set_unity_matrix

  !-----------------------------------------------------------------------------------
  ! if we dont parallelize in y, this can be largely simplified

   Subroutine transpose_cmplx(n1, n2, in_matrix, ssoff, transposed_matrix, ddoff)

    Integer, Intent(in):: n1, n2, ssoff, ddoff
    Complex, Intent(in):: in_matrix(0:, 0:)
    Complex, Intent(out):: transposed_matrix(0:, 0:)

    Complex, Dimension(n2/n_procs_y, n1/n_procs_y, 0:n_procs_y-1) :: sbuf, rbufc
    Integer:: n1l, n2l, i1
    Integer:: pp, ierr

    n1l = n1/n_procs_y
    n2l = n2/n_procs_y

    if (n_procs_y.eq.1) then
       transposed_matrix(ddoff:ddoff+n2-1,:)=Transpose(in_matrix(ssoff:ssoff+n1-1,:))
    else
       Do pp = 0, n_procs_y-1
          i1 = ssoff+pp*n1l
          sbuf(:,:,pp) = Transpose(in_matrix(i1:i1+n1l-1, :))
       Enddo
       Call mpi_alltoall(&
            sbuf, n1l*n2l, MPI_COMPLEX_TYPE,&
            rbufc, n1l*n2l, MPI_COMPLEX_TYPE,&
            mpi_comm_y, ierr)
       Do pp = 0, n_procs_y-1
          i1 = ddoff + pp*n2l
          transposed_matrix(i1:i1+n2l-1,:) = rbufc(:,:,pp)
       Enddo
    endif

  End Subroutine transpose_cmplx


  Subroutine transpose_real(n1, n2, in_matrix, ssoff, transposed_matrix, ddoff)
    Integer, Intent(in):: n1, n2, ssoff, ddoff
    real, Intent(in):: in_matrix(0:, 0:)
    real, Intent(out):: transposed_matrix(0:, 0:)
    real, Dimension(n2/n_procs_y, n1/n_procs_y, 0:n_procs_y-1) :: sbuf, rbufc
    Integer:: n1l, n2l, i1
    Integer:: pp, ierr

    n1l = n1/n_procs_y
    n2l = n2/n_procs_y

    if (n_procs_y.eq.1) then
       transposed_matrix(ddoff:ddoff+n2-1,:)=Transpose(in_matrix(ssoff:ssoff+n1-1,:))
    else
        Do pp = 0, n_procs_y-1
          i1 = ssoff+pp*n1l
          sbuf(:,:,pp) = Transpose(in_matrix(i1:i1+n1l-1, :))
        Enddo
        Call mpi_alltoall(&
              sbuf, n1l*n2l, MPI_REAL_TYPE,&
              rbufc, n1l*n2l, MPI_REAL_TYPE,&
              mpi_comm_y, ierr)
         Do pp = 0, n_procs_y-1
          i1 = ddoff + pp*n2l
          transposed_matrix(i1:i1+n2l-1,:) = rbufc(:,:,pp)
       Enddo
    endif

  End Subroutine transpose_real


  !-----------------------------------------------------------------------------------
  !fft related subroutines XGC grid
   subroutine gto_real_n(field_in,field_out,myli0)
     integer,intent(in) :: myli0
     complex, dimension(0:myli0-1,lj1:lj2), intent(IN) :: field_in
     real, dimension(0:myli0-1,mylj1:mylj2), intent(OUT) :: field_out

     complex, dimension(0:nj0-1,0:myli0/n_procs_y-1):: tmp_cmplx
     real, dimension(0:2*nj0*res_fact-1,0:myli0/n_procs_y-1) :: tmp_re

     call transpose_cmplx(myli0,lj0,field_in,0,tmp_cmplx,0)

     if (nj0.eq.1) then
        call to_real_y_arb_lin(tmp_cmplx,tmp_re)
     else
        call to_real_y_arb_nl(tmp_cmplx,tmp_re)
     endif

     call transpose_real(mylj0,myli0,tmp_re,0,field_out,0)

   end subroutine gto_real_n


   subroutine gto_fourier_n(field_in,field_out,myli0)
     integer,intent(in) :: myli0
     real, dimension(0:myli0-1,mylj1:mylj2) :: field_in
     complex,dimension(0:myli0-1,lj1:lj2) :: field_out

     complex, dimension(0:nj0-1,0:myli0/n_procs_y-1):: tmp_cmplx
     real, dimension(0:mylj0-1,0:myli0/n_procs_y-1) :: tmp_re

        call transpose_real(myli0,mylj0,field_in,0,tmp_re,0)

        if (nj0.eq.1) then
           call to_fourier_y_arb_lin(tmp_re,tmp_cmplx)
        else
           call to_fourier_y_arb_nl(tmp_re,tmp_cmplx)
        endif

        call transpose_cmplx(lj0,myli0,tmp_cmplx,0,field_out,0)

   end subroutine gto_fourier_n

end module poloidal_planes_aux
