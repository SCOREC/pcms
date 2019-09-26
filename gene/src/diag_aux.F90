#include "switches.h"
#include "redef.h"
#include "intrinsic_sizes.h"
!> This auxiliary module contains arrays and subroutines
!! that are shared by various diagnostics (modules)
Module diagnostics_auxiliary
  use aux_fields
  use geometry, only: flux_geomfac, geom
  use par_mod
  Use vel_space, only: mat_00, mat_10, gyro_op_wrapper
  use x_derivatives, only: x_deriv_exc
  use axpy

  Implicit None

  PUBLIC :: mem_est_diag_aux, initialize_mats_diag_aux, &
       &finalize_mats_diag_aux, calc_dfields_dy

  PUBLIC :: mat_20, mat_01, nmat_00, nmat_20, nmat_01, mat_30, mat_11
  PUBLIC :: pTj_qjB0
  PUBLIC :: calc_moments_trap_pass

  PRIVATE

  REAL, DIMENSION(:,:,:,:,:), ALLOCATABLE:: mat_20, mat_01, nmat_00, nmat_20, nmat_01
  REAL, DIMENSION(:,:,:,:,:,:), ALLOCATABLE:: mat_30, mat_11
  REAL, DIMENSION(:,:), ALLOCATABLE :: pTj_qjB0

CONTAINS

!!!******************************************************************!!!
!!!******************************************************************!!!
  !>Give an estimate of the memory requirements of this module
  Real Function mem_est_diag_aux(mem_req_in)
    real:: mem_req_in
    real:: mem_loc=0

    mem_loc = 0.0

    !module variables
    !mat_20, mat_01
    mem_loc=mem_loc+2*SIZE_OF_REAL_MB*pi0*pj0*lk0*ll0*lm0
    !mat_30, mat_11
    mem_loc=mem_loc+2*SIZE_OF_REAL_MB*pi0*pj0*lklmn0

    mem_est_diag_aux=mem_req_in+mem_loc
  End Function mem_est_diag_aux

!!!******************************************************************!!!
  SUBROUTINE initialize_mats_diag_aux
    INTEGER :: k, l, m, n

    ALLOCATE(mat_20(pi1:pi2, pj1:pj2, lk1:lk2, ll1:ll2, lm1:lm2))
    ALLOCATE(mat_30(pi1:pi2, pj1:pj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
    ALLOCATE(mat_01(pi1:pi2, pj1:pj2, lk1:lk2, ll1:ll2, lm1:lm2))
    ALLOCATE(mat_11(pi1:pi2, pj1:pj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
    IF (n_fields .GT. 2 .or. all_moms) THEN
       ALLOCATE(nmat_00(pi1:pi2,pj1:pj2,lk1:lk2,ll1:ll2,lm1:lm2))
       ALLOCATE(nmat_20(pi1:pi2,pj1:pj2,lk1:lk2,ll1:ll2,lm1:lm2))
       ALLOCATE(nmat_01(pi1:pi2,pj1:pj2,lk1:lk2,ll1:ll2,lm1:lm2))

       ALLOCATE(pTj_qjB0(lk1:lk2,ln1:ln2))
       DO n = ln1, ln2
          DO k = lk1, lk2
             pTj_qjB0(k,n) = spec(n)%temp / (spec(n)%charge * geom%Bfield(pi1,pj1,k))
          END DO
       END DO
    ENDIF

    DO m=lm1,lm2
       DO l=ll1,ll2
          DO k=lk1,lk2
             ! temperatures
             mat_20(:,:,k,l,m) = vp(l)*vp(l) * mat_00(:,:,k,l,m)
             mat_01(:,:,k,l,m) = mu(m)*geom%Bfield(pi1:pi2,:,k) * mat_00(:,:,k,l,m)
             DO n=ln1,ln2
                mat_30(:,:,k,l,m,n) = vp(l)*vp(l) * mat_10(:,:,k,l,m,n) !q_||,||
                mat_11(:,:,k,l,m,n) = mu(m)*geom%Bfield(pi1:pi2,:,k) *&
                     mat_10(:,:,k,l,m,n) !q_||,perp
             END DO

             IF (n_fields .GT. 2 .or. all_moms) THEN
               nmat_00(:,:,k,l,m) = & ! densI1
                 mu(m) * geom%Bfield(pi1:pi2,:,k) * mat_00(:,:,k,l,m)
               nmat_20(:,:,k,l,m) = & ! TparI1
                 mu(m) * geom%Bfield(pi1:pi2,:,k) * mat_20(:,:,k,l,m)
               nmat_01(:,:,k,l,m) = & ! TperpI1
                 mu(m) * geom%Bfield(pi1:pi2,:,k) * mat_01(:,:,k,l,m)
             END IF
          END DO
       END DO
    END DO

  END SUBROUTINE initialize_mats_diag_aux

!!!******************************************************************!!!

  SUBROUTINE finalize_mats_diag_aux
    DEALLOCATE(mat_20,mat_30,mat_01,mat_11)
    IF (n_fields .GT. 2 .or. all_moms) DEALLOCATE(nmat_00,nmat_20,nmat_01,pTj_qjB0)
  END SUBROUTINE finalize_mats_diag_aux


!!!******************************************************************!!!
!!!***************** additional subroutines  **** *******************!!!
!!!******************************************************************!!!


  !> Calculates the y-derivative of the electromagnetic fields
  !! for computing the ExB drift velocity that enters the fluxes
  Subroutine calc_dfields_dy(dfields_dy,k,o)
    complex, dimension(li1:li2,lj1:lj2),intent(out):: dfields_dy
    integer,intent(in):: o, k
    complex, dimension(lbi:ubi,lj1:lj2):: field_wb
    integer:: j

      if (yx_order) then
         if (y_local) then
            do j=lj1,lj2
               dfields_dy(:,j)=imag*ki(:)*emfields(:,j,k,o)*flux_geomfac(pi1,pj1,k)
            end do
         else
            !for .not. y_local: compute the y derivative, then multiply geom_fac
            field_wb(li1:li2,:) = emfields(:,:,k,o)
            call x_deriv_exc(field_wb,dfields_dy)
            do j=lj1,lj2
               dfields_dy(:,j)=dfields_dy(:,j)*flux_geomfac(:,pj1,k)
            enddo
         end if
      else
         do j=lj1,lj2
            dfields_dy(:,j)=imag*kj(j)*emfields(:,j,k,o)
            if (x_local) then
               dfields_dy(:,j) = dfields_dy(:,j)*flux_geomfac(pi1,pj1,k)
            else
               dfields_dy(:,j) = dfields_dy(:,j)*flux_geomfac(:,pj1,k)
            endif
         end do
      end if

  End Subroutine calc_dfields_dy


  SUBROUTINE calc_moments_trap_pass(n_mom,p_dist,p_mom,n,diag_trap_levels,diag_Blev)
    Implicit none
    INTEGER, INTENT(in):: n_mom
    !LOGICAL, INTENT(in):: withbounds
    integer, intent(in) :: n, diag_trap_levels
    real,dimension(0:9),intent(in) :: diag_Blev
    COMPLEX, DIMENSION(:,:,:,:,:),INTENT(in):: p_dist
    COMPLEX, DIMENSION(li1:li2,lj1:lj2,lk1:lk2,1:n_mom,0:diag_trap_levels),INTENT(out):: p_mom
    !INTEGER, DIMENSION(1:n_mom), optional :: p_gy_op

    ! Local variables
    COMPLEX, DIMENSION(li1:li2,lj1:lj2) :: g_av
    real ::energy
    INTEGER:: k,l,m,i_mom,joff,koff,loff,moff,noff,trap_level

    PERFON('calc_moments')

    p_mom = cmplx(0.0,0.0)

    ! the passed array p_dist is running locally from 1 to the end
    ! dependent on if we have the boundaries included, the calculation of
    ! the inner point indices is different
    ! the first two dimensions are always without boundaries, that means
    ! running over 1:li0 1:lj0
    !
    joff = 1-lj1
    koff = 1-lk1+nzb
    loff = 1-ll1+nvb
    moff = 1-lm1+nwb
    noff=1-ln1


    DO m=lm1,lm2
       DO l=ll1,ll2
          DO k=lk1,lk2

             IF (xy_local) THEN
                energy = vp(l)*vp(l)+mu(m)*geom%Bfield(pi1,pj1,k)

                CALL gyro_op_wrapper(p_dist(:,:,k+koff,l+loff,m+moff),g_av,k,m,n,0)

                do trap_level=0,diag_trap_levels
                   if ((energy.Lt.mu(m)*diag_Blev(trap_level)).or.(trap_level.eq.diag_trap_levels)) then
                      CALL axpy_ij(lij0,mat_00(pi1,pj1,k,l,m),g_av,p_mom(:,:,k,1,trap_level))
                      CALL axpy_ij(lij0,mat_20(pi1,pj1,k,l,m),g_av,p_mom(:,:,k,2,trap_level))
                      CALL axpy_ij(lij0,mat_01(pi1,pj1,k,l,m),g_av,p_mom(:,:,k,3,trap_level))
                      CALL axpy_ij(lij0,mat_30(pi1,pj1,k,l,m,n),g_av,p_mom(:,:,k,4,trap_level))
                      CALL axpy_ij(lij0,mat_11(pi1,pj1,k,l,m,n),g_av,p_mom(:,:,k,5,trap_level))
                      CALL axpy_ij(lij0,mat_10(pi1,pj1,k,l,m,n),g_av,p_mom(:,:,k,6,trap_level))
                      exit
                   endif
                enddo

                IF (n_fields .GT. 2 .or. all_moms) THEN
                   CALL gyro_op_wrapper(p_dist(:,:,k+koff,l+loff,m+moff),g_av,k,m,n,1)

                   do trap_level=0,diag_trap_levels
                      if ((energy.Lt.mu(m)*diag_Blev(trap_level)).or.(trap_level.eq.diag_trap_levels)) then
                         do i_mom = 7,8
                            CALL axpy_ij(lij0,nmat_00(pi1,pj1,k,l,m),g_av,p_mom(:,:,k,i_mom,trap_level))
                            CALL axpy_ij(lij0,nmat_20(pi1,pj1,k,l,m),g_av,p_mom(:,:,k,i_mom,trap_level))
                            CALL axpy_ij(lij0,nmat_01(pi1,pj1,k,l,m),g_av,p_mom(:,:,k,i_mom,trap_level))
                            exit
                         enddo
                      endif
                   enddo
                 ENDIF

             ELSE
                ! in the following call, we assume that nyb = 0. If that is not the
                ! case, one has to rewrite the gyro_average method
                p_mom(:,:,k,1,0:) = p_mom(:,:,k,1,0:) + diag_a_mom(p_dist,mat_00(:,:,k,l,m),&
                                  & k,l,m,n,joff,koff,loff,moff, &
                                  & diag_trap_levels,diag_Blev)
                p_mom(:,:,k,2,0:) = p_mom(:,:,k,2,0:) + diag_a_mom(p_dist,mat_20(:,:,k,l,m),&
                                 & k,l,m,n,joff,koff,loff,moff, &
                                 & diag_trap_levels,diag_Blev)
                p_mom(:,:,k,3,0:) = p_mom(:,:,k,3,0:) + diag_a_mom(p_dist,mat_01(:,:,k,l,m),&
                                 & k,l,m,n,joff,koff,loff,moff, &
                                 & diag_trap_levels,diag_Blev)
                p_mom(:,:,k,4,0:) = p_mom(:,:,k,4,0:) + diag_a_mom(p_dist,mat_30(:,:,k,l,m,n),&
                                 & k,l,m,n,joff,koff,loff,moff, &
                                 & diag_trap_levels,diag_Blev)
                p_mom(:,:,k,5,0:) = p_mom(:,:,k,5,0:) + diag_a_mom(p_dist,mat_11(:,:,k,l,m,n),&
                                 & k,l,m,n,joff,koff,loff,moff, &
                                 & diag_trap_levels,diag_Blev)
                p_mom(:,:,k,6,0:) = p_mom(:,:,k,6,0:) + diag_a_mom(p_dist,mat_10(:,:,k,l,m,n),&
                                 & k,l,m,n,joff,koff,loff,moff, &
                                 & diag_trap_levels,diag_Blev)
             END IF
          ENDDO
       ENDDO
    ENDDO

    ! In p_mom we have now in each process the moments over
    ! the process local velocity space
    ! Note: Sum over all MPI processes needs to be done

    PERFOFF

  END SUBROUTINE calc_moments_trap_pass


  function diag_a_mom(p_dist,p_mat,k,l,m,n,joff,koff,loff,moff,diag_trap_levels,diag_Blev) result(p_mom)
    integer, intent(in) :: diag_trap_levels
    COMPLEX, DIMENSION(:,:,:,:,:),INTENT(in):: p_dist
    COMPLEX, DIMENSION(li1:li2,lj1:lj2,0:diag_trap_levels) :: p_mom
    REAL, DIMENSION(pi1:pi2,pj1:pj2),INTENT(in):: p_mat
    real, dimension(0:9),intent(in) :: diag_Blev
    integer, intent(in) :: k,l,m,n
    integer, intent(in) :: joff,koff,loff,moff

    !local variables
    COMPLEX, DIMENSION(li1:li2,lj1:lj2) :: dist_mat
    COMPLEX, DIMENSION(li1:li2,lj1:lj2) :: g_av
    integer :: j

#ifndef GY_WITH_MAT

    CALL gyro_op_wrapper(p_dist(:,:,k+koff,l+loff,m+moff),g_av,k,m,n,0)

    do j=lj1,lj2
       g_av(:,j)= g_av(:,j)* p_mat(:,pj1)
    enddo

#else

    do j=lj1,lj2
       dist_mat(:,j)= p_dist(:,j+joff,k+koff,l+loff,m+moff) &
                  & * p_mat(:,pj1)
    enddo

    CALL gyro_op_wrapper(dist_mat,g_av,k,m,n,0)

#endif

    p_mom = split_trap_pass(g_av,k,l,m,diag_trap_levels,diag_Blev)

  end function diag_a_mom


  function split_trap_pass(g_av,k,l,m,diag_trap_levels,diag_Blev) result(p_mom)
    integer, intent(in) :: diag_trap_levels
    COMPLEX, DIMENSION(li1:li2,lj1:lj2),intent(in) :: g_av
    real, dimension(0:9),intent(in) :: diag_Blev
    integer, intent(in) :: k,l,m
    COMPLEX, DIMENSION(li1:li2,lj1:lj2,0:diag_trap_levels):: p_mom

    ! Local variables
    real :: energy
    integer :: trap_level
    integer :: i

    p_mom(:,:,:)=cmplx(0.0,0.0)

    do trap_level=0,diag_trap_levels
       do i=li1,li2
          !\todo: check the following energy x-dependence (temperature missing?)
          energy = vp(l)*vp(l)+mu(m)*geom%Bfield(i,pj1,k)
          if ((energy.Lt.mu(m)*diag_Blev(trap_level)).or.(trap_level.eq.diag_trap_levels)) then
             p_mom(i,:,trap_level) =  g_av(i,:)
             !to have this exit doing the work we need to swith the do loops, but this will jump in memory
             !exit
          endif
       end do
    end do

  end function split_trap_pass


end Module diagnostics_auxiliary
