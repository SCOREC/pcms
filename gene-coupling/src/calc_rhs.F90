#include "redef.h"
#include "intrinsic_sizes.h"
#include "switches.h"

!>Computes the full right hand side of the gyrokinetic Vlasov equation
!!from a given modified distribution function g_1
!!
!!This module calls the primary implementations of the RHS of the gyrokinetic
!!equation; it is also used by the derived methods, e.g. to compute an explicit
!!representation of the linear operator or for the matrix-free methods used in
!!PETSc/SLEPc
Module calc_rhs
  Use aux_fields
  Use antenna, only: antenna_type, add_dApar_dt_antenna
  Use par_mod
  Use prefactors
  use boundaries
  !USE nonlinearity
  use parallel_nonlin
  use f0_term
  use dchidxy_terms
  use dgdxy_terms
  use dfdxy_terms
  use dchidz_term
  use dfdzv_terms
  use dzv_terms
  Use collisions
  Use aux_fields
  use gyro_average_df_mod, only: gyro_average_df_wb
  use blockindex
  use sources_mod
  use RK_coefficients
  use axpy
  use numerical_damping
  use communications, only: reduce_per_thread, threadlocal_mpi_comm_y, my_secnds
  use nonlinear_term_mod
  use df_nonlinear_term_mod
  use df_arakawa_nonlinear_term_mod
  use ff_nonlinear_term_mod
  use ff_nonlinear_term_h_mod
#ifdef WITH_C_NONLIN
  use df_nonlinear_c_mod
#endif
#ifdef WITH_CUDA_NONLIN
  use df_nonlinear_cuda_mod
#endif
#ifdef WITH_SPLIT_NONLIN
  use df_split_nonlinear_term_mod
#endif
#ifdef WITHOMP_NONLIN
  use df_omp_nonlinear_term_mod
#endif
#ifdef WITHOMPXY_NONLIN
  use df_ompxy_nonlinear_term_mod
#endif
#ifdef WITH_MIC_NONLIN
  use df_mic_nonlinear_term_mod
#endif
  ! x_derivative module is necessary to determine correct ranges
  ! for chi_block computations
  use x_derivatives, only: radscheme
  USE all_rhs_terms, only: this_nonlinear_term
#ifdef COUPLE_XGC
  use diagnostics, only: exec_all_diags
#endif
  use compute_f, only: f_, h_
  Use antenna, only: antenna_type, add_dApar_dt_antenna, antenna_contrib

  implicit None

  Public:: initialize_CalFullRhs, CalFullRhs, calc_rhs_only, &
       & finalize_CalFullRhs, mem_est_calc_rhs

  public:: rhs_nl, rhs_f0, f0_heat_src_active

  !for operator splitting of collisions and vlasov part in the time scheme
  public:: rhs_only_coll, rhs_with_coll

  !the following pointers are also used by get_energy_terms
  public :: bar_emfields
  public :: ptr_barchi, ptr_dbarchidxy, ptr_dgdxy, chi_block
  public :: this_nonlinear_term

  private

  integer :: init_status = 0

  logical:: rhs_nl=.false., rhs_f0=.false.
  logical:: rhs_only_coll=.false.
  logical:: rhs_with_coll=.true.

  logical:: with_dfdperp, with_g_update

  complex, dimension(:,:,:,:,:,:),allocatable, target :: bar_emfields
  complex, dimension(:,:,:),allocatable, target :: p_barchi
  complex, dimension(:,:,:,:),allocatable, target :: p_dbarchidxy, p_dgdxy
  complex, dimension(:,:,:,:),allocatable,target :: rhs_big
  complex, dimension(:,:,:),allocatable,target:: rhs_small
  complex, dimension(:,:,:), allocatable, target :: g_block_storage
  complex, dimension(:,:,:), allocatable :: chi_block

  !$DIR ATTRIBUTES ALIGN:64 :: rhs_small, chi_block

  !pointers are employed in this module in order to allow for 'optional'
  !arguments in subroutine calls
  !complex, dimension(:,:,:,:,:,:),pointer :: ptr_bar_emfields
  complex, dimension(:,:,:),pointer :: ptr_barchi
  complex, dimension(:,:,:,:),pointer :: ptr_dbarchidxy, ptr_dgdxy

Contains

  function mem_est_calc_rhs(mem_req_in)
    real:: mem_req_in, mem_est_calc_rhs
    real:: mem_loc
    logical :: to_deallocate

    mem_loc=0.
    mem_loc=mem_est_klmn_conv_arrs(mem_req_in)
    mem_loc=mem_est_dfdzv(mem_loc)
    if((hyp_x.gt.0).or.(hyp_y.gt.0)) &
         & mem_loc=mem_est_dfdxy(mem_loc)
    mem_loc=mem_est_dchidxy(mem_loc)
    mem_loc=mem_est_dchidz(mem_loc)
    mem_loc=mem_est_dgdxy(mem_loc)
    !    mem_loc=mem_est_coll(mem_loc)
    !if (nonlinear) mem_loc = mem_est_nonlinearity(mem_loc)
    if (nonlinear) then
       if (.not.associated(this_nonlinear_term)) then
          if (xy_local) then
             !flux tube: only yx_order = F
             if (.not.nonlin_h) then
                allocate(ff_nonlinear_term_t::this_nonlinear_term)
             else
                allocate(ff_nonlinear_term_h_t::this_nonlinear_term)
             endif
             to_deallocate = .true.
          else !.not.xy_local
                if (arakawa) then
                   allocate(df_arakawa_nonlinear_term_t::this_nonlinear_term)
                else
#ifdef WITH_C_NONLIN
                   allocate(df_nonlinear_c_t::this_nonlinear_term)
#elif defined(WITH_SPLIT_NONLIN)
                   allocate(df_split_nonlinear_term_t::this_nonlinear_term)
#elif defined(WITH_CUDA_NONLIN)
                   allocate(df_nonlinear_cuda_t::this_nonlinear_term)
#elif defined(WITH_MIC_NONLIN)
                   allocate(df_mic_nonlinear_term_t::this_nonlinear_term)
#elif defined(WITHOMP_NONLIN)
                   allocate(df_omp_nonlinear_term_t::this_nonlinear_term)
#elif defined(WITHOMPXY_NONLIN)
                   allocate(df_ompxy_nonlinear_term_t::this_nonlinear_term)
#else
                   allocate(df_nonlinear_term_t::this_nonlinear_term)
#endif
                end if
                to_deallocate = .true.
             !end if
          end if
          call this_nonlinear_term%construct(equil_par_curr)
       else
          ! already associated
          to_deallocate = .false.
       end if

       if (associated(this_nonlinear_term)) then
          mem_loc = this_nonlinear_term%mem_est(mem_loc)
          if (to_deallocate) then
             call this_nonlinear_term%destruct()
             deallocate(this_nonlinear_term)
          end if
       else
          WRITE(*,"(A)") "this_nonlinear_term is not associated in calc_rhs. This should not happen!"
          STOP
       end if
    end if

    !p_barchi
    if (nonlinear) mem_loc=mem_loc+&
         &lij0*lklmn0/nblocks*SIZE_OF_COMPLEX/(1024.)**2

    if (parallel_nl) mem_loc = mem_est_parallel_nonlin(mem_loc)

    if (.not.xy_local) then
       mem_loc = mem_est_f1_sources(mem_loc)
       !bar_emfields
       mem_loc = mem_loc + 1.*lij0*lz0*lm0*ln0*n_fields*SIZE_OF_COMPLEX/(1024.)**2
       !p_dbarchidxy, p_dgdxy
       if (nonlinear) mem_loc=mem_loc+4.*lij0*lklmn0/nblocks*SIZE_OF_COMPLEX/(1024.)**2
    endif

    mem_loc = mem_est_f0_term(mem_loc)

    mem_est_calc_rhs=mem_req_in+mem_loc

  end function mem_est_calc_rhs


  !>Wrapper for the initialization routines of the different implementations
  subroutine initialize_CalFullRhs(for_petsc)
    implicit none
    logical, optional:: for_petsc
    logical:: dfdzv_replace_rhs
#ifdef WITHOMP
    integer :: omp_get_thread_num
#endif
    integer :: my_thread

    character(len=MAX_TYPENAME_LENGTH) :: type_of_nonlin

    rhs_with_coll = ((collision_op.ne.'none').and.(.not.coll_split))

    if(present(for_petsc)) then
       with_g_update=.false.
    else
       with_g_update=low_mem
    endif

    if (rhs_nl) then
       my_thread = 0
       if (.not.associated(this_nonlinear_term)) then
          if (xy_local) then
             !flux-tub: only yx_order = F
             if (.not.nonlin_h) then
                allocate(ff_nonlinear_term_t::this_nonlinear_term)
             else
                allocate(ff_nonlinear_term_h_t::this_nonlinear_term)
             endif
          else
             IF (arakawa) THEN
                allocate(df_arakawa_nonlinear_term_t::this_nonlinear_term)
             else
#ifdef WITH_C_NONLIN
                allocate(df_nonlinear_c_t::this_nonlinear_term)
#elif defined(WITH_SPLIT_NONLIN)
                allocate(df_split_nonlinear_term_t::this_nonlinear_term)
#elif defined(WITH_CUDA_NONLIN)
                allocate(df_nonlinear_cuda_t::this_nonlinear_term)
#elif defined(WITH_MIC_NONLIN)
                allocate(df_mic_nonlinear_term_t::this_nonlinear_term)
#elif defined(WITHOMP_NONLIN)
                allocate(df_omp_nonlinear_term_t::this_nonlinear_term)
#elif defined(WITHOMPXY_NONLIN)
                allocate(df_ompxy_nonlinear_term_t::this_nonlinear_term)
#else
                allocate(df_nonlinear_term_t::this_nonlinear_term)
#endif
             END IF
          end if
          call this_nonlinear_term%construct(equil_par_curr)
       end if

       select type (tmp=>this_nonlinear_term)
#ifdef WITH_SPLIT_NONLIN
       class is (df_split_nonlinear_term_t)
          if (gpu_cpu_ratio.lt.0.0) then
             call this_nonlinear_term%autotune(lbg0)
          else
             call this_nonlinear_term%set_gpu_cpu_ratio(gpu_cpu_ratio)
          end if
#endif
#ifdef WITH_CUDA_NONLIN
       class is (df_nonlinear_cuda_t)
          call this_nonlinear_term%SetCudaDevice(1)
#endif
       end select
       call this_nonlinear_term%set_blocksize(lbg0)
       call this_nonlinear_term%initialize()
       if (parallel_nl) call initialize_parallel_nonlin
       if ((mype.eq.0).and.(my_thread.eq.0).and.print_ini_msg) then
          type_of_nonlin = this_nonlinear_term%getType()
          write(*,"(3A,I8)") "nonlinear type is ",type_of_nonlin, &
               &" with a blocksize of ",this_nonlinear_term%lbg0
       end if
    endif

    rhs_with_coll = ((collision_op.ne.'none').and.(.not.coll_split))

    ! following statement has been moved to discretization.F90
    !lbg0 = lklmn0/nblocks
    call initialize_klmn_conv_arrs

    if (init_status.eq.0) then
       if(with_g_update) then
          if(nblocks.eq.1) then
             allocate(rhs_big(li1:li2, lj1:lj2, 1:lbg0, 1:nblocks))
          else
             !with cache blocking
             allocate(rhs_small(li1:li2,lj1:lj2,1:lbg0))
             if (rhs_with_coll) then
                !collisions are performed without cache blocking
                allocate(rhs_big(li1:li2, lj1:lj2, 1:lbg0, 1:nblocks))
             endif
          end if
       endif

       if (.not.xy_local) then
          !allocate(bar_emfields(li1:li2,lj1:lj2,lbz:ubz,lm1:lm2,ln1:ln2,1:n_fields))
          !ptr_bar_emfields => bar_emfields
          allocate(g_block_storage(lbi:ubi,lj1:lj2,1:lbg0))
          allocate(chi_block(lbi:ubi,lj1:lj2,1:lbg0))
          if (rhs_nl) then
             call this_nonlinear_term%allocate_arrays(p_dbarchidxy,p_dgdxy)
             ptr_dbarchidxy => p_dbarchidxy
             ptr_dgdxy => p_dgdxy
          else
             NULLIFY(ptr_dbarchidxy)
             NULLIFY(ptr_dgdxy)
          endif
       endif
       if  (.not.allocated(bar_emfields)) allocate(bar_emfields(li1:li2,lj1:lj2,lbz:ubz,lm1:lm2,ln1:ln2,1:n_fields))
       init_status = 1
    endif

    if ((perf_vec(2).eq.1).and.rhs_nl) then
       if (.not.allocated(p_barchi)) allocate(p_barchi(li1:li2,lj1:lj2,1:lbg0))
       ptr_barchi => p_barchi
    else
       if (allocated(p_barchi)) deallocate(p_barchi)
       nullify(ptr_barchi)
    endif

    if (collision_op.ne.'none') call initialize_add_coll

    dfdzv_replace_rhs=.true.
    if (rhs_with_coll.and.(.not.with_g_update)) dfdzv_replace_rhs=.false.
    !because collisions initialize a_rhs and rhs_block points on a block of a_rhs.
    !(with g_update, an additional rhs_big contains collisions so that rhs_block
    !must be initialized by the dfdzv_terms)

    if (.not.x_local.and.explicit_buffer) call initialize_krookBuffer_operator

    if(.not.arakawa_zv) then
       call initialize_dfdzv(dfdzv_replace_rhs)
       call initialize_dchidz
    else
       call initialize_dzv(dfdzv_replace_rhs)
    endif
    if((hyp_x.gt.0).or.(hyp_y.gt.0).or.(hyp_perp.gt.0).or.Erad.ne.0.or.(GyroLES)) then
       with_dfdperp=.true.
    else
       with_dfdperp=.false.
    end if

    if (with_dfdperp) call initialize_dfdxy

    call initialize_dchidxy

    call initialize_dgdxy

    if (rhs_f0) call initialize_f0_term

    if (with_sources) call initialize_f1_sources

  end subroutine initialize_CalFullRhs


  !>Computes the full right hand side of the gyrokinetic equation for a given g_1. First the
  !!fields f_ and emfields are computed, then the right hand side of the gyrokinetic equation
  !!(see thesis tbg) is computed with the routine chosen implementation
  Subroutine CalFullRhs(p_g_,rhs,stage)
    ! INPUT
    ! p_g_ : distribution function g
    Complex, Dimension(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2),Intent(INOUT):: p_g_
    !****
    ! rhs is the right hand side, which is to be calculated from
    ! the given quantities
    Complex, Dimension(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2),Intent(INOUT) :: rhs
    integer:: stage
    real :: local_sum,global_sum
    logical,parameter :: OUTPUT=.false.

    LOGICAL :: overflow_exit, underflow_exit
    logical :: reset=.false.

    if (rhs_only_coll) then
       call calc_aux_fields(p_g_,emfields,f_,comp_h=coll_on_h,h_out=h_)
       if (coll_on_h) then
          call equ_collisions(h_,rhs,replace_rhs=.true.)
       else
          call equ_collisions(f_,rhs,replace_rhs=.true.)
       endif
    else
       !call calc_aux_fields(p_g_,emfields,f_,comp_h=arakawa_zv,h_out=h_)
       call calc_aux_fields(p_g_,emfields,f_,comp_h=(arakawa_zv.or.coll_on_h),h_out=h_,stage=stage)
       call calc_rhs_only(f_, p_g_, emfields, rhs, stage, h_)
    endif

#ifdef COUPLE_XGC
    if (stage.eq.1) call exec_all_diags(itime,time,overflow_exit,underflow_exit,reset)
#endif

    IF (OUTPUT) THEN
       CALL calculate_test_sum(f_(:,:,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2),local_sum, global_sum)
       IF (mype.EQ.0) WRITE(*,"(A,ES20.12)") "f_ = ",global_sum
    END IF

  end Subroutine CalFullRhs


  !>RHS without field solve for implicit shape conversion of p_g_, rhs, etc.
  !!shall be removed once strip mining has been applied to all parts of the
  !!Vlasov equation
  Subroutine calc_rhs_only (p_f_, p_g_, p_emfields, a_rhs, stage, p_h_)
    ! INPUT
    ! p_f_ : distribution function f1
    ! p_g_ : distribution function g
    ! p_emfields : fields
    ! a_rhs: right hand side for the old RK schemes / a_ vector for the low memory schemes
    ! stage: RK stage
    ! p_h_ : distribution function h1 (relevant for the arakawa_zv scheme)
    Complex, Dimension(li1:li2, lj1:lj2, 1:lzvwn0),Intent(INOUT):: p_f_
    Complex, Dimension(li1:li2, lj1:lj2, 1:lbg0, 1:nblocks),Intent(INOUT):: p_g_
    Complex, Dimension(li1:li2, lj1:lj2, lbz:ubz, 1:n_fields),Intent(INOUT):: p_emfields
    integer:: stage
    complex, dimension(li1:li2, lj1:lj2, 1:lbg0, 1:nblocks),target,intent(inout) :: a_rhs
    Complex, Dimension(li1:li2, lj1:lj2, 1:lzvwn0),optional,Intent(INOUT):: p_h_
    complex, dimension(:,:,:,:),pointer:: rhs, p_a_
    complex, dimension(:,:,:),pointer:: rhs_block, g_block
    integer :: iblock, lbg1,lbg2, my_thread,klmn

    DEBUG(1,"Start RHS")
    PERFON('CalFRhs0')

    if(stage.eq.1) then
       ve_max = 0.0
       ve_pnl_max = 0.0
    endif

    if (.not.xy_local) then
       call set_bar_emfields(p_emfields, bar_emfields)
       if ((p_has_0_mode).and.(heatsrc_on.or.partsrc_on)) then
          call compute_f_vabs_allspec(p_f_)
       end if
    endif

    if(with_g_update) p_a_=>a_rhs

    !collisions do not work with strip mining because of velocity space/species integrals
    if (rhs_with_coll) then
       if(with_g_update) then
          rhs=>rhs_big
       else
          rhs=>a_rhs
       end if
       if (coll_on_h) then
          call equ_collisions(p_h_,rhs,replace_rhs=.true.)
       else
          call equ_collisions(p_f_,rhs,replace_rhs=.true.)
       endif
    endif

    !computation is done in sub-blocks of g-like arrays to optimize cache usage

    my_thread = 0

    if(with_g_update) then
       if(nblocks.eq.1) then
          rhs_block=>rhs_big(:,:,:,1)
       else
          rhs_block=>rhs_small
       endif
    endif

    if (.not.xy_local) g_block => g_block_storage

    do iblock=1,nblocks
       lbg1 =(iblock-1)*lbg0+1
       lbg2 =lbg1+lbg0-1

       ! do the boundary exchange in x direction once for the whole block in advance
       ! This improves performance due to lower latency impact and larger messages.
       ! But we need a copy of g and chi which contains a whole block.
       if (.not.xy_local) then
          PERFON('blex_chi')
          if (n_fields.gt.1) then
             do klmn=lbg1,lbg2
                chi_block(li1:li2,lj1:lj2,klmn-lbg1+1) = &
                     &bar_emfields(li1:li2,lj1:lj2,sk(klmn),sm(klmn),sn(klmn),1) &
                     &- vTvpar(sl(klmn),sn(klmn)) * &
                     & bar_emfields(li1:li2,lj1:lj2,sk(klmn),sm(klmn),sn(klmn),2)
             end do
             call exchange_x(x_boundary_block,chi_block)
          else
             ! here chi is independent of v_par, therefore one exchange for each j,k,m,n
             ! is enough
             do klmn=lbg1,lbg2
                chi_block(li1:li2,lj1:lj2,klmn-lbg1+1) = &
                     &bar_emfields(li1:li2,lj1:lj2,sk(klmn),sm(klmn),sn(klmn),1)
             end do
             call exchange_x(x_boundary_block,chi_block)

          end if
          PERFOFF
       end if
       if (.not.xy_local) then
          PERFON('blex_g')
          g_block(li1:li2,lj1:lj2,1:lbg0) = p_g_(li1:li2,lj1:lj2,1:lbg0, iblock)
          call exchange_x(x_boundary_block,g_block)
          PERFOFF
       end if

       if(.not.with_g_update) rhs_block=>a_rhs(:,:,:,iblock)

       if(arakawa_zv) then
          PERFON('dzv_ak')
          call equ_dzv(p_h_,rhs_block,lbg1,lbg2)
          PERFOFF
          if (hyp_on_h) then
             if (hypz_compensation) call equ_comp_hypz(p_emfields,bar_emfields,rhs_block,lbg1,lbg2)
          else
             PERFON('hyp_zv_ak')
             call add_hypz_ak(p_f_,rhs_block,lbg1,lbg2)
             call add_hypv_ak(p_f_,rhs_block,lbg1,lbg2)
             PERFOFF
          endif
       else
          PERFON('eq_dfdzv')
          ! standard case (should be threadsafe)
          call equ_dfdzv(p_f_,rhs_block,lbg1,lbg2)
          PERFOFF
       endif

       ! standard explicit_buffer=.false.
       if (.not.x_local.and.explicit_buffer) call add_kBuffer_explicit(p_f_,rhs_block,lbg1,lbg2)
       !>todo: check if perpendicular hyperdiffusion should act on f or h!
       if (with_dfdperp) then
          if (hyp_on_h) then
             call add_dfdxy(p_h_,rhs_block,lbg1,lbg2) !We need to add hyp_x and hyp_y on h
          else
             call add_dfdxy(p_f_,rhs_block,lbg1,lbg2) !We need to add hyp_x and hyp_y on f
          endif
       endif
       if (rhs_f0) call add_f0_term(rhs_block,lbg1,lbg2,time)


       PERFON('dchidxy')
       if (xy_local) then
          call add_dchidxy(p_emfields,rhs_block,ptr_barchi,lbg1,lbg2,.true.)
       else
          call add_dchidxy(chi_block,rhs_block,ptr_dbarchidxy,lbg1,lbg2,.true.)
       end if
       PERFOFF

       PERFON('dchidz')
       if (.not.arakawa_zv) call add_dchidz(p_emfields, bar_emfields, &
            &rhs_block, lbg1, lbg2)
       PERFOFF



       PERFON('dgdxy')
       if (xy_local) then
          call add_dgdxy(p_g_, rhs_block, ptr_dgdxy, pdg1di, pdg1dj, lbg1,lbg2)
       else
          call add_dgdxy(g_block, rhs_block, ptr_dgdxy, pdg1di, pdg1dj, lbg1,lbg2)
       end if
       PERFOFF


       if (rhs_nl) then
          PERFON('add_nl')
          if (nonlinear) then !perpendicular nonlinearity
             if (xy_local) then
                if (.not.nonlin_h) then
                   call this_nonlinear_term%add(p_g_, ptr_dgdxy, p_emfields,&
                        &ptr_barchi, ptr_dbarchidxy, rhs_block, lbg1, lbg2, stage)
                else
                   call this_nonlinear_term%add(p_h_, ptr_dgdxy, p_emfields,&
                        &ptr_barchi, ptr_dbarchidxy, rhs_block, lbg1, lbg2, stage)
                endif
             else
                select type(tmp=>this_nonlinear_term)
                class is (df_arakawa_nonlinear_term_t)
                   ! this copy operation is only necessary for arakawa scheme

                   ptr_barchi(:,:,:) = chi_block(li1:li2,lj1:lj2,:)
                end select
                call this_nonlinear_term%add(g_block, ptr_dgdxy, p_emfields,&
                     &ptr_barchi, ptr_dbarchidxy, rhs_block, lbg1, lbg2, stage)
             end if
          endif
          PERFOFF
          if (parallel_nl) then
             PERFON('paral_nl')
             call add_parallel_nonlin(p_f_,p_emfields,&
                  &bar_emfields, ptr_dbarchidxy, rhs_block, lbg1, lbg2, stage)
             PERFOFF
          end if
       endif

       PERFON('dApardt')
       if (antenna_type.eq.3.and.antenna_contrib) call add_dApar_dt_antenna(rhs_block,lbg1,lbg2)
       PERFOFF



       if (with_sources.and.(.not.precond_approx)) then
          PERFON('sources')
          call add_f1_sources(rhs_block,lbg1,lbg2,stage)
          PERFOFF
       end if

       if(with_g_update) then
          PERFON('up_ag')
          !add the collision contribution if rhs_block.ne.rhs

          if ((nblocks.gt.1).and.(rhs_with_coll)) then
             call axpy_ij(lij0*lbg0,1.0,rhs(:,:,:,iblock),rhs_block)
             ! rhs_block=rhs_block+rhs(:,:,:,iblock)
          end if
          p_g_(:,:,:,iblock)=p_a_(:,:,:,iblock)+a_rk(stage)*dt*rhs_block
          if(stage.lt.rkstages) then
             call axpy_ij(lij0*lbg0,b_rk(stage)*dt,rhs_block,&
                  &p_a_(:,:,:,iblock))
             !p_a_(:,:,:,iblock)=p_a_(:,:,:,iblock)+b_rk(stage)*dt*rhs_block
          else
             call ccopy(lij0*lbg0,p_g_(:,:,:,iblock),1,p_a_(:,:,:,iblock),1)
!             p_a_(:,:,:,iblock)=p_g_(:,:,:,iblock)
          endif
          PERFOFF
       endif
    enddo

    PERFOFF
    DEBUG(1,"End RHS")

  End Subroutine calc_rhs_only


  !>Wrapper for the cleanup routines of the various implementations
  subroutine finalize_CalFullRhs

    if (rhs_nl) then
       if (allocated(p_barchi)) deallocate(p_barchi)
       nullify(ptr_barchi)
    endif

    if (.not.x_local.and.explicit_buffer) call finalize_krookBuffer_operator

    if (with_sources) call finalize_f1_sources

    if (.not.xy_local) then
       !deallocate(bar_emfields)
       deallocate(g_block_storage)
       deallocate(chi_block)
       if (rhs_nl) then
          call this_nonlinear_term%free_arrays(p_dbarchidxy,p_dgdxy)
          !deallocate(p_dbarchidxy,p_dgdxy)
          nullify(ptr_dbarchidxy,ptr_dgdxy)
       endif
    endif
    if (allocated(bar_emfields)) deallocate(bar_emfields)

    if(arakawa_zv) then
       call finalize_dzv
    else
       call finalize_dfdzv
       call finalize_dchidz
    end if

    if (collision_op.ne.'none') call finalize_add_coll
    if (with_dfdperp) call finalize_dfdxy

    if(rhs_f0) call finalize_f0_term

    call finalize_dchidxy


    call finalize_dgdxy
    if (rhs_nl) then
       if (parallel_nl) call finalize_parallel_nonlin
       call this_nonlinear_term%finalize()
       call this_nonlinear_term%destruct()
       deallocate(this_nonlinear_term)
    endif

    call finalize_klmn_conv_arrs

    if(with_g_update) then
       if(nblocks.eq.1) then
          deallocate(rhs_big)
       else
          deallocate(rhs_small)
          if (rhs_with_coll) deallocate(rhs_big)
       end if
    endif

    init_status = 0

  end subroutine finalize_CalFullRhs


  subroutine set_bar_emfields(p_emfields,p_bar_emfields)
    Complex, Dimension(li1:li2, lj1:lj2, lbz:ubz, 1:n_fields),Intent(IN):: p_emfields
    Complex, Dimension(li1:li2, lj1:lj2, lbz:ubz, lm1:lm2, ln1:ln2, 1:n_fields),Intent(OUT):: p_bar_emfields

    integer :: m,n,o

    do o=1,n_fields
       do n=ln1,ln2
          do m=lm1,lm2
             call gyro_average_df_wb(p_emfields(:,:,:,o),&
                  &p_bar_emfields(:,:,:,m,n,o),m,n)

             ! boundary exchange only for phi, as we only need derivatives
             ! of phi in x,y and z direction
             if (o==1) call exchange_z(p_bar_emfields(:,:,:,m,n,1))
          enddo
       enddo
    enddo

  end subroutine set_bar_emfields

End Module calc_rhs
