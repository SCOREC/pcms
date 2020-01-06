#include "redef.h"
#include "intrinsic_sizes.h"
#include "switches.h"
!>Computes the contribution of collisions to the RHS
!!
!!Implements the linearized Landau-Boltzmann collision operator for an arbitrary number of species
!!(the reduction to pitch angle scattering is also implemented for testing purposes).
!!First the test particle contribution is computed, the resulting violation of energy and
!!parallel momentum conservation is then used in the model operator for the field particle
!!contribution.
Module collisions
  Use par_mod
  Use communications
  Use vel_space
  Use geometry, only: geom, trpeps, q0, major_R!, magn_geometry
  use boundaries
  Use aux_func, only: j0, j1
  use profiles_mod, only: zeff, zeff0, zeff_prof
 ! use coordinates
  use GaussQuadrature, only: mu_grid_type
  use gyro_average_ff_mod, only: jfac
  use gyro_average_df_mod, only: gyro_average_df
  Implicit None
  public:: check_par_collisions, initialize_collisions, initialize_add_coll, equ_collisions, &
       finalize_add_coll, finalize_collisions, coll, nu_ei, nustar_i, nustar_e, coll_cons_model, &
       mem_est_collisions, spacediff, spacediff_off, set_coll_defaults, adiabatic_mass, ev_coll_est
  private

  integer :: init_status = 0
  integer :: init_status_fp = 0

  Real:: coll=0.0, nu_ei=0.0, nustar_i=0.0, nustar_e=0.0
  character(len=16):: coll_cons_model = 'xu_rosenbluth'
  Logical:: em_conserve=.true.
  integer:: em_cons=1
  real:: adiabatic_mass = 2.0
  real, dimension(:,:),allocatable :: zeff_spec

  logical:: spacediff=.false., spacediff_off=.true.

  Real, Dimension(:,:,:,:,:,:), Allocatable:: cd11,cd12,cd21,cd22,cr1,cr2
  Real, Dimension(:), Allocatable :: dmuin, dvinv, cvf_mu, cvf_vp
  Real:: dvin
  Real:: ev_coll_est

  !for sugama implementation
  Real, Dimension(:,:,:,:,:), Allocatable:: T_coeff1,T_coeff2,T_coeff3
  Real, Dimension(:,:,:,:,:), Allocatable:: T_coeff4,T_coeff5,T_coeff6

  Real, Dimension(:,:,:,:,:,:), Allocatable:: prefacbuff1,prefacbuff2,prefacbuff3
  Real, Dimension(:,:,:,:,:,:), Allocatable:: prefacbuff4,prefacbuff5,prefacbuff6
  Real, Dimension(:), Allocatable:: vpp, muu
  Real, Dimension(:), Allocatable:: stencil_first,stencil_second
  Real, Dimension(:,:,:,:,:,:), Allocatable:: particle_diffpart
  Complex, Dimension(:,:), Allocatable:: particle_buff

  Real, Dimension(:), Allocatable:: dvinverse,dmuinverse
  Real, Dimension(:,:,:), Allocatable:: vperpcoord,dvperpinverse

  Real, Dimension(:,:,:,:,:), Allocatable:: particle_coeff,sigma_bar
  Real, Dimension(:,:,:), Allocatable:: sigma_bar_int,sigma_int_num,sigma_int_den
  Complex, Dimension(:,:,:,:,:), Allocatable:: particle_fac

  Real, Dimension(:,:,:,:), Allocatable:: mat_00_sub

  !for energy/momentum conservation
  complex, dimension(:,:,:,:,:), allocatable:: mom_cmat, en_cmat
  Real, Dimension(:,:,:,:,:), Allocatable:: loc_mat_10, mat_01_20
  complex,dimension(:,:,:,:),allocatable :: m_fac_sum, e_fac_sum

  !em_cons_2
  complex, dimension(:,:,:,:,:,:), allocatable:: cons_mom_cmat, cons_en_cmat
  complex,dimension(:,:,:,:,:),allocatable:: rec_buff

  !for add_testpart_2
  complex,dimension(:,:,:,:,:),allocatable:: vfcoll_2

  !for add_testpart_3
  real,dimension(:,:,:,:,:,:,:),allocatable:: collsten_big
  complex,dimension(:,:,:,:,:,:),allocatable:: vfcoll_big
  integer,dimension(9):: shiftvec
  integer,dimension(9):: stenlb, stenub
#ifdef EXTERNAL_ERF
  EXTERNAL erf
#endif

  !variables for spatial diffusion part
  complex,dimension(:,:,:),allocatable:: kpf
  real, dimension(:,:,:),allocatable:: kperp2
  real, dimension(:,:,:,:,:,:),allocatable:: sp_pref

  !f/f0-implementation
  complex,dimension(:,:,:,:,:,:), target, Allocatable:: loc_f_fm

  integer::n_spec_coll  !size of collision matrix

  !Variable needed for adding the collisional moments
  real,dimension(:,:,:,:,:,:), allocatable:: moments_colla, moments_collb, moments_collc
  real,dimension(:,:,:,:,:,:), allocatable:: moments_coll1, moments_coll2, moments_coll3
  real,dimension(:,:,:,:,:), allocatable:: moments_coll4, moments_coll5, moments_coll6
  real,dimension(:,:,:,:,:,:), allocatable:: C_ab_TO_mom, C_ab_TO_en
  real,dimension(:,:,:,:,:), allocatable:: momentum_input, energy_input
  real,dimension(:,:,:,:,:,:), allocatable:: prefacterm1, prefacterm2, prefacterm3
  real,dimension(:,:,:,:,:,:), allocatable:: prefacterm4, prefacterm5, prefacterm6
  real,dimension(:,:,:,:,:,:), allocatable :: x_coeff1, x_coeff2, x_coeff3
  real,dimension(:,:,:,:,:), allocatable :: x_coeff4, x_coeff5, x_coeff6
  real,dimension(:,:,:,:,:,:), allocatable :: y_coeff1, y_coeff2, y_coeff3
  real,dimension(:,:,:,:,:,:), allocatable :: y_coeff4, y_coeff5, y_coeff6
  real,dimension(:,:,:,:,:,:), allocatable :: y_coeffa, y_coeffb, y_coeffc
  real,dimension(:,:,:,:), allocatable :: denominator1, denominator2, denominatora, denominatorb
  real,dimension(:,:,:,:), allocatable :: qintegral, maxwellianint, sintegral
  real,dimension(:,:,:,:,:,:), allocatable :: qfactor, pfactor, sfactor
  real,dimension(:,:,:,:,:), allocatable :: j1fac
  complex, dimension(:,:,:), allocatable :: coll_buff
  complex,dimension(:,:), allocatable :: gyro_in_a,gyro_in_b,gyro_out_a,gyro_out_b
  real,dimension(:,:,:,:), allocatable :: x_coeff3_int
  real,dimension(:,:,:), allocatable :: x_coeff6_int
  real,dimension(:,:,:,:), allocatable :: R_int,momentum_int1,momentum_int2,quartic_int
  complex,dimension(:,:,:,:,:,:,:), allocatable :: collision_mom
  complex,dimension(:,:,:,:,:), allocatable :: collision_mom2

Contains

  !>Give an estimate of the memory requirements of this module
  Real Function mem_est_collisions(mem_req_in)
    real:: mem_req_in
    real:: mem_loc=0

  if ((collision_op) .eq. 'sugama') then

    !T_coeffs and particle_diffpart
    mem_loc= 6*SIZE_OF_REAL_MB*pi0*pj0*lk0*ll0*lm0*ln0
    mem_loc= SIZE_OF_REAL_MB*pi0*pj0*lk0*ll0*lm0*ln0*n_spec
    !particle_coeff
    mem_loc= mem_loc+SIZE_OF_REAL_MB*pi0*pj0*lk0*ll0*lm0*ln0
    !particle_fac
    mem_loc= mem_loc+SIZE_OF_COMPLEX_MB*lijk0*ln0*n_spec
    !particle_buff
    mem_loc= mem_loc+SIZE_OF_COMPLEX_MB*lij0

  else
    !cd11 etc.
    mem_loc=3*SIZE_OF_REAL_MB*pi0*pj0*lk0*(ll0+1)*lm0*n_spec_coll*ln0
    mem_loc=mem_loc+3*SIZE_OF_REAL_MB*pi0*pj0*lk0*ll0*(lm0+1)*n_spec_coll*ln0
    !mom_cmat, en_cmat
    mem_loc=mem_loc+2*SIZE_OF_COMPLEX_MB*pi0*lklmn0
    !loc_mat_10, mat_01_20
    mem_loc=mem_loc+2*SIZE_OF_REAL_MB*pi0*lklmn0

    !mom_fac, en_fac
    mem_loc=mem_loc+2*SIZE_OF_COMPLEX_MB*lijk0*ln0*n_spec

    !vfcoll
    if (perf_vec(4)==1) mem_loc=mem_loc+SIZE_OF_COMPLEX_MB*(ll0+2)*(lm0+2)*n_spec_coll

    if (perf_vec(4)==2) then
       !vfcoll_2
       mem_loc=mem_loc+SIZE_OF_COMPLEX_MB*lij0*(ll0+2)*(lm0+2)*n_spec_coll
       !fcd11 etc., fluvp, flumu
       mem_loc=mem_loc+8*SIZE_OF_COMPLEX_MB*lij0*(ll0+2)*(lm0+2)*n_spec_coll
    endif

    if (perf_vec(4)==3) then
       !collsten_big
       mem_loc=mem_loc+SIZE_OF_REAL_MB*li0*lj0*lz0*lv0*lw0*9*ln0*n_spec_coll
       !vfcoll_big
       mem_loc=mem_loc+SIZE_OF_COMPLEX_MB*li0*lj0*lz0*lv0*lw0*n_spec_coll
       !collsten (is only temporary....but anyway)
       if(.not.xy_local) then
          mem_loc=mem_loc+li0*SIZE_OF_REAL_MB*lz0*lv0*lw0*ln0*n_spec_coll*9
       endif
    endif
  endif

    select case(em_cons)
    case (1)
       !mom_cmat, en_cmat
       mem_loc=mem_loc+2*SIZE_OF_COMPLEX_MB*pi0*lklmn0
       !mom_coeff, en_sum
       mem_loc = mem_loc + SIZE_OF_COMPLEX_MB*5*pi0*pj0*lk0*ln0
       !n_cmat
       mem_loc = mem_loc + SIZE_OF_COMPLEX_MB*pi0*pj0*lklmn0
       !e_fac_sum, m_fac_sum
       mem_loc=mem_loc+2*SIZE_OF_COMPLEX_MB*lijk0*n_spec
    case (2)
       !cons_mom_cmat, cons_en_cmat
       mem_loc=mem_loc+2*SIZE_OF_COMPLEX_MB*pi0*lklmn0*n_spec
       !cons_mom_coeff, cons_en_sum
       mem_loc = mem_loc + SIZE_OF_COMPLEX_MB*5*pi0*pj0*lk0*ln0*n_spec
       !cons_n_cmat
       mem_loc = mem_loc + SIZE_OF_COMPLEX_MB*pi0*pj0*lklmn0*n_spec
       !rec_buff
       mem_loc=mem_loc+SIZE_OF_COMPLEX_MB*lijk0*n_spec*n_spec
    case (3)
       !moments_coll
       mem_loc=mem_loc+3*SIZE_OF_REAL_MB*pi0*pj0*lk0*ll0*lm0*ln0*(1+2*n_spec)
       !x_coeff's
       mem_loc=mem_loc+3*SIZE_OF_REAL_MB*pi0*pj0*lk0*ll0*lm0*ln0*(1+n_spec)
       !y_coeff's
       mem_loc=mem_loc+9*SIZE_OF_REAL_MB*pi0*pj0*lk0*ll0*lm0*ln0*n_spec
       !collision_mom
       mem_loc=mem_loc+SIZE_OF_COMPLEX_MB*3*lijk0*n_spec*n_spec*2
       !collision_mom2
       mem_loc=mem_loc+SIZE_OF_COMPLEX_MB*3*lijk0*n_spec
       !qfactor
       mem_loc=mem_loc+SIZE_OF_REAL_MB*pi0*pj0*lk0*ll0*lm0*ln0*n_spec
       !Bessel functions
       mem_loc=mem_loc+SIZE_OF_REAL_MB*lijk0*lm0*ln0
       !coll_buff
       mem_loc=mem_loc+SIZE_OF_COMPLEX_MB*lij0*3
    endselect

    if (spacediff) then
      !kperp2, kpf
      mem_loc = mem_loc+(SIZE_OF_COMPLEX_MB+SIZE_OF_REAL_MB)*lijk0
      !sp_pref
      mem_loc = mem_loc+SIZE_OF_REAL_MB*pi0*pj0*lklmn0
    endif

    !loc_f_fm
    if (coll_f_fm_on) mem_loc=mem_loc+SIZE_OF_COMPLEX_MB*lij0*lz0*lv0*lw0*ln0

    mem_est_collisions=mem_req_in+mem_loc
  End Function mem_est_collisions

  !> sets the defaults for the collision input parameters.
  subroutine set_coll_defaults

    coll=0.0
    coll_cons_model = 'xu_rosenbluth'
    coll_order = 'second'
    em_conserve=.true.
    spacediff=.false.
    spacediff_off=.true.
    !collision_op = 'none'  !set in par_in.F90, because needed in vel_space at the moment
    !coll_f_fm_on = .false. !set in par_in.F90, "
    adiabatic_mass = 2.0  !mass of adiabatic ion species (if any) in units of m_proton

  end subroutine set_coll_defaults

  subroutine check_par_collisions
    logical :: write_pe
    integer :: n_ions
    real:: Ti, Te, me

    write_pe = (mype.eq.0).and.(print_ini_msg)
    if (collision_op.ne.'none') then
       if (coll.eq.0.) then
          if (write_pe) write(*,"(A)") "collision frequency is 0, deactivating collision term "
          collision_op = 'none'
       else
          if (write_pe) write(*,"(A)") "with collisions"
       endif
    endif

    if (collision_op.ne.'none') then
       if ((spacediff).and.(.not.xy_local)) then
          if (write_pe) write(*,"(A)") "nonlocal simulations do not support the spacediff term - exit"
          stop
       endif
       if (mod(nv0,2).ne.0) write(*,"(A)") 'WARNING: if nv0 is odd, dt_max might be very small!'

       if (collision_op=='krook') coll_cons_model = 'none'
       !note: benchmarks for pure pitch-angle collisions might require to set coll_cons_model='none'
       !pitch_all is pitch-angle collisions between all species:
       !conservation terms are allowed, they take into account momentum transfer
       !and correct for energy-transfer that may be introduced by the numerical scheme

       !backwards compatibility: set defaults for coll_cons_model
       if (.not. em_conserve) coll_cons_model = 'none'
       if (coll_cons_model .ne. 'none') em_conserve = .true.

       !set logicals for different fieldpart implementations
       select case(coll_cons_model)
       case('xu_rosenbluth','default','on')
          !current default:
          coll_cons_model='xu_rosenbluth'
          !implementation following flm thesis (Xu, Rosenbluth, Phys.Fluids B 1991)
          em_cons=1
       case('self_adj')
          !implementation following hkd thesis, (Stephan Brunner notes 2010)
          !+modification to conserve (gyrocenter) density, parallel momentum and energy
          ! to machine precision
          em_cons=2
       case('nakata')
          !implementation which follows the papers (Sugama, Phys. Plasmas 2009) and
          ! (Nakata, Computer Physics Communications, 2015). This operator modifies
          !the term which breaks self-adjointness symmetry, conserves perpendicular
          !momentum, and incorporates finite Larmour radius effects.
          em_cons=3
       case('none')
          !no field-particle model
          em_conserve = .false.
          em_cons=0
       case default
          if (write_pe) write(*,"(A)") "no valid choice for the parameter coll_cons_model in general namelist - exit"
          stop
       end select

       call get_nions(n_ions, Ti, Te, me)

       !stop also when Zeff is set (also -1) and n_ions>1
       if ((Zeff.ne.1.0).and.(n_ions.gt.1)) then
          if (write_pe) write(*,"(A)") "Zeff.ne.1.0 is not allowed with multiple ion species -- exit"
          stop
       endif

       n_spec_coll = n_spec
       !in case of adiabatic ions, we add e-i collisions, since they are same order as e-e collisions.
       if (n_ions==0) then
          !if (mype==0) write(*,"(A)") 'collisions include extra adiabatic ion species'
          select case (collision_op)
          case ('landau','pitch-all','on')
             !extra adiabatic ion species
             n_spec_coll = n_spec+1
          case default
            if (mype.eq.0) write(*,"(A)") 'ERROR: invalid collision_op for adiabatic ions -- exit!'
            stop
          endselect
       endif
    endif

    !collison operator splitting is only possible in initial value runs
    if (comp_type.ne.'IV') then
       coll_split=.false.
       coll_split_scheme='none'
    endif

    if (.not.coll_split) coll_split_scheme = 'none'

    select case (coll_split_scheme)
    case('EE1','RKC1','RKC2','RKC3','RKC4','RKCa','RK3','RK4','RK2')
       coll_split = .true.
    case('none')
       coll_split = .false.
    case default
       if (mype.eq.0) write(*,"(3A)") "coll_split_scheme = '",coll_split_scheme,"' is unknown, please choose a valid scheme"
       stop
    end select

    if (collision_op.eq.'none') then
       !deactivate collision-related switches
       coll_f_fm_on = .false.
       coll_split = .false.
       coll_split_scheme = 'none'
    endif

    if (coll_split.and.(coll.gt.0.05).and.write_pe) &
         write(*,"(A)") 'WARNING: coll_split is used at high collisionality: accuracy might decrease!'

  end subroutine check_par_collisions

  subroutine compute_nuei(nu_ei,nustar_i,nustar_e)
    real,intent(out)::nu_ei, nustar_i, nustar_e
    integer :: n, n_ions, firstion
    real:: Ti, Te, me

    call get_nions(n_ions, Ti, Te, me)

    !compute coll to nu_ei conversion
    ![Hinton Hazeltine definition, see gene documentation]
    !in normalized units c_ref/L_ref
    !here we sum over all ion species and multiply Zeff in case of one ion species
    nu_ei = 0.0
    do n=0,n_spec-1
       if (spec(n)%charge>0.and.(.not.spec(n)%passive)) then
          nu_ei=nu_ei+4*spec(n)%dens*spec(n)%charge**2/Te**2/me**0.5*coll
       endif
    enddo
    if (n_ions==1) nu_ei=nu_ei*zeff0
    !adiabatic ions: use zeff0=sum_i n_i q_i**2   (at x0)
    if (n_ions==0) nu_ei = 4*zeff0/Te**2/me**0.5*coll

    !trpeps available from geometry?
    if (trpeps==0.0) then
       !write(*,"(4A)") 'WARNING: no nustar output for ',magn_geometry, ' geometry:',&
       !     &' trpeps is not set'
       return
    endif

    !compute nustar
    !for multiple ion speciesn nustar_i is only approximately calculated
    !for the first ion species (colliding with itself and all others)
    nustar_i = 0.0
    nustar_e = 0.0
    firstion=-1

    do n=0,n_spec-1
       if (spec(n)%charge>0.and.(.not.spec(n)%passive)) then
          !n is the species scattered off
          if (firstion.lt.0) firstion=n
          nustar_i=nustar_i+8.0/3.0/pi**0.5*q0/trpeps**1.5*major_R&
               &*spec(n)%dens*spec(firstion)%charge**2*spec(n)%charge**2&
               &/Ti**2*coll
          nustar_e=nustar_e+16.0/3.0/pi**0.5*q0/trpeps**1.5*major_R&
               &*spec(n)%dens*spec(n)%charge**2/Te**2*coll
       endif
    enddo
    if (n_ions==1) then
       nustar_i = nustar_i*zeff0
       nustar_e = nustar_e*zeff0
    endif
    !adiabatic ions: use Zeff=sum_i n_i q_i**2
    if (n_ions==0) &
         nustar_e=16.0/3.0/pi**0.5*q0/trpeps**1.5*major_R*zeff0/Te**2*coll

    !no electrons..
    if (n_ions==n_spec) nustar_e=0.0

  end subroutine compute_nuei

  subroutine get_nions(n_ions, Ti, Te, me)
    integer,intent(out):: n_ions
    real,intent(out):: Ti, Te, me
    integer:: n, firstion

    Te=1.0
    Ti=1.0
    me=0.0002723
    n_ions = 0
    firstion = -1
    do n=0,n_spec-1
       if (spec(n)%charge.gt.0.and.(.not.spec(n)%passive)) then
          n_ions = n_ions+1
          if (firstion.lt.0) then
             firstion=n
             Ti=spec(firstion)%temp
          endif
       elseif (spec(n)%charge .lt.0.and.(.not.spec(n)%passive)) then
          Te=spec(n)%temp
          me=spec(n)%mass
       endif
    enddo
  end subroutine get_nions



!-----------------------------------------
!WRAPPER FOR THE DIFFERENT IMPLEMENTATIONS

  !>Initializes the alternative implementations of the collision operator
  subroutine initialize_add_coll

    if ((init_status.ne.perf_vec(4)).and.(init_status.gt.0).and.(collision_op.ne.'sugama')) &
         call finalize_add_coll

    if ((init_status_fp.ne.em_cons).and.(init_status_fp.gt.0)) &
         call finalize_fieldpart

    if ((init_status==0).and.(collision_op.ne.'sugama')) then
      ! if (collision_op .eq. 'sugama') then
      !    call initialize_testdifferentialpart
      ! else
       select case(perf_vec(4))
       case(2)
          call initialize_testpart_2
       case(3)
          call initialize_testpart_3
       end select
      ! endif
    endif

    init_status = perf_vec(4)

    if ((init_status_fp==0).and.(em_conserve)) &
        call initialize_fieldpart

  end subroutine initialize_add_coll


  !>Computes the contribution of collisions from f_ and adds it to a g_1 type array containing the
  !!contributions of the other terms
  subroutine equ_collisions(loc_f,loc_k,replace_rhs)
    !>f_ type array with exchanged ghost cells in vpar direction
    complex,dimension(li1:li2, lj1:lj2, lbz:ubz, lbv:ubv, lbw:ubw, ln1:ln2),target,intent(inout) :: loc_f
    !>g_1 type array containing the result
    complex,dimension(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2),intent(inout) :: loc_k
    !>replace loc_k with result or add result to loc_k
    complex, dimension(:,:,:,:,:,:),pointer :: ptr_f
    !>pointer for input array (f_over_fm or f_)
    logical::replace_rhs
    complex,dimension(li1:li2, lj1:lj2, lk1:lk2, ln1:ln2, 0:n_spec-1):: mom_fac, en_fac

    PERFON('equ_coll')

    if (replace_rhs) loc_k=(0.,0.)

    if (coll_f_fm_on) then
       !collisions require input f_fm
       call divide_by_fm(loc_f,loc_f_fm)
       call exchange_mu(loc_f_fm)
       ptr_f=>loc_f_fm
    else
       !collisions require input f_
       call exchange_mu(loc_f)
       ptr_f=>loc_f
    endif

    if (collision_op.ne.'krook') then
       if (collision_op .eq. 'sugama') then
          select case (perf_vec(4))
          case (1)
             if (n_procs_v .gt. 1) then
                call add_testdifferentialpart(ptr_f,loc_k)
             else
                call add_testdifferentialpart_alt(ptr_f,loc_k)
             endif
          case default
             stop "this perf_vec(4) option is not implemented"
          end select
       else
          select case (perf_vec(4))
          case (1)
             call add_testpart_1(ptr_f,loc_k,mom_fac,en_fac)
          case (2)
             call add_testpart_2(ptr_f,loc_k,mom_fac,en_fac)
          case (3)
             call add_testpart_3(ptr_f,loc_k,mom_fac,en_fac)
          end select
       endif
    else
       loc_k=-coll*loc_f(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
    endif

    if(.not.precond_approx) then
       select case(em_cons)
       case (1)
          call add_fieldpart_1(mom_fac,en_fac,loc_k)
       case (2)
          call add_fieldpart_2(mom_fac,en_fac,loc_k)
       case (3)
          if (xy_local) then
             call add_momentpart_local(ptr_f,loc_k)
          else
             call add_momentpart_global(ptr_f,loc_k)
          endif
       end select

       !add k_perp term
       if (spacediff) then
          call add_spacepart(ptr_f,loc_k)
       end if
    endif

    PERFOFF

  end subroutine equ_collisions

  !>Deletes arrays used by the alternative implementations of the collision operator
  subroutine finalize_add_coll

  if (collision_op .ne. 'sugama') then
    select case(init_status)
    case(2)
       call finalize_testpart_2
    case(3)
       call finalize_testpart_3
    end select
  endif

    if (em_conserve) call finalize_fieldpart

    init_status = 0

  end subroutine finalize_add_coll

  !>temporary subroutine: divides input distribution by fm and stores the result in loc_f_fm
  !loc_f_fm is then input to the collision routines if coll_f_fm_on is set true
  subroutine divide_by_fm(loc_f,loc_f_fm)
    !>f_ type array with exchanged ghost cells in vpar and mu directions
    complex,dimension(li1:li2, lj1:lj2, lbz:ubz, lbv:ubv, lbw:ubw, ln1:ln2) :: loc_f
    !>f_ type array with exchanged ghost cells in vpar and mu directions
    complex,dimension(li1:li2, lj1:lj2, lbz:ubz, lbv:ubv, lbw:ubw, ln1:ln2),intent(inout) :: loc_f_fm

    integer:: j,k,l,m,n,pni,lv,uv

    loc_f_fm = (0.0,0.0)

    !exclude unnessecary outer v points (otherwise, fm=0 -> division by 0)
    lv=lbv
    if(lv.lt.0) lv=0
    uv=ubv
    if(uv.gt.nv0-1) uv=nv0-1

    PERFON('div_fm')

    pni=pn1
    do n=ln1,ln2
       if (pn0.gt.1) pni=n
       do m=lm1,lm2
          do l=lv,uv !with inner v ghost cells!
             do k=lk1,lk2
                do j=lj1,lj2
                   if (xy_local) then
                      loc_f_fm(:,j,k,l,m,n) = loc_f(:,j,k,l,m,n)/fm(pi1,pj1,k,l,m,pni)
                   else
                      loc_f_fm(:,j,k,l,m,n) = loc_f(:,j,k,l,m,n)/fm(:,pj1,k,l,m,pni)
                   end if
                end do
             end do
          end do
       end do
    end do
    !if (mype.eq.0) print*,sum(abs(loc_f_fm(:,:,lk1:lk2,ll1:ll2,lm1:lm2,:)))
    PERFOFF
  end subroutine divide_by_fm


!---------
!TESTPART1

  !>Adds the test particle contribution and computes the violations of energy and parallel momentum
  !!conservation
  !!
  !!Alternative implementation 1: based on loops, easy to understand but bad performance on
  !!some machines
  subroutine add_testpart_1(loc_f,loc_k,mom_fac,en_fac)
    !>f_ type array with exchanged ghost cells in vpar and mu directions
    complex,dimension(li1:li2, lj1:lj2, lbz:ubz, lbv:ubv, lbw:ubw, ln1:ln2),intent(in):: loc_f
    !>g_1 type array to which the result is added
    complex,dimension(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2),intent(inout):: loc_k
    !>parallel momentum violation induced by test particle contribution
    complex,dimension(li1:li2, lj1:lj2, lk1:lk2, ln1:ln2, 0:n_spec-1),intent(out):: mom_fac
    !>energy violation induced by test particle contribution
    complex,dimension(li1:li2, lj1:lj2, lk1:lk2, ln1:ln2, 0:n_spec-1),intent(out):: en_fac

    complex,dimension(ll1-1:ll2+1,lm1-1:lm2+1, 0:n_spec_coll-1) :: vfcoll
    complex ::fluvp, flumu, fcd11,fcd12,fcd21,fcd22,fcr1,fcr2
    integer:: i,j,ii,k,l,m,n,iv,sp,ll,ul,lm,um

    PERFON('collt1')
    mom_fac=(0.,0.)
    en_fac=(0.,0.)
    do n=ln1,ln2
#ifdef WITHOMP_COLLISIONS
       !$omp parallel do default(NONE) &
       !$omp shared(my_pew,n_procs_w,my_pev,n_procs_v,n_spec,lijk0,lk0,lj0,cd11,cd12,cd21,cd22,cr1,cr2,loc_f,loc_k) &
       !$omp shared(n,li1,li2,pi1,lj1,lj2,pj1,lk1,ll1,ll2,dvin,dvinv,dmuin,cvf_mu,cvf_vp,lm1,lm2,xy_local) &
       !$omp shared(n_spec_coll,loc_mat_10,mat_01_20,en_fac,mom_fac,em_conserve) &
       !$omp private(i,j,ii,k,l,m,ll,ul,lm,um,flumu,fluvp,iv,vfcoll,fcd11,fcd12,fcd21,fcd22,fcr1,fcr2,sp)
#endif
       do iv=0,lijk0-1
          k=lk1+Modulo(iv,lk0)
          j=lj1+Modulo(iv/lk0,lj0)
          i=li1+iv/(lk0*lj0)
          vfcoll=(0.,0.)
          if (xy_local) then
             ii=pi1
          else
             ii=i
          end if
          do m=lm1,lm2
             if ((my_pew.eq.0).and.(m.eq.lm1)) then
                lm=m+1
             else
                lm=m-1
             end if
             if ((my_pew.eq.(n_procs_w-1)).and.(m.eq.lm2)) then
                um=m-1
             else
                um=m+1
             end if
             do l=ll1-1,ll2
                if ((my_pev.eq.0).and.(l.eq.ll1-1)) then
                   ll=l+1
                else
                   ll=l
                end if
                if ((my_pev.eq.(n_procs_v-1)).and.(l.eq.ll2)) then
                   ul=ll2
                else
                   ul=l+1
                end if
                fcd11=(loc_f(i,j,k,ul,m,n)-loc_f(i,j,k,ll,m,n))*dvin
                fcd12=((loc_f(i,j,k,ll,um,n)-loc_f(i,j,k,ll,m,n))*dmuin(m)+ &
                     (loc_f(i,j,k,ll,m,n)-loc_f(i,j,k,ll,lm,n))*dmuin(m-1)+ &
                     (loc_f(i,j,k,ul,um,n)-loc_f(i,j,k,ul,m,n))*dmuin(m)+&
                     (loc_f(i,j,k,ul,m,n)-loc_f(i,j,k,ul,lm,n))*dmuin(m-1))*0.25
                fcr1=(loc_f(i,j,k,ll,m,n)+loc_f(i,j,k,ul,m,n))*0.5
                do sp=0,n_spec_coll-1
                   fluvp=cd11(ii,k,l,m,sp,n)*fcd11+  &
                        cd12(ii,k,l,m,sp,n)*fcd12+ &
                        cr1(ii,k,l,m,sp,n)*fcr1
                   vfcoll(l,m,sp)=vfcoll(l,m,sp)+fluvp*cvf_vp(l)
                   vfcoll(l+1,m,sp)=vfcoll(l+1,m,sp)-fluvp*cvf_vp(l+1)
                enddo
             enddo
          enddo
          do m=lm1-1,lm2
             if ((my_pew.eq.0).and.(m.eq.(lm1-1))) then
                lm=lm1+1
             else
                lm=m
             end if
             if ((my_pew.eq.(n_procs_w-1)).and.(m.eq.lm2)) then
                um=lm2-1
             else
                um=m+1
             end if
             do l=ll1,ll2
                if ((my_pev.eq.0).and.(l.eq.ll1)) then
                   ll=ll1
                else
                   ll=l-1
                end if
                if ((my_pev.eq.(n_procs_v-1)).and.(l.eq.ll2)) then
                   ul=ll2
                else
                   ul=l+1
                end if
                fcd22=(loc_f(i,j,k,l,um,n)-loc_f(i,j,k,l,lm,n))*dmuin(m)
                fcd21=(loc_f(i,j,k,ul,lm,n)-loc_f(i,j,k,ll,lm,n)+ &
                     loc_f(i,j,k,ul,um,n)-loc_f(i,j,k,ll,um,n))*0.25*dvinv(l)
                fcr2=(loc_f(i,j,k,l,lm,n)+loc_f(i,j,k,l,um,n))*0.5
                do sp=0,n_spec_coll-1
                   flumu=cd22(ii,k,l,m,sp,n)*fcd22 +&
                        cd21(ii,k,l,m,sp,n)*fcd21 +&
                        cr2(ii,k,l,m,sp,n)*fcr2
                   vfcoll(l,m,sp)=vfcoll(l,m,sp)+flumu*cvf_mu(m)
                   vfcoll(l,m+1,sp)=vfcoll(l,m+1,sp)-flumu*cvf_mu(m+1)
                enddo
             enddo
          enddo
          !contributions from different s are added
          if (em_conserve) then
             do sp=0,n_spec-1 !conservation terms only for kinetic species
                do m=lm1,lm2
                   do l=ll1,ll2
                      mom_fac(i,j,k,n,sp)=mom_fac(i,j,k,n,sp)-vfcoll(l,m,sp)*loc_mat_10(ii,k,l,m,n)
                      en_fac(i,j,k,n,sp)=en_fac(i,j,k,n,sp)-vfcoll(l,m,sp)*mat_01_20(ii,k,l,m,n)
                   enddo
                enddo
             enddo
          endif
          do sp=0,n_spec_coll-1 !test-particle collisions also for extra (adiabatic ion) species (if any)
             do m=lm1,lm2
                do l=ll1,ll2
                   loc_k(i,j,k,l,m,n)=loc_k(i,j,k,l,m,n)+vfcoll(l,m,sp)
                enddo
             enddo
          enddo
       enddo
#ifdef WITHOMP_COLLISIONS
       !$omp end parallel do
#endif
    enddo
    PERFOFF

  end subroutine add_testpart_1

!---------
!TESTPART2

  subroutine initialize_testpart_2
    allocate(vfcoll_2(li1:li2,lj1:lj2,ll1-1:ll2+1,lm1-1:lm2+1, 0:n_spec_coll-1))
  end subroutine initialize_testpart_2

  !>Adds the test particle contribution and computes the violations of energy and parallel momentum
  !!conservation
  !!
  !!Alternative implementation 2: data is mostly handled in x-y blocks, efficient for bigger problems
  subroutine add_testpart_2(loc_f,loc_k, mom_fac, en_fac)
    !>f_ type array with exchanged ghost cells in vpar and mu directions
    complex,dimension(li1:li2, lj1:lj2, lbz:ubz, lbv:ubv, lbw:ubw, ln1:ln2),intent(in) :: loc_f
    !>g_1 type array to which the result is added
    complex,dimension(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2),intent(inout) :: loc_k
    !>parallel momentum violation induced by test particle contribution
    complex,dimension(li1:li2, lj1:lj2, lk1:lk2, ln1:ln2, 0:n_spec-1),intent(out) :: mom_fac
    !>energy violation induced by test particle contribution
    complex,dimension(li1:li2, lj1:lj2, lk1:lk2, ln1:ln2, 0:n_spec-1),intent(out) :: en_fac

    complex,dimension(li1:li2,lj1:lj2) :: fcd11,fcd12,fcd21,fcd22,fcr1,fcr2
    complex,dimension(li1:li2,lj1:lj2) :: fluvp, flumu

    integer:: i,j,k,l,m,n,sp,ll,ul,lm,um

    PERFON('collt2')
    mom_fac=(0.,0.)
    en_fac=(0.,0.)

    do n=ln1,ln2
       do k=lk1,lk2
          vfcoll_2=(0.,0.)
          PERFON_I('ct2_l1')
          do m=lm1,lm2
             if (m.eq.0) then
                !if ((my_pew.eq.0).and.(m.eq.lm1)) then
                !lm=m+1
                lm=1
             else
                lm=m-1
             end if
             if (m.eq.nw0-1) then
             !if ((my_pew.eq.(n_procs_w-1)).and.(m.eq.lm2)) then
                um=m-1
             else
                um=m+1
             end if
             do l=ll1-1,ll2
                if (l.eq.-1) then
                   !if ((my_pev.eq.0).and.(l.eq.ll1-1)) then
                   !ll=l+1
                   ll=0
                else
                   ll=l
                end if
                if (l.eq.nv0-1) then
                   !if ((my_pev.eq.(n_procs_v-1)).and.(l.eq.ll2)) then
                   ul=ll2
                else
                   ul=l+1
                end if
                fcd11=(loc_f(:,:,k,ul,m,n)-loc_f(:,:,k,ll,m,n))*dvin
                fcr1=(loc_f(:,:,k,ll,m,n)+loc_f(:,:,k,ul,m,n))*0.5
                fcd12=(dmuin(m-1)-dmuin(m))*0.5*fcr1&
                     +loc_f(:,:,k,ll,um,n)*dmuin(m)*0.25&
                     +loc_f(:,:,k,ul,um,n)*dmuin(m)*0.25&
                     -loc_f(:,:,k,ll,lm,n)*dmuin(m-1)*0.25&
                     -loc_f(:,:,k,ul,lm,n)*dmuin(m-1)*0.25
                do sp=0,n_spec_coll-1
                   if (xy_local) then
                      fluvp=cd11(pi1,k,l,m,sp,n)*fcd11+  &
                           cd12(pi1,k,l,m,sp,n)*fcd12+ &
                           cr1(pi1,k,l,m,sp,n)*fcr1
                   else
                      do j=lj1,lj2
                         do i=li1,li2
                            fluvp(i,j)=cd11(i,k,l,m,sp,n)*fcd11(i,j)+  &
                                 cd12(i,k,l,m,sp,n)*fcd12(i,j)+ &
                                 cr1(i,k,l,m,sp,n)*fcr1(i,j)
                         enddo
                      enddo
                   endif
                   vfcoll_2(:,:,l,m,sp)=vfcoll_2(:,:,l,m,sp)+fluvp(li1:li2,lj1:lj2)*cvf_vp(l)
                   vfcoll_2(:,:,l+1,m,sp)=vfcoll_2(:,:,l+1,m,sp)-fluvp(li1:li2,lj1:lj2)*cvf_vp(l+1)
                enddo
             enddo
          enddo
          PERFOFF_I
          PERFON_I('ct2_l2')
          do m=lm1-1,lm2
             if ((my_pew.eq.0).and.(m.eq.(lm1-1))) then
                lm=lm1+1
             else
                lm=m
             end if
             if ((my_pew.eq.(n_procs_w-1)).and.(m.eq.lm2)) then
                um=lm2-1
             else
                um=m+1
             end if
             do l=ll1,ll2
                if ((my_pev.eq.0).and.(l.eq.ll1)) then
                   ll=ll1
                else
                   ll=l-1
                end if
                if ((my_pev.eq.(n_procs_v-1)).and.(l.eq.ll2)) then
                   ul=ll2
                else
                   ul=l+1
                end if
                fcd22=(loc_f(:,:,k,l,um,n)-loc_f(:,:,k,l,lm,n))*dmuin(m)
                fcr2=(loc_f(:,:,k,l,lm,n)+loc_f(:,:,k,l,um,n))*0.5
                fcd21=(loc_f(:,:,k,ul,lm,n)-loc_f(:,:,k,ll,lm,n)+ &
                     loc_f(:,:,k,ul,um,n)-loc_f(:,:,k,ll,um,n))*0.25*dvinv(l)
                do sp=0,n_spec_coll-1
                   if (xy_local) then
                      flumu=cd22(pi1,k,l,m,sp,n)*fcd22 +&
                           cd21(pi1,k,l,m,sp,n)*fcd21 +&
                           cr2(pi1,k,l,m,sp,n)*fcr2
                   else
                      do j=lj1,lj2
                         do i=li1,li2
                            flumu(i,j)=cd22(i,k,l,m,sp,n)*fcd22(i,j) +&
                                 cd21(i,k,l,m,sp,n)*fcd21(i,j) +&
                                 cr2(i,k,l,m,sp,n)*fcr2(i,j)
                         enddo
                      enddo
                   endif
                   vfcoll_2(:,:,l,m,sp)=vfcoll_2(:,:,l,m,sp)+flumu(li1:li2,lj1:lj2)*cvf_mu(m)
                   vfcoll_2(:,:,l,m+1,sp)=vfcoll_2(:,:,l,m+1,sp)-flumu(li1:li2,lj1:lj2)*cvf_mu(m+1)
                enddo
             enddo
          enddo
          PERFOFF_I
          PERFON_I('ct2_l3')
          do sp=0,n_spec_coll-1
             do m=lm1,lm2
                do l=ll1,ll2
                   !contributions from different s are added
                   if (em_conserve.and.(sp.lt.n_spec)) then !conservation for kinetic species only
                      if (xy_local) then
                         mom_fac(:,:,k,n,sp)=mom_fac(:,:,k,n,sp)-vfcoll_2(:,:,l,m,sp)*loc_mat_10(pi1,k,l,m,n)
                         en_fac(:,:,k,n,sp)=en_fac(:,:,k,n,sp)-vfcoll_2(:,:,l,m,sp)*mat_01_20(pi1,k,l,m,n)
                      else
                         do j=lj1,lj2
                            mom_fac(:,j,k,n,sp)=mom_fac(:,j,k,n,sp)-vfcoll_2(:,j,l,m,sp)*loc_mat_10(:,k,l,m,n)
                            en_fac(:,j,k,n,sp)=en_fac(:,j,k,n,sp)-vfcoll_2(:,j,l,m,sp)*mat_01_20(:,k,l,m,n)
                         enddo
                      endif
                   endif
                   !collisional contribution by all species including extra adiabatic ions (if any)
                   loc_k(:,:,k,l,m,n)=loc_k(:,:,k,l,m,n)+vfcoll_2(:,:,l,m,sp)
                enddo
             enddo
          enddo
          PERFOFF_I
       enddo ! k
    enddo

    PERFOFF

  end subroutine add_testpart_2

  subroutine finalize_testpart_2
    deallocate(vfcoll_2)
  end subroutine finalize_testpart_2

!---------
!TESTPART3

  !> Initializes the big arrays needed for the add_testpart_3 implementation
  subroutine initialize_testpart_3
    real,dimension(:,:,:,:,:,:,:),allocatable:: collsten
    integer:: i,j,ij,ii,k,l,m,n,sp,sten, xyvar

    xyvar=1
    if(.not.xy_local) xyvar=xyvar*li0

    allocate(collsten(0:xyvar-1,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2,0:n_spec_coll-1,9))

    shiftvec=(/-lv0-1, -lv0, -lv0+1, &
         -1, 0, +1, &
         +lv0-1, +lv0, +lv0+1/)

    stenlb=lv0*nwb+nvb
    stenub=lv0*lw0-stenlb-1

    !boundaries in mu and v_\parallel
    if(my_pev.eq.0) stenlb=stenlb+(/1,0,0,1,0,0,1,0,0/)
    if(my_pev.eq.(n_procs_v-1)) stenub=stenub+(/0,0,-1,0,0,-1,0,0,-1/)

    if(my_pew.eq.0) stenlb=stenlb+(/lv0,lv0,lv0,0,0,0,0,0,0/)
    if(my_pew.eq.(n_procs_w-1)) stenub=stenub+(/0,0,0,0,0,0,-lv0,-lv0,-lv0/)

    collsten=0.

    do k=lk1,lk2
       do m=lm1,lm2
          do l=ll1,ll2
             do n=ln1,ln2
                do sp=0,n_spec_coll-1
                   do j=pj1,pj2
                      do i=pi1,pi2
                         ij=li0*(j-pj1)+(i-pi1)
                         !this part is for the computation of the contribution due to differences in the v parallel flux
                         collsten(ij,k,l,m,n,sp,1)=cvf_vp(l)*cd12(i,k,l-1,m,sp,n)*0.25*dmuin(m-1)
                         collsten(ij,k,l,m,n,sp,2)=cvf_vp(l)*(cd12(i,k,l-1,m,sp,n)*0.25*dmuin(m-1)&
                              -cd12(i,k,l,m,sp,n)*0.25*dmuin(m-1))
                         collsten(ij,k,l,m,n,sp,3)=-cvf_vp(l)*cd12(i,k,l,m,sp,n)*0.25*dmuin(m-1)
                         collsten(ij,k,l,m,n,sp,4)=-cvf_vp(l)*(cr1(i,k,l-1,m,sp,n)*0.5 - cd11(i,k,l-1,m,sp,n)*dvin &
                              +cd12(i,k,l-1,m,sp,n)*0.25*dmuin(m-1) - cd12(i,k,l-1,m,sp,n)*0.25*dmuin(m))
                         collsten(ij,k,l,m,n,sp,5)=cvf_vp(l)*(cr1(i,k,l,m,sp,n)*0.5 - cd11(i,k,l,m,sp,n)*dvin &
                              + cd12(i,k,l,m,sp,n)*0.25*dmuin(m-1) - cd12(i,k,l,m,sp,n)*0.25*dmuin(m) &
                              - cd11(i,k,l-1,m,sp,n)*dvin + cd12(i,k,l-1,m,sp,n)*0.25*dmuin(m) &
                              - cd12(i,k,l-1,m,sp,n)*0.25*dmuin(m-1) - cr1(i,k,l-1,m,sp,n)*0.5)
                         collsten(ij,k,l,m,n,sp,6)=cvf_vp(l)*(cd11(i,k,l,m,sp,n)*dvin - cd12(i,k,l,m,sp,n)*0.25*dmuin(m) &
                              + cd12(i,k,l,m,sp,n)*0.25*dmuin(m-1) + cr1(i,k,l,m,sp,n)*0.5)
                         collsten(ij,k,l,m,n,sp,7)=-cvf_vp(l)*cd12(i,k,l-1,m,sp,n)*0.25*dmuin(m)
                         collsten(ij,k,l,m,n,sp,8)=cvf_vp(l)*(cd12(i,k,l,m,sp,n)*0.25*dmuin(m)&
                              - cd12(i,k,l-1,m,sp,n)*0.25*dmuin(m))
                         collsten(ij,k,l,m,n,sp,9)=cvf_vp(l)*cd12(i,k,l,m,sp,n)*0.25*dmuin(m)
                         !this part is for the computation of the contribution due to differences in the mu flux
                         collsten(ij,k,l,m,n,sp,1)= collsten(ij,k,l,m,n,sp,1)+ cvf_mu(m)*cd21(i,k,l,m-1,sp,n)*0.25*dvinv(l)
                         collsten(ij,k,l,m,n,sp,2)= collsten(ij,k,l,m,n,sp,2)+ cvf_mu(m)*(cd22(i,k,l,m-1,sp,n)*dmuin(m-1) &
                              - cr2(i,k,l,m-1,sp,n)*0.5)
                         collsten(ij,k,l,m,n,sp,3)= collsten(ij,k,l,m,n,sp,3)- cvf_mu(m)*cd21(i,k,l,m-1,sp,n)*0.25*dvinv(l)
                         collsten(ij,k,l,m,n,sp,4)= collsten(ij,k,l,m,n,sp,4)+ cvf_mu(m)*(cd21(i,k,l,m-1,sp,n)*0.25*dvinv(l) &
                              - cd21(i,k,l,m,sp,n)*0.25*dvinv(l))
                         collsten(ij,k,l,m,n,sp,5)= collsten(ij,k,l,m,n,sp,5)+ cvf_mu(m)*(- cd22(i,k,l,m,sp,n)*dmuin(m) &
                              + cr2(i,k,l,m,sp,n)*0.5 &
                              - cd22(i,k,l,m-1,sp,n)*dmuin(m-1) - cr2(i,k,l,m-1,sp,n)*0.5)
                         collsten(ij,k,l,m,n,sp,6)= collsten(ij,k,l,m,n,sp,6)+ cvf_mu(m)*(cd21(i,k,l,m,sp,n)*0.25*dvinv(l) &
                              - cd21(i,k,l,m-1,sp,n)*0.25*dvinv(l))
                         collsten(ij,k,l,m,n,sp,7)= collsten(ij,k,l,m,n,sp,7)- cvf_mu(m)*cd21(i,k,l,m,sp,n)*0.25*dvinv(l)
                         collsten(ij,k,l,m,n,sp,8)= collsten(ij,k,l,m,n,sp,8)+ cvf_mu(m)*(cd22(i,k,l,m,sp,n)*dmuin(m) &
                              + cr2(i,k,l,m,sp,n)*0.5)
                         collsten(ij,k,l,m,n,sp,9)= collsten(ij,k,l,m,n,sp,9)+ cvf_mu(m)*cd21(i,k,l,m,sp,n)*0.25*dvinv(l)
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo

    !include v-space boundary conditions

#if 1
    do sten=1,9
       do sp=0,n_spec_coll-1
          do n=ln1,ln2
             select case (sten)
             case (1:3)
                if(my_pew.eq.0) then
                   collsten(:,:,:,lm1,n,sp,sten+6)=collsten(:,:,:,lm1,n,sp,sten+6)+collsten(:,:,:,lm1,n,sp,sten)
                   collsten(:,:,:,lm1,n,sp,sten)=0.
                end if
             case(7:9)
                if(my_pew.eq.(n_procs_w-1)) then
                   collsten(:,:,:,lm2,n,sp,sten-6)=collsten(:,:,:,lm2,n,sp,sten-6)+collsten(:,:,:,lm2,n,sp,sten)
                   collsten(:,:,:,lm2,n,sp,sten)=0.
                end if
             end select
             select case (sten)
             case (1,4,7)
                do m=lm1,lm2
                   if(my_pev.eq.0) then
                      collsten(:,:,ll1,m,n,sp,sten+1)=collsten(:,:,ll1,m,n,sp,sten+1)+collsten(:,:,ll1,m,n,sp,sten)
                      collsten(:,:,ll1,m,n,sp,sten)=0.
                   end if
                end do
             case(3,6,9)
                do m=lm1,lm2
                   if(my_pev.eq.(n_procs_v-1)) then
                      collsten(:,:,ll2,m,n,sp,sten-1)=collsten(:,:,ll2,m,n,sp,sten-1)+collsten(:,:,ll2,m,n,sp,sten)
                      collsten(:,:,ll2,m,n,sp,sten)=0.
                   end if
                end do
             end select
          enddo
       enddo
    enddo
#endif

    allocate(collsten_big(0:li0*lj0-1,lbz:ubz,lbv:ubv,lbw:ubw,9,ln1:ln2,0:n_spec_coll-1))
    collsten_big=0.
    do sp=0,n_spec_coll-1
       do n=ln1,ln2
          do sten=1,9
             do m=lbw,ubw
                do l=lbv,ubv
                   do k=lk1,lk2
                      do ij=0,lij0-1
                         if (xy_local) then
                            ii=0
                         else
                            ii=Modulo(ij,li0)
                         endif
                         collsten_big(ij,k,l,m,sten,n,sp)=collsten(ii,k,l,m,n,sp,sten)
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo

    deallocate(collsten)

    allocate(vfcoll_big(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw, 0:n_spec_coll-1))


  end subroutine initialize_testpart_3

  !>Adds the test particle contribution and computes the violations of energy and parallel momentum
  !!conservation
  !!
  !!Alternative implementation 3: a large array is precomputed, and only big blocks of data are
  !!multiplied - very efficient for a small problem size to cache ratio
  subroutine add_testpart_3(loc_f,loc_k,mom_fac,en_fac)
    !>f_ type array with exchanged ghost cells in vpar and mu directions
    complex,dimension(li1:li2, lj1:lj2, lbz:ubz, lbv:ubv, lbw:ubw, ln1:ln2),intent(in) :: loc_f
    !>g_1 type array to which the result is added
    complex,dimension(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2),intent(inout) :: loc_k
    !>parallel momentum violation induced by test particle contribution
    complex,dimension(li1:li2, lj1:lj2, lk1:lk2, ln1:ln2, 0:n_spec-1),intent(out) :: mom_fac
    !>energy violation induced by test particle contribution
    complex,dimension(li1:li2, lj1:lj2, lk1:lk2, ln1:ln2, 0:n_spec-1),intent(out) :: en_fac
    integer:: i,j,k,l,m,n,sp

    PERFON('collt3')
    mom_fac=(0.,0.)
    en_fac=(0.,0.)
    do n=ln1,ln2
       call teststen2(n,loc_f,collsten_big,vfcoll_big)
       if (em_conserve) then !contributions from different s are added
          PERFON('cal_vio')
          do sp=0,n_spec-1
             do k=lk1,lk2
                do m=lm1,lm2
                   do l=ll1,ll2
                      if (xy_local) then
                         mom_fac(:,:,k,n,sp)=mom_fac(:,:,k,n,sp)-vfcoll_big(:,:,k,l,m,sp)*loc_mat_10(pi1,k,l,m,n)
                         en_fac(:,:,k,n,sp)=en_fac(:,:,k,n,sp)-vfcoll_big(:,:,k,l,m,sp)*mat_01_20(pi1,k,l,m,n)
                      else
                         do j=lj1,lj2
                            mom_fac(:,j,k,n,sp)=mom_fac(:,j,k,n,sp)-vfcoll_big(:,j,k,l,m,sp)*loc_mat_10(:,k,l,m,n)
                            en_fac(:,j,k,n,sp)=en_fac(:,j,k,n,sp)-vfcoll_big(:,j,k,l,m,sp)*mat_01_20(:,k,l,m,n)
                         enddo
                      endif
                   enddo
                enddo
             enddo
          enddo
          PERFOFF
       endif
       do sp=0,n_spec_coll-1
          loc_k(:,:,:,:,:,n)=loc_k(:,:,:,:,:,n)+vfcoll_big(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,sp)
       enddo
    enddo
    PERFOFF

  end subroutine add_testpart_3

  !>Computes the test particle contribution in add_testpart_3 using a big precomputed matrix
  subroutine teststen2(n,loc_f_,collsten_gr,vfcoll_)
    complex,dimension(0:lij0*lz0-1, 0:lv0*lw0-1, ln1:ln2),intent(in) :: loc_f_ !<f_ type array with exchanged ghost cells
    real, dimension(0:lij0*lz0-1,0:lv0*lw0-1,9,ln1:ln2,0:n_spec_coll-1),intent(in) ::collsten_gr !<precomputed array (stencil)
    integer,intent(in):: n !<species index
    complex,dimension(0:lij0*lz0-1,0:lv0*lw0-1, 0:n_spec_coll-1),intent(out) :: vfcoll_ !<f_1 type output array
    integer:: sp

#if 1
    integer:: shift,sten

    vfcoll_=0.
    do sp=0,n_spec_coll-1
       !stencil index 5 corresponds to zero shift; the result has therefore the same dimensions as vfcoll_
       !and can be used to initialize vfcoll_
       sten=5
       shift=shiftvec(sten)
       vfcoll_(:,stenlb(sten):stenub(sten),sp)=&
            collsten_gr(:,stenlb(sten):stenub(sten),sten,n,sp)*loc_f_(:,stenlb(sten)+shift:stenub(sten)+shift,n)
       do sten=1,9
          if (sten.ne.5) then !already taken into account above
             shift=shiftvec(sten)
             vfcoll_(:,stenlb(sten):stenub(sten),sp)=vfcoll_(:,stenlb(sten):stenub(sten),sp)&
                  +collsten_gr(:,stenlb(sten):stenub(sten),sten,n,sp)*loc_f_(:,stenlb(sten)+shift:stenub(sten)+shift,n)
          endif
       enddo
    enddo
#else
    !only works with v-space boundaries
    do sp=0,n_spec_coll-1
       vfcoll_(:,stenlb:stenub,sp)=&
            collsten_gr(:,stenlb:stenub,1,n,sp)*loc_f_(:,stenlb-lv0-1:stenub-lv0-1,n)+&
            collsten_gr(:,stenlb:stenub,2,n,sp)*loc_f_(:,stenlb-lv0  :stenub-lv0  ,n)+&
            collsten_gr(:,stenlb:stenub,3,n,sp)*loc_f_(:,stenlb-lv0+1:stenub-lv0+1,n)+&
            collsten_gr(:,stenlb:stenub,4,n,sp)*loc_f_(:,stenlb-1    :stenub-1    ,n)+&
            collsten_gr(:,stenlb:stenub,5,n,sp)*loc_f_(:,stenlb      :stenub      ,n)+&
            collsten_gr(:,stenlb:stenub,6,n,sp)*loc_f_(:,stenlb+1    :stenub+1    ,n)+&
            collsten_gr(:,stenlb:stenub,7,n,sp)*loc_f_(:,stenlb+lv0-1:stenub+lv0-1,n)+&
            collsten_gr(:,stenlb:stenub,8,n,sp)*loc_f_(:,stenlb+lv0  :stenub+lv0  ,n)+&
            collsten_gr(:,stenlb:stenub,9,n,sp)*loc_f_(:,stenlb+lv0+1:stenub+lv0+1,n)
    enddo
#endif

  end subroutine teststen2

  subroutine finalize_testpart_3
    deallocate(collsten_big,vfcoll_big)
  end subroutine finalize_testpart_3


  !>Add some description here
  subroutine add_testdifferentialpart(loc_f,loc_k)
    !>f_ type array with exchanged ghost cells in vpar and mu directions
    complex,dimension(li1:li2, lj1:lj2, lbz:ubz, lbv:ubv, lbw:ubw, ln1:ln2),intent(in) :: loc_f
    !>g_1 type array to which the result is added
    complex,dimension(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2),intent(inout) :: loc_k
    integer:: l,m,n,k,j,sp

    if (xy_local) then

      ! if (coll_order .eq. 'fourth') then
#ifndef COLLSIX

          do n=ln1,ln2
             do m=lm1,lm2
                do l=ll1,ll2
                   do k=lk1,lk2

                      loc_k(:,:,k,l,m,n)=loc_k(:,:,k,l,m,n)+T_coeff1(pi1,k,l,m,n)&
                           *(-loc_f(:,:,k,l+2,m,n)+16*(loc_f(:,:,k,l+1,m,n)+loc_f(:,:,k,l-1,m,n))&
                           -loc_f(:,:,k,l-2,m,n))&
                           +T_coeff2(pi1,k,l,m,n)*(-loc_f(:,:,k,l,m+2,n)+16*(loc_f(:,:,k,l,m+1,n)&
                           +loc_f(:,:,k,l,m-1,n))-loc_f(:,:,k,l,m-2,n))&
                           +T_coeff3(pi1,k,l,m,n)*(loc_f(:,:,k,l+2,m+2,n)+loc_f(:,:,k,l-2,m-2,n)&
                           -loc_f(:,:,k,l+2,m-2,n)-loc_f(:,:,k,l-2,m+2,n)&
                           +8*(loc_f(:,:,k,l+2,m-1,n)&
                           -loc_f(:,:,k,l+2,m+1,n)+loc_f(:,:,k,l+1,m-2,n)-loc_f(:,:,k,l+1,m+2,n)&
                           +loc_f(:,:,k,l-1,m+2,n)-loc_f(:,:,k,l-1,m-2,n)+loc_f(:,:,k,l-2,m+1,n)&
                           -loc_f(:,:,k,l-2,m-1,n)&
                           +8*(loc_f(:,:,k,l+1,m+1,n)+loc_f(:,:,k,l-1,m-1,n)-loc_f(:,:,k,l+1,m-1,n)&
                           -loc_f(:,:,k,l-1,m+1,n))))&
                           +T_coeff4(pi1,k,l,m,n)*(-loc_f(:,:,k,l+2,m,n)+8*(loc_f(:,:,k,l+1,m,n)&
                           -loc_f(:,:,k,l-1,m,n))+loc_f(:,:,k,l-2,m,n))&
                           +T_coeff5(pi1,k,l,m,n)*(-loc_f(:,:,k,l,m+2,n)+8*(loc_f(:,:,k,l,m+1,n)&
                           -loc_f(:,:,k,l,m-1,n))+loc_f(:,:,k,l,m-2,n))&
                           +T_coeff6(pi1,k,l,m,n)*loc_f(:,:,k,l,m,n)

                   enddo
                enddo
             enddo
          enddo

     !  elseif (coll_order .eq. 'order6') then
#else

          do n=ln1,ln2
             do m=lm1,lm2
                do l=ll1,ll2
                   do k=lk1,lk2

                      loc_k(:,:,k,l,m,n)=loc_k(:,:,k,l,m,n)+T_coeff1(pi1,k,l,m,n)&
                           *(loc_f(:,:,k,l+3,m,n)-13.5*(loc_f(:,:,k,l+2,m,n)+loc_f(:,:,k,l-2,m,n))&
                           +135*(loc_f(:,:,k,l+1,m,n)+loc_f(:,:,k,l-1,m,n))&
                           +loc_f(:,:,k,l-3,m,n))&
                           +T_coeff2(pi1,k,l,m,n)&
                           *(loc_f(:,:,k,l,m+3,n)-13.5*(loc_f(:,:,k,l,m+2,n)+loc_f(:,:,k,l,m-2,n))&
                           +135*(loc_f(:,:,k,l,m+1,n)+loc_f(:,:,k,l,m-1,n))&
                           +loc_f(:,:,k,l,m-3,n))&
                           +T_coeff3(pi1,k,l,m,n)&
                           *(loc_f(:,:,k,l+3,m+3,n)-loc_f(:,:,k,l+3,m-3,n)-loc_f(:,:,k,l-3,m+3,n)&
                           +loc_f(:,:,k,l-3,m-3,n)&
                           +9*(loc_f(:,:,k,l+2,m-3,n)+loc_f(:,:,k,l+3,m-2,n)-loc_f(:,:,k,l+3,m+2,n)&
                           -loc_f(:,:,k,l+2,m+3,n)+loc_f(:,:,k,l-2,m+3,n)+loc_f(:,:,k,l-3,m+2,n)&
                           -loc_f(:,:,k,l-3,m-2,n)-loc_f(:,:,k,l-2,m-3,n)&
                           +9*(loc_f(:,:,k,l+2,m+2,n)+loc_f(:,:,k,l-2,m-2,n)-loc_f(:,:,k,l+2,m-2,n)&
                           -loc_f(:,:,k,l-2,m+2,n)))&
                           +45*(loc_f(:,:,k,l+3,m+1,n)-loc_f(:,:,k,l+3,m-1,n)-loc_f(:,:,k,l+1,m-3,n)&
                           +loc_f(:,:,k,l-1,m-3,n)+loc_f(:,:,k,l-3,m-1,n)-loc_f(:,:,k,l-3,m+1,n)&
                           -loc_f(:,:,k,l-1,m+3,n)+loc_f(:,:,k,l+1,m+3,n)&
                           +45*(loc_f(:,:,k,l+1,m+1,n)+loc_f(:,:,k,l-1,m-1,n)-loc_f(:,:,k,l+1,m-1,n)&
                           -loc_f(:,:,k,l-1,m+1,n)))&
                           +405*(loc_f(:,:,k,l+2,m-1,n)+loc_f(:,:,k,l+1,m-2,n)-loc_f(:,:,k,l-1,m-2,n)&
                           -loc_f(:,:,k,l-2,m-1,n)+loc_f(:,:,k,l-2,m+1,n)+loc_f(:,:,k,l-1,m+2,n)&
                           -loc_f(:,:,k,l+1,m+2,n)-loc_f(:,:,k,l+2,m+1,n)))&
                           +T_coeff4(pi1,k,l,m,n)&
                           *(loc_f(:,:,k,l+3,m,n)-9*(loc_f(:,:,k,l+2,m,n)-loc_f(:,:,k,l-2,m,n))&
                           +45*(loc_f(:,:,k,l+1,m,n)-loc_f(:,:,k,l-1,m,n))&
                           -loc_f(:,:,k,l-3,m,n))&
                           +T_coeff5(pi1,k,l,m,n)&
                           *(loc_f(:,:,k,l,m+3,n)-9*(loc_f(:,:,k,l,m+2,n)-loc_f(:,:,k,l,m-2,n))&
                           +45*(loc_f(:,:,k,l,m+1,n)-loc_f(:,:,k,l,m-1,n))&
                           -loc_f(:,:,k,l,m-3,n))&
                           +T_coeff6(pi1,k,l,m,n)*loc_f(:,:,k,l,m,n)


                   enddo
                enddo
             enddo
          enddo

    !   endif
#endif

       particle_fac = (0.0,0.0)

       do n=ln1,ln2
          do sp=0,n_spec-1
             do m=lm1,lm2
                do l=ll1,ll2
                   do k=lk1,lk2

                        particle_fac(:,:,k,sp,n)=particle_fac(:,:,k,sp,n)+particle_diffpart(pi1,k,l,m,sp,n)*loc_f(:,:,k,l,m,n)

                   enddo
                enddo
             enddo
          enddo
       enddo

       call my_complex_sum_vw(particle_fac,size(particle_fac))

       do n=ln1,ln2
          do m=lm1,lm2
             do l=ll1,ll2
                do k=lk1,lk2
                   particle_buff=particle_fac(:,:,k,0,n)
                   do sp=1,n_spec-1
                      particle_buff(:,:)=particle_buff(:,:)+particle_fac(:,:,k,sp,n)
                   enddo
                   loc_k(:,:,k,l,m,n)=loc_k(:,:,k,l,m,n)&
                        -particle_coeff(pi1,k,l,m,n)*particle_buff(:,:)
                enddo
             enddo
          enddo
       enddo

    else !.not.xy_local

     !  if (coll_order .eq. 'fourth') then
#ifndef COLLSIX

          do n=ln1,ln2
             do m=lm1,lm2
                do l=ll1,ll2
                   do k=lk1,lk2
                      do j=lj1,lj2


                         loc_k(:,j,k,l,m,n)=loc_k(:,j,k,l,m,n)+T_coeff1(:,k,l,m,n)&
                              *(-loc_f(:,j,k,l+2,m,n)+16*(loc_f(:,j,k,l+1,m,n)+loc_f(:,j,k,l-1,m,n))&
                              -loc_f(:,j,k,l-2,m,n))&
                              +T_coeff2(:,k,l,m,n)*(-loc_f(:,j,k,l,m+2,n)+16*(loc_f(:,j,k,l,m+1,n)&
                              +loc_f(:,j,k,l,m-1,n))-loc_f(:,j,k,l,m-2,n))&
                              +T_coeff3(:,k,l,m,n)*(loc_f(:,j,k,l+2,m+2,n)+loc_f(:,j,k,l-2,m-2,n)&
                              -loc_f(:,j,k,l+2,m-2,n)-loc_f(:,j,k,l-2,m+2,n)&
                              +8*(loc_f(:,j,k,l+2,m-1,n)&
                              -loc_f(:,j,k,l+2,m+1,n)+loc_f(:,j,k,l+1,m-2,n)-loc_f(:,j,k,l+1,m+2,n)&
                              +loc_f(:,j,k,l-1,m+2,n)-loc_f(:,j,k,l-1,m-2,n)+loc_f(:,j,k,l-2,m+1,n)&
                              -loc_f(:,j,k,l-2,m-1,n)&
                              +8*(loc_f(:,j,k,l+1,m+1,n)+loc_f(:,j,k,l-1,m-1,n)-loc_f(:,j,k,l+1,m-1,n)&
                              -loc_f(:,j,k,l-1,m+1,n))))&
                              +T_coeff4(:,k,l,m,n)*(-loc_f(:,j,k,l+2,m,n)+8*(loc_f(:,j,k,l+1,m,n)&
                              -loc_f(:,j,k,l-1,m,n))+loc_f(:,j,k,l-2,m,n))&
                              +T_coeff5(:,k,l,m,n)*(-loc_f(:,j,k,l,m+2,n)+8*(loc_f(:,j,k,l,m+1,n)&
                              -loc_f(:,j,k,l,m-1,n))+loc_f(:,j,k,l,m-2,n))&
                              +T_coeff6(:,k,l,m,n)*loc_f(:,j,k,l,m,n)


                      enddo
                   enddo
                enddo
             enddo
          enddo


     !  elseif (coll_order .eq. 'order6') then
#else

          do n=ln1,ln2
             do m=lm1,lm2
                do l=ll1,ll2
                   do k=lk1,lk2
                      do j=lj1,lj2


                         loc_k(:,j,k,l,m,n)=loc_k(:,j,k,l,m,n)+T_coeff1(:,k,l,m,n)&
                              *(loc_f(:,j,k,l+3,m,n)-13.5*(loc_f(:,j,k,l+2,m,n)+loc_f(:,j,k,l-2,m,n))&
                              +135*(loc_f(:,j,k,l+1,m,n)+loc_f(:,j,k,l-1,m,n))&
                              +loc_f(:,j,k,l-3,m,n))&
                              +T_coeff2(:,k,l,m,n)&
                              *(loc_f(:,j,k,l,m+3,n)-13.5*(loc_f(:,j,k,l,m+2,n)+loc_f(:,j,k,l,m-2,n))&
                              +135*(loc_f(:,j,k,l,m+1,n)+loc_f(:,j,k,l,m-1,n))&
                              +loc_f(:,j,k,l,m-3,n))&
                              +T_coeff3(:,k,l,m,n)&
                              *(loc_f(:,j,k,l+3,m+3,n)-loc_f(:,j,k,l+3,m-3,n)-loc_f(:,j,k,l-3,m+3,n)&
                              +loc_f(:,j,k,l-3,m-3,n)&
                              +9*(loc_f(:,j,k,l+2,m-3,n)+loc_f(:,j,k,l+3,m-2,n)-loc_f(:,j,k,l+3,m+2,n)&
                              -loc_f(:,j,k,l+2,m+3,n)+loc_f(:,j,k,l-2,m+3,n)+loc_f(:,j,k,l-3,m+2,n)&
                              -loc_f(:,j,k,l-3,m-2,n)-loc_f(:,j,k,l-2,m-3,n)&
                              +9*(loc_f(:,j,k,l+2,m+2,n)+loc_f(:,j,k,l-2,m-2,n)-loc_f(:,j,k,l+2,m-2,n)&
                              -loc_f(:,j,k,l-2,m+2,n)))&
                              +45*(loc_f(:,j,k,l+3,m+1,n)-loc_f(:,j,k,l+3,m-1,n)-loc_f(:,j,k,l+1,m-3,n)&
                              +loc_f(:,j,k,l-1,m-3,n)+loc_f(:,j,k,l-3,m-1,n)-loc_f(:,j,k,l-3,m+1,n)&
                              -loc_f(:,j,k,l-1,m+3,n)+loc_f(:,j,k,l+1,m+3,n)&
                              +45*(loc_f(:,j,k,l+1,m+1,n)+loc_f(:,j,k,l-1,m-1,n)-loc_f(:,j,k,l+1,m-1,n)&
                              -loc_f(:,j,k,l-1,m+1,n)))&
                              +405*(loc_f(:,j,k,l+2,m-1,n)+loc_f(:,j,k,l+1,m-2,n)-loc_f(:,j,k,l-1,m-2,n)&
                              -loc_f(:,j,k,l-2,m-1,n)+loc_f(:,j,k,l-2,m+1,n)+loc_f(:,j,k,l-1,m+2,n)&
                              -loc_f(:,j,k,l+1,m+2,n)-loc_f(:,j,k,l+2,m+1,n)))&
                              +T_coeff4(:,k,l,m,n)&
                              *(loc_f(:,j,k,l+3,m,n)-9*(loc_f(:,j,k,l+2,m,n)-loc_f(:,j,k,l-2,m,n))&
                              +45*(loc_f(:,j,k,l+1,m,n)-loc_f(:,j,k,l-1,m,n))&
                              -loc_f(:,j,k,l-3,m,n))&
                              +T_coeff5(:,k,l,m,n)&
                              *(loc_f(:,j,k,l,m+3,n)-9*(loc_f(:,j,k,l,m+2,n)-loc_f(:,j,k,l,m-2,n))&
                              +45*(loc_f(:,j,k,l,m+1,n)-loc_f(:,j,k,l,m-1,n))&
                              -loc_f(:,j,k,l,m-3,n))&
                              +T_coeff6(:,k,l,m,n)*loc_f(:,j,k,l,m,n)

                      enddo
                   enddo
                enddo
             enddo
          enddo

    !   endif
#endif

       particle_fac = (0.0,0.0)

       do n=ln1,ln2
          do sp=0,n_spec-1
             do m=lm1,lm2
                do l=ll1,ll2
                   do k=lk1,lk2
                      do j=lj1,lj2

                         particle_fac(:,j,k,sp,n)=particle_fac(:,j,k,sp,n)+particle_diffpart(:,k,l,m,sp,n)*loc_f(:,j,k,l,m,n)

                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo

       call my_complex_sum_vw(particle_fac,size(particle_fac))

       do n=ln1,ln2
          do m=lm1,lm2
             do l=ll1,ll2
                do k=lk1,lk2
                   do j=lj1,lj2
                      particle_buff(:,j)=particle_fac(:,j,k,0,n)
                      do sp=1,n_spec-1
                         particle_buff(:,j)=particle_buff(:,j)+particle_fac(:,j,k,sp,n)
                      enddo

                      loc_k(:,j,k,l,m,n)=loc_k(:,j,k,l,m,n)&
                           -particle_coeff(:,k,l,m,n)*particle_buff(:,j)

                   enddo
                enddo
             enddo
          enddo
       enddo

    endif

  end subroutine add_testdifferentialpart

  !>Same as testdifferentialpart_1 but for n_procs_v=1 without vpar ghost cells
  subroutine add_testdifferentialpart_alt(loc_f,loc_k)
    !>f_ type array with exchanged ghost cells in vpar and mu directions
    complex,dimension(li1:li2, lj1:lj2, lbz:ubz, lbv:ubv, lbw:ubw, ln1:ln2),intent(in) :: loc_f
    !>g_1 type array to which the result is added
    complex,dimension(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2),intent(inout) :: loc_k
    integer:: l,m,n,k,j,sp

    if (xy_local) then

     !  if (coll_order .eq. 'fourth') then
#ifndef COLLSIX


          do n=ln1,ln2
             do m=lm1,lm2
                do l=ll1,ll2
                   do k=lk1,lk2

                      loc_k(:,:,k,l,m,n)=loc_k(:,:,k,l,m,n)+T_coeff2(pi1,k,l,m,n)&
                           *(-loc_f(:,:,k,l,m+2,n)+16*(loc_f(:,:,k,l,m+1,n)&
                           +loc_f(:,:,k,l,m-1,n))-loc_f(:,:,k,l,m-2,n))&
                           +T_coeff5(pi1,k,l,m,n)*(-loc_f(:,:,k,l,m+2,n)+8*(loc_f(:,:,k,l,m+1,n)&
                           -loc_f(:,:,k,l,m-1,n))+loc_f(:,:,k,l,m-2,n))&
                           +T_coeff6(pi1,k,l,m,n)*loc_f(:,:,k,l,m,n)

                   enddo
                enddo
             enddo
          enddo


          do n=ln1,ln2
             do m=lm1,lm2
                do k=lk1,lk2

                   loc_k(:,:,k,0,m,n)=loc_k(:,:,k,0,m,n)+T_coeff1(pi1,k,0,m,n)&
                        *(-loc_f(:,:,k,2,m,n)+16*loc_f(:,:,k,1,m,n))&
                        +T_coeff3(pi1,k,0,m,n)*(loc_f(:,:,k,2,m+2,n)&
                        -loc_f(:,:,k,2,m-2,n)&
                        +8*(loc_f(:,:,k,2,m-1,n)&
                        -loc_f(:,:,k,2,m+1,n)+loc_f(:,:,k,1,m-2,n)-loc_f(:,:,k,1,m+2,n)&
                        +8*(loc_f(:,:,k,1,m+1,n)-loc_f(:,:,k,1,m-1,n))))&
                        +T_coeff4(pi1,k,0,m,n)*(-loc_f(:,:,k,2,m,n)+8*loc_f(:,:,k,1,m,n))

                   loc_k(:,:,k,1,m,n)=loc_k(:,:,k,1,m,n)+T_coeff1(pi1,k,1,m,n)&
                        *(-loc_f(:,:,k,3,m,n)+16*(loc_f(:,:,k,2,m,n)+loc_f(:,:,k,0,m,n)))&
                        +T_coeff3(pi1,k,1,m,n)*(loc_f(:,:,k,3,m+2,n)&
                        -loc_f(:,:,k,3,m-2,n)&
                        +8*(loc_f(:,:,k,3,m-1,n)&
                        -loc_f(:,:,k,3,m+1,n)+loc_f(:,:,k,2,m-2,n)-loc_f(:,:,k,2,m+2,n)&
                        +loc_f(:,:,k,0,m+2,n)-loc_f(:,:,k,0,m-2,n)&
                        +8*(loc_f(:,:,k,2,m+1,n)+loc_f(:,:,k,0,m-1,n)-loc_f(:,:,k,2,m-1,n)&
                        -loc_f(:,:,k,0,m+1,n))))&
                        +T_coeff4(pi1,k,1,m,n)*(-loc_f(:,:,k,3,m,n)+8*(loc_f(:,:,k,2,m,n)&
                        -loc_f(:,:,k,0,m,n)))

                enddo
                do l=2,nv0-3
                   do k=lk1,lk2

                      loc_k(:,:,k,l,m,n)=loc_k(:,:,k,l,m,n)+T_coeff1(pi1,k,l,m,n)&
                           *(-loc_f(:,:,k,l+2,m,n)+16*(loc_f(:,:,k,l+1,m,n)+loc_f(:,:,k,l-1,m,n))&
                           -loc_f(:,:,k,l-2,m,n))&
                           +T_coeff3(pi1,k,l,m,n)*(loc_f(:,:,k,l+2,m+2,n)+loc_f(:,:,k,l-2,m-2,n)&
                           -loc_f(:,:,k,l+2,m-2,n)-loc_f(:,:,k,l-2,m+2,n)&
                           +8*(loc_f(:,:,k,l+2,m-1,n)&
                           -loc_f(:,:,k,l+2,m+1,n)+loc_f(:,:,k,l+1,m-2,n)-loc_f(:,:,k,l+1,m+2,n)&
                           +loc_f(:,:,k,l-1,m+2,n)-loc_f(:,:,k,l-1,m-2,n)+loc_f(:,:,k,l-2,m+1,n)&
                           -loc_f(:,:,k,l-2,m-1,n)&
                           +8*(loc_f(:,:,k,l+1,m+1,n)+loc_f(:,:,k,l-1,m-1,n)-loc_f(:,:,k,l+1,m-1,n)&
                           -loc_f(:,:,k,l-1,m+1,n))))&
                           +T_coeff4(pi1,k,l,m,n)*(-loc_f(:,:,k,l+2,m,n)+8*(loc_f(:,:,k,l+1,m,n)&
                           -loc_f(:,:,k,l-1,m,n))+loc_f(:,:,k,l-2,m,n))

                   enddo
                enddo
                do k=lk1,lk2

                   loc_k(:,:,k,nv0-2,m,n)=loc_k(:,:,k,nv0-2,m,n)+T_coeff1(pi1,k,nv0-2,m,n)&
                        *(16*(loc_f(:,:,k,nv0-1,m,n)+loc_f(:,:,k,nv0-3,m,n))&
                        -loc_f(:,:,k,nv0-4,m,n))&
                        +T_coeff3(pi1,k,nv0-2,m,n)*(loc_f(:,:,k,nv0-4,m-2,n)&
                        -loc_f(:,:,k,nv0-4,m+2,n)&
                        +8*(loc_f(:,:,k,nv0-1,m-2,n)-loc_f(:,:,k,nv0-1,m+2,n)&
                        +loc_f(:,:,k,nv0-3,m+2,n)-loc_f(:,:,k,nv0-3,m-2,n)+loc_f(:,:,k,nv0-4,m+1,n)&
                        -loc_f(:,:,k,nv0-4,m-1,n)&
                        +8*(loc_f(:,:,k,nv0-1,m+1,n)+loc_f(:,:,k,nv0-3,m-1,n)-loc_f(:,:,k,nv0-1,m-1,n)&
                        -loc_f(:,:,k,nv0-3,m+1,n))))&
                        +T_coeff4(pi1,k,nv0-2,m,n)*(8*(loc_f(:,:,k,nv0-1,m,n)&
                        -loc_f(:,:,k,nv0-3,m,n))+loc_f(:,:,k,nv0-4,m,n))


                   loc_k(:,:,k,nv0-1,m,n)=loc_k(:,:,k,nv0-1,m,n)+T_coeff1(pi1,k,nv0-1,m,n)&
                        *(16*loc_f(:,:,k,nv0-2,m,n)&
                        -loc_f(:,:,k,nv0-3,m,n))&
                        +T_coeff3(pi1,k,nv0-1,m,n)*(loc_f(:,:,k,nv0-3,m-2,n)&
                        -loc_f(:,:,k,nv0-3,m+2,n)&
                        +8*(loc_f(:,:,k,nv0-2,m+2,n)-loc_f(:,:,k,nv0-2,m-2,n)+loc_f(:,:,k,nv0-3,m+1,n)&
                        -loc_f(:,:,k,nv0-3,m-1,n)&
                        +8*(loc_f(:,:,k,nv0-2,m-1,n)-loc_f(:,:,k,nv0-2,m+1,n))))&
                        +T_coeff4(pi1,k,nv0-1,m,n)*(-8*loc_f(:,:,k,nv0-2,m,n)+loc_f(:,:,k,nv0-3,m,n))


                enddo
             enddo
          enddo

     !  elseif (coll_order .eq. 'order6') then
#else


          do n=ln1,ln2
             do m=lm1,lm2
                do l=ll1,ll2
                   do k=lk1,lk2

                      loc_k(:,:,k,l,m,n)=loc_k(:,:,k,l,m,n)+T_coeff2(pi1,k,l,m,n)&
                           *(loc_f(:,:,k,l,m+3,n)-13.5*(loc_f(:,:,k,l,m+2,n)+loc_f(:,:,k,l,m-2,n))&
                           +135*(loc_f(:,:,k,l,m+1,n)+loc_f(:,:,k,l,m-1,n))&
                           +loc_f(:,:,k,l,m-3,n))&
                           +T_coeff5(pi1,k,l,m,n)&
                           *(loc_f(:,:,k,l,m+3,n)-9*(loc_f(:,:,k,l,m+2,n)-loc_f(:,:,k,l,m-2,n))&
                           +45*(loc_f(:,:,k,l,m+1,n)-loc_f(:,:,k,l,m-1,n))&
                           -loc_f(:,:,k,l,m-3,n))&
                           +T_coeff6(pi1,k,l,m,n)*loc_f(:,:,k,l,m,n)


                   enddo
                enddo
             enddo
          enddo



          do n=ln1,ln2
             do m=lm1,lm2
                do k=lk1,lk2

                   loc_k(:,:,k,0,m,n)=loc_k(:,:,k,0,m,n)+T_coeff1(pi1,k,0,m,n)&
                        *(loc_f(:,:,k,3,m,n)-13.5*loc_f(:,:,k,2,m,n)&
                        +135*loc_f(:,:,k,1,m,n))&
                        +T_coeff3(pi1,k,0,m,n)&
                        *(loc_f(:,:,k,3,m+3,n)-loc_f(:,:,k,3,m-3,n)&
                        +9*(loc_f(:,:,k,2,m-3,n)+loc_f(:,:,k,3,m-2,n)-loc_f(:,:,k,3,m+2,n)&
                        -loc_f(:,:,k,2,m+3,n)&
                        +9*(loc_f(:,:,k,2,m+2,n)-loc_f(:,:,k,2,m-2,n)))&
                        +45*(loc_f(:,:,k,3,m+1,n)-loc_f(:,:,k,3,m-1,n)-loc_f(:,:,k,1,m-3,n)&
                        +loc_f(:,:,k,1,m+3,n)&
                        +45*(loc_f(:,:,k,1,m+1,n)-loc_f(:,:,k,1,m-1,n)))&
                        +405*(loc_f(:,:,k,2,m-1,n)+loc_f(:,:,k,1,m-2,n)&
                        -loc_f(:,:,k,1,m+2,n)-loc_f(:,:,k,2,m+1,n)))&
                        +T_coeff4(pi1,k,0,m,n)&
                        *(loc_f(:,:,k,3,m,n)-9*loc_f(:,:,k,2,m,n)&
                        +45*loc_f(:,:,k,1,m,n))


                   loc_k(:,:,k,1,m,n)=loc_k(:,:,k,1,m,n)+T_coeff1(pi1,k,1,m,n)&
                        *(loc_f(:,:,k,4,m,n)-13.5*loc_f(:,:,k,3,m,n)&
                        +135*(loc_f(:,:,k,2,m,n)+loc_f(:,:,k,0,m,n)))&
                        +T_coeff3(pi1,k,1,m,n)&
                        *(loc_f(:,:,k,4,m+3,n)-loc_f(:,:,k,4,m-3,n)&
                        +9*(loc_f(:,:,k,3,m-3,n)+loc_f(:,:,k,4,m-2,n)-loc_f(:,:,k,4,m+2,n)&
                        -loc_f(:,:,k,3,m+3,n)&
                        +9*(loc_f(:,:,k,3,m+2,n)-loc_f(:,:,k,3,m-2,n)))&
                        +45*(loc_f(:,:,k,4,m+1,n)-loc_f(:,:,k,4,m-1,n)-loc_f(:,:,k,2,m-3,n)&
                        +loc_f(:,:,k,0,m-3,n)&
                        -loc_f(:,:,k,0,m+3,n)+loc_f(:,:,k,2,m+3,n)&
                        +45*(loc_f(:,:,k,2,m+1,n)+loc_f(:,:,k,0,m-1,n)-loc_f(:,:,k,2,m-1,n)&
                        -loc_f(:,:,k,0,m+1,n)))&
                        +405*(loc_f(:,:,k,3,m-1,n)+loc_f(:,:,k,2,m-2,n)-loc_f(:,:,k,0,m-2,n)&
                        +loc_f(:,:,k,0,m+2,n)&
                        -loc_f(:,:,k,2,m+2,n)-loc_f(:,:,k,3,m+1,n)))&
                        +T_coeff4(pi1,k,1,m,n)&
                        *(loc_f(:,:,k,4,m,n)-9*loc_f(:,:,k,3,m,n)&
                        +45*(loc_f(:,:,k,2,m,n)-loc_f(:,:,k,0,m,n)))


                   loc_k(:,:,k,2,m,n)=loc_k(:,:,k,2,m,n)+T_coeff1(pi1,k,2,m,n)&
                        *(loc_f(:,:,k,5,m,n)-13.5*(loc_f(:,:,k,4,m,n)+loc_f(:,:,k,0,m,n))&
                        +135*(loc_f(:,:,k,3,m,n)+loc_f(:,:,k,1,m,n)))&
                        +T_coeff3(pi1,k,2,m,n)&
                        *(loc_f(:,:,k,5,m+3,n)-loc_f(:,:,k,5,m-3,n)&
                        +9*(loc_f(:,:,k,4,m-3,n)+loc_f(:,:,k,5,m-2,n)-loc_f(:,:,k,5,m+2,n)&
                        -loc_f(:,:,k,4,m+3,n)+loc_f(:,:,k,0,m+3,n)&
                        -loc_f(:,:,k,0,m-3,n)&
                        +9*(loc_f(:,:,k,4,m+2,n)+loc_f(:,:,k,0,m-2,n)-loc_f(:,:,k,4,m-2,n)&
                        -loc_f(:,:,k,0,m+2,n)))&
                        +45*(loc_f(:,:,k,5,m+1,n)-loc_f(:,:,k,5,m-1,n)-loc_f(:,:,k,3,m-3,n)&
                        +loc_f(:,:,k,1,m-3,n)&
                        -loc_f(:,:,k,1,m+3,n)+loc_f(:,:,k,3,m+3,n)&
                        +45*(loc_f(:,:,k,3,m+1,n)+loc_f(:,:,k,1,m-1,n)-loc_f(:,:,k,3,m-1,n)&
                        -loc_f(:,:,k,1,m+1,n)))&
                        +405*(loc_f(:,:,k,4,m-1,n)+loc_f(:,:,k,3,m-2,n)-loc_f(:,:,k,1,m-2,n)&
                        -loc_f(:,:,k,0,m-1,n)+loc_f(:,:,k,0,m+1,n)+loc_f(:,:,k,1,m+2,n)&
                        -loc_f(:,:,k,3,m+2,n)-loc_f(:,:,k,4,m+1,n)))&
                        +T_coeff4(pi1,k,2,m,n)&
                        *(loc_f(:,:,k,5,m,n)-9*(loc_f(:,:,k,4,m,n)-loc_f(:,:,k,0,m,n))&
                        +45*(loc_f(:,:,k,3,m,n)-loc_f(:,:,k,1,m,n)))

                enddo
                do l=3,nv0-4
                   do k=lk1,lk2

                      loc_k(:,:,k,l,m,n)=loc_k(:,:,k,l,m,n)+T_coeff1(pi1,k,l,m,n)&
                           *(loc_f(:,:,k,l+3,m,n)-13.5*(loc_f(:,:,k,l+2,m,n)+loc_f(:,:,k,l-2,m,n))&
                           +135*(loc_f(:,:,k,l+1,m,n)+loc_f(:,:,k,l-1,m,n))&
                           +loc_f(:,:,k,l-3,m,n))&
                           +T_coeff3(pi1,k,l,m,n)&
                           *(loc_f(:,:,k,l+3,m+3,n)-loc_f(:,:,k,l+3,m-3,n)-loc_f(:,:,k,l-3,m+3,n)&
                           +loc_f(:,:,k,l-3,m-3,n)&
                           +9*(loc_f(:,:,k,l+2,m-3,n)+loc_f(:,:,k,l+3,m-2,n)-loc_f(:,:,k,l+3,m+2,n)&
                           -loc_f(:,:,k,l+2,m+3,n)+loc_f(:,:,k,l-2,m+3,n)+loc_f(:,:,k,l-3,m+2,n)&
                           -loc_f(:,:,k,l-3,m-2,n)-loc_f(:,:,k,l-2,m-3,n)&
                           +9*(loc_f(:,:,k,l+2,m+2,n)+loc_f(:,:,k,l-2,m-2,n)-loc_f(:,:,k,l+2,m-2,n)&
                           -loc_f(:,:,k,l-2,m+2,n)))&
                           +45*(loc_f(:,:,k,l+3,m+1,n)-loc_f(:,:,k,l+3,m-1,n)-loc_f(:,:,k,l+1,m-3,n)&
                           +loc_f(:,:,k,l-1,m-3,n)+loc_f(:,:,k,l-3,m-1,n)-loc_f(:,:,k,l-3,m+1,n)&
                           -loc_f(:,:,k,l-1,m+3,n)+loc_f(:,:,k,l+1,m+3,n)&
                           +45*(loc_f(:,:,k,l+1,m+1,n)+loc_f(:,:,k,l-1,m-1,n)-loc_f(:,:,k,l+1,m-1,n)&
                           -loc_f(:,:,k,l-1,m+1,n)))&
                           +405*(loc_f(:,:,k,l+2,m-1,n)+loc_f(:,:,k,l+1,m-2,n)-loc_f(:,:,k,l-1,m-2,n)&
                           -loc_f(:,:,k,l-2,m-1,n)+loc_f(:,:,k,l-2,m+1,n)+loc_f(:,:,k,l-1,m+2,n)&
                           -loc_f(:,:,k,l+1,m+2,n)-loc_f(:,:,k,l+2,m+1,n)))&
                           +T_coeff4(pi1,k,l,m,n)&
                           *(loc_f(:,:,k,l+3,m,n)-9*(loc_f(:,:,k,l+2,m,n)-loc_f(:,:,k,l-2,m,n))&
                           +45*(loc_f(:,:,k,l+1,m,n)-loc_f(:,:,k,l-1,m,n))&
                           -loc_f(:,:,k,l-3,m,n))

                   enddo
                enddo
                do k=lk1,lk2

                   loc_k(:,:,k,nv0-3,m,n)=loc_k(:,:,k,nv0-3,m,n)+T_coeff1(pi1,k,nv0-3,m,n)&
                        *(-13.5*(loc_f(:,:,k,nv0-1,m,n)+loc_f(:,:,k,nv0-5,m,n))&
                        +135*(loc_f(:,:,k,nv0-2,m,n)+loc_f(:,:,k,nv0-4,m,n))&
                        +loc_f(:,:,k,nv0-6,m,n))&
                        +T_coeff3(pi1,k,nv0-3,m,n)&
                        *(-loc_f(:,:,k,nv0-6,m+3,n)&
                        +loc_f(:,:,k,nv0-6,m-3,n)&
                        +9*(loc_f(:,:,k,nv0-1,m-3,n)&
                        -loc_f(:,:,k,nv0-1,m+3,n)+loc_f(:,:,k,nv0-5,m+3,n)+loc_f(:,:,k,nv0-6,m+2,n)&
                        -loc_f(:,:,k,nv0-6,m-2,n)-loc_f(:,:,k,nv0-5,m-3,n)&
                        +9*(loc_f(:,:,k,nv0-1,m+2,n)+loc_f(:,:,k,nv0-5,m-2,n)-loc_f(:,:,k,nv0-1,m-2,n)&
                        -loc_f(:,:,k,nv0-5,m+2,n)))&
                        +45*(-loc_f(:,:,k,nv0-2,m-3,n)&
                        +loc_f(:,:,k,nv0-4,m-3,n)+loc_f(:,:,k,nv0-6,m-1,n)-loc_f(:,:,k,nv0-6,m+1,n)&
                        -loc_f(:,:,k,nv0-4,m+3,n)+loc_f(:,:,k,nv0-2,m+3,n)&
                        +45*(loc_f(:,:,k,nv0-2,m+1,n)+loc_f(:,:,k,nv0-4,m-1,n)-loc_f(:,:,k,nv0-2,m-1,n)&
                        -loc_f(:,:,k,nv0-4,m+1,n)))&
                        +405*(loc_f(:,:,k,nv0-1,m-1,n)+loc_f(:,:,k,nv0-2,m-2,n)-loc_f(:,:,k,nv0-4,m-2,n)&
                        -loc_f(:,:,k,nv0-5,m-1,n)+loc_f(:,:,k,nv0-5,m+1,n)+loc_f(:,:,k,nv0-4,m+2,n)&
                        -loc_f(:,:,k,nv0-2,m+2,n)-loc_f(:,:,k,nv0-1,m+1,n)))&
                        +T_coeff4(pi1,k,nv0-3,m,n)&
                        *(-9*(loc_f(:,:,k,nv0-1,m,n)-loc_f(:,:,k,nv0-5,m,n))&
                        +45*(loc_f(:,:,k,nv0-2,m,n)-loc_f(:,:,k,nv0-4,m,n))&
                        -loc_f(:,:,k,nv0-6,m,n))


                   loc_k(:,:,k,nv0-2,m,n)=loc_k(:,:,k,nv0-2,m,n)+T_coeff1(pi1,k,nv0-2,m,n)&
                        *(-13.5*loc_f(:,:,k,nv0-4,m,n)&
                        +135*(loc_f(:,:,k,nv0-1,m,n)+loc_f(:,:,k,nv0-3,m,n))&
                        +loc_f(:,:,k,nv0-5,m,n))&
                        +T_coeff3(pi1,k,nv0-2,m,n)&
                        *(-loc_f(:,:,k,nv0-5,m+3,n)&
                        +loc_f(:,:,k,nv0-5,m-3,n)&
                        +9*(loc_f(:,:,k,nv0-4,m+3,n)+loc_f(:,:,k,nv0-5,m+2,n)&
                        -loc_f(:,:,k,nv0-5,m-2,n)-loc_f(:,:,k,nv0-4,m-3,n)&
                        +9*(loc_f(:,:,k,nv0-4,m-2,n)&
                        -loc_f(:,:,k,nv0-4,m+2,n)))&
                        +45*(-loc_f(:,:,k,nv0-1,m-3,n)&
                        +loc_f(:,:,k,nv0-3,m-3,n)+loc_f(:,:,k,nv0-5,m-1,n)-loc_f(:,:,k,nv0-5,m+1,n)&
                        -loc_f(:,:,k,nv0-3,m+3,n)+loc_f(:,:,k,nv0-1,m+3,n)&
                        +45*(loc_f(:,:,k,nv0-1,m+1,n)+loc_f(:,:,k,nv0-3,m-1,n)-loc_f(:,:,k,nv0-1,m-1,n)&
                        -loc_f(:,:,k,nv0-3,m+1,n)))&
                        +405*(loc_f(:,:,k,nv0-1,m-2,n)-loc_f(:,:,k,nv0-3,m-2,n)&
                        -loc_f(:,:,k,nv0-4,m-1,n)+loc_f(:,:,k,nv0-4,m+1,n)+loc_f(:,:,k,nv0-3,m+2,n)&
                        -loc_f(:,:,k,nv0-1,m+2,n)))&
                        +T_coeff4(pi1,k,nv0-2,m,n)&
                        *(9*loc_f(:,:,k,nv0-4,m,n)&
                        +45*(loc_f(:,:,k,nv0-1,m,n)-loc_f(:,:,k,nv0-3,m,n))&
                        -loc_f(:,:,k,nv0-5,m,n))


                   loc_k(:,:,k,nv0-1,m,n)=loc_k(:,:,k,nv0-1,m,n)+T_coeff1(pi1,k,nv0-1,m,n)&
                        *(-13.5*loc_f(:,:,k,nv0-3,m,n)&
                        +135*loc_f(:,:,k,nv0-2,m,n)&
                        +loc_f(:,:,k,nv0-4,m,n))&
                        +T_coeff3(pi1,k,nv0-1,m,n)&
                        *(-loc_f(:,:,k,nv0-4,m+3,n)&
                        +loc_f(:,:,k,nv0-4,m-3,n)&
                        +9*(loc_f(:,:,k,nv0-3,m+3,n)+loc_f(:,:,k,nv0-4,m+2,n)&
                        -loc_f(:,:,k,nv0-4,m-2,n)-loc_f(:,:,k,nv0-3,m-3,n)&
                        +9*(loc_f(:,:,k,nv0-3,m-2,n)&
                        -loc_f(:,:,k,nv0-3,m+2,n)))&
                        +45*(loc_f(:,:,k,nv0-2,m-3,n)+loc_f(:,:,k,nv0-4,m-1,n)-loc_f(:,:,k,nv0-4,m+1,n)&
                        -loc_f(:,:,k,nv0-2,m+3,n)&
                        +45*(loc_f(:,:,k,nv0-2,m-1,n)&
                        -loc_f(:,:,k,nv0-2,m+1,n)))&
                        +405*(-loc_f(:,:,k,nv0-2,m-2,n)&
                        -loc_f(:,:,k,nv0-3,m-1,n)+loc_f(:,:,k,nv0-3,m+1,n)+loc_f(:,:,k,nv0-2,m+2,n)))&
                        +T_coeff4(pi1,k,nv0-1,m,n)&
                        *(9*loc_f(:,:,k,nv0-3,m,n)&
                        +45*(-loc_f(:,:,k,nv0-2,m,n))&
                        -loc_f(:,:,k,nv0-4,m,n))


                enddo
             enddo
          enddo

    !   endif
#endif

       particle_fac = (0.0,0.0)

       do n=ln1,ln2
          do sp=0,n_spec-1
             do m=lm1,lm2
                do l=ll1,ll2
                   do k=lk1,lk2

                      particle_fac(:,:,k,sp,n)=particle_fac(:,:,k,sp,n)+particle_diffpart(pi1,k,l,m,sp,n)*loc_f(:,:,k,l,m,n)

                   enddo
                enddo
             enddo
          enddo
       enddo

       call my_complex_sum_vw(particle_fac,size(particle_fac))

       do n=ln1,ln2
          do m=lm1,lm2
             do l=ll1,ll2
                do k=lk1,lk2
                   particle_buff=particle_fac(:,:,k,0,n)
                   do sp=1,n_spec-1
                      particle_buff(:,:)=particle_buff(:,:)+particle_fac(:,:,k,sp,n)
                   enddo
                   loc_k(:,:,k,l,m,n)=loc_k(:,:,k,l,m,n)&
                        -particle_coeff(pi1,k,l,m,n)*particle_buff(:,:)
                enddo
             enddo
          enddo
       enddo

    else

     !  if (coll_order .eq. 'fourth') then
#ifndef COLLSIX


          do n=ln1,ln2
             do m=lm1,lm2
                do l=ll1,ll2
                   do k=lk1,lk2
                      do j=lj1,lj2


                         loc_k(:,j,k,l,m,n)=loc_k(:,j,k,l,m,n)+T_coeff2(:,k,l,m,n)&
                              *(-loc_f(:,j,k,l,m+2,n)+16*(loc_f(:,j,k,l,m+1,n)&
                              +loc_f(:,j,k,l,m-1,n))-loc_f(:,j,k,l,m-2,n))&
                              +T_coeff5(:,k,l,m,n)*(-loc_f(:,j,k,l,m+2,n)+8*(loc_f(:,j,k,l,m+1,n)&
                              -loc_f(:,j,k,l,m-1,n))+loc_f(:,j,k,l,m-2,n))&
                              +T_coeff6(:,k,l,m,n)*loc_f(:,j,k,l,m,n)


                      enddo
                   enddo
                enddo
             enddo
          enddo


          do n=ln1,ln2
             do m=lm1,lm2
                do k=lk1,lk2
                   do j=lj1,lj2

                      loc_k(:,j,k,0,m,n)=loc_k(:,j,k,0,m,n)+T_coeff1(:,k,0,m,n)&
                           *(-loc_f(:,j,k,2,m,n)+16*loc_f(:,j,k,1,m,n))&
                           +T_coeff3(:,k,0,m,n)*(loc_f(:,j,k,2,m+2,n)&
                           -loc_f(:,j,k,2,m-2,n)&
                           +8*(loc_f(:,j,k,2,m-1,n)&
                           -loc_f(:,j,k,2,m+1,n)+loc_f(:,j,k,1,m-2,n)-loc_f(:,j,k,1,m+2,n)&
                           +8*(loc_f(:,j,k,1,m+1,n)-loc_f(:,j,k,1,m-1,n))))&
                           +T_coeff4(:,k,0,m,n)*(-loc_f(:,j,k,2,m,n)+8*loc_f(:,j,k,1,m,n))


                      loc_k(:,j,k,1,m,n)=loc_k(:,j,k,1,m,n)+T_coeff1(:,k,1,m,n)&
                           *(-loc_f(:,j,k,3,m,n)+16*(loc_f(:,j,k,2,m,n)+loc_f(:,j,k,0,m,n)))&
                           +T_coeff3(:,k,1,m,n)*(loc_f(:,j,k,3,m+2,n)&
                           -loc_f(:,j,k,3,m-2,n)&
                           +8*(loc_f(:,j,k,3,m-1,n)&
                           -loc_f(:,j,k,3,m+1,n)+loc_f(:,j,k,2,m-2,n)-loc_f(:,j,k,2,m+2,n)&
                           +loc_f(:,j,k,0,m+2,n)-loc_f(:,j,k,0,m-2,n)&
                           +8*(loc_f(:,j,k,2,m+1,n)+loc_f(:,j,k,0,m-1,n)-loc_f(:,j,k,2,m-1,n)&
                           -loc_f(:,j,k,0,m+1,n))))&
                           +T_coeff4(:,k,1,m,n)*(-loc_f(:,j,k,3,m,n)+8*(loc_f(:,j,k,2,m,n)&
                           -loc_f(:,j,k,0,m,n)))

                   enddo
                enddo
                do l=2,nv0-3
                   do k=lk1,lk2
                      do j=lj1,lj2

                         loc_k(:,j,k,l,m,n)=loc_k(:,j,k,l,m,n)+T_coeff1(:,k,l,m,n)&
                              *(-loc_f(:,j,k,l+2,m,n)+16*(loc_f(:,j,k,l+1,m,n)+loc_f(:,j,k,l-1,m,n))&
                              -loc_f(:,j,k,l-2,m,n))&
                              +T_coeff3(:,k,l,m,n)*(loc_f(:,j,k,l+2,m+2,n)+loc_f(:,j,k,l-2,m-2,n)&
                              -loc_f(:,j,k,l+2,m-2,n)-loc_f(:,j,k,l-2,m+2,n)&
                              +8*(loc_f(:,j,k,l+2,m-1,n)&
                              -loc_f(:,j,k,l+2,m+1,n)+loc_f(:,j,k,l+1,m-2,n)-loc_f(:,j,k,l+1,m+2,n)&
                              +loc_f(:,j,k,l-1,m+2,n)-loc_f(:,j,k,l-1,m-2,n)+loc_f(:,j,k,l-2,m+1,n)&
                              -loc_f(:,j,k,l-2,m-1,n)&
                              +8*(loc_f(:,j,k,l+1,m+1,n)+loc_f(:,j,k,l-1,m-1,n)-loc_f(:,j,k,l+1,m-1,n)&
                              -loc_f(:,j,k,l-1,m+1,n))))&
                              +T_coeff4(:,k,l,m,n)*(-loc_f(:,j,k,l+2,m,n)+8*(loc_f(:,j,k,l+1,m,n)&
                              -loc_f(:,j,k,l-1,m,n))+loc_f(:,j,k,l-2,m,n))

                      enddo
                   enddo
                enddo
                do k=lk1,lk2
                   do j=lj1,lj2

                      loc_k(:,j,k,nv0-2,m,n)=loc_k(:,j,k,nv0-2,m,n)+T_coeff1(:,k,nv0-2,m,n)&
                           *(16*(loc_f(:,j,k,nv0-1,m,n)+loc_f(:,j,k,nv0-3,m,n))&
                           -loc_f(:,j,k,nv0-4,m,n))&
                           +T_coeff3(:,k,nv0-2,m,n)*(loc_f(:,j,k,nv0-4,m-2,n)&
                           -loc_f(:,j,k,nv0-4,m+2,n)&
                           +8*(loc_f(:,j,k,nv0-1,m-2,n)-loc_f(:,j,k,nv0-1,m+2,n)&
                           +loc_f(:,j,k,nv0-3,m+2,n)-loc_f(:,j,k,nv0-3,m-2,n)+loc_f(:,j,k,nv0-4,m+1,n)&
                           -loc_f(:,j,k,nv0-4,m-1,n)&
                           +8*(loc_f(:,j,k,nv0-1,m+1,n)+loc_f(:,j,k,nv0-3,m-1,n)-loc_f(:,j,k,nv0-1,m-1,n)&
                           -loc_f(:,j,k,nv0-3,m+1,n))))&
                           +T_coeff4(:,k,nv0-2,m,n)*(8*(loc_f(:,j,k,nv0-1,m,n)&
                           -loc_f(:,j,k,nv0-3,m,n))+loc_f(:,j,k,nv0-4,m,n))


                      loc_k(:,j,k,nv0-1,m,n)=loc_k(:,j,k,nv0-1,m,n)+T_coeff1(:,k,nv0-1,m,n)&
                           *(16*loc_f(:,j,k,nv0-2,m,n)&
                           -loc_f(:,j,k,nv0-3,m,n))&
                           +T_coeff3(:,k,nv0-1,m,n)*(loc_f(:,j,k,nv0-3,m-2,n)&
                           -loc_f(:,j,k,nv0-3,m+2,n)&
                           +8*(loc_f(:,j,k,nv0-2,m+2,n)-loc_f(:,j,k,nv0-2,m-2,n)+loc_f(:,j,k,nv0-3,m+1,n)&
                           -loc_f(:,j,k,nv0-3,m-1,n)&
                           +8*(loc_f(:,j,k,nv0-2,m-1,n)-loc_f(:,j,k,nv0-2,m+1,n))))&
                           +T_coeff4(:,k,nv0-1,m,n)*(-8*loc_f(:,j,k,nv0-2,m,n)+loc_f(:,j,k,nv0-3,m,n))

                   enddo
                enddo
             enddo
          enddo

      ! elseif (coll_order .eq. 'order6') then
#else

          do n=ln1,ln2
             do m=lm1,lm2
                do l=ll1,ll2
                   do k=lk1,lk2
                      do j=lj1,lj2


                         loc_k(:,j,k,l,m,n)=loc_k(:,j,k,l,m,n)+T_coeff2(:,k,l,m,n)&
                              *(loc_f(:,j,k,l,m+3,n)-13.5*(loc_f(:,j,k,l,m+2,n)+loc_f(:,j,k,l,m-2,n))&
                              +135*(loc_f(:,j,k,l,m+1,n)+loc_f(:,j,k,l,m-1,n))&
                              +loc_f(:,j,k,l,m-3,n))&
                              +T_coeff5(:,k,l,m,n)&
                              *(loc_f(:,j,k,l,m+3,n)-9*(loc_f(:,j,k,l,m+2,n)-loc_f(:,j,k,l,m-2,n))&
                              +45*(loc_f(:,j,k,l,m+1,n)-loc_f(:,j,k,l,m-1,n))&
                              -loc_f(:,j,k,l,m-3,n))&
                              +T_coeff6(:,k,l,m,n)*loc_f(:,j,k,l,m,n)

                      enddo
                   enddo
                enddo
             enddo
          enddo


          do n=ln1,ln2
             do m=lm1,lm2
                do k=lk1,lk2
                   do j=lj1,lj2

                      loc_k(:,j,k,0,m,n)=loc_k(:,j,k,0,m,n)+T_coeff1(:,k,0,m,n)&
                           *(loc_f(:,j,k,3,m,n)-13.5*loc_f(:,j,k,2,m,n)&
                           +135*loc_f(:,j,k,1,m,n))&
                           +T_coeff3(:,k,0,m,n)&
                           *(loc_f(:,j,k,3,m+3,n)-loc_f(:,j,k,3,m-3,n)&
                           +9*(loc_f(:,j,k,2,m-3,n)+loc_f(:,j,k,3,m-2,n)-loc_f(:,j,k,3,m+2,n)&
                           -loc_f(:,j,k,2,m+3,n)&
                           +9*(loc_f(:,j,k,2,m+2,n)-loc_f(:,j,k,2,m-2,n)))&
                           +45*(loc_f(:,j,k,3,m+1,n)-loc_f(:,j,k,3,m-1,n)-loc_f(:,j,k,1,m-3,n)&
                           +loc_f(:,j,k,1,m+3,n)&
                           +45*(loc_f(:,j,k,1,m+1,n)-loc_f(:,j,k,1,m-1,n)))&
                           +405*(loc_f(:,j,k,2,m-1,n)+loc_f(:,j,k,1,m-2,n)&
                           -loc_f(:,j,k,1,m+2,n)-loc_f(:,j,k,2,m+1,n)))&
                           +T_coeff4(:,k,0,m,n)&
                           *(loc_f(:,j,k,3,m,n)-9*loc_f(:,j,k,2,m,n)&
                           +45*loc_f(:,j,k,1,m,n))


                      loc_k(:,j,k,1,m,n)=loc_k(:,j,k,1,m,n)+T_coeff1(:,k,1,m,n)&
                           *(loc_f(:,j,k,4,m,n)-13.5*loc_f(:,j,k,3,m,n)&
                           +135*(loc_f(:,j,k,2,m,n)+loc_f(:,j,k,0,m,n)))&
                           +T_coeff3(:,k,1,m,n)&
                           *(loc_f(:,j,k,4,m+3,n)-loc_f(:,j,k,4,m-3,n)&
                           +9*(loc_f(:,j,k,3,m-3,n)+loc_f(:,j,k,4,m-2,n)-loc_f(:,j,k,4,m+2,n)&
                           -loc_f(:,j,k,3,m+3,n)&
                           +9*(loc_f(:,j,k,3,m+2,n)-loc_f(:,j,k,3,m-2,n)))&
                           +45*(loc_f(:,j,k,4,m+1,n)-loc_f(:,j,k,4,m-1,n)-loc_f(:,j,k,2,m-3,n)&
                           +loc_f(:,j,k,0,m-3,n)&
                           -loc_f(:,j,k,0,m+3,n)+loc_f(:,j,k,2,m+3,n)&
                           +45*(loc_f(:,j,k,2,m+1,n)+loc_f(:,j,k,0,m-1,n)-loc_f(:,j,k,2,m-1,n)&
                           -loc_f(:,j,k,0,m+1,n)))&
                           +405*(loc_f(:,j,k,3,m-1,n)+loc_f(:,j,k,2,m-2,n)-loc_f(:,j,k,0,m-2,n)&
                           +loc_f(:,j,k,0,m+2,n)&
                           -loc_f(:,j,k,2,m+2,n)-loc_f(:,j,k,3,m+1,n)))&
                           +T_coeff4(:,k,1,m,n)&
                           *(loc_f(:,j,k,4,m,n)-9*loc_f(:,j,k,3,m,n)&
                           +45*(loc_f(:,j,k,2,m,n)-loc_f(:,j,k,0,m,n)))


                      loc_k(:,j,k,2,m,n)=loc_k(:,j,k,2,m,n)+T_coeff1(:,k,2,m,n)&
                           *(loc_f(:,j,k,5,m,n)-13.5*(loc_f(:,j,k,4,m,n)+loc_f(:,j,k,0,m,n))&
                           +135*(loc_f(:,j,k,3,m,n)+loc_f(:,j,k,1,m,n)))&
                           +T_coeff3(:,k,2,m,n)&
                           *(loc_f(:,j,k,5,m+3,n)-loc_f(:,j,k,5,m-3,n)&
                           +9*(loc_f(:,j,k,4,m-3,n)+loc_f(:,j,k,5,m-2,n)-loc_f(:,j,k,5,m+2,n)&
                           -loc_f(:,j,k,4,m+3,n)+loc_f(:,j,k,0,m+3,n)&
                           -loc_f(:,j,k,0,m-3,n)&
                           +9*(loc_f(:,j,k,4,m+2,n)+loc_f(:,j,k,0,m-2,n)-loc_f(:,j,k,4,m-2,n)&
                           -loc_f(:,j,k,0,m+2,n)))&
                           +45*(loc_f(:,j,k,5,m+1,n)-loc_f(:,j,k,5,m-1,n)-loc_f(:,j,k,3,m-3,n)&
                           +loc_f(:,j,k,1,m-3,n)&
                           -loc_f(:,j,k,1,m+3,n)+loc_f(:,j,k,3,m+3,n)&
                           +45*(loc_f(:,j,k,3,m+1,n)+loc_f(:,j,k,1,m-1,n)-loc_f(:,j,k,3,m-1,n)&
                           -loc_f(:,j,k,1,m+1,n)))&
                           +405*(loc_f(:,j,k,4,m-1,n)+loc_f(:,j,k,3,m-2,n)-loc_f(:,j,k,1,m-2,n)&
                           -loc_f(:,j,k,0,m-1,n)+loc_f(:,j,k,0,m+1,n)+loc_f(:,j,k,1,m+2,n)&
                           -loc_f(:,j,k,3,m+2,n)-loc_f(:,j,k,4,m+1,n)))&
                           +T_coeff4(:,k,2,m,n)&
                           *(loc_f(:,j,k,5,m,n)-9*(loc_f(:,j,k,4,m,n)-loc_f(:,j,k,0,m,n))&
                           +45*(loc_f(:,j,k,3,m,n)-loc_f(:,j,k,1,m,n)))

                   enddo
                enddo
                do l=3,nv0-4
                   do k=lk1,lk2
                      do j=lj1,lj2

                         loc_k(:,j,k,l,m,n)=loc_k(:,j,k,l,m,n)+T_coeff1(:,k,l,m,n)&
                              *(loc_f(:,j,k,l+3,m,n)-13.5*(loc_f(:,j,k,l+2,m,n)+loc_f(:,j,k,l-2,m,n))&
                              +135*(loc_f(:,j,k,l+1,m,n)+loc_f(:,j,k,l-1,m,n))&
                              +loc_f(:,j,k,l-3,m,n))&
                              +T_coeff3(:,k,l,m,n)&
                              *(loc_f(:,j,k,l+3,m+3,n)-loc_f(:,j,k,l+3,m-3,n)-loc_f(:,j,k,l-3,m+3,n)&
                              +loc_f(:,j,k,l-3,m-3,n)&
                              +9*(loc_f(:,j,k,l+2,m-3,n)+loc_f(:,j,k,l+3,m-2,n)-loc_f(:,j,k,l+3,m+2,n)&
                              -loc_f(:,j,k,l+2,m+3,n)+loc_f(:,j,k,l-2,m+3,n)+loc_f(:,j,k,l-3,m+2,n)&
                              -loc_f(:,j,k,l-3,m-2,n)-loc_f(:,j,k,l-2,m-3,n)&
                              +9*(loc_f(:,j,k,l+2,m+2,n)+loc_f(:,j,k,l-2,m-2,n)-loc_f(:,j,k,l+2,m-2,n)&
                              -loc_f(:,j,k,l-2,m+2,n)))&
                              +45*(loc_f(:,j,k,l+3,m+1,n)-loc_f(:,j,k,l+3,m-1,n)-loc_f(:,j,k,l+1,m-3,n)&
                              +loc_f(:,j,k,l-1,m-3,n)+loc_f(:,j,k,l-3,m-1,n)-loc_f(:,j,k,l-3,m+1,n)&
                              -loc_f(:,j,k,l-1,m+3,n)+loc_f(:,j,k,l+1,m+3,n)&
                              +45*(loc_f(:,j,k,l+1,m+1,n)+loc_f(:,j,k,l-1,m-1,n)-loc_f(:,j,k,l+1,m-1,n)&
                              -loc_f(:,j,k,l-1,m+1,n)))&
                              +405*(loc_f(:,j,k,l+2,m-1,n)+loc_f(:,j,k,l+1,m-2,n)-loc_f(:,j,k,l-1,m-2,n)&
                              -loc_f(:,j,k,l-2,m-1,n)+loc_f(:,j,k,l-2,m+1,n)+loc_f(:,j,k,l-1,m+2,n)&
                              -loc_f(:,j,k,l+1,m+2,n)-loc_f(:,j,k,l+2,m+1,n)))&
                              +T_coeff4(:,k,l,m,n)&
                              *(loc_f(:,j,k,l+3,m,n)-9*(loc_f(:,j,k,l+2,m,n)-loc_f(:,j,k,l-2,m,n))&
                              +45*(loc_f(:,j,k,l+1,m,n)-loc_f(:,j,k,l-1,m,n))&
                              -loc_f(:,j,k,l-3,m,n))

                      enddo
                   enddo
                enddo
                do k=lk1,lk2
                   do j=lj1,lj2

                      loc_k(:,j,k,nv0-3,m,n)=loc_k(:,j,k,nv0-3,m,n)+T_coeff1(:,k,nv0-3,m,n)&
                           *(-13.5*(loc_f(:,j,k,nv0-1,m,n)+loc_f(:,j,k,nv0-5,m,n))&
                           +135*(loc_f(:,j,k,nv0-2,m,n)+loc_f(:,j,k,nv0-4,m,n))&
                           +loc_f(:,j,k,nv0-6,m,n))&
                           +T_coeff3(:,k,nv0-3,m,n)&
                           *(-loc_f(:,j,k,nv0-6,m+3,n)&
                           +loc_f(:,j,k,nv0-6,m-3,n)&
                           +9*(loc_f(:,j,k,nv0-1,m-3,n)&
                           -loc_f(:,j,k,nv0-1,m+3,n)+loc_f(:,j,k,nv0-5,m+3,n)+loc_f(:,j,k,nv0-6,m+2,n)&
                           -loc_f(:,j,k,nv0-6,m-2,n)-loc_f(:,j,k,nv0-5,m-3,n)&
                           +9*(loc_f(:,j,k,nv0-1,m+2,n)+loc_f(:,j,k,nv0-5,m-2,n)-loc_f(:,j,k,nv0-1,m-2,n)&
                           -loc_f(:,j,k,nv0-5,m+2,n)))&
                           +45*(-loc_f(:,j,k,nv0-2,m-3,n)&
                           +loc_f(:,j,k,nv0-4,m-3,n)+loc_f(:,j,k,nv0-6,m-1,n)-loc_f(:,j,k,nv0-6,m+1,n)&
                           -loc_f(:,j,k,nv0-4,m+3,n)+loc_f(:,j,k,nv0-2,m+3,n)&
                           +45*(loc_f(:,j,k,nv0-2,m+1,n)+loc_f(:,j,k,nv0-4,m-1,n)-loc_f(:,j,k,nv0-2,m-1,n)&
                           -loc_f(:,j,k,nv0-4,m+1,n)))&
                           +405*(loc_f(:,j,k,nv0-1,m-1,n)+loc_f(:,j,k,nv0-2,m-2,n)-loc_f(:,j,k,nv0-4,m-2,n)&
                           -loc_f(:,j,k,nv0-5,m-1,n)+loc_f(:,j,k,nv0-5,m+1,n)+loc_f(:,j,k,nv0-4,m+2,n)&
                           -loc_f(:,j,k,nv0-2,m+2,n)-loc_f(:,j,k,nv0-1,m+1,n)))&
                           +T_coeff4(:,k,nv0-3,m,n)&
                           *(-9*(loc_f(:,j,k,nv0-1,m,n)-loc_f(:,j,k,nv0-5,m,n))&
                           +45*(loc_f(:,j,k,nv0-2,m,n)-loc_f(:,j,k,nv0-4,m,n))&
                           -loc_f(:,j,k,nv0-6,m,n))


                      loc_k(:,j,k,nv0-2,m,n)=loc_k(:,j,k,nv0-2,m,n)+T_coeff1(:,k,nv0-2,m,n)&
                           *(-13.5*loc_f(:,j,k,nv0-4,m,n)&
                           +135*(loc_f(:,j,k,nv0-1,m,n)+loc_f(:,j,k,nv0-3,m,n))&
                           +loc_f(:,j,k,nv0-5,m,n))&
                           +T_coeff3(:,k,nv0-2,m,n)&
                           *(-loc_f(:,j,k,nv0-5,m+3,n)&
                           +loc_f(:,j,k,nv0-5,m-3,n)&
                           +9*(loc_f(:,j,k,nv0-4,m+3,n)+loc_f(:,j,k,nv0-5,m+2,n)&
                           -loc_f(:,j,k,nv0-5,m-2,n)-loc_f(:,j,k,nv0-4,m-3,n)&
                           +9*(loc_f(:,j,k,nv0-4,m-2,n)&
                           -loc_f(:,j,k,nv0-4,m+2,n)))&
                           +45*(-loc_f(:,j,k,nv0-1,m-3,n)&
                           +loc_f(:,j,k,nv0-3,m-3,n)+loc_f(:,j,k,nv0-5,m-1,n)-loc_f(:,j,k,nv0-5,m+1,n)&
                           -loc_f(:,j,k,nv0-3,m+3,n)+loc_f(:,j,k,nv0-1,m+3,n)&
                           +45*(loc_f(:,j,k,nv0-1,m+1,n)+loc_f(:,j,k,nv0-3,m-1,n)-loc_f(:,j,k,nv0-1,m-1,n)&
                           -loc_f(:,j,k,nv0-3,m+1,n)))&
                           +405*(loc_f(:,j,k,nv0-1,m-2,n)-loc_f(:,j,k,nv0-3,m-2,n)&
                           -loc_f(:,j,k,nv0-4,m-1,n)+loc_f(:,j,k,nv0-4,m+1,n)+loc_f(:,j,k,nv0-3,m+2,n)&
                           -loc_f(:,j,k,nv0-1,m+2,n)))&
                           +T_coeff4(:,k,nv0-2,m,n)&
                           *(9*loc_f(:,j,k,nv0-4,m,n)&
                           +45*(loc_f(:,j,k,nv0-1,m,n)-loc_f(:,j,k,nv0-3,m,n))&
                           -loc_f(:,j,k,nv0-5,m,n))


                      loc_k(:,j,k,nv0-1,m,n)=loc_k(:,j,k,nv0-1,m,n)+T_coeff1(:,k,nv0-1,m,n)&
                           *(-13.5*loc_f(:,j,k,nv0-3,m,n)&
                           +135*loc_f(:,j,k,nv0-2,m,n)&
                           +loc_f(:,j,k,nv0-4,m,n))&
                           +T_coeff3(:,k,nv0-1,m,n)&
                           *(-loc_f(:,j,k,nv0-4,m+3,n)&
                           +loc_f(:,j,k,nv0-4,m-3,n)&
                           +9*(loc_f(:,j,k,nv0-3,m+3,n)+loc_f(:,j,k,nv0-4,m+2,n)&
                           -loc_f(:,j,k,nv0-4,m-2,n)-loc_f(:,j,k,nv0-3,m-3,n)&
                           +9*(loc_f(:,j,k,nv0-3,m-2,n)&
                           -loc_f(:,j,k,nv0-3,m+2,n)))&
                           +45*(loc_f(:,j,k,nv0-2,m-3,n)+loc_f(:,j,k,nv0-4,m-1,n)-loc_f(:,j,k,nv0-4,m+1,n)&
                           -loc_f(:,j,k,nv0-2,m+3,n)&
                           +45*(loc_f(:,j,k,nv0-2,m-1,n)&
                           -loc_f(:,j,k,nv0-2,m+1,n)))&
                           +405*(-loc_f(:,j,k,nv0-2,m-2,n)&
                           -loc_f(:,j,k,nv0-3,m-1,n)+loc_f(:,j,k,nv0-3,m+1,n)+loc_f(:,j,k,nv0-2,m+2,n)))&
                           +T_coeff4(:,k,nv0-1,m,n)&
                           *(9*loc_f(:,j,k,nv0-3,m,n)&
                           +45*(-loc_f(:,j,k,nv0-2,m,n))&
                           -loc_f(:,j,k,nv0-4,m,n))

                   enddo
                enddo
             enddo
          enddo

     !  endif
#endif

       particle_fac = (0.0,0.0)

       do n=ln1,ln2
          do sp=0,n_spec-1
             do m=lm1,lm2
                do l=ll1,ll2
                   do k=lk1,lk2
                      do j=lj1,lj2

                         particle_fac(:,j,k,sp,n)=particle_fac(:,j,k,sp,n)+particle_diffpart(:,k,l,m,sp,n)*loc_f(:,j,k,l,m,n)

                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo

       call my_complex_sum_vw(particle_fac,size(particle_fac))

       do n=ln1,ln2
          do m=lm1,lm2
             do l=ll1,ll2
                do k=lk1,lk2
                   do j=lj1,lj2
                         particle_buff(:,j)=particle_fac(:,j,k,0,n)
                         do sp=1,n_spec-1
                            particle_buff(:,j)=particle_buff(:,j)+particle_fac(:,j,k,sp,n)
                         enddo
                         loc_k(:,j,k,l,m,n)=loc_k(:,j,k,l,m,n)&
                              -particle_coeff(:,k,l,m,n)*particle_buff(:,j)
                   enddo
                enddo
             enddo
          enddo
       enddo

    endif

  end subroutine add_testdifferentialpart_alt



!---------
!FIELDPART

  !>Initializes the matrices to compute energy and parallel momentum and the conservation terms of
  !!the field particle contribution
  subroutine initialize_fieldpart
    integer:: i,j,k,l,m,n,pn,sp,pni
    complex,dimension(pi1:pi2):: en_coeff, n_coeff
    complex,dimension(:,:,:), allocatable:: mom_coeff
    complex,dimension(:,:,:,:), allocatable:: en_sum
    complex,dimension(:,:,:,:,:), allocatable:: n_cmat

    real:: x_a, x_b, alpha, kappa

    complex,dimension(:,:,:,:), allocatable:: cons_mom_coeff
    complex,dimension(:,:,:,:,:), allocatable:: cons_en_sum
    complex,dimension(:,:,:,:,:,:), allocatable:: cons_n_cmat

    real:: gamma_star, vTa_vTb, theta_ab, theta_ba, rho, postest
    real:: kperpendicular2

    PERFON('ini_fp')

    select case(em_cons)
    case (1)

       allocate(loc_mat_10(pi1:pi2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2), &
            mat_01_20(pi1:pi2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))

       allocate(mom_cmat(pi1:pi2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2), &
            en_cmat(pi1:pi2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2), &
            m_fac_sum(li1:li2, lj1:lj2, lk1:lk2, 0:n_spec-1), &
            e_fac_sum(li1:li2, lj1:lj2, lk1:lk2, 0:n_spec-1), &
            mom_coeff(pi1:pi2,lk1:lk2,ln1:ln2), &
            en_sum(pi1:pi2,lk1:lk2,1:4,ln1:ln2), &
            n_cmat(pi1:pi2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2))
    case (2)

       allocate(loc_mat_10(pi1:pi2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2), &
            mat_01_20(pi1:pi2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))

       allocate(rec_buff(li1:li2, lj1:lj2, lk1:lk2, ln1:ln2, 0:n_spec-1),&
            cons_mom_cmat(pi1:pi2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2, 0:n_spec-1), &
            cons_en_cmat(pi1:pi2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2, 0:n_spec-1),&
            cons_mom_coeff(pi1:pi2,lk1:lk2,ln1:ln2,0:n_spec-1), &
            cons_en_sum(pi1:pi2,lk1:lk2,1:4,ln1:ln2,0:n_spec-1), &
            cons_n_cmat(pi1:pi2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2,0:n_spec-1))
    case (3)
       allocate(moments_coll1(pi1:pi2,ll1:ll2,lk1:lk2,lm1:lm2,0:n_spec-1,ln1:ln2),&
            moments_coll2(pi1:pi2,ll1:ll2,lk1:lk2,lm1:lm2,0:n_spec-1,ln1:ln2),&
            moments_coll3(pi1:pi2,ll1:ll2,lk1:lk2,lm1:lm2,0:n_spec-1,ln1:ln2),&
            moments_coll4(pi1:pi2,ll1:ll2,lk1:lk2,lm1:lm2,ln1:ln2),&
            moments_coll5(pi1:pi2,ll1:ll2,lk1:lk2,lm1:lm2,ln1:ln2),&
            moments_coll6(pi1:pi2,ll1:ll2,lk1:lk2,lm1:lm2,ln1:ln2),&
            x_coeff1(pi1:pi2,ll1:ll2,lk1:lk2,lm1:lm2,0:n_spec-1,ln1:ln2),&
            x_coeff2(pi1:pi2,ll1:ll2,lk1:lk2,lm1:lm2,0:n_spec-1,ln1:ln2),&
            x_coeff3(pi1:pi2,ll1:ll2,lk1:lk2,lm1:lm2,0:n_spec-1,ln1:ln2),&
            x_coeff4(pi1:pi2,ll1:ll2,lk1:lk2,lm1:lm2,ln1:ln2),&
            x_coeff5(pi1:pi2,ll1:ll2,lk1:lk2,lm1:lm2,ln1:ln2),&
            x_coeff6(pi1:pi2,ll1:ll2,lk1:lk2,lm1:lm2,ln1:ln2),&
            y_coeff1(pi1:pi2,ll1:ll2,lk1:lk2,lm1:lm2,0:n_spec-1,ln1:ln2),&
            y_coeff2(pi1:pi2,ll1:ll2,lk1:lk2,lm1:lm2,0:n_spec-1,ln1:ln2),&
            y_coeff3(pi1:pi2,ll1:ll2,lk1:lk2,lm1:lm2,0:n_spec-1,ln1:ln2),&
            y_coeff4(pi1:pi2,ll1:ll2,lk1:lk2,lm1:lm2,0:n_spec-1,ln1:ln2),&
            y_coeff5(pi1:pi2,ll1:ll2,lk1:lk2,lm1:lm2,0:n_spec-1,ln1:ln2),&
            y_coeff6(pi1:pi2,ll1:ll2,lk1:lk2,lm1:lm2,0:n_spec-1,ln1:ln2),&
            y_coeffa(pi1:pi2,ll1:ll2,lk1:lk2,lm1:lm2,0:n_spec-1,ln1:ln2),&
            y_coeffb(pi1:pi2,ll1:ll2,lk1:lk2,lm1:lm2,0:n_spec-1,ln1:ln2),&
            y_coeffc(pi1:pi2,ll1:ll2,lk1:lk2,lm1:lm2,0:n_spec-1,ln1:ln2),&
            denominator1(pi1:pi2,lk1:lk2,ln1:ln2,0:(n_spec-1)),&
            denominator2(pi1:pi2,lk1:lk2,ln1:ln2,0:(n_spec-1)),&
            denominatora(pi1:pi2,lk1:lk2,ln1:ln2,0:(n_spec-1)),&
            denominatorb(pi1:pi2,lk1:lk2,ln1:ln2,0:(n_spec-1)),&
            qfactor(pi1:pi2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2,0:(n_spec-1)),&
            sfactor(pi1:pi2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2,0:(n_spec-1)),&
            pfactor(pi1:pi2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2,0:(n_spec-1)),&
            qintegral(pi1:pi2,lk1:lk2,ln1:ln2,0:(n_spec-1)),&
            sintegral(pi1:pi2,lk1:lk2,ln1:ln2,0:(n_spec-1)),&
            maxwellianint(pi1:pi2,lk1:lk2,ln1:ln2,0:(n_spec-1)),&
            quartic_int(pi1:pi2,lk1:lk2,ln1:ln2,0:n_spec-1),&
            R_int(pi1:pi2,lk1:lk2,ln1:ln2,0:n_spec-1),&
            momentum_int1(pi1:pi2,lk1:lk2,ln1:ln2,0:n_spec-1),&
            momentum_int2(pi1:pi2,lk1:lk2,ln1:ln2,0:n_spec-1))

    endselect

    ! normalised matrices for compensation of energy and parallel momentum loss due to collisions
    select case (em_cons)
    case (1)

    ! matrices for calculation of energy and parallel momentum
    loc_mat_10=0.0; mat_01_20=0.0
    do n=ln1,ln2
       do m=lm1,lm2
          do l=ll1,ll2
             do k=lk1,lk2
                do i=pi1,pi2
                   loc_mat_10(i,k,l,m,n)= spec(n)%dens*sqrt(spec(n)%temp*spec(n)%mass) &
                        *mat_00(i,pj1,k,l,m)*vp(l)
                   mat_01_20(i,k,l,m,n)= spec(n)%dens*spec(n)%temp &
                        *mat_00(i,pj1,k,l,m)*(mu(m)*geom%Bfield(i,pj1,k)+vp(l)**2)
                enddo
             enddo
          enddo
       enddo
    enddo

       mom_cmat=(0.,0.)
       en_cmat=(0.,0.)
       n_cmat=(0.,0.)
       mom_coeff=0.
       en_sum=0.
       do n=ln1,ln2
          if (pn0==1) then
             pn=pn1
          else
             pn=n
          endif
          do m=lm1,lm2
             do l=ll1,ll2
                do k=lk1,lk2
                   do i=pi1,pi2
                      mom_cmat(i,k,l,m,n)=fm(i,pj1,k,l,m,pn)*vp(l)
                      en_cmat(i,k,l,m,n)=fm(i,pj1,k,l,m,pn)*(vp(l)**2+mu(m)*geom%Bfield(i,pj1,k))
                      n_cmat(i,k,l,m,n)=fm(i,pj1,k,l,m,pn)
                      mom_coeff(i,k,n)=mom_coeff(i,k,n)+ mom_cmat(i,k,l,m,n)*loc_mat_10(i,k,l,m,n)
                      en_sum(i,k,1,n)= en_sum(i,k,1,n) + en_cmat(i,k,l,m,n)*mat_01_20(i,k,l,m,n)
                      en_sum(i,k,2,n)= en_sum(i,k,2,n) + n_cmat(i,k,l,m,n)*mat_01_20(i,k,l,m,n)
                      en_sum(i,k,3,n)= en_sum(i,k,3,n) + en_cmat(i,k,l,m,n)*mat_00(i,pj1,k,l,m)
                      en_sum(i,k,4,n)= en_sum(i,k,4,n) + n_cmat(i,k,l,m,n)*mat_00(i,pj1,k,l,m)
                   enddo
                enddo
             enddo
          enddo
       enddo

       call my_complex_sum_vw(mom_coeff,size(mom_coeff))
       call my_complex_sum_vw(en_sum,size(en_sum))

       do n=ln1,ln2
          do k=lk1,lk2
             en_coeff=1/(en_sum(:,k,1,n)-en_sum(:,k,2,n)*en_sum(:,k,3,n)/en_sum(:,k,4,n))
             n_coeff=-en_sum(:,k,3,n)/en_sum(:,k,4,n)
             do m=lm1,lm2
                do l=ll1,ll2
                   mom_cmat(:,k,l,m,n)= mom_cmat(:,k,l,m,n)/mom_coeff(:,k,n)
                   en_cmat(:,k,l,m,n)= en_coeff*(en_cmat(:,k,l,m,n)+n_coeff*n_cmat(:,k,l,m,n))
                enddo
             enddo
          enddo
       enddo

       deallocate(mom_coeff, en_sum, n_cmat)

    case (2)

       ! matrices for calculation of energy and parallel momentum
       loc_mat_10=0.0; mat_01_20=0.0
       do n=ln1,ln2
          do m=lm1,lm2
             do l=ll1,ll2
                do k=lk1,lk2
                   loc_mat_10(:,k,l,m,n)= spec(n)%dens*sqrt(spec(n)%temp*spec(n)%mass) &
                        *mat_00(:,pj1,k,l,m)*vp(l)
                   mat_01_20(:,k,l,m,n)= spec(n)%dens*spec(n)%temp &
                        *mat_00(:,pj1,k,l,m)*(mu(m)*geom%Bfield(pi1:pi2,pj1,k)+vp(l)**2)
                enddo
             enddo
          enddo
       enddo

       cons_mom_cmat=(0.,0.)
       cons_en_cmat=(0.,0.)
       cons_n_cmat=(0.,0.)
       cons_mom_coeff=0.
       cons_en_sum=0.
       do sp=0,n_spec-1
          do n=ln1,ln2
             if (pn0==1) then
                pn=pn1
             else
                pn=n
             endif
             do m=lm1,lm2
                do l=ll1,ll2
                   do k=lk1,lk2
                      do i=pi1,pi2
                         vTa_vTb = sqrt(spec(n)%temp*spec(n)%temp_prof(i)/spec(n)%mass&
                              /(spec(sp)%temp*spec(sp)%temp_prof(i)/spec(sp)%mass))
                         x_b = sqrt(vp(l)**2+geom%Bfield(i,pj1,k)*mu(m))*sqrt(2.)*vTa_vTb
                         !global: profile x_b(x)
                         cons_mom_cmat(i,k,l,m,n,sp)=fm(i,pj1,k,l,m,pn)*vp(l)&
                              *sqrt(2.*pi)/2.*(1+spec(sp)%mass/spec(n)%mass)**1.5*6.0*Hfunc(x_b)
                         cons_en_cmat(i,k,l,m,n,sp)=fm(i,pj1,k,l,m,pn)&
                              *sqrt(2.*pi)/2.*(1+spec(sp)%mass/spec(n)%mass)**1.5&
                              *(((1+spec(n)%mass/spec(sp)%mass)*x_b**2-1)*Hfunc(x_b)-Kfunc(x_b))
                         cons_n_cmat(i,k,l,m,n,sp)=fm(i,pj1,k,l,m,pn)
                         cons_mom_coeff(i,k,n,sp)=cons_mom_coeff(i,k,n,sp) + cons_mom_cmat(i,k,l,m,n,sp)*loc_mat_10(i,k,l,m,n)
                         cons_en_sum(i,k,1,n,sp)= cons_en_sum(i,k,1,n,sp) + cons_en_cmat(i,k,l,m,n,sp)*mat_01_20(i,k,l,m,n)
                         cons_en_sum(i,k,2,n,sp)= cons_en_sum(i,k,2,n,sp) + cons_n_cmat(i,k,l,m,n,sp)*mat_01_20(i,k,l,m,n)
                         cons_en_sum(i,k,3,n,sp)= cons_en_sum(i,k,3,n,sp) + cons_en_cmat(i,k,l,m,n,sp)*mat_00(i,pj1,k,l,m)
                         cons_en_sum(i,k,4,n,sp)= cons_en_sum(i,k,4,n,sp) + cons_n_cmat(i,k,l,m,n,sp)*mat_00(i,pj1,k,l,m)
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo

       call my_complex_sum_vw(cons_mom_coeff,size(cons_mom_coeff))
       call my_complex_sum_vw(cons_en_sum,size(cons_en_sum))

       do sp=0,n_spec-1
          do n=ln1,ln2
             do k=lk1,lk2
                en_coeff=1/(cons_en_sum(:,k,1,n,sp)-cons_en_sum(:,k,2,n,sp)*cons_en_sum(:,k,3,n,sp)/cons_en_sum(:,k,4,n,sp))
                n_coeff=-cons_en_sum(:,k,3,n,sp)/cons_en_sum(:,k,4,n,sp)
                do m=lm1,lm2
                   do l=ll1,ll2
                      cons_mom_cmat(:,k,l,m,n,sp)= cons_mom_cmat(:,k,l,m,n,sp)/cons_mom_coeff(:,k,n,sp)
                      cons_en_cmat(:,k,l,m,n,sp)= en_coeff*(cons_en_cmat(:,k,l,m,n,sp)+n_coeff*cons_n_cmat(:,k,l,m,n,sp))
                   enddo
                enddo
             enddo
          enddo
       enddo

       deallocate(cons_mom_coeff, cons_en_sum, cons_n_cmat)


    case(3)

       maxwellianint=0.0
       qintegral=0.0
       sintegral=0.0
       denominator1=0.0
       denominator2=0.0
       denominatora=0.0
       denominatorb=0.0
       quartic_int=0.0
       R_int=0.0
       momentum_int1=0.0
       momentum_int2=0.0

       do sp=0,n_spec-1
          pni=pn1
          do n=ln1,ln2
             if (pn0.gt.1) pni=n
             do l=ll1,ll2
                do m=lm1,lm2
                   do k=lk1,lk2
                         do i=pi1,pi2

                            postest=geom%Bfield(i,pj1,k)
                            if (postest < 0) then
                               Write (*,*) 'Bfield amplitude must be positive!'
                               STOP
                            endif

                            kappa = spec(sp)%dens*spec(sp)%dens_prof(i)*coll*zeff_spec(i,sp)*((2*spec(n)%charge*spec(sp)%charge)&
                                 **2)/sqrt(spec(n)%mass*((spec(n)%temp*spec(n)%temp_prof(i))**3))

                            theta_ab=sqrt((spec(n)%mass+spec(sp)%mass)*spec(n)%temp*spec(n)%temp_prof(i)&
                                 /(spec(n)%mass*spec(sp)%temp*spec(sp)%temp_prof(i)&
                                 +spec(sp)%mass*spec(n)%temp*spec(n)%temp_prof(i)))

                            theta_ba=sqrt((spec(sp)%mass+spec(n)%mass)*spec(sp)%temp*spec(sp)%temp_prof(i)&
                                 /(spec(sp)%mass*spec(n)%temp*spec(n)%temp_prof(i)&
                                 +spec(n)%mass*spec(sp)%temp*spec(sp)%temp_prof(i)))

                            alpha = sqrt(spec(n)%temp*spec(n)%temp_prof(i)/spec(n)%mass&
                                 /(spec(sp)%temp*spec(sp)%temp_prof(i)/spec(sp)%mass))

                            x_a=sqrt((vp(l)**2+mu(m)*geom%Bfield(i,pj1,k))/spec(n)%temp_prof(i))
                            x_b=x_a*alpha


                            pfactor(i,k,l,m,n,sp)= theta_ab*(C_ab_TO_mom(i,k,l,m,n,sp)-(16.0/(3.0*sqrt(pi)))*coll&
                                 *spec(sp)%dens*spec(sp)%dens_prof(i)*(((spec(n)%charge*spec(sp)%charge)**2)&
                                 /sqrt(spec(n)%mass*((spec(n)%temp*spec(n)%temp_prof(i))**3)))*fm(i,pj1,k,l,m,pni)&
                                 *vp(l)*(spec(n)%mass/spec(n)%temp)*alpha*(theta_ab-1)/sqrt(1+alpha**2))

                            qfactor(i,k,l,m,n,sp)= fm(i,pj1,k,l,m,pni)*((3*sqrt(pi)/2.0)&
                                 *Rfunc(alpha,x_b)+2*alpha*(theta_ab-1)*(x_a**2-1.5)/((1+alpha**2)**1.5))

                            sfactor(i,k,l,m,n,sp)= theta_ab*(C_ab_TO_en(i,k,l,m,n,sp)-(16.0/(3.0*sqrt(pi)))*coll&
                                 *spec(sp)%dens*spec(sp)%dens_prof(i)*(((spec(n)%charge*spec(sp)%charge)**2)&
                                 /sqrt(spec(n)%mass*((spec(n)%temp*spec(n)%temp_prof(i))**3)))*fm(i,pj1,k,l,m,pni)&
                                 *(x_a**2-1.5)*(theta_ab-1)*(2*alpha/((1+alpha**2)**1.5)))

                            maxwellianint(i,k,n,sp)=maxwellianint(i,k,n,sp)&
                                 +fm(i,pj1,k,l,m,pni)*mat_00(i,pj1,k,l,m)

                            quartic_int(i,k,n,sp)= quartic_int(i,k,n,sp)&
                                 +fm(i,pj1,k,l,m,pni)*(x_a**2-1.5)*(x_a**2)*mat_00(i,pj1,k,l,m)

                            R_int(i,k,n,sp)=R_int(i,k,n,sp)&
                                 +fm(i,pj1,k,l,m,pni)*Rfunc(alpha,x_b)*(x_a**2)*mat_00(i,pj1,k,l,m)

                            momentum_int1(i,k,n,sp)=momentum_int1(i,k,n,sp)&
                                 +fm(i,pj1,k,l,m,pni)*(vp(l)**2)*mat_00(i,pj1,k,l,m)

                            momentum_int2(i,k,n,sp)=momentum_int2(i,k,n,sp)&
                                 +fm(i,pj1,k,l,m,pni)*(vp(l)**2)*Hfunc(sqrt(2.0)*x_b)*sqrt(2.0)&
                                 *mat_00(i,pj1,k,l,m)

                         enddo
                   enddo
                enddo
             enddo
          enddo
       enddo

       call my_real_sum_vw(maxwellianint,size(maxwellianint))
       call my_real_sum_vw(quartic_int,size(quartic_int))
       call my_real_sum_vw(R_int,size(R_int))
       call my_real_sum_vw(momentum_int1,size(momentum_int1))
       call my_real_sum_vw(momentum_int2,size(momentum_int2))

       do sp=0,n_spec-1
          do n=ln1,ln2
             do l=ll1,ll2
                do m=lm1,lm2
                   do k=lk1,lk2

                   qintegral(:,k,n,sp)=qintegral(:,k,n,sp)&
                        +qfactor(:,k,l,m,n,sp)*mat_00(:,pj1,k,l,m)

                   sintegral(:,k,n,sp)=sintegral(:,k,n,sp)&
                        +sfactor(:,k,l,m,n,sp)*mat_00(:,pj1,k,l,m)

                   enddo
                enddo
             enddo
          enddo
       enddo

       call my_real_sum_vw(qintegral,size(qintegral))
       call my_real_sum_vw(sintegral,size(sintegral))

       do sp=0,n_spec-1
          pni=pn1
          do n=ln1,ln2
             if (pn0.gt.1) pni=n
             do l=ll1,ll2
                do m=lm1,lm2
                   do k=lk1,lk2
                         do i=pi1,pi2

                            theta_ab=sqrt((spec(n)%mass+spec(sp)%mass)*spec(n)%temp*spec(n)%temp_prof(i)&
                                 /(spec(n)%mass*spec(sp)%temp*spec(sp)%temp_prof(i)&
                                 +spec(sp)%mass*spec(n)%temp*spec(n)%temp_prof(i)))

                            theta_ba=sqrt((spec(sp)%mass+spec(n)%mass)*spec(sp)%temp*spec(sp)%temp_prof(i)&
                                 /(spec(sp)%mass*spec(n)%temp*spec(n)%temp_prof(i)&
                                 +spec(n)%mass*spec(sp)%temp*spec(sp)%temp_prof(i)))

                            alpha = sqrt(spec(n)%temp*spec(n)%temp_prof(i)/spec(n)%mass&
                                 /(spec(sp)%temp*spec(sp)%temp_prof(i)/spec(sp)%mass))

                            x_a=sqrt((vp(l)**2+mu(m)*geom%Bfield(i,pj1,k))/spec(n)%temp_prof(i))
                            x_b=x_a*alpha

                            denominator1(i,k,n,sp)=denominator1(i,k,n,sp)+fm(i,pj1,k,l,m,pni)&
                                 *(vp(l)**2)*((3*sqrt(pi)/2.0)*Hfunc(sqrt(2.0)*x_b)*sqrt(2.0)&
                                 +(theta_ab-1)/((1+alpha**2)**1.5))*mat_00(i,pj1,k,l,m)

                            denominatora(i,k,n,sp)=denominatora(i,k,n,sp)+vp(l)*pfactor(i,k,l,m,n,sp)&
                                 *mat_00(i,pj1,k,l,m)

                            denominator2(i,k,n,sp)=denominator2(i,k,n,sp)+(x_a**2)&
                                 *(qfactor(i,k,l,m,n,sp)-fm(i,pj1,k,l,m,pni)*qintegral(i,k,n,sp)&
                                 /maxwellianint(i,k,n,sp))*mat_00(i,pj1,k,l,m)

                            denominatorb(i,k,n,sp)=denominatorb(i,k,n,sp)+(x_a**2)&
                                 *(sfactor(i,k,l,m,n,sp)-fm(i,pj1,k,l,m,pni)*sintegral(i,k,n,sp)&
                                 /maxwellianint(i,k,n,sp))*mat_00(i,pj1,k,l,m)


                         enddo
                   enddo
                enddo
             enddo
          enddo
       enddo

       call my_real_sum_vw(denominator1,size(denominator1))
       call my_real_sum_vw(denominator2,size(denominator2))
       call my_real_sum_vw(denominatora,size(denominatora))
       call my_real_sum_vw(denominatorb,size(denominatorb))


       x_coeff4 = 0.0
       x_coeff5 = 0.0
       x_coeff6 = 0.0

       do sp=0,n_spec-1
          pni=pn1
          do n=ln1,ln2
             if (pn0.gt.1) pni=n
             do l=ll1,ll2
                do m=lm1,lm2
                   do k=lk1,lk2
                         do i=pi1,pi2

                            rho = sqrt(2.0*spec(n)%mass*spec(n)%temp*mu(m)&
                                 /(geom%Bfield(i,pj1,k)))/spec(n)%charge

                            theta_ab=sqrt((spec(n)%mass+spec(sp)%mass)*spec(n)%temp*spec(n)%temp_prof(i)&
                                 /(spec(n)%mass*spec(sp)%temp*spec(sp)%temp_prof(i)&
                                 +spec(sp)%mass*spec(n)%temp*spec(n)%temp_prof(i)))

                            theta_ba=sqrt((spec(sp)%mass+spec(n)%mass)*spec(sp)%temp*spec(sp)%temp_prof(i)&
                                 /(spec(sp)%mass*spec(n)%temp*spec(n)%temp_prof(i)&
                                 +spec(n)%mass*spec(sp)%temp*spec(sp)%temp_prof(i)))

                            alpha = sqrt(spec(n)%temp*spec(n)%temp_prof(i)/spec(n)%mass&
                                 /(spec(sp)%temp*spec(sp)%temp_prof(i)/spec(sp)%mass))

                            gamma_star=(16*coll/(3.0*sqrt(pi)))*((spec(n)%charge*spec(sp)%charge)**2)&
                                 *fm(i,pj1,k,l,m,pni)*spec(sp)%dens*spec(sp)%dens_prof(i)&
                                 *zeff_spec(i,sp)/sqrt(spec(n)%mass*((spec(n)%temp&
                                 *spec(n)%temp_prof(i))**3))

                            x_a=sqrt((vp(l)**2+mu(m)*geom%Bfield(i,pj1,k))/spec(n)%temp_prof(i))
                            x_b=x_a*alpha

                            x_coeff1(i,l,k,m,sp,n)= (theta_ab-1)*fm(i,pj1,k,l,m,pni)*vp(l)&
                                 /momentum_int1(i,k,n,sp)

                            x_coeff2(i,l,k,m,sp,n)= (theta_ab-1)*fm(i,pj1,k,l,m,pni)&
                                 *sqrt(mu(m)*geom%Bfield(i,pj1,k))/momentum_int1(i,k,n,sp)

                            x_coeff3(i,l,k,m,sp,n)= (theta_ab-1)*fm(i,pj1,k,l,m,pni)&
                                 *(x_a**2-3/2.0)/quartic_int(i,k,n,sp)

                            x_coeff4(i,l,k,m,n)=x_coeff4(i,l,k,m,n)-gamma_star*theta_ab*(alpha*(theta_ab-1)/sqrt(1+alpha**2))&
                                 *vp(l)*((3.0*sqrt(pi)/2.0)*Hfunc(sqrt(2.0)*x_b)*sqrt(2.0)+(theta_ab-1)/((1+alpha**2)**1.5))&
                                 /denominator1(i,k,n,sp)

                            x_coeff5(i,l,k,m,n)=x_coeff5(i,l,k,m,n)-gamma_star*theta_ab*(alpha*(theta_ab-1)/sqrt(1+alpha**2))&
                                 *sqrt(mu(m)*geom%Bfield(i,pj1,k))&
                                 *((3.0*sqrt(pi)/2.0)*Hfunc(sqrt(2.0)*x_b)*sqrt(2.0)+(theta_ab-1)/((1+alpha**2)**1.5))&
                                 /denominator1(i,k,n,sp)

                            x_coeff6(i,l,k,m,n)=x_coeff6(i,l,k,m,n)-gamma_star*(theta_ab-1)*theta_ab&
                                 *(2*alpha/((1+alpha**2)**1.5))*(qfactor(i,k,l,m,n,sp)/fm(i,pj1,k,l,m,pni)&
                                 -qintegral(i,k,n,sp)/maxwellianint(i,k,n,sp))/denominator2(i,k,n,sp)

                            !Field-part coefficients

                            y_coeffa(i,l,k,m,sp,n)= -(spec(sp)%dens/spec(n)%dens)&
                                 *sqrt(spec(sp)%mass*spec(sp)%temp/(spec(n)%temp*spec(n)%mass))&
                                 *pfactor(i,k,l,m,n,sp)/denominatora(i,k,n,sp)

                            y_coeff1(i,l,k,m,sp,n)= -(theta_ba-1)*(spec(sp)%dens/spec(n)%dens)&
                                 *sqrt(spec(sp)%mass*spec(sp)%temp/(spec(n)%temp*spec(n)%mass))*((3*sqrt(pi)/2.0)&
                                 *Hfunc(sqrt(2.0)*x_b)*sqrt(2.0)+(theta_ab-1)/((1+alpha**2)**1.5))&
                                 *vp(l)*fm(i,pj1,k,l,m,pni)/denominator1(i,k,n,sp)

                            y_coeffb(i,l,k,m,sp,n)= -(spec(sp)%dens/spec(n)%dens)&
                                 *sqrt(spec(sp)%mass*spec(sp)%temp/(spec(n)%temp*spec(n)%mass))&
                                 *pfactor(i,k,l,m,n,sp)*(sqrt(mu(m)*geom%Bfield(i,pj1,k))/vp(l))/denominatora(i,k,n,sp)

                            y_coeff2(i,l,k,m,sp,n)= -(theta_ba-1)*(spec(sp)%dens/spec(n)%dens)&
                                 *sqrt(spec(sp)%mass*spec(sp)%temp/(spec(n)%temp*spec(n)%mass))*((3*sqrt(pi)/2.0)&
                                 *Hfunc(sqrt(2.0)*x_b)*sqrt(2.0)+(theta_ab-1)/((1+alpha**2)**1.5))&
                                 *fm(i,pj1,k,l,m,pni)*sqrt(mu(m)*geom%Bfield(i,pj1,k))&
                                 /denominator1(i,k,n,sp)

                            y_coeffc(i,l,k,m,sp,n)= -(sfactor(i,k,l,m,n,sp)-(sintegral(i,k,n,sp)&
                                 /maxwellianint(i,k,n,sp))*fm(i,pj1,k,l,m,pni))&
                                 *(spec(sp)%temp*spec(sp)%temp_prof(i)*spec(sp)%dens&
                                 /(spec(n)%temp*spec(n)%temp_prof(i)*spec(n)%dens))&
                                 /denominatorb(i,k,n,sp)

                            y_coeff3(i,l,k,m,sp,n)= -(theta_ba-1)&
                                 *(qfactor(i,k,l,m,n,sp)-(qintegral(i,k,n,sp)&
                                 /maxwellianint(i,k,n,sp))*fm(i,pj1,k,l,m,pni))&
                                 *(spec(sp)%temp*spec(sp)%temp_prof(i)*spec(sp)%dens&
                                 /(spec(n)%temp*spec(n)%temp_prof(i)*spec(n)%dens))&
                                 /denominator2(i,k,n,sp)

                            y_coeff4(i,l,k,m,sp,n)= gamma_star*((theta_ba-1)&
                                 *(spec(n)%dens_prof(i)/spec(sp)%dens_prof(i))&
                                 *sqrt((spec(n)%dens_prof(i)/spec(sp)%dens_prof(i))**3)&
                                 /((1+1/alpha**2)**1.5))*theta_ba*(1+1/alpha**2)*(1/alpha)&
                                 *((3*sqrt(pi)/2.0)*Hfunc(sqrt(2.0)*x_b)*sqrt(2.0)+(theta_ab-1)&
                                 /((1+alpha**2)**1.5))*(spec(n)%temp/spec(sp)%temp)*vp(l)&
                                 /denominator1(i,k,n,sp)

                            y_coeff5(i,l,k,m,sp,n)= gamma_star*((theta_ba-1)&
                                 *(spec(n)%dens_prof(i)/spec(sp)%dens_prof(i))&
                                 *sqrt((spec(n)%dens_prof(i)/spec(sp)%dens_prof(i))**3)&
                                 /((1+1/alpha**2)**1.5))*theta_ba*(1+1/alpha**2)*(1/alpha)&
                                 *((3*sqrt(pi)/2.0)*Hfunc(sqrt(2.0)*x_b)*sqrt(2.0)+(theta_ab-1)&
                                 /((1+alpha**2)**1.5))*(spec(n)%temp/spec(sp)%temp)*sqrt(mu(m)&
                                 *geom%Bfield(i,pj1,k))/denominator1(i,k,n,sp)

                            y_coeff6(i,l,k,m,sp,n)= 2*gamma_star*(theta_ba-1)*(1/alpha)&
                                 *(spec(n)%dens_prof(i)/spec(sp)%dens_prof(i))&
                                 *sqrt(spec(n)%temp_prof(i)/spec(sp)%temp_prof(i))*(1/(1+1/alpha**2)**1.5)&
                                 *theta_ba*(qfactor(i,k,l,m,n,sp)/fm(i,pj1,k,l,m,pni)-qintegral(i,k,n,sp)&
                                 /maxwellianint(i,k,n,sp))*sqrt(spec(n)%mass*spec(n)%temp/(spec(sp)%mass&
                                 *spec(sp)%temp))/denominator2(i,k,n,sp)

                            if (coll_f_fm_on) then

                               moments_coll1(i,l,k,m,sp,n)=-mat_00(i,pj1,k,l,m)*coll*(1+alpha**2)*alpha&
                                    *(8.0*spec(sp)%dens*spec(sp)%dens_prof(i)&
                                    *(((spec(n)%charge*spec(sp)%charge)**2)&
                                    /sqrt(spec(n)%mass*((spec(n)%temp*spec(n)%temp_prof(i))**3)))&
                                    *Hfunc(sqrt(2.0)*x_b)*sqrt(2.0)*vp(l))*fm(i,pj1,k,l,m,pni)
                               moments_coll2(i,l,k,m,sp,n)=-mat_00(i,pj1,k,l,m)*coll*(1+alpha**2)&
                                    *alpha*sqrt(mu(m)*geom%Bfield(i,pj1,k))&
                                    *(8.0*spec(sp)%dens*spec(sp)%dens_prof(i)&
                                    *(((spec(n)%charge*spec(sp)%charge)**2)&
                                    /sqrt(spec(n)%mass*((spec(n)%temp*spec(n)%temp_prof(i))**3)))&
                                    *Hfunc(sqrt(2.0)*x_b)*sqrt(2.0))*fm(i,pj1,k,l,m,pni)
                               moments_coll3(i,l,k,m,sp,n)=-mat_00(i,pj1,k,l,m)*coll&
                                    *Rfunc(alpha,x_b)*(8.0*spec(sp)%dens*spec(sp)%dens_prof(i)&
                                    *((spec(n)%charge*spec(sp)%charge)**2)&
                                    /sqrt(spec(n)%mass*((spec(n)%temp*spec(n)%temp_prof(i))**3)))&
                                    *fm(i,pj1,k,l,m,pni)
                               moments_coll4(i,l,k,m,n)=mat_00(i,pj1,k,l,m)*vp(l)*fm(i,pj1,k,l,m,pni)
                               moments_coll5(i,l,k,m,n)=mat_00(i,pj1,k,l,m)&
                                    *sqrt(mu(m)*geom%Bfield(i,pj1,k))*fm(i,pj1,k,l,m,pni)
                               moments_coll6(i,l,k,m,n)=mat_00(i,pj1,k,l,m)*(x_a**2-3/(2.0))*fm(i,pj1,k,l,m,pni)
                            else
                               moments_coll1(i,l,k,m,sp,n)=-mat_00(i,pj1,k,l,m)*(1+alpha**2)*alpha*coll&
                                    *(8.0*spec(sp)%dens*spec(sp)%dens_prof(i)*(((spec(n)%charge*spec(sp)%charge)**2)&
                                    /sqrt(spec(n)%mass*((spec(n)%temp*spec(n)%temp_prof(i))**3)))&
                                    *Hfunc(sqrt(2.0)*x_b)*sqrt(2.0)*vp(l))
                               moments_coll2(i,l,k,m,sp,n)=-mat_00(i,pj1,k,l,m)*coll&
                                    *(1+alpha**2)*alpha*sqrt(mu(m)*geom%Bfield(i,pj1,k))&
                                    *(8.0*spec(sp)%dens*spec(sp)%dens_prof(i)*(((spec(n)%charge*spec(sp)%charge)**2)&
                                    /sqrt(spec(n)%mass*((spec(n)%temp*spec(n)%temp_prof(i))**3)))&
                                    *Hfunc(sqrt(2.0)*x_b)*sqrt(2.0))
                               moments_coll3(i,l,k,m,sp,n)=-mat_00(i,pj1,k,l,m)*Rfunc(alpha,x_b)*coll&
                                    *((spec(n)%charge*spec(sp)%charge)**2)&
                                    *(8.0*spec(sp)%dens*spec(sp)%dens_prof(i)&
                                    /sqrt(spec(n)%mass*((spec(n)%temp*spec(n)%temp_prof(i))**3)))
                               moments_coll4(i,l,k,m,n)=mat_00(i,pj1,k,l,m)*vp(l)
                               moments_coll5(i,l,k,m,n)=mat_00(i,pj1,k,l,m)*sqrt(mu(m)*geom%Bfield(i,pj1,k))
                               moments_coll6(i,l,k,m,n)=mat_00(i,pj1,k,l,m)*(x_a**2-3/(2.0))
                            endif
                         enddo
                   enddo
                enddo
             enddo
          enddo
       enddo

       deallocate(qfactor,qintegral,maxwellianint,denominator1,denominator2,denominatora,denominatorb,sfactor,sintegral,pfactor)
       deallocate(R_int,quartic_int,momentum_int1,momentum_int2)

       if(allocated(C_ab_TO_mom)) deallocate(C_ab_TO_mom)
       if(allocated(C_ab_TO_en)) deallocate(C_ab_TO_en)

       !!! Now we correct for the particle conservation error in x_coeff3 and x_coeff6


       allocate(x_coeff3_int(pi1:pi2,lk1:lk2,ln1:ln2,0:n_spec-1),&
            x_coeff6_int(pi1:pi2,lk1:lk2,ln1:ln2))

       x_coeff3_int=0.0
       x_coeff6_int=0.0

       do n=ln1,ln2
          do m=lm1,lm2
             do k=lk1,lk2
                do l=ll1,ll2
                   do sp=0,n_spec-1
                      x_coeff3_int(:,k,n,sp)=x_coeff3_int(:,k,n,sp)&
                           +x_coeff3(:,l,k,m,sp,n)*mat_00(:,pj1,k,l,m)
                   enddo
                   x_coeff6_int(:,k,n)=x_coeff6_int(:,k,n)&
                        +x_coeff6(:,l,k,m,n)*mat_00(:,pj1,k,l,m)
                enddo
             enddo
          enddo
       enddo


       call my_real_sum_vw(x_coeff3_int,size(x_coeff3_int))
       call my_real_sum_vw(x_coeff6_int,size(x_coeff6_int))

       do n=ln1,ln2
          do m=lm1,lm2
             do k=lk1,lk2
                do l=ll1,ll2
                   do sp=0,n_spec-1
                      x_coeff3(:,l,k,m,sp,n)=x_coeff3(:,l,k,m,sp,n)&
                           -particle_coeff(:,k,l,m,n)*x_coeff3_int(:,k,n,sp)
                   enddo
                   !    x_coeff6(:,:,k,l,m,n)=x_coeff6(:,:,k,l,m,n)&
                   !         -particle_coeff(:,:,k,l,m,n)*x_coeff6_int(:,:,k,n)
                enddo
             enddo
          enddo
       enddo

       deallocate(x_coeff3_int,x_coeff6_int)

       if (xy_local) then
          allocate(j1fac(li1:li2,lj1:lj2,lk1:lk2,lm1:lm2,ln1:ln2))
          allocate(coll_buff(li1:li2,lj1:lj2,1:3))
          coll_buff = (0.0,0.0)

          do n=ln1,ln2
             do m=lm1,lm2
                do k=lk1,lk2
                   do j=lj1,lj2
                      do i=li1,li2

                         rho = sqrt(2.0*spec(n)%mass*spec(n)%temp*mu(m)&
                              /(geom%Bfield(pi1,pj1,k)))/spec(n)%charge

                         kperpendicular2=geom%gii(pi1,pj1,k)*ki(i)**2+2*geom%gij(pi1,pj1,k)&
                              *ki(i)*kj(j)+geom%gjj(pi1,pj1,k)*kj(j)**2

                         j1fac(i,j,k,m,n)=j1(sqrt(kperpendicular2)*rho)

                      enddo
                   enddo
                enddo
             enddo
          enddo
       endif

       if (.not.xy_local) then
          deallocate(moments_coll2)
          deallocate(moments_coll5)
          deallocate(x_coeff2,x_coeff5)
          deallocate(y_coeff2,y_coeff5,y_coeffb)
          allocate(gyro_in_a(li1:li2,lj1:lj2))
          allocate(gyro_in_b(li1:li2,lj1:lj2))
          allocate(gyro_out_a(li1:li2,lj1:lj2))
          allocate(gyro_out_b(li1:li2,lj1:lj2))
          gyro_in_a = (0.0,0.0)
          gyro_in_b = (0.0,0.0)
          gyro_out_a = (0.0,0.0)
          gyro_out_b = (0.0,0.0)
       endif

       if (xy_local) then
          allocate(collision_mom(li1:li2, lj1:lj2, 1:3, lk1:lk2, 0:n_spec-1, 0:n_spec-1, 1:2))
          allocate(collision_mom2(li1:li2, lj1:lj2, 1:3, lk1:lk2, 0:n_spec-1))
       else
          allocate(collision_mom(1:2, li1:li2, lj1:lj2, lk1:lk2, 0:n_spec-1, 0:n_spec-1, 1:2))
          allocate(collision_mom2(1:2, li1:li2, lj1:lj2, lk1:lk2, 0:n_spec-1))
       endif

    endselect

    init_status_fp = em_cons
    PERFOFF

  end subroutine initialize_fieldpart

  !>Adds the field particle part of the collision operator
  !!
  !!Exchanges the momentum and energy violations and adds conservation terms to the RHS
  subroutine add_fieldpart_1(m_fac,e_fac,l_k)
    implicit none
    !> g_1 type array to which the contribution is added
    complex,dimension(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2),intent(inout) :: l_k
    complex,dimension(li1:li2, lj1:lj2, lk1:lk2, ln1:ln2, 0:n_spec-1),intent(inout) :: m_fac !<parallel momentum that is transferred
    complex,dimension(li1:li2, lj1:lj2, lk1:lk2, ln1:ln2, 0:n_spec-1),intent(inout) :: e_fac !<energy that is transferred

    integer:: i,j,k,l,m,n !,ii

    PERFON('collf1')

    !overall energy/momentum lost on background species
    !PERFON('momfac')
    !sum over second-last (local species) index
    m_fac_sum = sum(m_fac,4)
    e_fac_sum = sum(e_fac,4)

    call my_complex_sum_vwspec(m_fac_sum,lijk0*n_spec)
    call my_complex_sum_vwspec(e_fac_sum,lijk0*n_spec)
    !PERFOFF

    !PERFON('compens')
    ! compensate for energy/momentum violation
#if 0
    do n=ln1,ln2
#ifdef WITHOMP_COLLISIONS
       !$omp parallel do default(none) &
       !$omp private(i,j,k) &
       !$omp shared(n,l_k,m_fac_sum,e_fac_sum,mom_cmat,en_cmat,xy_local)&
       !$omp shared(pi1,pj1,lk1,lj1,li1,ll1,lm1,lijk0,llm0,lj0,lk0,li0)
#endif
       do iv=0,lijk0-1
          k=lk1+Modulo(iv,lk0)
          j=lj1+Modulo(iv/lk0,lj0)
          i=li1+iv/(lk0*lj0)
          if (xy_local) then
             ii=pi1
          else
             ii=i
          endif
          call caxpy(llm0,m_fac_sum(i,j,k,n),mom_cmat(ii,k,ll1,lm1,n),lk0,l_k(i,j,k,ll1,lm1,n),lijk0)
          call caxpy(llm0,e_fac_sum(i,j,k,n),en_cmat(ii,k,ll1,lm1,n),lk0,l_k(i,j,k,ll1,lm1,n),lijk0)
       enddo
#ifdef WITHOMP_COLLISIONS
       !$omp end parallel do
#endif
    enddo
#else
    do n=ln1,ln2
       do m=lm1,lm2
          do l=ll1,ll2
             do k=lk1,lk2
                if (xy_local) then
                   l_k(:,:,k,l,m,n)=l_k(:,:,k,l,m,n)+mom_cmat(pi1,k,l,m,n)*m_fac_sum(:,:,k,n) &
                        +en_cmat(pi1,k,l,m,n)*e_fac_sum(:,:,k,n)
                else
                   do j=lj1,lj2
                      do i=li1,li2
                         l_k(i,j,k,l,m,n)=l_k(i,j,k,l,m,n)+mom_cmat(i,k,l,m,n)*m_fac_sum(i,j,k,n) &
                              +en_cmat(i,k,l,m,n)*e_fac_sum(i,j,k,n)
                      enddo
                   enddo
                endif
             enddo
          enddo
       enddo
    enddo
#endif
    !PERFOFF
    PERFOFF

  end subroutine add_fieldpart_1

  !>Adds the field particle part of the collision operator - self-adjoint implementation
  !!
  !!Exchanges the momentum and energy violations and adds conservation terms to the RHS
  subroutine add_fieldpart_2(m_fac,e_fac,l_k)
    !> g_1 type array to which the contribution is added
    complex,dimension(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2),intent(inout) :: l_k
    complex,dimension(li1:li2, lj1:lj2, lk1:lk2, ln1:ln2, 0:n_spec-1),intent(inout) :: m_fac !<parallel momentum that is transferred
    complex,dimension(li1:li2, lj1:lj2, lk1:lk2, ln1:ln2, 0:n_spec-1),intent(inout) :: e_fac !<energy that is transferred

    integer:: i,j,ii,k,l,m,n,sp,ierr
    integer:: sb1,sblock

    PERFON('collf2')

    !overall energy/momentum lost on background species
    !PERFON('momfac')
    call my_complex_sum_vw(m_fac,lijk0*n_spec*ln0)
    call my_complex_sum_vw(e_fac,lijk0*n_spec*ln0)

    !transpose species submatrix for en amd mom transfer
    !exchange blocks between processes (global transpose):
    call mpi_alltoall(m_fac,lijk0*ln0*ln0,MPI_COMPLEX_TYPE,&
         rec_buff,lijk0*ln0*ln0,MPI_COMPLEX_TYPE,mpi_comm_spec,ierr)
    do sblock=0,n_procs_s-1
       sb1=ln0*sblock
       do n=0,ln0-1
          do sp=0,ln0-1
             !locally transpose the blocks
             m_fac(:,:,:,ln1+n,sb1+sp)=rec_buff(:,:,:,ln1+sp,sb1+n)
          enddo
       enddo
    enddo

    call mpi_alltoall(e_fac,lijk0*ln0*ln0,MPI_COMPLEX_TYPE,&
         rec_buff,lijk0*ln0*ln0,MPI_COMPLEX_TYPE,mpi_comm_spec,ierr)
    do sblock=0,n_procs_s-1
       sb1=ln0*sblock
       do n=0,ln0-1
          do sp=0,ln0-1
             e_fac(:,:,:,ln1+n,sb1+sp)=rec_buff(:,:,:,ln1+sp,sb1+n)
          enddo
       enddo
    enddo

    !PERFOFF

    !PERFON('compens')
    do sp=0,n_spec-1
       do n=ln1,ln2
          do m=lm1,lm2
             do l=ll1,ll2
                do k=lk1,lk2
                   if (xy_local) then
                      !terms for different species are added
                      l_k(:,:,k,l,m,n)=l_k(:,:,k,l,m,n)&
                           +cons_mom_cmat(pi1,k,l,m,n,sp)*m_fac(:,:,k,n,sp) &
                           +cons_en_cmat(pi1,k,l,m,n,sp)*e_fac(:,:,k,n,sp)
                   else
                      do j=lj1,lj2
                         do i=li1,li2
                            if (xy_local) then
                               ii=pi1
                            else
                               ii=i
                            endif
                            l_k(i,j,k,l,m,n)=l_k(i,j,k,l,m,n)&
                                 +cons_mom_cmat(ii,k,l,m,n,sp)*m_fac(i,j,k,n,sp) &
                                 +cons_en_cmat(ii,k,l,m,n,sp)*e_fac(i,j,k,n,sp)
                         enddo
                      enddo
                   endif
                enddo
             enddo
          enddo
       enddo
    enddo

    !PERFOFF
    PERFOFF

  end subroutine add_fieldpart_2

  Subroutine add_momentpart_local(loc_f,loc_k)
    Implicit None
    complex,dimension(li1:li2, lj1:lj2, lbz:ubz, lbv:ubv, lbw:ubw, ln1:ln2),intent(in):: loc_f
    complex,dimension(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2),intent(inout):: loc_k
    Integer :: n,sp,m,l,k

    !PERFON('collf3')


    collision_mom = (0.0,0.0)
    collision_mom2 = (0.0,0.0)

    do n=ln1,ln2
       do sp=0,(n_spec-1)
          if (spec(n)%temp /= spec(sp)%temp) then
             do m=lm1,lm2
                do k=lk1,lk2
                   coll_buff = (0.0,0.0)
                   do l=ll1,ll2
                      coll_buff(:,:,1)=coll_buff(:,:,1)+moments_coll1(pi1,l,k,m,sp,n)&
                           *loc_f(:,:,k,l,m,n)
                      coll_buff(:,:,2)=coll_buff(:,:,2)+moments_coll2(pi1,l,k,m,sp,n)&
                           *loc_f(:,:,k,l,m,n)
                      coll_buff(:,:,3)=coll_buff(:,:,3)+moments_coll3(pi1,l,k,m,sp,n)&
                           *loc_f(:,:,k,l,m,n)
                   enddo
                   collision_mom(:,:,1,k,n,sp,1)=collision_mom(:,:,1,k,n,sp,1)+jfac(:,:,k,m,n)*coll_buff(:,:,1)
                   collision_mom(:,:,2,k,n,sp,1)=collision_mom(:,:,2,k,n,sp,1)+j1fac(:,:,k,m,n)*coll_buff(:,:,2)
                   collision_mom(:,:,3,k,n,sp,1)=collision_mom(:,:,3,k,n,sp,1)+jfac(:,:,k,m,n)*coll_buff(:,:,3)
                enddo
             enddo
          endif
       enddo
    enddo


    do n=ln1,ln2
       do m=lm1,lm2
          do k=lk1,lk2
             coll_buff = (0.0,0.0)
             do l=ll1,ll2
                coll_buff(:,:,1)=coll_buff(:,:,1)+moments_coll4(pi1,l,k,m,n)&
                     *loc_f(:,:,k,l,m,n)
                coll_buff(:,:,2)=coll_buff(:,:,2)+moments_coll5(pi1,l,k,m,n)&
                     *loc_f(:,:,k,l,m,n)
                coll_buff(:,:,3)=coll_buff(:,:,3)+moments_coll6(pi1,l,k,m,n)&
                     *loc_f(:,:,k,l,m,n)
             enddo
             collision_mom2(:,:,1,k,n)=collision_mom2(:,:,1,k,n)+jfac(:,:,k,m,n)*coll_buff(:,:,1)
             collision_mom2(:,:,2,k,n)=collision_mom2(:,:,2,k,n)+j1fac(:,:,k,m,n)*coll_buff(:,:,2)
             collision_mom2(:,:,3,k,n)=collision_mom2(:,:,3,k,n)+jfac(:,:,k,m,n)*coll_buff(:,:,3)
          enddo
       enddo
    enddo

    do n=ln1,ln2
       do sp=0,(n_spec-1)
          do m=lm1,lm2
             do k=lk1,lk2
                coll_buff = (0.0,0.0)
                do l=ll1,ll2
                   coll_buff(:,:,1)=coll_buff(:,:,1)+moments_colla(pi1,l,k,m,sp,n)&
                        *loc_f(:,:,k,l,m,n)
                   coll_buff(:,:,2)=coll_buff(:,:,2)+moments_collb(pi1,l,k,m,sp,n)&
                        *loc_f(:,:,k,l,m,n)
                   coll_buff(:,:,3)=coll_buff(:,:,3)+moments_collc(pi1,l,k,m,sp,n)&
                        *loc_f(:,:,k,l,m,n)
                enddo
                collision_mom(:,:,1,k,n,sp,2)=collision_mom(:,:,1,k,n,sp,2)+jfac(:,:,k,m,n)*coll_buff(:,:,1)
                collision_mom(:,:,2,k,n,sp,2)=collision_mom(:,:,2,k,n,sp,2)+j1fac(:,:,k,m,n)*coll_buff(:,:,2)
                collision_mom(:,:,3,k,n,sp,2)=collision_mom(:,:,3,k,n,sp,2)+jfac(:,:,k,m,n)*coll_buff(:,:,3)
             enddo
          enddo
       enddo
    enddo

    call my_complex_sum_vwspec(collision_mom,size(collision_mom))
    call my_complex_sum_vwspec(collision_mom2,size(collision_mom2))

    do n=ln1,ln2
       do sp=0,(n_spec-1)
          if (spec(n)%temp /= spec(sp)%temp) then
             do m=lm1,lm2
                do k=lk1,lk2
                   coll_buff(:,:,1)=jfac(:,:,k,m,n)*collision_mom(:,:,1,k,n,sp,1)
                   coll_buff(:,:,2)=j1fac(:,:,k,m,n)*collision_mom(:,:,2,k,n,sp,1)
                   coll_buff(:,:,3)=jfac(:,:,k,m,n)*collision_mom(:,:,3,k,n,sp,1)
                   do l=ll1,ll2
                      loc_k(:,:,k,l,m,n)=loc_k(:,:,k,l,m,n)&
                           +x_coeff1(pi1,l,k,m,sp,n)*coll_buff(:,:,1)&
                           +x_coeff2(pi1,l,k,m,sp,n)*coll_buff(:,:,2)&
                           +x_coeff3(pi1,l,k,m,sp,n)*coll_buff(:,:,3)
                   enddo
                   coll_buff(:,:,1)=jfac(:,:,k,m,n)*collision_mom(:,:,1,k,sp,n,1)
                   coll_buff(:,:,2)=j1fac(:,:,k,m,n)*collision_mom(:,:,2,k,sp,n,1)
                   coll_buff(:,:,3)=jfac(:,:,k,m,n)*collision_mom(:,:,3,k,sp,n,1)
                   do l=ll1,ll2
                      loc_k(:,:,k,l,m,n)=loc_k(:,:,k,l,m,n)&
                           +y_coeff1(pi1,l,k,m,sp,n)*coll_buff(:,:,1)&
                           +y_coeff2(pi1,l,k,m,sp,n)*coll_buff(:,:,2)&
                           +y_coeff3(pi1,l,k,m,sp,n)*coll_buff(:,:,3)
                   enddo
                   coll_buff(:,:,1)=jfac(:,:,k,m,n)*collision_mom2(:,:,1,k,sp)
                   coll_buff(:,:,2)=j1fac(:,:,k,m,n)*collision_mom2(:,:,2,k,sp)
                   coll_buff(:,:,3)=jfac(:,:,k,m,n)*collision_mom2(:,:,3,k,sp)
                   do l=ll1,ll2
                      loc_k(:,:,k,l,m,n)=loc_k(:,:,k,l,m,n)&
                           +y_coeff4(pi1,l,k,m,sp,n)*coll_buff(:,:,1)&
                           +y_coeff5(pi1,l,k,m,sp,n)*coll_buff(:,:,2)&
                           +y_coeff6(pi1,l,k,m,sp,n)*coll_buff(:,:,3)
                   enddo
                enddo
             enddo
          endif
       enddo
    enddo

    do n=ln1,ln2
       do m=lm1,lm2
          do k=lk1,lk2
             coll_buff(:,:,1)=jfac(:,:,k,m,n)*collision_mom2(:,:,1,k,n)
             coll_buff(:,:,2)=j1fac(:,:,k,m,n)*collision_mom2(:,:,2,k,n)
             coll_buff(:,:,3)=jfac(:,:,k,m,n)*collision_mom2(:,:,3,k,n)
             do l=ll1,ll2
                loc_k(:,:,k,l,m,n)=loc_k(:,:,k,l,m,n)&
                     +x_coeff4(pi1,l,k,m,n)*coll_buff(:,:,1)&
                     +x_coeff5(pi1,l,k,m,n)*coll_buff(:,:,2)&
                     +x_coeff6(pi1,l,k,m,n)*coll_buff(:,:,3)
             enddo
          enddo
       enddo
    enddo

    do n=ln1,ln2
       do sp=0,(n_spec-1)
          do m=lm1,lm2
             do k=lk1,lk2
                coll_buff(:,:,1)=jfac(:,:,k,m,n)*collision_mom(:,:,1,k,sp,n,2)
                coll_buff(:,:,2)=j1fac(:,:,k,m,n)*collision_mom(:,:,2,k,sp,n,2)
                coll_buff(:,:,3)=jfac(:,:,k,m,n)*collision_mom(:,:,3,k,sp,n,2)
                do l=ll1,ll2
                   loc_k(:,:,k,l,m,n)=loc_k(:,:,k,l,m,n)&
                        +y_coeffa(pi1,l,k,m,sp,n)*coll_buff(:,:,1)&
                        +y_coeffb(pi1,l,k,m,sp,n)*coll_buff(:,:,2)&
                        +y_coeffc(pi1,l,k,m,sp,n)*coll_buff(:,:,3)
                enddo
             enddo
          enddo
       enddo
    enddo


  ! PERFOFF

  End subroutine add_momentpart_local


  Subroutine add_momentpart_global(loc_f,loc_k)
    Implicit None
    complex,dimension(li1:li2, lj1:lj2, lbz:ubz, lbv:ubv, lbw:ubw, ln1:ln2),intent(in):: loc_f
    complex,dimension(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2),intent(inout):: loc_k
    Integer :: n,sp,m,l,k,j,i

    !PERFON('collf3')

    collision_mom = (0.0,0.0)
    collision_mom2 = (0.0,0.0)

    do n=ln1,ln2
       do sp=0,(n_spec-1)
          if (spec(n)%temp /= spec(sp)%temp) then
             do m=lm1,lm2
                do k=lk1,lk2
                   do j=lj1,lj2
                      do i=li1,li2
                         gyro_in_a(i,j) = (0.0,0.0)
                         gyro_in_b(i,j) = (0.0,0.0)
                         do l=ll1,ll2
                            gyro_in_a(i,j)=gyro_in_a(i,j)+moments_coll1(i,l,k,m,sp,n)&
                                 *loc_f(i,j,k,l,m,n)
                            gyro_in_b(i,j)=gyro_in_b(i,j)+moments_coll3(i,l,k,m,sp,n)&
                                 *loc_f(i,j,k,l,m,n)
                         enddo
                         !call gyro_average_df(gyro_in_a,gyro_out_a,k,m,n,.true.,'N')
                         !call gyro_average_df(gyro_in_b,gyro_out_b,k,m,n,.true.,'N')
                         !collision_mom(1,i,j,k,n,sp,1)=collision_mom(1,i,j,k,n,sp,1)+gyro_out_a(i,j)
                         !collision_mom(2,i,j,k,n,sp,1)=collision_mom(2,i,j,k,n,sp,1)+gyro_out_b(i,j)
                         collision_mom(1,i,j,k,n,sp,1)=collision_mom(1,i,j,k,n,sp,1)+gyro_in_a(i,j)
                         collision_mom(2,i,j,k,n,sp,1)=collision_mom(2,i,j,k,n,sp,1)+gyro_in_b(i,j)

                      enddo
                   enddo
                enddo
             enddo
          endif
       enddo
    enddo


    do n=ln1,ln2
       do m=lm1,lm2
          do k=lk1,lk2
             do j=lj1,lj2
                do i=li1,li2
                   gyro_in_a(i,j)=(0.0,0.0)
                   gyro_in_b(i,j)=(0.0,0.0)
                   do l=ll1,ll2
                      gyro_in_a(i,j)=gyro_in_a(i,j)+loc_f(i,j,k,l,m,n)*moments_coll4(i,l,k,m,n)
                      gyro_in_b(i,j)=gyro_in_b(i,j)+loc_f(i,j,k,l,m,n)*moments_coll6(i,l,k,m,n)
                   enddo
                   !call gyro_average_df(gyro_in_a,gyro_out_a,k,m,n,.true.,'N')
                   !call gyro_average_df(gyro_in_b,gyro_out_b,k,m,n,.true.,'N')
                   !collision_mom2(1,i,j,k,n)=collision_mom2(1,i,j,k,n)+gyro_out_a(i,j)
                   !collision_mom2(2,i,j,k,n)=collision_mom2(2,i,j,k,n)+gyro_out_b(i,j)
                   collision_mom2(1,i,j,k,n)=collision_mom2(1,i,j,k,n)+gyro_in_a(i,j)
                   collision_mom2(2,i,j,k,n)=collision_mom2(2,i,j,k,n)+gyro_in_b(i,j)
                enddo
             enddo
          enddo
       enddo
    enddo

    do n=ln1,ln2
       do sp=0,(n_spec-1)
          do m=lm1,lm2
             do k=lk1,lk2
                do j=lj1,lj2
                   do i=li1,li2
                      gyro_in_a(i,j) = (0.0,0.0)
                      gyro_in_b(i,j) = (0.0,0.0)
                      do l=ll1,ll2
                         gyro_in_a(i,j)=gyro_in_a(i,j)+moments_colla(i,l,k,m,sp,n)&
                              *loc_f(i,j,k,l,m,n)
                         gyro_in_b(i,j)=gyro_in_b(i,j)+moments_collc(i,l,k,m,sp,n)&
                              *loc_f(i,j,k,l,m,n)
                      enddo
                      !call gyro_average_df(gyro_in_a,gyro_out_a,k,m,n,.true.,'N')
                      !call gyro_average_df(gyro_in_b,gyro_out_b,k,m,n,.true.,'N')
                      !collision_mom(1,i,j,k,n,sp,2)=collision_mom(1,i,j,k,n,sp,2)+gyro_out_a(i,j)
                      !collision_mom(2,i,j,k,n,sp,2)=collision_mom(2,i,j,k,n,sp,2)+gyro_out_b(i,j)
                      collision_mom(1,i,j,k,n,sp,2)=collision_mom(1,i,j,k,n,sp,2)+gyro_in_a(i,j)
                      collision_mom(2,i,j,k,n,sp,2)=collision_mom(2,i,j,k,n,sp,2)+gyro_in_b(i,j)

                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo

    call my_complex_sum_vwspec(collision_mom,size(collision_mom))
    call my_complex_sum_vwspec(collision_mom2,size(collision_mom2))

    do n=ln1,ln2
       do sp=0,(n_spec-1)
          if (spec(n)%temp /= spec(sp)%temp) then
             do m=lm1,lm2
                do k=lk1,lk2
                   gyro_in_a=(0.0,0.0)
                   gyro_in_b=(0.0,0.0)
                   do j=lj1,lj2
                      do i=li1,li2
                         gyro_in_a(i,j)=collision_mom(1,i,j,k,n,sp,1)
                         gyro_in_b(i,j)=collision_mom(2,i,j,k,n,sp,1)
                      enddo
                   enddo
                   gyro_out_a=(0.0,0.0)
                   gyro_out_b=(0.0,0.0)
                   !call gyro_average_df(gyro_in_a,gyro_out_a,k,m,n,.true.,'C')
                   !call gyro_average_df(gyro_in_b,gyro_out_b,k,m,n,.true.,'C')
                   do l=ll1,ll2
                      do j=lj1,lj2
                         do i=li1,li2
                            !loc_k(i,j,k,l,m,n)=loc_k(i,j,k,l,m,n)+x_coeff1(i,l,k,m,sp,n)*gyro_out_a(i,j)&
                            !     +x_coeff3(i,l,k,m,sp,n)*gyro_out_b(i,j)
                            loc_k(i,j,k,l,m,n)=loc_k(i,j,k,l,m,n)+x_coeff1(i,l,k,m,sp,n)*gyro_in_a(i,j)&
                                 +x_coeff3(i,l,k,m,sp,n)*gyro_in_b(i,j)
                         enddo
                      enddo
                   enddo
                   gyro_in_a=(0.0,0.0)
                   gyro_in_b=(0.0,0.0)
                   do j=lj1,lj2
                      do i=li1,li2
                         gyro_in_a(i,j)=collision_mom(1,i,j,k,sp,n,1)
                         gyro_in_b(i,j)=collision_mom(2,i,j,k,sp,n,1)
                      enddo
                   enddo
                   gyro_out_a=(0.0,0.0)
                   gyro_out_b=(0.0,0.0)
                   !call gyro_average_df(gyro_in_a,gyro_out_a,k,m,n,.true.,'C')
                   !call gyro_average_df(gyro_in_b,gyro_out_b,k,m,n,.true.,'C')
                   do l=ll1,ll2
                      do j=lj1,lj2
                         do i=li1,li2
                            !loc_k(i,j,k,l,m,n)=loc_k(i,j,k,l,m,n)+gyro_out_a(i,j)*y_coeff1(i,l,k,m,sp,n)&
                            !     +gyro_out_b(i,j)*y_coeff3(i,l,k,m,sp,n)
                            loc_k(i,j,k,l,m,n)=loc_k(i,j,k,l,m,n)+gyro_in_a(i,j)*y_coeff1(i,l,k,m,sp,n)&
                                 +gyro_in_b(i,j)*y_coeff3(i,l,k,m,sp,n)
                         enddo
                      enddo
                   enddo
                   gyro_in_a=(0.0,0.0)
                   gyro_in_b=(0.0,0.0)
                   do j=lj1,lj2
                      do i=li1,li2
                         gyro_in_a(i,j)=collision_mom2(1,i,j,k,sp)
                         gyro_in_b(i,j)=collision_mom2(2,i,j,k,sp)
                      enddo
                   enddo
                   gyro_out_a=(0.0,0.0)
                   gyro_out_b=(0.0,0.0)
                   !call gyro_average_df(gyro_in_a,gyro_out_a,k,m,n,.true.,'C')
                   !call gyro_average_df(gyro_in_b,gyro_out_b,k,m,n,.true.,'C')
                   do l=ll1,ll2
                      do j=lj1,lj2
                         do i=li1,li2
                            !loc_k(i,j,k,l,m,n)=loc_k(i,j,k,l,m,n)+gyro_out_a(i,j)*y_coeff4(i,l,k,m,sp,n)&
                            !     +gyro_out_b(i,j)*y_coeff6(i,l,k,m,sp,n)
                            loc_k(i,j,k,l,m,n)=loc_k(i,j,k,l,m,n)+gyro_in_a(i,j)*y_coeff4(i,l,k,m,sp,n)&
                                 +gyro_in_b(i,j)*y_coeff6(i,l,k,m,sp,n)
                         enddo
                      enddo
                   enddo
                enddo
             enddo
          endif
       enddo
    enddo


    do n=ln1,ln2
       do m=lm1,lm2
          do k=lk1,lk2
             gyro_in_a=(0.0,0.0)
             gyro_in_b=(0.0,0.0)
             do j=lj1,lj2
                do i=li1,li2
                   gyro_in_a(i,j)=collision_mom2(1,i,j,k,n)
                   gyro_in_b(i,j)=collision_mom2(2,i,j,k,n)
                enddo
             enddo
             gyro_out_a=(0.0,0.0)
             gyro_out_b=(0.0,0.0)
             !call gyro_average_df(gyro_in_a,gyro_out_a,k,m,n,.true.,'C')
             !call gyro_average_df(gyro_in_b,gyro_out_b,k,m,n,.true.,'C')
             do l=ll1,ll2
                do j=lj1,lj2
                   do i=li1,li2
                      !loc_k(i,j,k,l,m,n)=loc_k(i,j,k,l,m,n)+x_coeff4(i,l,k,m,n)*gyro_out_a(i,j)&
                      !     +x_coeff6(i,l,k,m,n)*gyro_out_b(i,j)
                      loc_k(i,j,k,l,m,n)=loc_k(i,j,k,l,m,n)+x_coeff4(i,l,k,m,n)*gyro_in_a(i,j)&
                           +x_coeff6(i,l,k,m,n)*gyro_in_b(i,j)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo


    do n=ln1,ln2
       do sp=0,(n_spec-1)
          do m=lm1,lm2
             do k=lk1,lk2
                gyro_in_a=(0.0,0.0)
                gyro_in_b=(0.0,0.0)
                do j=lj1,lj2
                   do i=li1,li2
                      gyro_in_a(i,j)=collision_mom(1,i,j,k,sp,n,2)
                      gyro_in_b(i,j)=collision_mom(2,i,j,k,sp,n,2)
                   enddo
                enddo
                gyro_out_a=(0.0,0.0)
                gyro_out_b=(0.0,0.0)
                !call gyro_average_df(gyro_in_a,gyro_out_a,k,m,n,.true.,'C')
                !call gyro_average_df(gyro_in_b,gyro_out_b,k,m,n,.true.,'C')
                do l=ll1,ll2
                   do j=lj1,lj2
                      do i=li1,li2
                         !loc_k(i,j,k,l,m,n)=loc_k(i,j,k,l,m,n)+gyro_out_a(i,j)*y_coeffa(i,l,k,m,sp,n)&
                         !     +gyro_out_b(i,j)*y_coeffc(i,l,k,m,sp,n)
                         loc_k(i,j,k,l,m,n)=loc_k(i,j,k,l,m,n)+gyro_in_a(i,j)*y_coeffa(i,l,k,m,sp,n)&
                              +gyro_in_b(i,j)*y_coeffc(i,l,k,m,sp,n)
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo


  ! PERFOFF

  End subroutine add_momentpart_global



  !>Deletes all arrays needed for the computation of the field particle contribution
  subroutine finalize_fieldpart
    select case (init_status_fp)
    case (1)
       deallocate(mom_cmat,en_cmat,e_fac_sum,m_fac_sum)
       deallocate(loc_mat_10,mat_01_20)
    case (2)
       deallocate(rec_buff,cons_mom_cmat,cons_en_cmat)
       deallocate(loc_mat_10,mat_01_20)
    case (3)
      if (xy_local) then
         deallocate(moments_coll1,moments_coll2,moments_coll3)
         deallocate(moments_coll4,moments_coll5,moments_coll6)
         deallocate(x_coeff1,x_coeff2,x_coeff3,x_coeff4,x_coeff5,x_coeff6)
         deallocate(y_coeff1,y_coeff2,y_coeff3,y_coeff4,y_coeff5,y_coeff6)
         deallocate(y_coeffa,y_coeffb,y_coeffc)
         deallocate(collision_mom)
         deallocate(collision_mom2)
         deallocate(j1fac)
         deallocate(coll_buff)
      else
         deallocate(moments_coll1,moments_coll3)
         deallocate(moments_coll4,moments_coll6)
         deallocate(x_coeff1,x_coeff3,x_coeff4,x_coeff6)
         deallocate(y_coeff1,y_coeff3,y_coeff4,y_coeff6)
         deallocate(y_coeffa,y_coeffc)
         deallocate(collision_mom)
         deallocate(collision_mom2)
         deallocate(gyro_in_a,gyro_in_b,gyro_out_a,gyro_out_b)
      endif

    endselect
    init_status_fp = 0
  end subroutine finalize_fieldpart


  !>Initializes the collision operator
  subroutine initialize_collisions
    real, dimension(-1:nw0):: mugl
    real, dimension(-1:nv0):: vpgl
    real:: vbetr, nuprime, velquo, fm_, Bloc
    integer:: i, j, k, l, m, n, sprime, ll, ul, lm, um, sp, pni
    real:: second_deriv_fac, mixed_deriv_fac, first_deriv_fac
    logical:: active_spec
    real:: Tprime, mi, ni, qi2
    real::tmp4,tmp5,tmp6
    real:: nu_parallel, x_a, x_b, nu_D, nu_0, alpha, kappa
    real:: erf
    integer:: mupupup,mupup,mup,mdown,mdowndown,mdowndowndown,lup,lupup,lupupup,&
              ldown,ldowndown,ldowndowndown
    integer:: highhighhigh_l,highhigh_l,high_l,low_l,lowlow_l,lowlowlow_l,&
              highhighhigh_m,highhigh_m,high_m,low_m,lowlow_m,lowlowlow_m
    integer:: lower_l,higher_l,lower_m,higher_m, l_prime, m_prime, mm

    PERFON('coll_ini')

    if (collision_op .eq. 'sugama') then


       !if (coll_order .eq. 'second') then
       !   Write (*,*) 'second order scheme is not implemented with sugama operator'
       !   STOP
       !elseif (coll_order .eq. 'fourth') then
#ifndef COLLSIX
          low_l=-2
          high_l=2
          low_m=-2
          high_m=2
          lower_l=ll1-2
          higher_l=ll2+2
          lower_m=lm1-2
          higher_m=lm2+2
          second_deriv_fac=1.0/12.0
          mixed_deriv_fac=1.0/144.0
          first_deriv_fac=1.0/12.0
       !elseif (coll_order .eq. 'order6') then
#else
          low_l=-3
          high_l=3
          low_m=-3
          high_m=3
          lower_l=ll1-3
          higher_l=ll2+3
          lower_m=lm1-3
          higher_m=lm2+3
          second_deriv_fac=1.0/90.0
          mixed_deriv_fac=1.0/3600.0
          first_deriv_fac=1.0/60.0
       !endif
#endif

       allocate(particle_diffpart(pi1:pi2,lk1:lk2,ll1:ll2,lm1:lm2,0:n_spec-1,ln1:ln2),&
            dvinverse(0:nv0-1),dmuinverse(0:nw0-1),&
            vperpcoord(pi1:pi2,lk1:lk2,0:nw0-1),&
            dvperpinverse(pi1:pi2,lk1:lk2,0:nw0-1),&
            moments_colla(pi1:pi2,ll1:ll2,lk1:lk2,lm1:lm2,0:n_spec-1,ln1:ln2),&
            moments_collc(pi1:pi2,ll1:ll2,lk1:lk2,lm1:lm2,0:n_spec-1,ln1:ln2),&
            C_ab_TO_mom(pi1:pi2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2,0:n_spec-1),&
            C_ab_TO_en(pi1:pi2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2,0:n_spec-1),&
            particle_buff(li1:li2,lj1:lj2),&
            prefacbuff1(pi1:pi2,lk1:lk2,lower_l:higher_l,lower_m:higher_m,ln1:ln2,0:n_spec-1),&
            prefacbuff2(pi1:pi2,lk1:lk2,lower_l:higher_l,lower_m:higher_m,ln1:ln2,0:n_spec-1),&
            prefacbuff3(pi1:pi2,lk1:lk2,lower_l:higher_l,lower_m:higher_m,ln1:ln2,0:n_spec-1),&
            prefacbuff4(pi1:pi2,lk1:lk2,lower_l:higher_l,lower_m:higher_m,ln1:ln2,0:n_spec-1),&
            prefacbuff5(pi1:pi2,lk1:lk2,lower_l:higher_l,lower_m:higher_m,ln1:ln2,0:n_spec-1),&
            prefacbuff6(pi1:pi2,lk1:lk2,lower_l:higher_l,lower_m:higher_m,ln1:ln2,0:n_spec-1),&
            mat_00_sub(pi1:pi2,lk1:lk2,lower_l:higher_l,lower_m:higher_m),&
            vpp((-3):(nv0+2)),&
            muu((-3):(nw0+2)),&
            stencil_second(low_l:high_l),&
            stencil_first(low_l:high_l))

       if (xy_local) then
          allocate(moments_collb(pi1:pi2,ll1:ll2,lk1:lk2,lm1:lm2,0:n_spec-1,ln1:ln2))
       endif

       do m=0,nw0-1
          if (m .eq. 0) then
             dmuinverse(m)= 1/(mu(m+1)-mu(m))
          elseif (m .eq. nw0-1) then
             dmuinverse(m)= 1/(mu(m)-mu(m-1))
          else
             dmuinverse(m)= 2/(mu(m+1)-mu(m-1))
          endif
       enddo

       do l=0,nv0-1
          if (l .eq. 0) then
             dvinverse(l)= 1/(vp(l+1)-vp(l))
          elseif (l .eq. nv0-1) then
             dvinverse(l)= 1/(vp(l)-vp(l-1))
          else
             dvinverse(l)= 2/(vp(l+1)-vp(l-1))
          endif
       enddo


       do m=0,nw0-1
          do k=lk1,lk2
             do i=pi1,pi2
                vperpcoord(i,k,m)=sqrt(mu(m)*geom%Bfield(i,pj1,k))
             enddo
          enddo
       enddo


       do m=0,nw0-1
          do k=lk1,lk2
             do i=pi1,pi2
                if (m .eq. 0) then
                   dvperpinverse(i,k,m)= 1/(vperpcoord(i,k,m+1)-vperpcoord(i,k,m))
                elseif (m .eq. nw0-1) then
                   dvperpinverse(i,k,m)= 1/(vperpcoord(i,k,m)-vperpcoord(i,k,m-1))
                else
                   dvperpinverse(i,k,m)= 2/(vperpcoord(i,k,m+1)-vperpcoord(i,k,m-1))
                endif
             enddo
          enddo
       enddo

       do m=lower_m,higher_m
          do l=lower_l,higher_l
             do k=lk1,lk2
                do i=pi1,pi2
                   if (((m<0).or.(m>nw0-1)).or.((l<0).or.(l>nv0-1))) then
                      mat_00_sub(i,k,l,m)=0.0
                   else
                      mat_00_sub(i,k,l,m) = pi * mu_weight(m)*&
                           &geom%Bfield(i,pj1,k)* vp_weight(l)
                   endif
                enddo
             enddo
          enddo
       enddo

       vpp = 0.0
       muu = 0.0
       do l=0,(nv0-1)
          vpp(l)=vp(l)
       enddo
       do m=0,(nw0-1)
          muu(m)=mu(m)
       enddo

       !quick check: could profile_mod set zeff0 correctly?
       if (zeff0.lt.0) then
          if (mype==0) write(*,"(A)") &
               &"zeff is negative and could not be set correctly by profile_mod -- exit"
          stop
       endif
       call compute_nuei(nu_ei,nustar_i,nustar_e)

       allocate(zeff_spec(pi1:pi2,0:n_spec_coll-1))
       !definition:  ne zeff = sum_i n_i (q_i/e)^2
       !in the collision frequency, zeff_spec occurs as a prefactor to n_iq_i**2
       !considering quasineutrality ne = sum_i (q_i/e)
       !zeff_spec is Zeff/(q_i/e)
       zeff_spec=1.0
       if (zeff.ne.1.0) then
          do n=0,n_spec-1
             if (spec(n)%charge.gt.0) then
                if (allocated(zeff_prof)) then
                   !ions get zeff profile
                   zeff_spec(pi1:pi2,n)=zeff_prof(pi1:pi2)/spec(n)%charge
                else
                   !ions get constant zeff0
                   zeff_spec(pi1:pi2,n)=zeff0/spec(n)%charge
                endif
             endif
          enddo
          !for possible extra adiabatic ion species, zeff_spec==zeff
          if (n_spec_coll.gt.n_spec) then
             if (allocated(zeff_prof)) then
                zeff_spec(pi1:pi2,n_spec_coll-1)=zeff_prof(pi1:pi2)
             else
                zeff_spec(pi1:pi2,n_spec_coll-1)=zeff0
             endif
          endif
       endif

       prefacbuff1= 0.0
       prefacbuff2= 0.0
       prefacbuff3= 0.0
       prefacbuff4= 0.0
       prefacbuff5= 0.0
       prefacbuff6= 0.0

       if (coll_f_fm_on) then
          if (mu_grid_type .eq. 'equidist') then
             Write (*,*) 'coll_f_fm_on is not implemented for this grid type in mu'
             STOP
          endif
       endif


       do sp=0,n_spec-1
          pni=pn1
          do n=ln1,ln2
             if (pn0.gt.1) pni=n
             do l=lower_l,higher_l
                do m=lower_m,higher_m
                   do k=lk1,lk2
                         do i=pi1,pi2

                            if (((m<0).or.(m>nw0-1)).or.((l<0).or.(l>nv0-1))) then

                               prefacbuff1(i,k,l,m,n,sp)=0.0
                               prefacbuff2(i,k,l,m,n,sp)=0.0
                               prefacbuff3(i,k,l,m,n,sp)=0.0
                               prefacbuff4(i,k,l,m,n,sp)=0.0
                               prefacbuff5(i,k,l,m,n,sp)=0.0
                               prefacbuff6(i,k,l,m,n,sp)=0.0

                            else

                               if (coll_f_fm_on) then
                                  kappa = spec(sp)%dens*spec(n)%dens_prof(i)&
                                       *((pi*spec(n)%temp_prof(i))**(-1.5))*exp(-(vp(l)**2+mu(m)*geom%Bfield(i,pj1,k)))&
                                       *spec(sp)%dens_prof(i)*coll*zeff_spec(i,sp)*((2*spec(n)%charge*spec(sp)%charge)**2)&
                                       /sqrt(spec(n)%mass*((spec(n)%temp*spec(n)%temp_prof(i))**3))
                               else
                                  kappa = spec(sp)%dens*spec(sp)%dens_prof(i)&
                                       *coll*zeff_spec(i,sp)*((2*spec(n)%charge*spec(sp)%charge)&
                                       **2)/sqrt(spec(n)%mass*((spec(n)%temp*spec(n)%temp_prof(i))**3))
                               endif

                               alpha = sqrt((spec(n)%temp*spec(n)%temp_prof(i)/spec(n)%mass)&
                                    /(spec(sp)%temp*spec(sp)%temp_prof(i)/spec(sp)%mass))

                               x_a=sqrt((vp(l)**2+mu(m)*geom%Bfield(i,pj1,k))/spec(n)%temp_prof(i))
                               x_b=x_a*alpha
                               nu_parallel= 2*alpha*(Hfunc(x_b*sqrt(2.0))*sqrt(2.0))/(x_a**2)
                               nu_D= (erf(x_b)-x_b*(Hfunc(x_b*sqrt(2.0))*sqrt(2.0)))/(x_a**3)
                               nu_0= alpha*(2/sqrt(pi))*exp(-x_b**2)/(x_a**2)

                               if (mu_grid_type .eq. 'equidist') then

                                  prefacbuff1(i,k,l,m,n,sp)=(kappa/2.0)*(nu_parallel*(vp(l)**2)+nu_D*mu(m)&
                                       *geom%Bfield(i,pj1,k))*(dvinverse(l)**2)*second_deriv_fac*mat_00_sub(i,k,l,m)

                                  prefacbuff2(i,k,l,m,n,sp)=kappa*(2*mu(m)/geom%Bfield(i,pj1,k))*(nu_D&
                                       *(vp(l)**2)+nu_parallel*mu(m)*geom%Bfield(i,pj1,k))*(dmuinverse(m)**2)&
                                       *second_deriv_fac*mat_00_sub(i,k,l,m)

                                  prefacbuff3(i,k,l,m,n,sp)=2*kappa*(nu_parallel-nu_D)*mu(m)*vp(l)*dvinverse(l)&
                                       *dmuinverse(m)*mixed_deriv_fac*mat_00_sub(i,k,l,m)

                                  prefacbuff4(i,k,l,m,n,sp)=kappa*nu_parallel*(x_a**2)*(1-alpha**2)*vp(l)&
                                       *dvinverse(l)*first_deriv_fac*mat_00_sub(i,k,l,m)

                                  prefacbuff5(i,k,l,m,n,sp)=kappa*(2*mu(m)*(nu_parallel*(x_a**2)*(1-alpha**2)&
                                       +(nu_D/2.0)*(1+(vp(l)**2)/(mu(m)*geom%Bfield(i,pj1,k))))+2*(nu_D&
                                       *(vp(l)**2)+nu_parallel*mu(m)*geom%Bfield(i,pj1,k))&
                                       /geom%Bfield(i,pj1,k))*dmuinverse(m)*first_deriv_fac*mat_00_sub(i,k,l,m)

                                  prefacbuff6(i,k,l,m,n,sp)=2*kappa*nu_0*(x_a**2)*mat_00_sub(i,k,l,m)

                               elseif (mu_grid_type .eq. 'eq_vperp') then

                                  prefacbuff1(i,k,l,m,n,sp)=(kappa/2.0)*(nu_parallel*(vp(l)**2)+nu_D&
                                       *(vperpcoord(i,k,m)**2))*(dvinverse(l)**2)*second_deriv_fac*mat_00_sub(i,k,l,m)

                                  prefacbuff2(i,k,l,m,n,sp)=(kappa/2.0)*(nu_D*(vp(l)**2)+nu_parallel&
                                       *(vperpcoord(i,k,m)**2))*(dvperpinverse(i,k,m)**2)*second_deriv_fac*mat_00_sub(i,k,l,m)

                                  prefacbuff3(i,k,l,m,n,sp)=kappa*(nu_parallel-nu_D)*vp(l)*vperpcoord(i,k,m)&
                                       *dvinverse(l)*dvperpinverse(i,k,m)*mixed_deriv_fac*mat_00_sub(i,k,l,m)

                                  if (coll_f_fm_on) then

                                     prefacbuff4(i,k,l,m,n,sp)=kappa*nu_parallel*(x_a**2)*(-(1+alpha**2))&
                                          *vp(l)*dvinverse(l)*first_deriv_fac*mat_00_sub(i,k,l,m)

                                     prefacbuff5(i,k,l,m,n,sp)=kappa*(nu_parallel*(x_a**2)*(-(1+alpha**2))+(nu_D/2.0)&
                                          *(1+(vp(l)/vperpcoord(i,k,m))**2))*vperpcoord(i,k,m)&
                                          *dvperpinverse(i,k,m)*first_deriv_fac*mat_00_sub(i,k,l,m)

                                     prefacbuff6(i,k,l,m,n,sp)=0.0

                                  else

                                     prefacbuff4(i,k,l,m,n,sp)=kappa*nu_parallel*(x_a**2)*(1-alpha**2)&
                                          *vp(l)*dvinverse(l)*first_deriv_fac*mat_00_sub(i,k,l,m)

                                     prefacbuff5(i,k,l,m,n,sp)=kappa*(nu_parallel*(x_a**2)*(1-alpha**2)+(nu_D/2.0)&
                                          *(1+(vp(l)/vperpcoord(i,k,m))**2))*vperpcoord(i,k,m)&
                                          *dvperpinverse(i,k,m)*first_deriv_fac*mat_00_sub(i,k,l,m)

                                     prefacbuff6(i,k,l,m,n,sp)=2*kappa*nu_0*(x_a**2)*mat_00_sub(i,k,l,m)

                                  endif
                               endif
                            endif
                         enddo
                   enddo
                enddo
             enddo
          enddo
       enddo

       particle_diffpart = 0.0
       moments_colla=0.0
       moments_collc=0.0

       !if (coll_order .eq. 'fourth') then
#ifndef COLLSIX
          stencil_second= (/ -1.0, 16.0, -30.0, 16.0, -1.0 /)
          stencil_first= (/ 1.0, -8.0, 0.0, 8.0, -1.0 /)
       !elseif (coll_order .eq. 'order6') then
#else
          stencil_second= (/ 1.0, -13.5, 135.0, -245.0, 135.0, -13.5, 1.0 /)
          stencil_first= (/ -1.0, 9.0, -45.0, 0.0, 45.0, -9.0, 1.0 /)
#endif
       !endif

       if (xy_local) then
          moments_collb=0.0
       endif

       do sp=0,n_spec-1
          do n=ln1,ln2
             do m=lm1,lm2
                do l=ll1,ll2
                   do k=lk1,lk2
                         do i=pi1,pi2

       Bloc = geom%Bfield(i,pj1,k) !for readability

       do l_prime=low_l,high_l
          ll=l+l_prime

          particle_diffpart(i,k,l,m,sp,n)=particle_diffpart(i,k,l,m,sp,n)&
               +stencil_second(l_prime)*prefacbuff1(i,k,ll,m,n,sp)&
               -stencil_first(l_prime)*prefacbuff4(i,k,ll,m,n,sp)


          moments_colla(i,l,k,m,sp,n)=moments_colla(i,l,k,m,sp,n)&
               +stencil_second(l_prime)*prefacbuff1(i,k,ll,m,n,sp)*vpp(ll)&
               -stencil_first(l_prime)*prefacbuff4(i,k,ll,m,n,sp)*vpp(ll)


          moments_collc(i,l,k,m,sp,n)=moments_collc(i,l,k,m,sp,n)&
               +stencil_second(l_prime)*prefacbuff1(i,k,ll,m,n,sp)*(vpp(ll)**2+Bloc*muu(m))/spec(n)%temp_prof(i)&
               -stencil_first(l_prime)*prefacbuff4(i,k,ll,m,n,sp)*(vpp(ll)**2+Bloc*muu(m))/spec(n)%temp_prof(i)

       enddo
       do m_prime=low_m,high_m
          mm=m+m_prime

          particle_diffpart(i,k,l,m,sp,n)=particle_diffpart(i,k,l,m,sp,n)&
               +stencil_second(m_prime)*prefacbuff2(i,k,l,mm,n,sp)&
               -stencil_first(m_prime)*prefacbuff5(i,k,l,mm,n,sp)


          moments_colla(i,l,k,m,sp,n)=moments_colla(i,l,k,m,sp,n)&
               +stencil_second(m_prime)*prefacbuff2(i,k,l,mm,n,sp)*vpp(l)&
               -stencil_first(m_prime)*prefacbuff5(i,k,l,mm,n,sp)*vpp(l)


          moments_collc(i,l,k,m,sp,n)=moments_collc(i,l,k,m,sp,n)&
               +stencil_second(m_prime)*prefacbuff2(i,k,l,mm,n,sp)*(vpp(l)**2+Bloc*muu(mm))/spec(n)%temp_prof(i)&
               -stencil_first(m_prime)*prefacbuff5(i,k,l,mm,n,sp)*(vpp(l)**2+Bloc*muu(mm))/spec(n)%temp_prof(i)

       enddo
       do l_prime=low_l,high_l
          do m_prime=low_m,high_m
             ll=l+l_prime
             mm=m+m_prime

             particle_diffpart(i,k,l,m,sp,n)=particle_diffpart(i,k,l,m,sp,n)&
                  +stencil_first(l_prime)*stencil_first(m_prime)&
                  *prefacbuff3(i,k,ll,mm,n,sp)

             moments_colla(i,l,k,m,sp,n)=moments_colla(i,l,k,m,sp,n)&
                  +stencil_first(l_prime)*stencil_first(m_prime)&
                  *prefacbuff3(i,k,ll,mm,n,sp)*vpp(ll)

             moments_collc(i,l,k,m,sp,n)=moments_collc(i,l,k,m,sp,n)&
                  +stencil_first(l_prime)*stencil_first(m_prime)&
                  *prefacbuff3(i,k,ll,mm,n,sp)*(vpp(ll)**2+Bloc*muu(mm))/spec(n)%temp_prof(i)

          enddo
       enddo

       particle_diffpart(i,k,l,m,sp,n)=particle_diffpart(i,k,l,m,sp,n)&
            +prefacbuff6(i,k,l,m,n,sp)

       moments_colla(i,l,k,m,sp,n)=moments_colla(i,l,k,m,sp,n)&
            +prefacbuff6(i,k,l,m,n,sp)*vpp(l)

       moments_collc(i,l,k,m,sp,n)=moments_collc(i,l,k,m,sp,n)&
            +prefacbuff6(i,k,l,m,n,sp)*(vpp(l)**2+Bloc*muu(m))/spec(n)%temp_prof(i)


                         enddo
                   enddo
                enddo
             enddo
          enddo
       enddo


    if (xy_local) then
       do sp=0,n_spec-1
          do n=ln1,ln2
             do m=lm1,lm2
                do l=ll1,ll2
                   do k=lk1,lk2
                         do i=pi1,pi2

      Bloc = geom%Bfield(i,pj1,k) !for readability

        moments_collb(i,l,k,m,sp,n)=moments_collb(i,l,k,m,sp,n)&
             +(sqrt(mu(m)*Bloc)/vp(l))*moments_colla(i,l,k,m,sp,n)


                         enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
   endif

       deallocate(prefacbuff1,prefacbuff2,prefacbuff3,prefacbuff4,prefacbuff5,prefacbuff6)
       deallocate(vpp,muu)

       allocate(particle_coeff(pi1:pi2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2),&
          sigma_bar(pi1:pi2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2),&
          sigma_bar_int(pi1:pi2,lk1:lk2,ln1:ln2),&
          sigma_int_num(pi1:pi2,lk1:lk2,ln1:ln2),&
          sigma_int_den(pi1:pi2,lk1:lk2,ln1:ln2))


       !!Particle Conservation part


       sigma_int_num= 0.0
       sigma_int_den= 0.0

       pni=pn1
       do n=ln1,ln2
          if (pn0.gt.1) pni=n
          do l=ll1,ll2
             do m=lm1,lm2
                do k=lk1,lk2
                      do i=pi1,pi2

                         x_a=sqrt((vp(l)**2+mu(m)*geom%Bfield(i,pj1,k))/spec(n)%temp_prof(i))

                         sigma_int_num(i,k,n)= sigma_int_num(i,k,n)&
                              +((x_a**2)*(x_a-8/(3*sqrt(pi)))*fm(i,pj1,k,l,m,pni))*mat_00_sub(i,k,l,m)

                         sigma_int_den(i,k,n)= sigma_int_den(i,k,n)&
                              +(x_a**2)*fm(i,pj1,k,l,m,pni)*mat_00_sub(i,k,l,m)

                      enddo
                enddo
             enddo
          enddo
       enddo

       call my_real_sum_vw(sigma_int_num,size(sigma_int_num))
       call my_real_sum_vw(sigma_int_den,size(sigma_int_den))

       pni=pn1
       do n=ln1,ln2
          if (pn0.gt.1) pni=n
          do l=ll1,ll2
             do m=lm1,lm2
                do k=lk1,lk2
                      do i=pi1,pi2

                         x_a=sqrt((vp(l)**2+mu(m)*geom%Bfield(i,pj1,k))/spec(n)%temp_prof(i))

                         sigma_bar(i,k,l,m,n)=(x_a-8/(3*sqrt(pi)))*fm(i,pj1,k,l,m,pni)&
                              -sigma_int_num(i,k,n)*fm(i,pj1,k,l,m,pni)&
                              /sigma_int_den(i,k,n)

                      enddo
                enddo
             enddo
          enddo
       enddo

       sigma_bar_int = 0.0

       do n=ln1,ln2
          do l=ll1,ll2
             do m=lm1,lm2
                do k=lk1,lk2
                   do i=pi1,pi2

                      sigma_bar_int(i,k,n)=sigma_bar_int(i,k,n)&
                           +sigma_bar(i,k,l,m,n)*mat_00_sub(i,k,l,m)

                   enddo
                enddo
             enddo
          enddo
       enddo

       call my_real_sum_vw(sigma_bar_int,size(sigma_bar_int))

       pni=pn1
       do n=ln1,ln2
          if (pn0.gt.1) pni=n
          do l=ll1,ll2
             do m=lm1,lm2
                do k=lk1,lk2
                   do i=pi1,pi2
                      particle_coeff(i,k,l,m,n)= (sigma_bar(i,k,l,m,n)/sigma_bar_int(i,k,n))
                   enddo
                enddo
             enddo
          enddo
       enddo

       deallocate(sigma_bar,sigma_bar_int,sigma_int_num,sigma_int_den)
       deallocate(mat_00_sub)


       allocate(prefacterm1(pi1:pi2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2,0:n_spec-1),&
            prefacterm2(pi1:pi2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2,0:n_spec-1),&
            prefacterm3(pi1:pi2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2,0:n_spec-1),&
            prefacterm4(pi1:pi2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2,0:n_spec-1),&
            prefacterm5(pi1:pi2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2,0:n_spec-1),&
            prefacterm6(pi1:pi2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2,0:n_spec-1),&
            momentum_input(pi1:pi2,lk1:lk2,(ll1-3):(ll2+3),(lm1-3):(lm2+3),ln1:ln2),&
            energy_input(pi1:pi2,lk1:lk2,(ll1-3):(ll2+3),(lm1-3):(lm2+3),ln1:ln2))


       prefacterm1=0.0
       prefacterm2=0.0
       prefacterm3=0.0
       prefacterm4=0.0
       prefacterm5=0.0
       prefacterm6=0.0

       do sp=0,n_spec-1
          pni=pn1
          do n=ln1,ln2
             if (pn0.gt.1) pni=n
             do l=ll1,ll2
                do m=lm1,lm2
                   do k=lk1,lk2
                         do i=pi1,pi2

                            alpha = sqrt((spec(n)%temp*spec(n)%temp_prof(i)/spec(n)%mass)&
                                    /(spec(sp)%temp*spec(sp)%temp_prof(i)/spec(sp)%mass))

                            if (coll_f_fm_on) then
                               kappa = spec(sp)%dens*spec(n)%dens_prof(i)&
                                    *((pi*spec(n)%temp_prof(i))**(-1.5))*exp(-(vp(l)**2+mu(m)*geom%Bfield(i,pj1,k)))&
                                    *spec(sp)%dens_prof(i)*coll*zeff_spec(i,sp)*((2*spec(n)%charge*spec(sp)%charge)**2)&
                                    /sqrt(spec(n)%mass*((spec(n)%temp*spec(n)%temp_prof(i))**3))
                            else
                               kappa = spec(sp)%dens*spec(sp)%dens_prof(i)&
                                    *coll*zeff_spec(i,sp)*((2.0*spec(n)%charge*spec(sp)%charge)&
                                    **2)/sqrt(spec(n)%mass*((spec(n)%temp*spec(n)%temp_prof(i))**3))
                            endif

                            x_a=sqrt((vp(l)**2+mu(m)*geom%Bfield(i,pj1,k))/spec(n)%temp_prof(i))
                            x_b=x_a*alpha
                            nu_parallel= 2*alpha*(Hfunc(x_b*sqrt(2.0))*sqrt(2.0))/(x_a**2)
                            nu_D= (erf(x_b)-x_b*(Hfunc(x_b*sqrt(2.0))*sqrt(2.0)))/(x_a**3)
                            nu_0= alpha*(2/sqrt(pi))*exp(-x_b**2)/(x_a**2)

                            if (mu_grid_type .eq. 'equidist') then

                               prefacterm1(i,k,l,m,n,sp)=prefacterm1(i,k,l,m,n,sp)+(kappa/2.0)*(nu_parallel*(vp(l)**2)&
                                    +nu_D*mu(m)*geom%Bfield(i,pj1,k))*(dvinverse(l)**2)*second_deriv_fac

                               prefacterm2(i,k,l,m,n,sp)=prefacterm2(i,k,l,m,n,sp)+kappa*(2*mu(m)/geom%Bfield(i,pj1,k))&
                                    *(nu_D*(vp(l)**2)+nu_parallel*mu(m)*geom%Bfield(i,pj1,k))*(dmuinverse(m)**2)*second_deriv_fac

                               prefacterm3(i,k,l,m,n,sp)=prefacterm3(i,k,l,m,n,sp)+2*kappa*(nu_parallel-nu_D)*mu(m)&
                                    *vp(l)*(dvinverse(l)*dmuinverse(m))*mixed_deriv_fac

                               prefacterm4(i,k,l,m,n,sp)=prefacterm4(i,k,l,m,n,sp)+kappa*nu_parallel*(x_a**2)*(1-alpha**2)&
                                    *vp(l)*(dvinverse(l))*first_deriv_fac

                               prefacterm5(i,k,l,m,n,sp)=prefacterm5(i,k,l,m,n,sp)+kappa*(2*mu(m)&
                                    *(nu_parallel*(x_a**2)*(1-alpha**2)&
                                    +(nu_D/2.0)*(1+(vp(l)**2)/(mu(m)*geom%Bfield(i,pj1,k))))+2*(nu_D&
                                    *(vp(l)**2)+nu_parallel*mu(m)*geom%Bfield(i,pj1,k))&
                                    /geom%Bfield(i,pj1,k))*(dmuinverse(m))*first_deriv_fac

                               prefacterm6(i,k,l,m,n,sp)=prefacterm6(i,k,l,m,n,sp)+2*kappa*nu_0*(x_a**2)

                            elseif (mu_grid_type .eq. 'eq_vperp') then

                               prefacterm1(i,k,l,m,n,sp)=prefacterm1(i,k,l,m,n,sp)+(kappa/2.0)*(nu_parallel&
                                    *(vp(l)**2)+nu_D*(vperpcoord(i,k,m)**2))*(dvinverse(l)**2)*second_deriv_fac

                               prefacterm2(i,k,l,m,n,sp)=prefacterm2(i,k,l,m,n,sp)+(kappa/2.0)*(nu_D*(vp(l)**2)+nu_parallel&
                                    *(vperpcoord(i,k,m)**2))*(dvperpinverse(i,k,m)**2)*second_deriv_fac

                               prefacterm3(i,k,l,m,n,sp)=prefacterm3(i,k,l,m,n,sp)&
                                    +kappa*(nu_parallel-nu_D)*vp(l)*vperpcoord(i,k,m)&
                                    *(dvinverse(l)*dvperpinverse(i,k,m))*mixed_deriv_fac

                               if (coll_f_fm_on) then

                                  prefacterm4(i,k,l,m,n,sp)=prefacterm4(i,k,l,m,n,sp)&
                                       +kappa*nu_parallel*(x_a**2)*(-(1+alpha**2))&
                                       *vp(l)*(dvinverse(l))*first_deriv_fac

                                  prefacterm5(i,k,l,m,n,sp)=prefacterm5(i,k,l,m,n,sp)&
                                       +kappa*(nu_parallel*(x_a**2)*(-(1+alpha**2))+(nu_D/2.0)&
                                       *(1+(vp(l)/vperpcoord(i,k,m))**2))*vperpcoord(i,k,m)&
                                       *(dvperpinverse(i,k,m))*first_deriv_fac

                                  prefacterm6(i,k,l,m,n,sp)=0.0

                               else

                                  prefacterm4(i,k,l,m,n,sp)=prefacterm4(i,k,l,m,n,sp)+kappa*nu_parallel*(x_a**2)*(1-alpha**2)&
                                       *vp(l)*(dvinverse(l))*first_deriv_fac

                                  prefacterm5(i,k,l,m,n,sp)=prefacterm5(i,k,l,m,n,sp)&
                                       +kappa*(nu_parallel*(x_a**2)*(1-alpha**2)+(nu_D/2.0)&
                                       *(1+(vp(l)/vperpcoord(i,k,m))**2))*vperpcoord(i,k,m)&
                                       *(dvperpinverse(i,k,m))*first_deriv_fac

                                  prefacterm6(i,k,l,m,n,sp)=prefacterm6(i,k,l,m,n,sp)+2*kappa*nu_0*(x_a**2)

                               endif
                            else
                               if (collision_op .eq. 'sugama') then
                                  Write (*,*) 'This collision scheme is not compatible with this grid type!!'
                                  STOP
                               endif
                            endif
                         enddo
                   enddo
                enddo
             enddo
          enddo
       enddo

       momentum_input = 0.0
       energy_input = 0.0

       pni=pn1
       do n=ln1,ln2
          if (pn0.gt.1) pni=n
          do m=lm1-3,lm2+3
             do l=ll1-3,ll2+3
                do k=lk1,lk2
                      do i=pi1,pi2

                         if ((l .lt. 0).or.(l .gt. nv0-1).or.(m .lt. 0).or.(m .gt. nw0-1)) then
                            momentum_input(i,k,l,m,n) = 0.0
                            energy_input(i,k,l,m,n) = 0.0
                         else

                            x_a=sqrt((vp(l)**2+geom%Bfield(i,pj1,k)*mu(m))/spec(n)%temp_prof(i))
                            if (coll_f_fm_on) then
                               momentum_input(i,k,l,m,n)=spec(n)%mass*vp(l)/spec(n)%temp
                               energy_input(i,k,l,m,n)=(x_a**2)
                            else
                               momentum_input(i,k,l,m,n)=spec(n)%dens_prof(i)*((pi*spec(n)%temp_prof(i))**(-1.5))&
                                    *exp(-x_a**2)*spec(n)%mass*vp(l)/spec(n)%temp

                               energy_input(i,k,l,m,n)=spec(n)%dens_prof(i)*((pi*spec(n)%temp_prof(i))**(-1.5))&
                                    *exp(-x_a**2)*(x_a**2)
                            endif
                         endif
                      enddo
                enddo
             enddo
          enddo
       enddo

       C_ab_TO_mom = 0.0
       C_ab_TO_en  = 0.0


       !if (coll_order .eq. 'fourth') then
#ifndef COLLSIX
          prefacterm6=prefacterm6-30*(prefacterm1+prefacterm2)
#else
       !elseif (coll_order .eq. 'order6') then
          prefacterm6=prefacterm6-245*(prefacterm1+prefacterm2)
#endif
       !endif

       !if (coll_order .eq. 'fourth') then
#ifndef COLLSIX
          stencil_second= (/ -1.0, 16.0, 0.0, 16.0, -1.0 /)
          stencil_first= (/ 1.0, -8.0, 0.0, 8.0, -1.0 /)
       !elseif (coll_order .eq. 'order6') then
#else
          stencil_second= (/ 1.0, -13.5, 135.0, 0.0, 135.0, -13.5, 1.0 /)
          stencil_first= (/ -1.0, 9.0, -45.0, 0.0, 45.0, -9.0, 1.0 /)
#endif
       !endif

       do sp=0,n_spec-1
          do n=ln1,ln2
             do m=lm1,lm2
                do l=ll1,ll2
                   do k=lk1,lk2

       do l_prime=low_l,high_l
          ll=l_prime+l

          C_ab_TO_mom(:,k,l,m,n,sp)=C_ab_TO_mom(:,k,l,m,n,sp)&
               +stencil_second(l_prime)*momentum_input(:,k,ll,m,n)*prefacterm1(:,k,l,m,n,sp)&
               +stencil_first(l_prime)*momentum_input(:,k,ll,m,n)*prefacterm4(:,k,l,m,n,sp)

          C_ab_TO_en(:,k,l,m,n,sp)=C_ab_TO_en(:,k,l,m,n,sp)&
               +stencil_second(l_prime)*energy_input(:,k,ll,m,n)*prefacterm1(:,k,l,m,n,sp)&
               +stencil_first(l_prime)*energy_input(:,k,ll,m,n)*prefacterm4(:,k,l,m,n,sp)

       enddo
       do m_prime=low_m,high_m
          mm=m_prime+m

          C_ab_TO_mom(:,k,l,m,n,sp)=C_ab_TO_mom(:,k,l,m,n,sp)&
               +stencil_second(m_prime)*momentum_input(:,k,l,mm,n)*prefacterm2(:,k,l,m,n,sp)&
               +stencil_first(m_prime)*momentum_input(:,k,l,mm,n)*prefacterm5(:,k,l,m,n,sp)

          C_ab_TO_en(:,k,l,m,n,sp)=C_ab_TO_en(:,k,l,m,n,sp)&
               +stencil_second(m_prime)*energy_input(:,k,l,mm,n)*prefacterm2(:,k,l,m,n,sp)&
               +stencil_first(m_prime)*energy_input(:,k,l,mm,n)*prefacterm5(:,k,l,m,n,sp)

       enddo
       do l_prime=low_l,high_l
          do m_prime=low_m,high_m
             ll=l_prime+l
             mm=m_prime+m

             C_ab_TO_mom(:,k,l,m,n,sp)=C_ab_TO_mom(:,k,l,m,n,sp)&
                  +stencil_first(m_prime)*stencil_first(l_prime)&
                  *momentum_input(:,k,ll,mm,n)*prefacterm3(:,k,l,m,n,sp)

             C_ab_TO_en(:,k,l,m,n,sp)=C_ab_TO_en(:,k,l,m,n,sp)&
                  +stencil_first(m_prime)*stencil_first(l_prime)&
                  *energy_input(:,k,ll,mm,n)*prefacterm3(:,k,l,m,n,sp)

          enddo
       enddo

       C_ab_TO_mom(:,k,l,m,n,sp)=C_ab_TO_mom(:,k,l,m,n,sp)&
            +momentum_input(:,k,l,m,n)*prefacterm6(:,k,l,m,n,sp)

       C_ab_TO_en(:,k,l,m,n,sp)=C_ab_TO_en(:,k,l,m,n,sp)&
            +energy_input(:,k,l,m,n)*prefacterm6(:,k,l,m,n,sp)


                   enddo
                enddo
             enddo
          enddo
       enddo

       deallocate(stencil_first,stencil_second)
       deallocate(momentum_input,energy_input)
       deallocate(dvinverse,dmuinverse,vperpcoord,dvperpinverse)

       allocate(T_coeff1(pi1:pi2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2),&
            T_coeff2(pi1:pi2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2),&
            T_coeff3(pi1:pi2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2),&
            T_coeff4(pi1:pi2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2),&
            T_coeff5(pi1:pi2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2),&
            T_coeff6(pi1:pi2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2),&
            particle_fac(li1:li2,lj1:lj2,lk1:lk2,0:n_spec-1,ln1:ln2))


       T_coeff1=prefacterm1(:,:,:,:,:,0)
       T_coeff2=prefacterm2(:,:,:,:,:,0)
       T_coeff3=prefacterm3(:,:,:,:,:,0)
       T_coeff4=prefacterm4(:,:,:,:,:,0)
       T_coeff5=prefacterm5(:,:,:,:,:,0)
       T_coeff6=prefacterm6(:,:,:,:,:,0)


       if (n_spec .gt. 1) then
          do sp=1,n_spec-1
             T_coeff1=T_coeff1+prefacterm1(:,:,:,:,:,sp)
             T_coeff2=T_coeff2+prefacterm2(:,:,:,:,:,sp)
             T_coeff3=T_coeff3+prefacterm3(:,:,:,:,:,sp)
             T_coeff4=T_coeff4+prefacterm4(:,:,:,:,:,sp)
             T_coeff5=T_coeff5+prefacterm5(:,:,:,:,:,sp)
             T_coeff6=T_coeff6+prefacterm6(:,:,:,:,:,sp)
          enddo
       endif

       deallocate(prefacterm1,prefacterm2,prefacterm3,prefacterm4,prefacterm5,prefacterm6)

       if (coll_f_fm_on)  then
         allocate(loc_f_fm(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2))
         loc_f_fm=0.
       end if

    else

       !quick check: could profile_mod set zeff0 correctly?
       if (zeff0.lt.0) then
          if (mype==0) write(*,"(A)") &
               &"zeff is negative and could not be set correctly by profile_mod -- exit"
          stop
       endif
       call compute_nuei(nu_ei,nustar_i,nustar_e)

       allocate(zeff_spec(pi1:pi2,0:n_spec_coll-1))
       !definition:  ne zeff = sum_i n_i (q_i/e)^2
       !in the collision frequency, zeff_spec occurs as a prefactor to n_iq_i**2
       !considering quasineutrality ne = sum_i (q_i/e)
       !zeff_spec is Zeff/(q_i/e)
       zeff_spec=1.0
       if (zeff.ne.1.0) then
          do n=0,n_spec-1
             if (spec(n)%charge.gt.0) then
                if (allocated(zeff_prof)) then
                   !ions get zeff profile
                   zeff_spec(pi1:pi2,n)=zeff_prof(pi1:pi2)/spec(n)%charge
                else
                   !ions get constant zeff0
                   zeff_spec(pi1:pi2,n)=zeff0/spec(n)%charge
                endif
             endif
          enddo
          !for possible extra adiabatic ion species, zeff_spec==zeff
          if (n_spec_coll.gt.n_spec) then
             if (allocated(zeff_prof)) then
                zeff_spec(pi1:pi2,n_spec_coll-1)=zeff_prof(pi1:pi2)
             else
                zeff_spec(pi1:pi2,n_spec_coll-1)=zeff0
             endif
          endif
       endif

       if (coll_f_fm_on)  then
         allocate(loc_f_fm(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2))
         loc_f_fm=0.
       end if

       allocate(cd11(pi1:pi2,lk1:lk2,ll1-1:ll2,lm1:lm2,0:n_spec_coll-1,ln1:ln2),&
            cd12(pi1:pi2,lk1:lk2,ll1-1:ll2,lm1:lm2,0:n_spec_coll-1,ln1:ln2),  &
            cd21(pi1:pi2,lk1:lk2,ll1:ll2,lm1-1:lm2,0:n_spec_coll-1,ln1:ln2),&
            cd22(pi1:pi2,lk1:lk2,ll1:ll2,lm1-1:lm2,0:n_spec_coll-1,ln1:ln2),  &
            cr1(pi1:pi2,lk1:lk2,ll1-1:ll2,lm1:lm2,0:n_spec_coll-1,ln1:ln2),&
            cr2(pi1:pi2,lk1:lk2,ll1:ll2,lm1-1:lm2,0:n_spec_coll-1,ln1:ln2), &
            cvf_mu(lm1-1:lm2+1),dmuin(lm1-1:lm2),cvf_vp(ll1-1:ll2+1),&
            dvinv(ll1-1:ll2+1))
       cd11=0.
       cd12=0.
       cd21=0.
       cd22=0.
       cr1=0.
       cr2=0.

       !construct global mu
       mugl(0:nw0-1)=mu
       mugl(-1)=0.
       mugl(nw0)=lw

       !construct global vp
       vpgl=0.
       vpgl(0:nv0-1)=vp

       !collision matrix
       do n=ln1,ln2
          do sprime=0,n_spec_coll-1
             if (sprime.lt.n_spec) active_spec=.not.spec(sprime)%passive
             if (sprime.ge.n_spec) active_spec=.true.
             if (active_spec) then
                do k=lk1,lk2
                   do j=pj1,pj2
                      do i=pi1,pi2
                         if (sprime.lt.n_spec) then !kinetic species
                            Tprime=spec(sprime)%temp*spec(sprime)%temp_prof(i)
                            nuprime=spec(n)%charge**2*zeff_spec(i,sprime)&
                                 *spec(sprime)%charge**2*sqrt(spec(n)%mass)&
                                 *spec(sprime)%dens*spec(sprime)%dens_prof(i)&
                                 *spec(sprime)%temp*spec(sprime)%temp_prof(i)&
                                 /((spec(n)%temp)**2.5*spec(sprime)%mass)
                            !x=v/v_j'(x) -> velquo has profile in sprime
                            velquo=sqrt((spec(n)%temp*spec(sprime)%mass)&
                                 /(spec(sprime)%temp*spec(sprime)%temp_prof(i)*spec(n)%mass))
                         else !extra species for adiabatic ions
                            !parameters relative to active (electron) species n
                            !Ti = Zeff*Te/tau  !including profile
                            !mi = adiabatic_mass*(m_p/m_e)*\hat{m}_e (adiabatic_mass is m_i/m_p)
                            !ni=ne              !actually tau (Zeff) has the density info in the adiabatic response
                            !qi2=qe**2*Zeff
                            !mi=adiabatic_mass*m_proton/m_electron*spec(n)%mass
                            mi=adiabatic_mass*m_proton/m_electron*spec(n)%mass
                            qi2=zeff_spec(i,sprime)*spec(n)%charge**2
                            Tprime=spec(n)%temp*spec(n)%temp_prof(i)*zeff_spec(i,sprime)/tau
                            ni=spec(n)%dens*spec(n)%dens_prof(i)
                            nuprime=spec(n)%charge**2*qi2*sqrt(spec(n)%mass)&
                                 *ni*Tprime/((spec(n)%temp)**2.5*mi)
                            velquo=sqrt((spec(n)%temp*mi)/(Tprime*spec(n)%mass))
                         endif
                         do m=lm1,lm2
                            do l=ll1-1,ll2
                               vbetr=sqrt(((vpgl(l)+vpgl(l+1))/2)**2+geom%Bfield(i,pj1,k)*mu(m))
                               cd11(i,k,l,m,sprime,n)= coll*nuprime*(mu(m)*geom%Bfield(i,pj1,k)*cfk1(velquo*vbetr)&
                                    +0.25*(vpgl(l)+vpgl(l+1))**2*cfk3(velquo*vbetr))/vbetr**5
                               cd12(i,k,l,m,sprime,n)=coll*nuprime&
                                    *6*(vpgl(l)+vpgl(l+1))/2*mu(m)*cfk2(velquo*vbetr)/vbetr**5
                               cr1(i,k,l,m,sprime,n)=coll*nuprime*2.&
                                    *spec(n)%temp/Tprime &
                                    *cfk3(velquo*vbetr)*(vpgl(l)+vpgl(l+1))/2/vbetr**3
                               if (coll_f_fm_on) then
                                  fm_  = spec(n)%dens_prof(i)*(pi*spec(n)%temp_prof(i))**(-1.5)&
                                       *exp(-vbetr**2/spec(n)%temp_prof(i))
                                  cd11(i,k,l,m,sprime,n) = cd11(i,k,l,m,sprime,n)*fm_
                                  cd12(i,k,l,m,sprime,n) = cd12(i,k,l,m,sprime,n)*fm_
                               If (em_cons==3) then
                               cr1(i,k,l,m,sprime,n)=0
                               Else
                                  cr1(i,k,l,m,sprime,n) =  cr1(i,k,l,m,sprime,n)*fm_&
                                    *(1-Tprime/(spec(n)%temp*spec(n)%temp_prof(i)))
                               end if
                               end if
                            enddo
                         enddo
                         do m=lm1-1,lm2
                            do l=ll1,ll2
                               vbetr=sqrt((mugl(m)+mugl(m+1))/2*geom%Bfield(i,pj1,k)+vp(l)**2)
                               cd21(i,k,l,m,sprime,n)=coll*nuprime&
                                    *6*vp(l)*(mugl(m)+mugl(m+1))/2*cfk2(velquo*vbetr)/vbetr**5
                               cd22(i,k,l,m,sprime,n)= coll*nuprime&
                                    *(4*vp(l)**2*(mugl(m)+mugl(m+1))/2/geom%Bfield(i,pj1,k)*cfk1(velquo*vbetr)+ &
                                    4*((mugl(m)+mugl(m+1))/2)**2*cfk3(velquo*vbetr))/vbetr**5
                               cr2(i,k,l,m,sprime,n)=coll*nuprime*4.&
                                    *spec(n)%temp/Tprime&
                                    *cfk3(velquo*vbetr)*(mugl(m)+mugl(m+1))/2/vbetr**3
                               if (coll_f_fm_on) then
                                  fm_  = spec(n)%dens_prof(i)*(pi*spec(n)%temp_prof(i))**(-1.5)&
                                       *exp(-vbetr**2/spec(n)%temp_prof(i))
                                  cd21(i,k,l,m,sprime,n) = cd21(i,k,l,m,sprime,n)*fm_
                                  cd22(i,k,l,m,sprime,n) = cd22(i,k,l,m,sprime,n)*fm_
                               If (em_cons==3) then
                                    cr2(i,k,l,m,sprime,n)=0
                               Else
                                  cr2 (i,k,l,m,sprime,n) = cr2 (i,k,l,m,sprime,n)*fm_&
                                      *(1-Tprime/(spec(n)%temp*spec(n)%temp_prof(i)))
                               end if
                               end if
                            enddo
                         enddo
                      enddo
                   enddo
                enddo
             endif
          enddo
       enddo

       !pitch-all:   modifies rosenbluth potentials, keeps i-e ,e-e, i-i terms
       !pitch-angle: modifies rosenbluth potentials, deletes i-i and i-e collisions, replaces e-e term by "hee"-model
       !pitch-ion:   modifies rosenbluth potentials, deletes i-i, i-e and e-e terms

       if ((collision_op.eq.'pitch-angle').or.(collision_op.eq.'pitch-ion')) then
          do n=ln1,ln2
             do sprime=0,n_spec-1
                if ((spec(n)%charge.lt.0).and.(spec(sprime)%charge.gt.0)) then
                   if (collision_op.eq.'pitch-angle') then !add ee-aproximation term
                      do m=lm1,lm2
                         do l=ll1-1,ll2
                            do k=lk1,lk2
                               do j=pj1,pj2
                                  do i=pi1,pi2
                                     vbetr=sqrt(((vpgl(l)+vpgl(l+1))/2)**2+geom%Bfield(i,pj1,k)*mu(m))
                                     cd11(i,k,l,m,sprime,n)= cd11(i,k,l,m,sprime,n)*(1+hee(vbetr)&
                                          /zeff_spec(i,sprime)/spec(sprime)%charge)
                                     cd12(i,k,l,m,sprime,n)= cd12(i,k,l,m,sprime,n)*(1+hee(vbetr)&
                                          /zeff_spec(i,sprime)/spec(sprime)%charge)
                                  enddo
                               enddo
                            enddo
                         enddo
                      enddo
                      do m=lm1-1,lm2
                         do l=ll1,ll2
                            do k=lk1,lk2
                               do i=pi1,pi2
                                  vbetr=sqrt((mugl(m)+mugl(m+1))/2*geom%Bfield(i,pj1,k)+vp(l)**2)
                                  cd21(i,k,l,m,sprime,n)= cd21(i,k,l,m,sprime,n)*(1+hee(vbetr)&
                                       /zeff_spec(i,sprime)/spec(sprime)%charge)
                                  cd22(i,k,l,m,sprime,n)= cd22(i,k,l,m,sprime,n)*(1+hee(vbetr)&
                                       /zeff_spec(i,sprime)/spec(sprime)%charge)
                               enddo
                            enddo
                         enddo
                      enddo
                   endif
                elseif (spec(n)%charge.lt.0.and.spec(sprime)%charge.lt.0) then
                   !delete e-e collisions that are "hee"-modeled in case of pitch-angle
                   cd11(:,:,:,:,sprime,n)=0.
                   cd12(:,:,:,:,sprime,n)=0.
                   cd21(:,:,:,:,sprime,n)=0.
                   cd22(:,:,:,:,sprime,n)=0.
                elseif (spec(n)%charge.gt.0.and.spec(sprime)%charge.gt.0) then
                   !delete i-i collisions
                   cd11(:,:,:,:,sprime,n)=0.
                   cd12(:,:,:,:,sprime,n)=0.
                   cd21(:,:,:,:,sprime,n)=0.
                   cd22(:,:,:,:,sprime,n)=0.
                elseif ((spec(n)%charge.gt.0).and.(spec(sprime)%charge.lt.0)) then
                   !delete the ions scattering off electrons part
                   cd11(:,:,:,:,sprime,n)=0.
                   cd12(:,:,:,:,sprime,n)=0.
                   cd21(:,:,:,:,sprime,n)=0.
                   cd22(:,:,:,:,sprime,n)=0.
                endif
             enddo
          enddo
       endif

      ! deallocate(zeff_spec)

       !only electrons scattering off ions are considered
       if (collision_op.eq.'only-elec') then
          do n=ln1,ln2
             do sprime=0,n_spec_coll-1
                !if ((spec(n)%name.eq.'electrons').and.(spec(sprime)%name.eq.'ions')) then
                if ((spec(n)%charge.lt.0).and.(spec(sprime)%charge.gt.0)) then
                   !keep landau e-i collisions
                else
                   cd11(:,:,:,:,sprime,n)=0.
                   cd12(:,:,:,:,sprime,n)=0.
                   cd21(:,:,:,:,sprime,n)=0.
                   cd22(:,:,:,:,sprime,n)=0.
                endif
             enddo
          enddo
       endif

       !length of control volume for flux calculation
       !do l=ll1,ll2
       !   cvf_vp(l)=1/(vpgl(l+1)-vpgl(l))
       !enddo
       !do m=lm1,lm2
       !   cvf_mu(m)=2/(mugl(m+1)-mugl(m-1))
       !end do

       !(incl.boundaries)
       ll=ll1-1
       if(my_pev.eq.0) ll=ll1
       ul=ll2+1
       if(my_pev.eq.n_procs_v-1) ul = ll2
       cvf_vp(ll:ul)=1./vp_weight(ll:ul)

       lm=lm1-1
       if(my_pew.eq.0) lm=lm1
       um=lm2+1
       if(my_pew.eq.n_procs_w-1) um = lm2
       cvf_mu(lm:um)=1./mu_weight(lm:um)

       !factors for central differences
       dvin=1/dv
       do l=ll1,ll2
          dvinv(l)=1./(vpgl(l+1)-vpgl(l))
       enddo
       do m=lm1-1,lm2
          dmuin(m)=1./(mugl(m+1)-mugl(m))
       end do

       !boundaries

       if (my_pev.eq.0) then
          cd11(:,:,ll1-1,:,:,:)=0.0
          cd12(:,:,ll1-1,:,:,:)=0.0
          cr1(:,:,ll1-1,:,:,:)=0.0
          dvinv(ll1)=2/dv
          cvf_vp(ll1-1)=0
       endif

       if (my_pev.eq.n_procs_v-1) then
          cd11(:,:,ll2,:,:,:)=0.0
          cd12(:,:,ll2,:,:,:)=0.0
          cr1(:,:,ll2,:,:,:)=0.0
          dvinv(ll2)=2/dv
          cvf_vp(ll2+1)=0
       endif

       if (my_pew.eq.0) then
          cd21(:,:,:,lm1-1,:,:)=0.0
          cd22(:,:,:,lm1-1,:,:)=0.0
          cr2(:,:,:,lm1-1,:,:)=0.0
          dmuin(lm1-1)=-dmuin(lm1)
          cvf_mu(lm1-1)=0.
       endif

       if (my_pew.eq.n_procs_w-1) then
          cd21(:,:,:,lm2,:,:)=0.0
          cd22(:,:,:,lm2,:,:)=0.0
          cr2(:,:,:,lm2,:,:)=0.0
          dmuin(lm2)=-dmuin(lm2-1)
          cvf_mu(lm2+1)=0.
       endif

    endif


    if (collision_op .ne. 'sugama') then

       !this computes maximum matrix element of the collision operator
       !which is in all tests found on the diagonal(sten=5) for n=electrons at z=0(outboard)
       !(copied from initialize_testpart_3, boundaries matter!)
        ev_coll_est = 0.0

       do k=lk1,lk2
          do m=lm1,lm2
             do l=ll1,ll2
                do n=ln1,ln2
                   do sprime=0,n_spec_coll-1
                      do i=pi1,pi2
                         !sten=5
                         !vp flux
                         tmp5 = cvf_vp(l)*(cr1(i,k,l,m,sprime,n)*0.5 - cd11(i,k,l,m,sprime,n)*dvin &
                              + cd12(i,k,l,m,sprime,n)*0.25*dmuin(m-1) - cd12(i,k,l,m,sprime,n)*0.25*dmuin(m) &
                              - cd11(i,k,l-1,m,sprime,n)*dvin + cd12(i,k,l-1,m,sprime,n)*0.25*dmuin(m) &
                              - cd12(i,k,l-1,m,sprime,n)*0.25*dmuin(m-1) - cr1(i,k,l-1,m,sprime,n)*0.5)
                         !mu flux
                         tmp5 = tmp5 + cvf_mu(m)*(- cd22(i,k,l,m,sprime,n)*dmuin(m) &
                              + cr2(i,k,l,m,sprime,n)*0.5 &
                              - cd22(i,k,l,m-1,sprime,n)*dmuin(m-1) - cr2(i,k,l,m-1,sprime,n)*0.5)
                         if(my_pev.eq.0.and.l==ll1) then
                            !sten=4
                            tmp4 =-cvf_vp(l)*(cr1(i,k,l-1,m,sprime,n)*0.5 - cd11(i,k,l-1,m,sprime,n)*dvin &
                                    +cd12(i,k,l-1,m,sprime,n)*0.25*dmuin(m-1) - cd12(i,k,l-1,m,sprime,n)*0.25*dmuin(m))
                            tmp4 = tmp4+ cvf_mu(m)*(cd21(i,k,l,m-1,sprime,n)*0.25*dvinv(l) &
                                 - cd21(i,k,l,m,sprime,n)*0.25*dvinv(l))
                            tmp5 = tmp5+tmp4
                         end if
                         if(my_pew.eq.(n_procs_w-1).and.(m==lm2)) then
                            !sten=5
                            tmp5 = 0
                         end if
                         if(my_pev.eq.(n_procs_v-1).and.(l==ll2)) then
                            !sten=6
                            tmp6 =cvf_vp(l)*(cd11(i,k,l,m,sprime,n)*dvin - cd12(i,k,l,m,sprime,n)*0.25*dmuin(m) &
                                 + cd12(i,k,l,m,sprime,n)*0.25*dmuin(m-1) + cr1(i,k,l,m,sprime,n)*0.5)
                            tmp6 = tmp6 + cvf_mu(m)*(cd21(i,k,l,m,sprime,n)*0.25*dvinv(l) &
                                 - cd21(i,k,l,m-1,sprime,n)*0.25*dvinv(l))
                            tmp5 = tmp5+tmp6
                         end if
                         if (abs(tmp5)>ev_coll_est) ev_coll_est=abs(tmp5)
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo

       !empiric scaling factor by comparing to a set of exact SLEPc coll. eigenvalues between 100 and 7500
       ev_coll_est = ev_coll_est*2.8
       Call my_real_max_to_all(ev_coll_est)

    endif

    if (spacediff) call initialize_spacepart

    PERFOFF

  end subroutine initialize_collisions

  !>Deletes all arrays related to the collision operator
  subroutine finalize_collisions

    if (collision_op .eq. 'sugama') then

       deallocate(T_coeff1,T_coeff2,T_coeff3,T_coeff4,T_coeff5,T_coeff6)
       deallocate(particle_coeff,particle_fac)
       deallocate(particle_diffpart)
       deallocate(moments_colla,moments_collc)
       deallocate(particle_buff)
       if(allocated(C_ab_TO_mom)) deallocate(C_ab_TO_mom)
       if(allocated(C_ab_TO_en)) deallocate(C_ab_TO_en)

       if (xy_local) then
          deallocate(moments_collb)
       endif
    endif

    if (coll_f_fm_on) deallocate(loc_f_fm)
    deallocate(zeff_spec)

    if (collision_op .ne. 'sugama') then
       deallocate(cd11,cd12,cd21,cd22)
       deallocate(cr1,cr2)
       deallocate(cvf_mu,dmuin,cvf_vp)
       deallocate(dvinv)
    endif
    if (spacediff) call finalize_spacepart


  end subroutine finalize_collisions


!-----------------
!SPATIAL DIFFUSION

  !> Initializes the arrays used for the spatial diffusion part of the collision operator
  subroutine initialize_spacepart
    real:: nupref, vabs, velquo, fm_
    integer:: i,j,k,l,m,n,sp
    real:: Tprime, mi, ni, qi2
    logical:: active_spec

    allocate(kperp2(li1:li2,lj1:lj2,lk1:lk2),kpf(li1:li2,lj1:lj2,lk1:lk2),&
         sp_pref(pi1:pi2,pj1:pj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2))

    do k=lk1,lk2
       do j=lj1,lj2
          do i=li1,li2
             kperp2(i,j,k)=geom%gii(pi1,pj1,k)*ki(i)**2+2*geom%gij(pi1,pj1,k)*ki(i)*kj(j)+geom%gjj(pi1,pj1,k)*kj(j)**2
          enddo
       enddo
    enddo
    sp_pref=0.
    do n=ln1,ln2
       do m=lm1,lm2
          do l=ll1,ll2
             do k=lk1,lk2
                do sp=0,n_spec_coll-1
                   if (sp.lt.n_spec) active_spec=.not.spec(sp)%passive
                   if (sp.ge.n_spec) active_spec=.true.
                   if (active_spec) then
                      do j=pj1,pj2
                         do i=pi1,pi2
                            if (sp.lt.n_spec) then
                               !regular active species sprime
                               nupref=-coll*zeff_spec(i,sp)*spec(sp)%charge**2*spec(sp)%dens*spec(sp)%dens_prof(i)&
                                    &*spec(sp)%temp*spec(sp)%temp_prof(i)*spec(n)%mass**1.5 &
                                    /(spec(sp)%mass*geom%Bfield(i,j,k)**2*spec(n)%temp**1.5)
                               velquo=sqrt((spec(n)%temp*spec(sp)%mass)&
                                    /(spec(sp)%temp*spec(sp)%temp_prof(i)*spec(n)%mass))
                            else
                               !extra ion species sprime for adiabatic ion
                               !simulation
                               !parameters relative to active (electron) species
                               !Ti = Zeff*Te/tau
                               !mi = adiabatic_mass*(m_p/m_e)*\hat{m}_e
                               !(adiabatic_mass is m_i/m_p)
                               !ni=ne              !this is anyway assumed in
                               !the adiabatic response
                               !qi2=qe**2*Zeff
                               mi=adiabatic_mass*m_proton/m_electron*spec(n)%mass
                               qi2=Zeff0*spec(n)%charge**2
                               Tprime=spec(n)%temp*spec(n)%temp_prof(i)*Zeff0/tau
                               ni=spec(n)%dens*spec(n)%dens_prof(i)
                               nupref=-coll*qi2*ni*Tprime*spec(n)%mass**1.5 &
                                    /(mi*geom%Bfield(i,j,k)**2*spec(n)%temp**1.5)
                               velquo=sqrt((spec(n)%temp*mi)/(Tprime*spec(n)%mass))
                            endif
                            !vabs=sqrt(vp(l)**2+geom%Bfield(i,j,k)*mu(m))
                            !sp_pref(i,j,k,l,m,n)=sp_pref(i,j,k,l,m,n)+nupref*(2*vabs**2*cfk1(velquo*vabs) &
                            !     +3*geom%Bfield(i,j,k)*mu(m)*cfk2(velquo*vabs))/vabs**5
                            if (coll_f_fm_on)   then
                               vabs=sqrt(vp(l)**2+geom%Bfield(i,j,k)*mu(m))
                               fm_ = spec(n)%dens_prof(i) * (pi*spec(n)%temp_prof(i))**(-1.5) * exp(-vabs**2)
                               !sp_pref(i,j,k,l,m,n) =  sp_pref(i,j,k,l,m,n) * fm_
                               sp_pref(i,j,k,l,m,n)=sp_pref(i,j,k,l,m,n)+fm_*nupref*(2*vabs**2*cfk1(velquo*vabs) &
                                    +3*geom%Bfield(i,j,k)*mu(m)*cfk2(velquo*vabs))/vabs**5
                            else
                               vabs=sqrt(vp(l)**2+geom%Bfield(i,j,k)*mu(m))
                               sp_pref(i,j,k,l,m,n)=sp_pref(i,j,k,l,m,n)+nupref*(2*vabs**2*cfk1(velquo*vabs) &
                                    +3*geom%Bfield(i,j,k)*mu(m)*cfk2(velquo*vabs))/vabs**5

                            end if
                         enddo
                      enddo
                   endif
                enddo
             enddo
          enddo
       enddo
    enddo

  end subroutine initialize_spacepart

  !>Add the spatial diffusion part of the collision operator
  subroutine add_spacepart(loc_f,loc_k)
    implicit none
    !>f_ type array with exchanged ghost cells in vpar and mu directions
    complex,dimension(li1:li2, lj1:lj2, lbz:ubz, lbv:ubv, lbw:ubw, ln1:ln2),intent(in):: loc_f
    !>g_1 type array to which the result is added
    complex,dimension(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2),intent(inout):: loc_k
    !complex,dimension(li1:li2, lj1:lj2, lk1:lk2, ln1:ln2, 0:n_spec-1),intent(inout):: mom_fac, en_fac
    !mom_fac and en_fac by spacepart appear to be small.

    integer:: k,l,m,n

    PERFON('coll_sp')
    !todo think about implementation for global code
    if (xy_local) then
       do n=ln1,ln2
          do m=lm1,lm2
             do l=ll1,ll2
                kpf=kperp2*loc_f(:,:,lk1:lk2,l,m,n)
                do k=lk1,lk2
                   loc_k(:,:,k,l,m,n)=loc_k(:,:,k,l,m,n)+kpf(:,:,k)*sp_pref(pi1,pj1,k,l,m,n)
                enddo
             enddo
          enddo
       enddo
    else
       !kperp**2 F_1
       !to be implemented
    endif
    PERFOFF

  end subroutine add_spacepart

  !> Deletes arrays used by the spatial diffusion part of the collision operator
  subroutine finalize_spacepart
    deallocate(kperp2,kpf,sp_pref)
  end subroutine finalize_spacepart

!-----------------
!special functions

  !>F_1 function (p. 23 of flm's thesis)
  function cfk1(vvar) result (erg)
    real:: vvar,erg
    real:: erf
    select case(collision_op)
    case('pitch-angle','pitch-ion','pitch-all')
       !!pitch-angle operator transformed to v_par, mu co-ordinates (andbu)
!       if (vvar.lt.1E-4) then !!small velocity approximation
!          erg = 8.0*vvar**3/(3.0*sqrt(pi))
!       elseif (vvar.gt.10) then
       if (vvar.gt.10) then !!avoid SIGFPE due to large arguments
          !!expand v_a>>v_tb   (electrons scattering off heavy ions, vvar>>1)
          erg=(2.0*vvar**2.0-1.0)
       else
          erg=2.0*vvar*exp(-vvar**2.0)/sqrt(pi)+(2.0*vvar**2.0-1.0)*erf(vvar)
       endif
    case default
       erg=2.0*vvar*exp(-vvar**2.0)/sqrt(pi)+(2.0*vvar**2.0-1.0)*erf(vvar)
    endselect
  end function cfk1

  !>F_2 function (p. 23 of flm's thesis)
  function cfk2(vvar) result (erg)
    real:: vvar,erg
    real:: erf
    select case(collision_op)
    case('pitch-angle','pitch-ion','pitch-all')
       !!pitch-angle operator transformed to v_par, mu co-ordinates (andbu)
       erg = -cfk1(vvar)/3.
       !!expand v_a>>v_tb   (electrons scattering off heavy ions, vvar>>1)
!      erg=-2./3.*vvar**2
    case default
       erg=(1-2/3.*vvar**2)*erf(vvar)-2*vvar*exp(-vvar**2)/sqrt(pi)
    endselect
  end function cfk2

  !>F_3 function (p. 23 of flm's thesis)
  function cfk3(vvar) result (erg)
    real:: vvar,erg
    real:: erf
    select case(collision_op)
    case('pitch-angle','pitch-ion','pitch-all')
       erg=0.
    case default
       erg=2*erf(vvar)-4*vvar*exp(-vvar**2)/sqrt(pi)
    endselect
  end function cfk3

  !>Function used in (old?) GS2 to modify the pitch-angle scattering operator to include
  !!electron-electron scattering, only used for testing purposes
  function hee(vvar) result (erg)
    real:: vvar,erg
    real:: erf

    erg=exp(-vvar**2)/sqrt(pi)/vvar+(1.-1./(2*vvar**2))*erf(vvar)

  end function hee

  !>H function (eq. 28 of POP 17,122301 (Vernay))
  function Hfunc(vvar) result (erg)
    real:: vvar,erg
    real:: erf
    erg=1./vvar**3*(erf(vvar/sqrt(2.0))-sqrt(2./pi)*vvar*exp(-vvar**2/2.))
  end function Hfunc

  !>K function (eq. 29 of POP 17,122301 (Vernay))
  function Kfunc(vvar) result (erg)
    real:: vvar,erg
    real:: erf
    erg=1./vvar**3*((vvar**2-1.)*erf(vvar/sqrt(2.0)) + sqrt(2./pi)*vvar*exp(-vvar**2/2.) )
  end function Kfunc

  !>R function (to be referenced in documentation later)
  Function Rfunc(alpha,x) result (erg)
  Implicit None
  real :: alpha, x, erg
  real :: erf
  erg=(erf(x)-x*(1+alpha**2)*(2.0/sqrt(pi))*exp(-x**2))/(alpha*x)
  End Function Rfunc
!
end Module collisions
