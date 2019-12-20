!---------------------- PETSC SOLVER ---------------------------------------------
#ifndef NO_PETSC
!
module xgc_solver_module
#include <petscversion.h>
#if PETSC_VERSION_LT(3,6,0) || defined(PETSC_OLD_360)
#include <finclude/petscdef.h>
#include <finclude/petscsysdef.h>
#else
#include <petsc/finclude/petscdef.h>
#include <petsc/finclude/petscsysdef.h>
#endif
  PetscLogStage::smstage
  PetscLogEvent::FormFunctionEvent
  type xgc_solver
     Mat::Amat
     Mat::rhs_mat
     Mat::rhs2_mat
     Vec::XVec
     Vec::BVec
     DM::da
     DM::daphi
     DM::dalam
     Mat::KKTmat,schur,AinvB
     Mat::Bmat
     Mat::CMat
     Vec::CConstBCVec
     Mat::Dmat
     Mat::A0Mat
     Vec::rVec2,xVec2,bVec2
     SNES::snes
     PC::pc
     !real (8),allocatable::n0(:),e_Te(:) ! FormFunction
     Vec::TeInv ! FormFunction
     Mat::BIdentMat
     !Mat::BMassMat
     Mat::FSAMass ! FormFunction
     real (8)::scale
     KSP::ksp
     VecScatter::from_petsc,to_petsc
     PetscInt,allocatable::xgc_petsc(:),petscloc_xgc(:)
     integer::n_rhs_mat
     integer::comm,mype,totalpe,prefix
  end type xgc_solver
  !
contains
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine delete_solver( this, ierr )
    use sml_module, only : sml_mype
    implicit none
    type(xgc_solver) :: this
    PetscErrorCode :: ierr
    if (this%KKTmat .ne. 0) then ! two field solver
      call MatDestroy(this%KKTmat,ierr)
      call MatDestroy(this%schur,ierr)
      call MatDestroy(this%AinvB,ierr)
      call MatDestroy(this%Dmat,ierr)
      call MatDestroy(this%Bmat,ierr)
      call MatDestroy(this%Cmat,ierr)
      call VecDestroy(this%CConstBCVec,ierr)
      call MatDestroy(this%A0Mat,ierr)
      call MatDestroy(this%BIdentMat,ierr)
      call MatDestroy(this%FSAMass,ierr)
      call VecDestroy(this%TeInv,ierr)
      call VecDestroy(this%rVec2,ierr)
      call VecDestroy(this%bVec2,ierr)
      call VecDestroy(this%xVec2,ierr)
      call SNESDestroy(this%snes,ierr)
      call DMDestroy(this%da,ierr)
      call DMDestroy(this%daphi,ierr)
      call DMDestroy(this%dalam,ierr)
      call PCDestroy(this%pc,ierr)
   end if
   if (allocated(this%xgc_petsc)) then
      deallocate(this%xgc_petsc)
      deallocate(this%petscloc_xgc)
   endif

   if (this%bVec .ne. 0 ) then
      call VecDestroy(this%xVec,ierr)
      call VecDestroy(this%bVec,ierr)
      call VecScatterDestroy(this%from_petsc,ierr)
      call VecScatterDestroy(this%to_petsc,ierr)
   endif

   if (this%ksp .ne. 0) then
      call KSPDestroy(this%ksp,ierr) 
   endif
   if (this%Amat .ne. 0) then
      call MatDestroy( this%Amat,ierr)
   endif
   
   if (this%rhs_mat .ne. 0) then
      call MatDestroy( this%rhs_mat,ierr)
   endif
   
   if (this%rhs2_mat .ne. 0) then
      call MatDestroy( this%rhs2_mat,ierr)
   endif
   
 end subroutine delete_solver
 
end module xgc_solver_module
#endif

#ifdef XGC1_EM
!
! time stepper class
!
module xgc_ts_module
#include <petscversion.h>
#if PETSC_VERSION_LT(3,6,0)
#include <finclude/petscsysdef.h>
#else
#include <petsc/finclude/petscsysdef.h>
#endif
  !use petscsysdef
  !use petsctsdef
  use boundary_class
  use grid_class

  PetscLogEvent::ADVMatrixInitEvent,ADVInitEvent,ADVSolveEvent,ADVFormIFuncEvent,ADVFormIJacEvent,ADVFormRHSEvent
  type xgc_ts
     TS::ts
     IS::iss(0:2)
     Mat::FJacobian2,FJacobian,muDel2 ! LHS
     Mat::mass ! mass matrix, somewhat depreciated for 3-field solver
     Vec::mass_inv_vec ! 3-field solver requires diagonal mass, lets keep this shadow
     VecScatter,dimension(0:2)::from_petsc,to_petsc,from_petsc_single_value
     PetscInt,ALLOCATABLE::xgc_petsc(:),petscloc_xgc(:),petsc_xgc(:)
     PetscInt::nnode ! size of XGC vectors in plane
     real (kind=8),pointer::dni_v(:),ui_v(:)
     PetscScalar::scales(0:2)

     type(grid_type), pointer :: grid
     type(boundary2_type), pointer :: bd
#ifdef OLD_INIT_FF
     integer, pointer :: tr(:,:)
     !rh note that with OLD_INIT_FF, dx=1/dl_parallel
     real (kind=8), pointer :: p(:,:,:), dx(:)
#else
     integer, pointer :: tr(:,:)
     real (kind=8), pointer :: p(:,:,:), dx(:,:)
#endif

  end type xgc_ts
  !
contains
  !
  subroutine ts_create(this)
    !use sml_module, only : sml_mype
    implicit none
    type(xgc_ts),intent(out)::this
    !
    this%ts = 0 ! PETSC_NULL_OBJECT
    this%mass = 0
    this%mass_inv_vec = 0
    this%FJacobian = 0
    this%FJacobian2 = 0    
    this%muDel2 = 0    
    this%scales(0) = 1d0
    this%scales(1) = 1d0
    this%scales(2) = 1d0
  end subroutine ts_create
  !
  subroutine ts_delete( this, ierr )
    use sml_module, only : sml_mype
    implicit none
#if PETSC_VERSION_LT(3,6,0)
#include <finclude/petsc.h>
#else
#include <petsc/finclude/petsc.h>
#endif
    type(xgc_ts)::this
    PetscErrorCode::ierr
    if (this%ts/=0) then
       if (this%mass/=0) call MatDestroy( this%mass, ierr );CHKERRQ(ierr)
       if (this%mass_inv_vec/=0) call VecDestroy( this%mass_inv_vec, ierr );CHKERRQ(ierr)
       call ISDestroy( this%iss(0), ierr );CHKERRQ(ierr)
       call ISDestroy( this%iss(1), ierr );CHKERRQ(ierr)
       call ISDestroy( this%iss(2), ierr );CHKERRQ(ierr)
       if (this%FJacobian/=0) call MatDestroy( this%FJacobian, ierr );CHKERRQ(ierr)
       if (this%FJacobian2/=0) call MatDestroy( this%FJacobian2, ierr );CHKERRQ(ierr)
       if (this%muDel2/=0) call MatDestroy( this%muDel2, ierr );CHKERRQ(ierr)
       call TSDestroy(this%ts,ierr);CHKERRQ(ierr)
       call VecScatterDestroy(this%to_petsc(0),ierr)
       call VecScatterDestroy(this%to_petsc(1),ierr)
       call VecScatterDestroy(this%to_petsc(2),ierr)
       call VecScatterDestroy(this%from_petsc(0),ierr)
       call VecScatterDestroy(this%from_petsc(1),ierr)
       call VecScatterDestroy(this%from_petsc(2),ierr)
       call VecScatterDestroy(this%from_petsc_single_value,ierr)       
       deallocate(this%xgc_petsc)
       deallocate(this%petsc_xgc)
    endif
  end subroutine ts_delete
end module xgc_ts_module
#endif


!
module psn_class
  use xgc_solver_module
#ifdef XGC1_EM
   use xgc_ts_module
#endif
  use boundary_class
  use mat_class
  type psn_type
     type(xgc_solver) :: solver00,solverH
#ifdef XGC1_EM
     type(xgc_ts)::ts
     KSP::kspmass ! solver for J in explicit solver
#endif
     real (8), allocatable :: iden_ff(:,:),idensity(:,:),idensity0(:),factor1(:),edensity0(:),edensity(:,:)  ! remove iden_ff
     real (8), allocatable :: iden_rho_ff(:,:,:)
     real (8), allocatable :: rhs1(:), rhs2(:)
     

     real (8), allocatable :: tempi_ev(:), B(:)
     real (8), allocatable :: tite(:),dpot(:,:),dden(:,:),rhoi2(:),tempe_ev(:), dpot_ff(:,:)   !! ion density
     real (8), allocatable :: pot0(:), add_pot0(:), pot0m(:)
     real (8), allocatable :: pot_rho_ff(:,:,:), E_rho_ff(:,:,:,:)

     real (8), allocatable :: E_para(:,:)
     real (8), allocatable :: E_perp_node(:,:,:), E_perp0_node(:,:)
     real (8), allocatable :: E_perp_tr(:,:,:), E_perp0_tr(:,:)

#ifdef OLD_INIT_FF
     real (8), allocatable :: bfollow_p(:,:,:)
     integer, allocatable :: bfollow_tr(:,:)
     real (8), allocatable :: bfollow_1_dx(:)
#endif

     real (8), allocatable :: ff_hdp_p(:,:,:) ! field following _ half dphi
     integer, allocatable :: ff_hdp_tr(:,:)
#ifndef OLD_INIT_FF
     real (8), allocatable :: ff_hdp_dx(:,:)
#endif
     real (8), allocatable :: ff_1dp_p(:,:,:) ! field following _ one dphi
     integer, allocatable :: ff_1dp_tr(:,:)
#ifndef OLD_INIT_FF
     real (8), allocatable :: ff_1dp_dx(:,:)
#endif

     real (8), allocatable :: iden00_1d(:), vol00(:), pot00_1d(:), eden00_1d(:), cden00_1d(:),n2ne00_1d(:)

     !f0_grid related
#ifndef F0_TOR_LINEAR
     real (8), allocatable :: iden_rho_f0(:,:)
     real (8), allocatable :: eden_f0(:)
#else
     real (8), allocatable :: iden_rho_f0(:,:,:)
     real (8), allocatable :: eden_f0(:,:)
#endif
     real (8), allocatable :: eden00_f0(:)
     real (8), allocatable :: idensity_f0(:,:), edensity_f0(:,:) !- real space density


     integer :: nwall
     real (8), allocatable :: sheath_pot(:)
     integer, allocatable :: wall_nodes(:)
     real (8), allocatable :: sheath_lost(:,:) , sheath_ilost(:,:)  ! ilost is ion only
     integer, allocatable :: node_to_wall(:)

     type(boundary2_type) :: cbd0_2, pbd0_2, cbdh_2, pbdh_2 ! charge / potential boundary of High-k solver / 0 solver
!#ifdef XGC_COUPLING_CORE_EDGE
!Seung-Hoe quick patch for BCs:
     type(boundary2_type) :: cbd0_tmp ! PATCH
!#endif

#ifdef XGC1_EM
     type(boundary2_type) :: febdh_2
#endif
     type (mat_type) :: gyro_avg_mat
     ! electron subcycling
     real (8), allocatable :: ddpotdt(:,:), dpotsave(:,:,:) 
     real (8), allocatable :: save_dpot(:,:,:), save_dpot0(:,:,:) 
     real (8), allocatable :: E00_ff(:,:,:)
     ! all plane information
#ifdef USE_CALC_GRADIENT
     real (8), allocatable :: pot_phi_real(:,:)
#endif
     real (8), allocatable :: pot_phi_ff(:,:,:), E_phi_ff(:,:,:,:), ddpotdt_phi(:,:,:)

     ! for simple00
     real (8), allocatable :: mn_eb2(:)


#ifdef XGC1_EM
     ! for hybrid scheme
     real (8), allocatable :: eden_hyb(:)
     real (8), allocatable :: eden_hyb0(:)
     real (8), allocatable :: A_par(:)
     real (8), allocatable :: A_par0(:)
     real (8), allocatable :: u_e(:)
     real (8), allocatable :: E_para_tot(:)
     real (8), allocatable :: ijpar_ff(:,:)
     real (8), allocatable :: Je(:)
#endif
  end type psn_type

contains
  subroutine psn_mem_alloc(psn,n,ntr,npsi,nrho,nhybrid)
    use sml_module
    implicit none
    type(psn_type) :: psn
    integer :: n, ntr, npsi, nrho, nhybrid ! ierr !node, nphi, ntriangle
    integer :: nphim1
    allocate( psn%iden_rho_ff(n,0:1,0:nrho), psn%iden_ff(n,0:1))
    allocate( psn%idensity(n,0:1), psn%dpot(n,-1:2), psn%dpot_ff(n,0:1))
    allocate( psn%edensity(n,0:1))  ! only when electron is on
    allocate( psn%rhs1(n), psn%rhs2(n) )

    allocate( psn%E_perp_tr(2,ntr,0:1),psn%E_perp0_tr(2,ntr))
    allocate( psn%E_perp_node(2,n,0:1),psn%E_perp0_node(2,n),psn%E_para(n,0:1) )

    allocate( psn%pot_rho_ff(0:1,0:nrho,n), psn%E_rho_ff(3,0:1,0:nrho,n))

    allocate( psn%idensity0(n), psn%edensity0(n) )
    allocate( psn%tempi_ev(n), psn%tite(n), psn%rhoi2(n), psn%tempe_ev(n), &
         psn%dden(n,1:1) , psn%pot0(n) )
#ifdef OLD_INIT_FF
    allocate( psn%bfollow_p(3,2,n), psn%bfollow_tr(2,n), psn%bfollow_1_dx(n) )
#endif
    allocate( psn%ff_hdp_p(3,n,0:1), psn%ff_hdp_tr(n,0:1) )
#ifndef OLD_INIT_FF
     allocate(psn%ff_hdp_dx(n,0:1))
#endif

    allocate( psn%ff_1dp_p(3,n,0:1), psn%ff_1dp_tr(n,0:1) )
#ifndef OLD_INIT_FF
     allocate(psn%ff_1dp_dx(n,0:1))
#endif

    ! 00
    allocate( psn%iden00_1d(npsi),psn%vol00(npsi), psn%pot00_1d(npsi), &
         &psn%eden00_1d(npsi), psn%cden00_1d(npsi),psn%n2ne00_1d(npsi) )

    if(sml_electron_on .or. sml_electron_hyb) then
       ! following two lines of memory allocation are for old delta-f scheme
       ! Those are not needed for code benchmark and linear calculations
       allocate(psn%ddpotdt(n,0:1),psn%dpotsave(n,0:1,4*nhybrid))  ! should be 4 * sml_nhybrid
       allocate(psn%E00_ff(2,0:1,n))
#ifdef PURE_RK4
       allocate(psn%save_dpot(n,0:1,4),psn%save_dpot0(n,0:1,4))
#else
       allocate(psn%save_dpot(n,0:1,2),psn%save_dpot0(n,0:1,2))
#endif

       !for all plane information
       nphim1 = sml_nphi_total-1
       allocate(psn%pot_phi_ff(0:1,n,0:nphim1),psn%ddpotdt_phi(n,0:1,0:nphim1))
#ifdef USE_CACL_GRADIENT
#else
       allocate( psn%E_phi_ff(3,0:1,n,0:nphim1) )
#endif

    endif
#ifdef USE_CALC_GRADIENT
       nphim1 = sml_nphi_total-1
       allocate(psn%pot_phi_real(n,0:nphim1))
#endif


    if(sml_f0_grid) then
       allocate(psn%B(n))
#ifndef F0_TOR_LINEAR
       allocate(psn%iden_rho_f0(n,0:nrho))
#else
       allocate(psn%iden_rho_f0(n,0:1,0:nrho))
#endif
       if(sml_electron_on) then
#ifndef F0_TOR_LINEAR
          allocate(psn%eden_f0(n))
#else
          allocate(psn%eden_f0(n,0:1))
#endif
          allocate(psn%eden00_f0(npsi))
       endif
#ifdef F0_CHARGE_N0
       allocate( psn%idensity_f0(n,0:1))
       if(sml_electron_on) then
          allocate( psn%edensity_f0(n,0:1))
       endif
#endif

    endif


    !if(sml_iter_solver .or. sml_poisson_solver_type/=0) 
    allocate(psn%pot0m(n))


    psn%idensity=0D0
    psn%idensity0=0D0
    psn%tempi_ev=0D0
    psn%tite=0D0
    psn%tempe_ev=0D0
    psn%dpot=0D0
    psn%edensity0=0D0
    psn%pot0=0D0
    psn%E_para=0D0
    psn%E_perp_tr=0D0
    psn%E_perp0_tr=0D0
    psn%E_perp_node=0D0
    psn%E_perp0_node=0D0
#ifdef OLD_INIT_FF
    psn%bfollow_p=0D0
    psn%bfollow_tr=0
    psn%bfollow_1_dx=0D0
#endif
    psn%rhoi2=0D0

    ! required
    if (allocated(psn%dpotsave)) then
       psn%dpotsave=0D0
    endif

    if(allocated(psn%save_dpot)) then
       psn%save_dpot=0D0
       psn%save_dpot0=0D0
    endif
    psn%eden00_1d=0D0
    psn%n2ne00_1d=0D0

#ifdef XGC1_EM
    if(sml_electron_hyb) then
      allocate( psn%eden_hyb(n))  ! only when hybrid model is on
      allocate( psn%eden_hyb0(n))  ! only when hybrid model is on
      allocate( psn%A_par(n))  ! only when hybrid model is on
      allocate( psn%A_par0(n))  ! only when hybrid model is on
      allocate( psn%u_e(n))  ! only when hybrid model is on
      allocate( psn%E_para_tot(n))  ! only when hybrid model is on
      allocate( psn%ijpar_ff(n,0:1))  ! only when hybrid model is on
      allocate( psn%Je(n) )

      psn%eden_hyb=0D0
      psn%eden_hyb0=0D0
      psn%A_par=0D0
      psn%A_par0=0D0
      psn%u_e=0D0
      psn%E_para_tot=0D0
      psn%ijpar_ff=0D0
      psn%Je=0D0
   endif
#endif

  end subroutine psn_mem_alloc

  subroutine psn_sheath_alloc(psn)
    use sml_module, only : sml_nthreads
    implicit none
    type(psn_type) :: psn
    allocate(psn%sheath_pot(psn%nwall))
    allocate(psn%wall_nodes(psn%nwall))
    allocate(psn%sheath_lost(psn%nwall,sml_nthreads))
    allocate(psn%sheath_ilost(psn%nwall, sml_nthreads))
    psn%sheath_lost=0D0
    psn%sheath_ilost=0D0
    ! allocation of node_to_wall is in sheath_init since grid%nnode is not visible
  end subroutine psn_sheath_alloc
end module psn_class
!
