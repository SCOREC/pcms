#include <petscversion.h>
#if PETSC_VERSION_LT(3,6,0)
#include <finclude/petscsnesdef.h>
#else
#include <petsc/finclude/petscsnesdef.h>
#endif
subroutine push_fluid(grid,psn,ipc)
  use petsc
  use sml_module
  use grid_class
  use psn_class
  use perf_monitor
  use omp_module, only: split_indices
  use eq_module, only: eq_den_v1
  implicit none
  type(grid_type),target :: grid
  type(psn_type),target :: psn
  integer, intent(in) :: ipc
  !
  real (8) :: dt_now, dt, norm
  real (kind=8),allocatable::dAdt(:),ddenedt(:)
  real (kind=8),allocatable::apar(:),J1(:) ! just for bechmark with implicit petsc solver
  integer::ith,i,i_beg(sml_nthreads),i_end(sml_nthreads)
  PetscErrorCode::ierr
  real (8),target::dum(grid%nnode)  !dummy variable for interface. It will be replaced with ion flow later.
  real (kind=8),allocatable::dentmp(:)

  !integer,parameter::debug=1
  integer :: debug

  Mat::mat
  Vec::vec,globalvec,globalvec2
  VecScatter::scat
  IS::is,is2
  PetscInt::low,npetscplane,npetsctot
  PetscInt,parameter::ione=1,targ_p=0
  PetscViewer::viewer
  MPI_Comm::w_comm

  debug=sml_hyb_debug

  if(sml_mype==targ_p) print*,sml_mype,') push_fluid: start'

  if(sml_gstep<2 .and. ipc==1) then
    call initial_condition(psn%A_par,psn,grid)
  endif
  !if(sml_gstep==2)psn%eden_hyb=0D0
  dt=sml_dt

  !save phase0 information
  select case(ipc)
  case(1)
     if (sml_nrk .gt. 1 .or. debug==2) then
       dt_now=0.5D0*dt
     else
       dt_now=dt
     endif
     !!!$OMP PARALLEL DO &
     !!!$OMP PRIVATE( ITH, I )
     !do ith=1,sml_nthreads
     !do i=i_beg(ith), i_end(ith)
     !      phase0(:,i)=ptl(i)%ph
     !enddo
     !enddo
     psn%A_par0=psn%A_par
     psn%eden_hyb0=psn%eden_hyb
  case(2)
     dt_now=dt
  end select

  if (sml_use_ts_solver) then
     ! debug
     if (debug==1) then
        if (sum(abs(psn%A_par))==targ_p) stop 'push_fluid:test vector == 0'
        if(sml_mype==targ_p) call write_vec(grid%nnode,psn%A_par,'petsc_lhs.m','petsclhs',ierr)        
        call VecCreateSeqWithArray(PETSC_COMM_SELF,ione,grid%nnode,psn%A_par,vec,ierr);CHKERRQ(ierr)
        call MatGetSubMatrix(psn%ts%FJacobian,psn%ts%iss(0),psn%ts%iss(3),MAT_INITIAL_MATRIX,mat,ierr);CHKERRQ(ierr)
        if (grid%nnode .lt. 200) then
           call  PetscObjectSetName(mat,'B0gradMat',ierr);CHKERRQ(ierr)
           call  PetscObjectGetComm(mat,w_comm,ierr);CHKERRQ(ierr)
           call  PetscViewerASCIIOpen(w_comm, 'B0gradMat.m', viewer,ierr);CHKERRQ(ierr)
           call  PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB,ierr);CHKERRQ(ierr)
           call  MatView(mat,viewer,ierr);CHKERRQ(ierr)
           call  PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
        end if
#if PETSC_VERSION_LT(3,6,0)
        call MatGetVecs(mat,globalvec2,globalvec,ierr);CHKERRQ(ierr)
#else
        call MatCreateVecs(mat,globalvec2,globalvec,ierr);CHKERRQ(ierr)
#endif
        ! get vector
        call MatGetSize(mat,npetsctot,PETSC_NULL_INTEGER,ierr);CHKERRQ(ierr)
        npetscplane = npetsctot/sml_nphi_total
        low = sml_plane_index*npetscplane ! global start my0 for plane
        call ISCreateGeneral(PETSC_COMM_SELF,npetscplane,psn%ts%petsc_xgc,PETSC_COPY_VALUES,is2,ierr);CHKERRQ(ierr)
        call ISCreateStride(PETSC_COMM_SELF,npetscplane,low,ione,is,ierr);CHKERRQ(ierr)
        call VecScatterCreate(globalvec,is,vec,is2,scat,ierr);CHKERRQ(ierr) 
        call ISDestroy(is,ierr);CHKERRQ(ierr)
        call ISDestroy(is2,ierr);CHKERRQ(ierr)

        ! put lhs vector in PETSc globalvec
        call VecScatterBegin(scat, vec, globalvec,INSERT_VALUES,SCATTER_REVERSE,ierr);CHKERRQ(ierr)
        call VecScatterEnd(  scat, vec, globalvec,INSERT_VALUES,SCATTER_REVERSE,ierr);CHKERRQ(ierr)
        call MatMult(mat, globalvec, globalvec2, ierr);CHKERRQ(ierr) ! globalvec2 has B0grad|| * lhs
        call VecScatterBegin(scat,globalvec2,vec,INSERT_VALUES,SCATTER_FORWARD,ierr)
        call VecScatterEnd(  scat,globalvec2,vec,INSERT_VALUES,SCATTER_FORWARD,ierr) ! 
        if(sml_mype==targ_p) call write_vec(grid%nnode,psn%A_par,'petsc_b0grad.m','petscb0grad',ierr)
        
        call VecScatterDestroy(scat,ierr)
        call VecDestroy(vec,ierr);CHKERRQ(ierr)
        call VecDestroy(globalvec,ierr);CHKERRQ(ierr)
        call VecDestroy(globalvec2,ierr);CHKERRQ(ierr)
        call MatDestroy(mat,ierr);CHKERRQ(ierr)
        call mpi_barrier(sml_comm,ierr)
        stop 'PETSc solver debugging'
     end if

     ! u_i will be added later and replace dum
     dum=0D0
     psn%ts%dni_v => psn%idensity(:,1) ! cache this, 
     psn%ts%ui_v => dum

     call ts_solve(psn%ts, dt_now, psn%eden_hyb, psn%A_par, psn%dpot(:,1), psn%Je, ierr)

     !rh Make Je electron current
     if (sml_hyb_ion_on) then
        psn%Je = psn%Je - psn%ijpar_ff(:,1)
     endif
     !rh Assuming constant density for now
     psn%u_e = psn%Je/sml_e_charge/eq_den_v1
     
  else
     ! debug
     if (debug==1) then
        if (sum(abs(psn%A_par))==targ_p) stop 'push_fluid:test vector == 0'
        if(sml_mype==targ_p) call write_vec(grid%nnode,psn%A_par,'xgc_lhs.m','xgclhs',ierr)
        allocate(apar(grid%nnode))
        allocate(J1(grid%nnode))
        call B0Grad_B0( apar, J1,  grid, psn, psn%A_par )
        if(sml_mype==targ_p) call write_vec(grid%nnode,apar,'xgc_b0grad.m','xgcb0grad',ierr)
        deallocate(apar, J1)
        call mpi_barrier(sml_comm,ierr)
        stop 'debugging'
     end if

     !call solve_poisson(grid,psn,1)    
     psn%pot0m = 0D0

     call split_indices(grid%nnode, sml_nthreads, i_beg, i_end)

     allocate(dentmp(grid%nnode))
     
     ! For all poloidal plane - without n=0 mode
     !$OMP PARALLEL DO &
     !$OMP PRIVATE( ITH, I )
     do ith=1,sml_nthreads     
        psn%edensity(:,1) = psn%eden_hyb
        do i=i_beg(ith),i_end(ith)
           dentmp(i)= (psn%idensity(i,1) - psn%edensity(i,1) )   !????? check later convert_001d_2_grid 
           !/psn%density0_full(:)  ! ni - ne - <ni - ne>
        enddo
     enddo

     ! smoothing of charge density -- RHS
     !call smooth_tr_connect(grid,dentmp3)
     !call smooth_pol(grid,dentmp3)
     call set_boundary2_values(dentmp,0D0,psn%pbd0_2)
     
     psn%dden(:,1)=dentmp  ! this can be optimized. redundant memory usage.
     call t_startf("POISSON_TURB_SOLVE")
     psn%dpot(:,1) = 0d0
     call petsc_solve(grid%nnode,dentmp,0,psn%dpot(:,1),psn%solver00%comm,psn%solver00,ierr)
     call t_stopf("POISSON_TURB_SOLVE")
     deallocate(dentmp)

     ! compute derivative for forward Euler
     allocate(ddenedt(grid%nnode), dAdt(grid%nnode))
     !call t_startf("PUSH_FLUID_DERIVATIVE") 
     call fluid_derivative(ddenedt, dAdt, psn, grid)
     !call t_stopf("PUSH_FLUID_DERIVATIVE")
     
     psn%A_par=psn%A_par0 + dAdt*dt_now
     psn%eden_hyb= psn%eden_hyb0 + ddenedt*dt_now

     deallocate(ddenedt, dAdt)
  end if ! solver switch

  ! debug
  if (debug==2) then
     if(sml_mype==targ_p) then
        call VecCreateSeqWithArray(PETSC_COMM_SELF,ione,grid%nnode,psn%eden_hyb,vec,ierr);CHKERRQ(ierr)
        if (sml_use_ts_solver) then
           call PetscObjectSetName(vec,'petscn1',ierr);CHKERRQ(ierr)
           call PetscViewerASCIIOpen(PETSC_COMM_SELF, 'petsc_n1.m', viewer,ierr);CHKERRQ(ierr)
        else
           call PetscObjectSetName(vec,'xgcn1',ierr);CHKERRQ(ierr)
           call PetscViewerASCIIOpen(PETSC_COMM_SELF, 'xgc_n1.m', viewer,ierr);CHKERRQ(ierr)
        end if
        call PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB,ierr);CHKERRQ(ierr)
        call VecView(vec,viewer,ierr);CHKERRQ(ierr)
        call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
        call VecDestroy(vec,ierr);CHKERRQ(ierr)
        
        call VecCreateSeqWithArray(PETSC_COMM_SELF,ione,grid%nnode,psn%A_par,vec,ierr);CHKERRQ(ierr)
        if (sml_use_ts_solver) then
           call PetscObjectSetName(vec,'petscapar',ierr);CHKERRQ(ierr)
           call PetscViewerASCIIOpen(PETSC_COMM_SELF, 'petsc_apar.m', viewer,ierr);CHKERRQ(ierr)
        else
           call PetscObjectSetName(vec,'xgcapar',ierr);CHKERRQ(ierr)
           call PetscViewerASCIIOpen(PETSC_COMM_SELF, 'xgc_apar.m', viewer,ierr);CHKERRQ(ierr)
        end if
        call PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB,ierr);CHKERRQ(ierr)
        call VecView(vec,viewer,ierr);CHKERRQ(ierr)
        call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
        call VecDestroy(vec,ierr);CHKERRQ(ierr)
        
        call VecCreateSeqWithArray(PETSC_COMM_SELF,ione,grid%nnode,psn%dpot(:,1),vec,ierr);CHKERRQ(ierr)
        if (sml_use_ts_solver) then
           call PetscObjectSetName(vec,'petscpot',ierr);CHKERRQ(ierr)
           call PetscViewerASCIIOpen(PETSC_COMM_SELF, 'petsc_pot.m', viewer,ierr);CHKERRQ(ierr)
        else
           call PetscObjectSetName(vec,'xgcpot',ierr);CHKERRQ(ierr)
           call PetscViewerASCIIOpen(PETSC_COMM_SELF, 'xgc_pot.m', viewer,ierr);CHKERRQ(ierr)
        end if
        call PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB,ierr);CHKERRQ(ierr)
        call VecView(vec,viewer,ierr);CHKERRQ(ierr)
        call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
        call VecDestroy(vec,ierr);CHKERRQ(ierr)
        
        call VecCreateSeqWithArray(PETSC_COMM_SELF,ione,grid%nnode,psn%Je,vec,ierr);CHKERRQ(ierr)
        if (sml_use_ts_solver) then
           call PetscObjectSetName(vec,'petscje',ierr);CHKERRQ(ierr)
           call PetscViewerASCIIOpen(PETSC_COMM_SELF, 'petsc_je.m', viewer,ierr);CHKERRQ(ierr)
        else
           call PetscObjectSetName(vec,'xgcje',ierr);CHKERRQ(ierr)
           call PetscViewerASCIIOpen(PETSC_COMM_SELF, 'xgc_je.m', viewer,ierr);CHKERRQ(ierr)
        end if
        call PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB,ierr);CHKERRQ(ierr)
        call VecView(vec,viewer,ierr);CHKERRQ(ierr)
        call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
        call VecDestroy(vec,ierr);CHKERRQ(ierr)
     end if
     call mpi_barrier(sml_comm,ierr)
     stop 'debugging'
  end if
     
  ! post process 
  call t_startf("PUSH_FLUID_TIME_ADVANCE")
  call set_boundary2_values(psn%a_par,0D0,psn%pbd0_2)
  if (sml_mode_select_on==1.and.(mod(sml_gstep,10)==0)) then
     call mode_selection_comb(sml_mode_select_n,grid,psn,psn%a_par)
  endif
  
  call set_boundary2_values(psn%eden_hyb,0D0,psn%febdh_2)
  if (sml_mode_select_on==1.and.(mod(sml_gstep,10)==0))then
     call mode_selection_comb(sml_mode_select_n,grid,psn,psn%eden_hyb)
  endif
  call t_stopf("PUSH_FLUID_TIME_ADVANCE")

end subroutine push_fluid

subroutine B0Grad_B0( apar, J1,  grid, psn, Apar_in ) ! apar has B_0 grad_par(n1)/B_0
  use grid_class
  use psn_class
  use eq_module
  use sml_module
  use perf_monitor
  implicit none
#include <petscversion.h>
#if PETSC_VERSION_LT(3,6,0)
#include <finclude/petsc.h90>
#else
#include <petsc/finclude/petsc.h90>
#endif
  type(grid_type),intent(in)::grid             !!! type(psn_type)::a_psn
  type(psn_type),intent(in)::psn
  real (kind=8),dimension( grid%nnode)::J1,apar, Apar_in !!! J1 is the perturbed current

  real (kind=8), allocatable :: u_par_over_b(:)
  real (kind=8), allocatable :: dparX(:)
  real (kind=8), allocatable :: lap_A_par(:)
  real (kind=8), allocatable :: bpar_grad_u(:)

  logical, external :: is_nan
  integer :: i, sgn
  Vec::vec
  PetscErrorCode::ierr
  PetscInt::ione

  allocate(u_par_over_b(grid%nnode))
  allocate(dparX(grid%nnode))
  allocate(bpar_grad_u(grid%nnode), lap_A_par(grid%nnode))

  !rh We do not need the Laplacian here, just calculate
  !rh B0.grad(Apar_in/(e B))
  !rh
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!! Get nabla^2(A_parallel)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call fem_del2( lap_A_par, Apar_in, grid, psn%pbd0_2 )
  !rh if (sml_mass_solve) then
     ione = 1
     call VecCreateSeqWithArray(PETSC_COMM_SELF,ione,grid%nnode,lap_A_par,vec,ierr);CHKERRQ(ierr)
     call VecScatterBegin(psn%solver00%to_petsc, vec, psn%solver00%bVec, INSERT_VALUES,SCATTER_FORWARD,ierr)
     CHKERRQ(ierr)
     call VecScatterEnd(  psn%solver00%to_petsc, vec, psn%solver00%bVec, INSERT_VALUES,SCATTER_FORWARD,ierr)
     CHKERRQ(ierr)

     call KSPSolve(psn%kspmass, psn%solver00%bVec, psn%solver00%xVec, ierr);CHKERRQ(ierr)

     call VecScatterBegin(psn%solver00%from_petsc, psn%solver00%xVec, vec, INSERT_VALUES,SCATTER_FORWARD,ierr)
     call VecScatterEnd(  psn%solver00%from_petsc, psn%solver00%xVec, vec, INSERT_VALUES,SCATTER_FORWARD,ierr)
     call VecDestroy(vec,ierr)
     lap_A_par = -lap_A_par
  !rh else
  !rh    lap_A_par= -lap_A_par/grid%node_area ! old solver
  !rh end if

  call set_boundary2_values(lap_A_par,0D0,psn%pbd0_2)

  if(sml_mode_select_on==1.and.(mod(sml_gstep,10)==0)) then
    call mode_selection_comb(sml_mode_select_n,grid,psn,lap_A_par)
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!! Evaluate B.grad(u_par,e/B) and db_perp.grad(u_par,e/B)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do i=1, grid%nnode
  !   J1(i) =  lap_A_par(i)/sml_mu0     ! perturbed current from laplace A_par
  !   if(is_nan(J1(i))) print *, "n1 nan", lap_A_par(i), Apar_in(i)
  !   u_par_over_b(i) = J1(i)/sml_e_charge/grid%bfield(4,i)
     u_par_over_b(i) = Apar_in(i)/sml_e_charge/grid%bfield(4,i)
  enddo

  ! Get the B grad_|| ne ue_||/B for each node
  call t_startf("GET_LEFT_RIGHT_PLANE")
  dparX(:)=0D0
#ifdef OLD_INIT_FF
  call GradParX(grid,psn%bfollow_tr,psn%bfollow_p,psn%bfollow_1_dx,psn%pbd0_2,u_par_over_b,dparX)
#else
  call GradParX(grid,psn%ff_1dp_tr,psn%ff_1dp_p,psn%ff_1dp_dx,psn%pbd0_2,u_par_over_b,dparX)
#endif
  call t_stopf("GET_LEFT_RIGHT_PLANE")

  do i=1, grid%nnode
    apar(i) =dparX(i)*grid%bfield(4,i)
  enddo

end subroutine B0Grad_B0

subroutine Ampere_eq(grid,      & ! (in)
     psn,                       & ! (inout)
     Pe_para,                   & ! (in)
     deltape,                   & ! (in)
     dAdt,                      & ! (out)
     v_e,                       & ! (out)
     delta_B_perp               & ! (out)
     )

  use grid_class
  use psn_class
  use sml_module
  use eq_module
  use perf_monitor
  implicit none
  !
  type(grid_type), intent(in) :: grid
  type(psn_type), intent(inout) :: psn
  real (kind=8), dimension(grid%nnode), intent(in) :: Pe_para, deltape
  real (kind=8), intent(out), dimension(grid%nnode,2) :: delta_B_perp, v_e
  real (kind=8), intent(out), dimension(grid%nnode)  :: dAdt
  !
  integer :: i, j, sgn
  real (kind=8) :: r, z, psi, psitmp, eta
  real (kind=8), allocatable :: grad_para_phi(:), dparX(:)
  real (kind=8), allocatable :: E_r_node(:), E_z_node(:)
  real (kind=8), allocatable :: B_perp_grad_pe(:), dum_v(:,:)
  real (kind=8) :: e_terms(5)
  real (kind=8) , external :: psi_interpol, I_interpol
  logical, external :: is_nan
  real (kind=8) :: hyb_ion

  allocate(grad_para_phi(grid%nnode), E_r_node(grid%nnode), E_z_node(grid%nnode))
  allocate(B_perp_grad_pe(grid%nnode))
  allocate(dparX(grid%nnode))

  if (sml_hyb_ion_on) then
     hyb_ion = 1D0
  else
     hyb_ion = 0D0
  endif

  !rh Coefficients for switching on/off different terms in E_par equation
  if (sml_hyb_alfven_test) then
    !rh keep only the eta*j term for consistency with implicit solver
    e_terms = (/0d0, 0D0, 0D0, 0D0, 1D0/)
  elseif (sml_hyb_tearing_test) then
    !rh e_terms(2) and e_terms(4) are not implemented yet!
    !rh Neglect grad_par(delta_pe) for now --->
    e_terms = (/0d0, 0D0, 0D0, 0D0, 1D0/)
  elseif (sml_hyb_linear) then
    e_terms = (/1d0, 0D0, 1D0, 0D0, 1D0/)
  else
    e_terms(:) = 1D0
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!! Get grad(dpot)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  grad_para_phi(:)=0D0

  call get_dpot_grad(grid,psn,E_r_node,E_z_node,grad_para_phi)
  !rh I choose to define this as b.grad(phi)
  grad_para_phi=-grad_para_phi
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!! Evaluate b.grad(dP_e)/(e n_0)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !rh This is now a subroutine (to be used in the implicit code)
  dparX(:)=0D0

  call t_startf("GET_LEFT_RIGHT_PLANE")
#ifdef OLD_INIT_FF
  call grad_par_delta_pe(grid,psn%bfollow_tr,psn%bfollow_p,psn%bfollow_1_dx,psn%pbd0_2,deltape,dparX)
#else
  call grad_par_delta_pe(grid,psn%ff_1dp_tr,psn%ff_1dp_p,psn%ff_1dp_dx,psn%pbd0_2,deltape,dparX)
#endif
  call t_stopf("GET_LEFT_RIGHT_PLANE")

  psn%E_para_tot = e_terms(1) * dparX

  !rh do i=1, grid%nnode
  !rh   if(is_nan(psn%E_para_tot(i))) print *, "psn%E_para_tot Nan", deltape_r(i),deltape_l(i),grid%dene(i)
  !rh enddo

  !if(sml_gstep>=2)psn%E_para_tot = 0D0


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!! Evaluate the perturbation of the magnetic field (dB_perp)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call get_dbperp(grid,psn%A_par,delta_B_perp)


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!! Evaluate dB_perp.grad(P_0e)/B
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   if(sml_deltaf_f0_mode==-2) then
     ! This is for constant temperature
     ! b <-- T_e*dB_perp/|B|
     ! u <-- n_0
     allocate(dum_v(grid%nnode,2))
     do i=1,grid%nnode
       dum_v(i,1)=grid%tempe(i)*delta_B_perp(i,1)/grid%bfield(4,i)
       dum_v(i,2)=grid%tempe(i)*delta_B_perp(i,2)/grid%bfield(4,i)
     enddo
     if (sml_old_grad_perp) then
       call b_dot_grad_u(B_perp_grad_pe,dum_v(:,1),dum_v(:,2),grid%dene,grid)
     else
       call v_dot_grad_perp(grid,dum_v,grid%dene,B_perp_grad_pe)
     endif
     deallocate(dum_v)
   elseif(sml_deltaf_f0_mode==-1) then
     do i =1, grid%nnode
       psitmp=1D0/grid%absgradpsi(i)  ! drdpsi
       B_perp_grad_pe(i)=-psitmp*Pe_para(i)*sml_f0_1_Ln*(delta_B_perp(i,2)*grid%gradpsi(i,1)+  &
                          delta_B_perp(i,2)*grid%gradpsi(i,2))/grid%bfield(4,i)
     enddo
   endif

   psn%E_para_tot = psn%E_para_tot - e_terms(3) * B_perp_grad_pe/sml_e_charge/grid%dene


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!! Evaluate final parallel electric field
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   do i=1, grid%nnode
      !if(eta>sml_eta)then
      !  psn%E_para_tot(i) =-eta*psn%u_e(i)*sml_e_charge*dene(i) !e n0 E||=eta*delta j||
      !else
        !rh psn%E_para_tot(i) =psn%E_para_tot(i)-sml_eta*psn%u_e(i)*sml_e_charge*grid%dene(i)
        psn%E_para_tot(i) =psn%E_para_tot(i) + e_terms(5) * sml_eta &
                                     * (hyb_ion * psn%ijpar_ff(i,1) - sml_e_charge * grid%dene(i) * psn%u_e(i))
      !endif
   enddo
   if(sml_mype==0.and.mod(sml_gstep,1000)==0) print *, 'delta j, j0=', sum((psn%E_para_tot/sml_eta)**2), sum((grid%j0_par)**2)


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!! Evaluate dA/dt=-b.grad(phi)-E_parallel
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Get A_par d A_par/dt=-grad_para phi-E||, use explicit method
   dAdt(:) = -grad_para_phi(:) - psn%E_para_tot(:)
   if(is_nan(dAdt)) print *, "dAdt nan"


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!! Evaluate ExB velocity
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   do i=1, grid%nnode
      v_e(i,1)=-E_z_node(i)*grid%bfield(3,i)/grid%bfield(4,i)**2  ! V_E_r
      v_e(i,2)=E_r_node(i)*grid%bfield(3,i)/grid%bfield(4,i)**2 ! V_E_z
   enddo


   deallocate(grad_para_phi)
   deallocate(E_r_node, E_z_node)
   deallocate(B_perp_grad_pe)
   deallocate(dparX)

   return

end subroutine Ampere_eq



subroutine fluid_derivative( ddenedt, dAdt, & ! (out)
     psn,           & ! (in and out) psn%eden_hyb
     grid           & ! (in) grid
     )
   
  use grid_class    
  use psn_class
  use ptl_module, only: ptl_charge
  use eq_module
  use sml_module
  use perf_monitor
  implicit none
#include <petscversion.h>
#if PETSC_VERSION_LT(3,6,0)
#include <finclude/petsc.h90>
#else
#include <petsc/finclude/petsc.h90>
#endif
  !
  type(grid_type), intent(in):: grid
  type(psn_type), intent(inout):: psn
  real (kind=8), intent(out), dimension(grid%nnode) :: ddenedt, dAdt
  !
  integer :: i, sgn
  real (kind=8) :: r, z, psi, tmp, psitmp
  real (kind=8), allocatable :: dpsi_dr(:), dpsi_dz(:)
  real (kind=8), allocatable :: u_par_over_b(:)
  real (kind=8), allocatable :: dparX(:)
  real (kind=8), allocatable :: cross_B_pe(:,:), cross_B_dpot(:,:)
  real (kind=8), allocatable :: u_e(:), u_i(:), lap_A_par(:), delta_pe(:)
  real (kind=8), allocatable :: v_e(:,:), delta_B_perp(:,:)
  real (kind=8), allocatable :: bpar_grad_u(:), bperp_grad_u(:), ve_dot_gradne(:), crossB_grad_pe(:)
  real (kind=8), allocatable :: crossB_grad_dpot(:), B_perp_grad_pe(:), Pe_para(:), ve_dot_graddne(:)
  real (kind=8), allocatable :: u0_par(:), dB_grad_u0(:)
  integer, parameter :: m1=3415, m2=3415
  real (kind=8) :: envelop, sml_f0_psi_c, sml_f0_1_psi_w
  real (kind=8) :: n_terms(7)
  real (kind=8) , external :: psi_interpol
  real (kind=8) :: hyb_ion
  logical, external :: is_nan
  Vec::vec
  PetscErrorCode::ierr
  PetscInt::ione
  
  allocate(u_par_over_b(grid%nnode))
  allocate(dparX(grid%nnode))
  allocate(cross_B_pe(grid%nnode,2), cross_B_dpot(grid%nnode,2))
  allocate(v_e(grid%nnode,2), delta_B_perp(grid%nnode,2))
  allocate(u_e(grid%nnode), u_i(grid%nnode), lap_A_par(grid%nnode), delta_pe(grid%nnode))
  allocate(bpar_grad_u(grid%nnode),bperp_grad_u(grid%nnode), ve_dot_gradne(grid%nnode), &
           crossB_grad_pe(grid%nnode), crossB_grad_dpot(grid%nnode), &
           ve_dot_graddne(grid%nnode))
  allocate(B_perp_grad_pe(grid%nnode), Pe_para(grid%nnode),dB_grad_u0(grid%nnode), u0_par(grid%nnode))
  allocate(dpsi_dr(grid%nnode),dpsi_dz(grid%nnode))

  if (sml_hyb_ion_on) then 
     hyb_ion = 1D0
  else
     hyb_ion = 0D0
  endif

  if (sml_hyb_alfven_test) then
    n_terms = (/1d0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0/)
  elseif (sml_hyb_tearing_test) then
    n_terms = (/1d0, 0D0, 1D0, 0D0, 0D0, 0D0, 0D0/)
  elseif (sml_hyb_linear) then
    n_terms = (/1d0, 0D0, 1D0, 1D0, 0D0, 1D0, 1D0/)
  else
    n_terms(:) = 1D0
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!! Set up basic values
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  sml_f0_psi_c=0.5*(sqrt(sml_inpsi)+sqrt(sml_outpsi))
  sml_f0_1_psi_w=1D0/( 0.4*(sqrt(sml_outpsi)-sqrt(sml_inpsi)) )

  do i=1, grid%nnode
     Pe_para(i)=grid%dene(i)*grid%tempe(i)
     delta_pe(i) = psn%eden_hyb(i)*grid%tempe(i)
  enddo
  
  do i=1, grid%nnode
    if(sml_mype==0.and.sml_gstep<10) write(1212,*) v_e(i,1), v_e(i,2)
  enddo


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!! Sanity check for tearing mode: print j_toroidal from Ampere's law
  !!!!!!!!!!!! (rh: Why is this called j0_par?)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! To test tearing mode
  if(sml_mype==0.and.sml_gstep==1) then
    open(unit=4321,file='j0.dat',position='append')
    do i=1, grid%nnOde
       write(4321,100) sqrt((grid%x(1,i)-eq_axis_r)**2+grid%x(2,i)**2), grid%j0_par(i)
    enddo
    close(4321)
  endif


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!! Get nabla^2(A_parallel)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call fem_del2( lap_A_par, psn%A_par, grid, psn%pbd0_2 )
  !rh Deactivated mass solve due to instability in tearing mode test
  !rh if (sml_mass_solve) then
     ione = 1
     call VecCreateSeqWithArray(PETSC_COMM_SELF,ione,grid%nnode,lap_A_par,vec,ierr);CHKERRQ(ierr)
     call VecScatterBegin(psn%solver00%to_petsc, vec, psn%solver00%bVec, INSERT_VALUES,SCATTER_FORWARD,ierr)
     CHKERRQ(ierr)
     call VecScatterEnd(  psn%solver00%to_petsc, vec, psn%solver00%bVec, INSERT_VALUES,SCATTER_FORWARD,ierr) 
     CHKERRQ(ierr)

     call KSPSolve(psn%kspmass, psn%solver00%bVec, psn%solver00%xVec, ierr);CHKERRQ(ierr)

     call VecScatterBegin(psn%solver00%from_petsc, psn%solver00%xVec, vec, INSERT_VALUES,SCATTER_FORWARD,ierr)
     call VecScatterEnd(  psn%solver00%from_petsc, psn%solver00%xVec, vec, INSERT_VALUES,SCATTER_FORWARD,ierr)  
     call VecDestroy(vec,ierr)
     !rh This is confusing, but fem_del2 calculates -grad_perp^2!!!
     lap_A_par = -lap_A_par
  !rh else
  !rh    lap_A_par= -lap_A_par/grid%node_area ! old solver
  !rh end if
 
  call set_boundary2_values(lap_A_par,0D0,psn%pbd0_2)

  if(sml_mode_select_on==1.and.(mod(sml_gstep,10)==0)) then
    call mode_selection_comb(sml_mode_select_n,grid,psn,lap_A_par)
  endif


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!! Evaluate parallel ion and electron flows
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !rh call ion_v_para(u_i,psn,grid)
  ! Get Ue_para
  do i=1, grid%nnode
     !rh ijpar is really the ion current
     psn%u_e(i) = hyb_ion * psn%ijpar_ff(i,1)/(grid%dene(i)*ptl_charge(1)) + lap_A_par(i)/sml_mu0/sml_e_charge/grid%dene(i)
     !rh J = -e n_0 u_e = - lap_A_par/sml_mu0 - e n_0 u_i
     psn%Je(i) = -lap_A_par(i)/sml_mu0 - hyb_ion * psn%ijpar_ff(i,1)
  enddo
  if(sml_mode_select_on==1.and.(mod(sml_gstep,10)==0)) then
    call mode_selection_comb(sml_mode_select_n,grid,psn,psn%u_e)
  endif


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!! Evaluate dA/dt, E_parallel, and db_perp
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! This needs updated psn%u_e
  call t_startf("FLUID_DERIV_AMPERE_EQ")
  call Ampere_eq(grid, psn, pe_para, delta_pe, dAdt, v_e, delta_B_perp)
  call t_stopf("FLUID_DERIV_AMPERE_EQ")


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!! Evaluate B.grad(u_par,e/B) and db_perp.grad(u_par,e/B)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do i=1, grid%nnode
     !rh u_par_over_b(i) = (u_i(i) * 0D0 + &   ! ion parallel velocity
     !rh     lap_A_par(i)/(sml_mu0*sml_e_charge)) / grid%bfield(4,i)      ! from laplace A_par
     !rh This avoids confusion -->
     u_par_over_b(i) = psn%u_e(i) * grid%dene(i) / grid%bfield(4,i)
     if(is_nan(u_par_over_b(i))) print *, "u_par_over_b nan", lap_A_par(i), psn%A_par(i)
  enddo

  ! Get the B grad_|| ne ue_||/B for each node
  call t_startf("GET_LEFT_RIGHT_PLANE")
  dparX(:)=0D0
#ifdef OLD_INIT_FF
  call GradParX(grid,psn%bfollow_tr,psn%bfollow_p,psn%bfollow_1_dx,psn%pbd0_2,u_par_over_b,dparX)
#else
  call GradParX(grid,psn%ff_1dp_tr,psn%ff_1dp_p,psn%ff_1dp_dx,psn%pbd0_2,u_par_over_b,dparX)
#endif
  call t_stopf("GET_LEFT_RIGHT_PLANE")

  do i=1, grid%nnode
    bpar_grad_u(i) =dparX(i)*grid%bfield(4,i)
  enddo

  ! Get the delta_B_perp dot grad ne ue_||/B
  if (sml_old_grad_perp) then
    call b_dot_grad_u(bperp_grad_u, delta_B_perp(:,1), delta_B_perp(:,2), u_par_over_b, grid)
  else
    call v_dot_grad_perp(grid,delta_B_perp,u_par_over_b,bperp_grad_u)
  endif


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!! Evaluate v_E.grad(n_0+dn_e)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Get v_e dot grad ne
  if(sml_deltaf_f0_mode==-2) then  !consistent grad ne with real profile
    if (sml_old_grad_perp) then
      call b_dot_grad_u(ve_dot_gradne, v_e(:,1), v_e(:,2), (grid%dene), grid)
      call b_dot_grad_u(ve_dot_graddne, v_e(:,1), v_e(:,2), (psn%eden_hyb), grid)
    else
      call v_dot_grad_perp(grid,v_e,(grid%dene),ve_dot_gradne)
      call v_dot_grad_perp(grid,v_e,(psn%eden_hyb),ve_dot_graddne)
    endif
  elseif(sml_deltaf_f0_mode==-1) then
  !  do i =1, grid%nnode
  !    psitmp=1D0/grid%absgradpsi(i)  ! drdpsi 
  !    dpsi_dr(i)= -sml_f0_1_Ln*grid%gradpsi(i,1)*psitmp*grid%dene(i)
  !    dpsi_dz(i)= -sml_f0_1_Ln*grid%gradpsi(i,2)*psitmp*grid%dene(i)
  !  enddo
  !  call b_dot_grad_u(ve_dot_gradne, dpsi_dr, dpsi_dz, -psn%dpot(:,1), grid)

    if (sml_old_grad_perp) then
      call b_dot_grad_u(ve_dot_graddne, v_e(:,1), v_e(:,2), psn%eden_hyb, grid)
    else
      call v_dot_grad_perp(grid,v_e,psn%eden_hyb,ve_dot_graddne)
    endif

    do i =1, grid%nnode
      psitmp = 1D0/grid%absgradpsi(i)  ! drdpsi
      envelop = exp( - ((sqrt(grid%psi(i)) - sml_f0_psi_c)*sml_f0_1_psi_w )**8 ) 
      ve_dot_gradne(i) = -envelop*sml_f0_1_Ln*psitmp*grid%dene(i)   &
                         * (v_e(i,1)*grid%gradpsi(i,1)+v_e(i,2)*grid%gradpsi(i,2))
    enddo
    !rh Allow for separation of linear and non-linear terms in the end
    !rh ve_dot_gradne = ve_dot_gradne + ve_dot_graddne
  endif


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!! Evaluate (B x grad(B)).grad(delta_pe)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Get B cross grad B dot grad delta pe /me /Omega_e/B^2
   do i=1, grid%nnode
      cross_B_pe(i,1) = grid%B_cross_gradB(1,i)/(sml_e_charge*grid%bfield(4,i)**3)
      cross_B_pe(i,2) = grid%B_cross_gradB(2,i)/(sml_e_charge*grid%bfield(4,i)**3)
   enddo

   if (sml_old_grad_perp) then
     call b_dot_grad_u(crossB_grad_pe, cross_B_pe(:,1), cross_B_pe(:,2), delta_pe, grid)
   else
     call v_dot_grad_perp(grid,cross_B_pe,delta_pe,crossB_grad_pe)
   endif

   ! add the artifical temperature gradient term in deltaf-mode
   !if(sml_deltaf_f0_mode==-1) then
   !  do i=1, grid%nnode
   !    psitmp=1D0/grid%absgradpsi(i)  ! drdpsi
   !    dpsi_dr(i)= -sml_f0_1_Lt_e*grid%gradpsi(i,1)*psitmp*grid%tempe(i)
   !    dpsi_dz(i)= -sml_f0_1_Lt_e*grid%gradpsi(i,2)*psitmp*grid%tempe(i)
   !    crossB_grad_pe(i) = crossB_grad_pe(i) +(cross_B_pe(i,1)*dpsi_dr(i)   &
   !                       +cross_B_pe(i,2)*dpsi_dz(i))*psn%eden_hyb(i)
   !  enddo
   !endif


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!! Evaluate (B x grad(B)).grad(phi)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   do i=1, grid%nnode
      !rh Does eden_hyb belong here???
      cross_B_dpot(i,1) = grid%B_cross_gradB(1,i)*2D0*(grid%dene(i)+psn%eden_hyb(i))/grid%bfield(4,i)**3
      cross_B_dpot(i,2) = grid%B_cross_gradB(2,i)*2D0*(grid%dene(i)+psn%eden_hyb(i))/grid%bfield(4,i)**3
   enddo
   !rh psn%dpot(:,1) is contiguous in memory, we can pass it as argument directly
   if (sml_old_grad_perp) then
     call b_dot_grad_u(crossB_grad_dpot, cross_B_dpot(:,1), cross_B_dpot(:,2), psn%dpot(:,1), grid)
   else
     call v_dot_grad_perp(grid,cross_B_dpot,psn%dpot(:,1),crossB_grad_dpot)
   endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!! Evaluate kink term dB_perp.grad(j_0/(e*B))
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !rh This is now a subroutine (to be used also in the implicit code)
  call kink_drive2(grid,psn%A_par,dB_grad_u0)
  if (sml_mype==0) then
    do i=1,grid%nnode
      write(9998,*),dB_grad_u0(i)
    enddo
  endif

  call kink_drive(grid,delta_B_perp,dB_grad_u0)
  if (sml_mype==0) then
    do i=1,grid%nnode
      write(9999,*),dB_grad_u0(i)
    enddo
  endif


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!! Combine all terms to evaluate d(dn_e)/dt
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !rh original -->
   !rh ddenedt = -(bpar_grad_u + bperp_grad_u*0D0 + ve_dot_gradne*0D0 - 2*crossB_grad_pe*0D0 + crossB_grad_dpot*0D0 &
   !rh             + dB_grad_u0)
   !rh Use n_terms to switch on/off desired terms
   ddenedt = -( n_terms(1) * bpar_grad_u + n_terms(2) * bperp_grad_u + n_terms(3) * dB_grad_u0  &
               +n_terms(4) * ve_dot_gradne + n_terms(5) * ve_dot_graddne  &
               -n_terms(6) * 2*crossB_grad_pe + n_terms(7) * crossB_grad_dpot )
                               !- 2*crossB_grad_pe, &   !rh Factor 2 --> dP_perp=dP_para??? Why "-", should be "+"???
                               !crossB_grad_dpot/) ,)    !rh Why "+", should be "-"???


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!! Some diagnostic output
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   if(sml_mype==0 .and. mod(sml_gstep,10)==0 ) then
     open(unit=1223,file='osc.dat',position='append')
     write(1223,100) sml_gstep*sml_dt, sqrt(sum(psn%dpot(:,1)**2)), sqrt(psn%dpot(m1,1)**2+psn%dpot(m2,1)**2), &
          sqrt(sum(psn%a_par**2)), sqrt(sum(psn%u_e**2)), sqrt(sum(crossB_grad_dpot**2)), psn%dpot(m1,1)
     100   FORMAT(7(D12.4, 2X))
     close(1223)
   endif
   if(sml_mype==0.and. mod(sml_gstep,10)==0) print *, 'A,dadt,n,dpot,dndt=', sml_gstep, sqrt(sum(psn%a_par**2)), sqrt(sum(dadt**2)), &
        sqrt(sum(psn%eden_hyb**2)), sqrt(sum(psn%dpot(:,1)**2)), sqrt(sum(ddenedt**2)),sqrt(sum(bpar_grad_u**2)),  &
        sqrt(sum(ve_dot_gradne**2)), sqrt(sum(crossB_grad_pe**2)), sqrt(sum(crossB_grad_dpot**2)), &
        sqrt(sum(dB_grad_u0**2))
  !if(sml_mype==0) print *, 'nu/B, u ,lapA =', u_par_over_b(m1), u_e(m1), lap_a_par(m1)


  deallocate(dpsi_dr, dpsi_dz, ve_dot_gradne, ve_dot_graddne)
  deallocate(bpar_grad_u, crossB_grad_dpot, dparX)
  deallocate(u_par_over_b, cross_B_pe, cross_B_dpot)
  deallocate(v_e, delta_B_perp, u_e, u_i, lap_A_par, delta_pe)
  deallocate(bperp_grad_u, crossB_grad_pe, u0_par)
  deallocate(B_perp_grad_pe, Pe_para, dB_grad_u0)

  return

end subroutine fluid_derivative



subroutine initial_condition(A_par,psn,grid)
  use sml_module
  use grid_class
  use psn_class
  use eq_module
  implicit none

  integer, parameter :: m=2, n=1
  real (8), parameter :: n0=1D-14
  type(grid_type) :: grid
  type(psn_type) :: psn
  integer :: i, j, ij
  real (8) :: r, z, theta, phi, psi, center, width, r_minor
  real (kind=8), intent(out) :: A_par(grid%nnode)


  center=0D0
  width=0.5
  A_par=0D0

  if(.true.)then
    do i = 1, grid%nnode
       r = grid%x(1,i)
       z = grid%x(2,i)
       phi = grid%phimin
       psi = grid%psi(i)
       r_minor=dsqrt((r-eq_axis_r)**2+(z-eq_axis_z)**2)
       theta=acos((r-eq_axis_r)/(r_minor+1D-23))
       if(z<eq_axis_z) then
          theta=sml_2pi-theta
       endif
       A_par(i) = n0*sin(sml_pi*(r_minor-center)/width)*cos(-n*(phi-sml_pi)+m*theta)*exp(-((psi/eq_x_psi-0.411D0)/0.05D0)**2)
       !A_par(i) = n0*cos(n*phi+m*theta)
       !rh Perturbation localized in radius and phi to excite a localized continuum of Alfven waves
       !rh A_par(i) = n0*exp(-((psi/eq_x_psi-0.45D0)/0.2)**2)*exp(-((phi-sml_pi)/(sml_pi/4D0))**2)
    enddo
  endif

  psn%eden_hyb=0D0
  if(.false.) then
    do i = 1, grid%nnode
       r = grid%x(1,i)
       z = grid%x(2,i)
       phi = grid%phimax
       psi = grid%psi(i)
       theta=acos((r-eq_axis_r)/(dsqrt((r-eq_axis_r)**2+(z-eq_axis_z)**2)+1D-23))
       if(z<eq_axis_z) then
          theta=sml_2pi-theta
       endif
       psn%eden_hyb(i) = 1D10*sin(sml_pi*(r_minor-0.05)/(0.50))*cos(-n*phi+m*theta)
       !psn%eden_hyb(i) = 5d7*2*(psi/eq_x_psi-0.411)/0.03*exp(-((psi/eq_x_psi-0.411)/0.03)**2)*cos(-n*phi+m*(theta-sml_pi/4D0)) 
       !rh psn%eden_hyb(i) = 1D10*exp(-((psi/eq_x_psi-0.45D0)/0.1)**2)*exp(-((phi-sml_pi)/(sml_pi/4D0))**2)
    enddo
  endif

  return

end subroutine initial_condition



subroutine ion_v_para(j_i,psn,grid)
  use grid_class
  use psn_class
  implicit none
  type(psn_type) :: psn
  type(grid_type) :: grid
  real (kind=8), intent(out), dimension(grid%nnode)  :: j_i

  j_i =psn%ijpar_ff(:,1)

  return
end subroutine ion_v_para


subroutine GradParX(grid,tr,p,dx,bd,input,output)
  use grid_class
  use boundary_class
  use sml_module
  implicit none
  include 'mpif.h'
  
  type(grid_type), intent(in) :: grid
  type(boundary2_type), intent(in) :: bd
#ifdef OLD_INIT_FF
  integer, intent(in) :: tr(2,grid%nnode)
  !rh note that with OLD_INIT_FF, dx=1/dl_parallel
  real (kind=8), intent(in) :: p(3,2,grid%nnode), dx(grid%nnode)
#else
  integer, intent(in) :: tr(grid%nnode,0:1)
  real (kind=8), intent(in) :: p(3,grid%nnode,0:1), dx(grid%nnode,0:1)
#endif
  real (kind=8), dimension(grid%nnode), intent(in) :: input
  real (kind=8), dimension(grid%nnode), intent(out) :: output

  real (kind=8), dimension(:), allocatable :: sendl, recvr, dum_r, dum_l, X_r, X_l
  real (kind=8) :: p_r(3), p_l(3), sgn, l_l, l_r, l_tot, mid_weight
  integer :: icount, isource, idest, isendtag, irecvtag, ierr
  integer :: istatus(mpi_status_size), idx1, idx2
  integer :: i, nd, nodes_r(3), nodes_l(3)
  logical, external :: is_nan

  allocate(sendl(grid%nnode), recvr(grid%nnode))
  allocate(dum_l(grid%nnode), dum_r(grid%nnode), X_r(grid%nnode), X_l(grid%nnode))

  sendl=input
  recvr(:)=0D0
  icount=grid%nnode
  isource= modulo(sml_intpl_mype+1,sml_intpl_totalpe)
  idest  = modulo(sml_intpl_mype-1+sml_intpl_totalpe,sml_intpl_totalpe)
  isendtag=sml_intpl_mype
  irecvtag=isource

  call mpi_sendrecv(sendl,icount,MPI_REAL8,idest,isendtag,&
       recvr,icount,MPI_REAL8,isource,irecvtag,sml_intpl_comm,istatus,ierr)
  dum_r = recvr

  sendl=input
  recvr(:)=0D0
  icount=grid%nnode
  isource= modulo(sml_intpl_mype-1+sml_intpl_totalpe,sml_intpl_totalpe)
  idest  = modulo(sml_intpl_mype+1,sml_intpl_totalpe)
  isendtag=sml_intpl_mype
  irecvtag=isource

  call mpi_sendrecv(sendl,icount,MPI_REAL8,idest,isendtag,&
       recvr,icount,MPI_REAL8,isource,irecvtag,sml_intpl_comm,istatus,ierr)
  dum_l = recvr

  deallocate(sendl, recvr)

#ifdef OLD_INIT_FF
  if(sml_bt_sign>0D0) then
     idx1=1
     idx2=2
  else
     idx1=2
     idx2=1
  endif
#else
  idx1=1
  idx2=0
#endif

  do i=1, grid%nnode
#ifdef OLD_INIT_FF
     nodes_r=grid%nd(:,tr(idx1,i))
     nodes_l=grid%nd(:,tr(idx2,i))
     p_r(:)=p(:,idx1,i)
     p_l(:)=p(:,idx2,i)
#else
     nodes_r=grid%nd(:,tr(i,idx1))
     nodes_l=grid%nd(:,tr(i,idx2))
     p_r(:)=p(:,i,idx1)
     p_l(:)=p(:,i,idx2)
#endif
     X_r(i) = 0D0
     X_l(i) = 0D0
     do nd=1,3
       !rh added is_inside check to make grad_par comparable to MakeGradParMat
       if (is_inside(nodes_r(nd),bd)) then
         X_r(i) = X_r(i) + dum_r(nodes_r(nd))*p_r(nd)
       endif
       if (is_inside(nodes_l(nd),bd)) then
         X_l(i) = X_l(i) + dum_l(nodes_l(nd))*p_l(nd)
       endif
     enddo
  enddo

  deallocate(dum_l, dum_r)

  if (sml_bt_sign .lt. 0D0) then
    sgn=-1D0
  else
    sgn=1D0
  endif

#ifdef OLD_INIT_FF
  !rh original derivative assuming equal distance between planes -1 and 0
  !rh and planes 0 and +1
  output(:) = sgn * (X_r(:)-X_l(:)) * dx(:)
#else
  do i=1, grid%nnode
    l_l   = dx(i,idx2)
    l_r   = dx(i,idx1)
    l_tot = l_l+l_r
    if (is_inside(i,bd)) then
      mid_weight=1D0
    else
      mid_weight=0D0
    endif
    !rh More accurate parallel derivative that takes assumes arbitrary
    !rh distances between planes -1, 0 and +1
    !rh d_parallel(X) = w(-1)*X(-1) + w(0)*X(0) + w(1)*X(1)
    !rh w(-1) = - dl(1)/(dl(-1)*(dl(-1)+dl(1)))
    !rh w(-1) =  (dl(1)-dl(-1))/(dl(-1)*dl(1))
    !rh w(1)  = + dl(-1)/(dl(1)*(dl(-1)+dl(1)))
    output(i) = sgn * ( - l_r/(l_l*l_tot)*X_l(i) &
                        + (l_r-l_l)/(l_l*l_r)*mid_weight*input(i) &
                        + l_l/(l_r*l_tot)*X_r(i) )
  enddo
#endif

  do i=1, grid%nnode
    if(is_nan(output(i))) print *, "GradParX Nan", X_r(i),X_l(i),input(i)
  enddo

  deallocate(X_r, X_l)

  return

end subroutine GradParX

subroutine GradPhiX(grid,input,output)
  use grid_class
  use boundary_class
  use sml_module
  use eq_module, only:eq_axis_r
  implicit none
  include 'mpif.h'

  type(grid_type), intent(in) :: grid
  !type(boundary2_type), intent(in) :: bd
  real (kind=8), dimension(grid%nnode), intent(in) :: input
  real (kind=8), dimension(grid%nnode), intent(out) :: output

  real (kind=8), dimension(:), allocatable :: sendl, recvr, dum_r, dum_l
  integer :: icount, isource, idest, isendtag, irecvtag, ierr
  integer :: istatus(mpi_status_size)
  integer :: i
  logical, external :: is_nan

  allocate(sendl(grid%nnode), recvr(grid%nnode))
  allocate(dum_l(grid%nnode), dum_r(grid%nnode))

  sendl=input
  recvr(:)=0D0
  icount=grid%nnode
  isource= modulo(sml_intpl_mype+1,sml_intpl_totalpe)
  idest  = modulo(sml_intpl_mype-1+sml_intpl_totalpe,sml_intpl_totalpe)
  isendtag=sml_intpl_mype
  irecvtag=isource

  call mpi_sendrecv(sendl,icount,MPI_REAL8,idest,isendtag,&
       recvr,icount,MPI_REAL8,isource,irecvtag,sml_intpl_comm,istatus,ierr)
  dum_r = recvr

  sendl=input
  recvr(:)=0D0
  icount=grid%nnode
  isource= modulo(sml_intpl_mype-1+sml_intpl_totalpe,sml_intpl_totalpe)
  idest  = modulo(sml_intpl_mype+1,sml_intpl_totalpe)
  isendtag=sml_intpl_mype
  irecvtag=isource

  call mpi_sendrecv(sendl,icount,MPI_REAL8,idest,isendtag,&
       recvr,icount,MPI_REAL8,isource,irecvtag,sml_intpl_comm,istatus,ierr)
  dum_l = recvr

  deallocate(sendl, recvr)

  if (sml_cylindrical) then
    do i=1,grid%nnode
      output(i) = (dum_r(i)-dum_l(i))/(eq_axis_r*grid%delta_phi)
    enddo
  else
    do i=1,grid%nnode
      output(i) = (dum_r(i)-dum_l(i))/(grid%x(1,i)*grid%delta_phi)
    enddo
  endif

  deallocate(dum_r, dum_l)

end subroutine GradPhiX

subroutine get_dpot_grad(grid,psn,E_r,E_z,E_para)
  use grid_class
  use psn_class
  use omp_module
  use perf_monitor
  use sml_module, only: sml_bt_sign, sml_nthreads, sml_turb_efield, sml_old_grad_perp
  implicit none

  type(grid_type), intent(in) :: grid
  type(psn_type), intent(in) :: psn
  real (kind=8), dimension(grid%nnode), intent(out) :: E_r, E_z, E_para
  real (kind=8), dimension(:), allocatable :: dpot, dparX
  real (kind=8), dimension(:,:), allocatable :: E_perp_tr
  real (kind=8) :: E(2), dp1, dp2, area_sum
  integer, dimension(sml_nthreads) :: itr_beg, itr_end, i_beg, i_end
  integer :: i, j, itr, ith, nodes(3), sgn

#ifdef EM_NO_PHI00

  E_para(:) = psn%E_para(:,1)
  E_r(:)    = psn%E_perp_node(1,:,1)
  E_z(:)    = psn%E_perp_node(2,:,1)

  return

#else

  E_r(:)    = 0D0
  E_z(:)    = 0D0
  E_para(:) = 0D0

  allocate(dpot(grid%nnode), dparX(grid%nnode))

  dparX(:)=0D0

  if(sml_turb_efield) then
    dpot(:)=psn%dpot(:,1)
  else
    dpot(:)=0D0
  endif

  !rh We may still be able to avoid this MPI communication
  !rh in case the poisson solver has done this already
  !rh Then we can use psn%dpot and transform to ff coordinates
  !rh to get dpot_r and dpot_l
  call t_startf("GET_LEFT_RIGHT_PLANE")
#ifdef OLD_INIT_FF
  call GradParX(grid,psn%bfollow_tr,psn%bfollow_p,psn%bfollow_1_dx,psn%pbd0_2,dpot,dparX)
#else
  call GradParX(grid,psn%ff_1dp_tr,psn%ff_1dp_p,psn%ff_1dp_dx,psn%pbd0_2,dpot,dparX)
#endif

  call t_stopf("GET_LEFT_RIGHT_PLANE")

  E_para = -dparX

  deallocate(dparX)

  if (sml_old_grad_perp) then
    !rh Old FE gradient operator

    allocate(E_perp_tr(2,grid%ntriangle))
    E_perp_tr(:,:)=0D0

    ! Divide triangles among OpenMP threads
    call split_indices(grid%ntriangle, sml_nthreads, itr_beg, itr_end)

    !obtain perpendicular E-field -- total E-field
!!$OMP PARALLEL DO &
!!$OMP PRIVATE( ITH, ITR, NODES, DP1, DP2)
    !do ith=1,sml_nthreads
       !do itr=itr_beg(ith), itr_end(ith)
       do itr=1, grid%ntriangle
          nodes(:)=grid%nd(:,itr)

          dp1=dpot(nodes(1))- dpot(nodes(3))
          dp2=dpot(nodes(2))- dpot(nodes(3))

          E_perp_tr(:,itr)= &
                  -( dp1*grid%mapping(1,1:2,itr) + dp2*grid%mapping(2,1:2,itr) )
       enddo
    !enddo

    ! Divide nodes among OpenMP threads
    call split_indices(grid%nnode, sml_nthreads, i_beg, i_end)

    ! Get area average of perp E-field

!!$OMP PARALLEL DO &
!!$OMP PRIVATE( ITH, I, AREA_SUM, E, J, ITR )
    !do ith=1,sml_nthreads
       !do i=i_beg(ith),i_end(ith)
       do i=1, grid%nnode
          area_sum=0D0
          E=0D0
          do j=1, grid%num_t_node(i)
             itr=grid%tr_node(j,i)
             area_sum=area_sum+grid%tr_area(itr)
             E(:)=E(:)+ E_perp_tr(:,itr) * grid%tr_area(itr)
          enddo
          E_r(i)=E(1)/area_sum
          E_z(i)=E(2)/area_sum
       enddo
    !enddo

    deallocate(E_perp_tr)

  else
    !rh 2nd order FD gradient
    call grid_deriv(grid,dpot,E_r,E_z)
    !rh we need -grad(dpot)
    E_r=-E_r
    E_z=-E_z
  endif

  deallocate(dpot)

  return

#endif

end subroutine get_dpot_grad


!rh Evaluate dB_perp = grad_perp(A_par) x b
subroutine get_dbperp(grid,A_par,delta_B_perp)
  use grid_class
  use omp_module
  use sml_module, only: sml_bt_sign, sml_nthreads, sml_old_grad_perp
  implicit none

  type(grid_type), intent(in) :: grid
  real (kind=8), intent(in) :: A_par(grid%nnode)
  real (kind=8), dimension(grid%nnode,2), intent(out) :: delta_B_perp

  real (kind=8), dimension(:,:), allocatable :: B_perp_tr
  real (kind=8) :: area_sum, B_perp(2), dp1, dp2, sgn
  integer, dimension(sml_nthreads) :: itr_beg, itr_end, i_beg, i_end
  integer :: nodes(3), i, j, itr, ith

  delta_B_perp(:,:) = 0D0
  sgn=sml_bt_sign

  if (sml_old_grad_perp) then
    !rh Old FE gradient

    allocate(B_perp_tr(2,grid%ntriangle))


    ! Divide triangles among OpenMP threads
    call split_indices(grid%ntriangle, sml_nthreads, itr_beg, itr_end)

    ! obtain grad(A_parallel)

!$OMP PARALLEL DO &
!$OMP PRIVATE( ITH, ITR, NODES, DP1, DP2)
    !do ith=1,sml_nthreads
       ! store grad-2d A_par in B_perp_tr
       !do itr=itr_beg(ith), itr_end(ith)
       do itr=1, grid%ntriangle
          nodes(:)=grid%nd(:,itr)

          dp1=A_par(nodes(1))- A_par(nodes(3))
          dp2=A_par(nodes(2))- A_par(nodes(3))
          ! To check whether this calculation is correct
          B_perp_tr(1,itr)= -sgn &        ! -sgn * d A_par/ dr
                 * ( dp1*grid%mapping(1,1,itr) + dp2*grid%mapping(2,1,itr) )
          B_perp_tr(2,itr)= sgn &       ! sgn * d A_par/ dz
                 * ( dp1*grid%mapping(1,2,itr) + dp2*grid%mapping(2,2,itr) )
       enddo
    !enddo

    !Second convert the B_perp value from each triangle to each node

    call split_indices(grid%nnode, sml_nthreads, i_beg, i_end)

    ! Get area average of perp B-field

!$OMP PARALLEL DO &
!$OMP PRIVATE( ITH, I, AREA_SUM, B_perp, J, ITR )
    !do ith=1,sml_nthreads

       !do i=i_beg(ith),i_end(ith)
       do i=1, grid%nnode
          area_sum=0D0
          B_perp=0D0
          do j=1, grid%num_t_node(i)
             itr=grid%tr_node(j,i)
             area_sum=area_sum+grid%tr_area(itr)
             B_perp(:)=B_perp(:)+ B_perp_tr(:,itr) * grid%tr_area(itr)
          enddo
          !rh Note the swap in R-Z -->
          delta_B_perp(i,1)=B_perp(2)/area_sum  ! delta_B_perp_r = -d(A_par)/dz
          delta_B_perp(i,2)=B_perp(1)/area_sum ! delta_B_perp_z = d(A_par)/dr
          !rh More accurate expression for delta_B_perp --->
          !rh delta_B_perp(i,1)=B_perp(2)/area_sum + grid%curl_nb(1,i)*A_par(i) ! delta_B_perp_r + curl_nb
          !rh delta_B_perp(i,2)=B_perp(1)/area_sum + grid%curl_nb(2,i)*A_par(i) ! delta_B_perp_z + culr_nb
       enddo

    !enddo

    deallocate(B_perp_tr)

  else
    !rh 2nd order FD gradient
    call grid_deriv(grid,A_par,delta_B_perp(:,2),delta_B_perp(:,1))
    !rh -d(A_par)/dr
    delta_B_perp(:,2)=-delta_B_perp(:,2)
    delta_B_perp(:,:)=sgn*delta_B_perp
    !rh delta_b_perp(:,1)=delta_B_perp(:,1)+grid%curl_nb(1,i)*A_par(i)
    !rh delta_b_perp(:,2)=delta_B_perp(:,2)+grid%curl_nb(2,i)*A_par(i)
  endif

  return

end subroutine get_dbperp


subroutine kink_drive(grid,delta_B_perp,dB_grad_u0)
  use grid_class
  use sml_module, only: sml_e_charge, sml_old_grad_perp
  implicit none
  type(grid_type), intent(in) :: grid
  real (kind=8), intent(in) :: delta_B_perp(grid%nnode,2)
  real (kind=8), intent(out) :: dB_grad_u0(grid%nnode)
  real (kind=8) :: u0_par(grid%nnode)
  integer :: i

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!! Evaluate kink term dB_perp.grad(j_0/(e*B))
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do i=1, grid%nnode
     u0_par(i) = -grid%j0_par(i)/(sml_e_charge*grid%bfield(4,i))
  enddo

  if (sml_old_grad_perp) then
    call b_dot_grad_u(dB_grad_u0, delta_B_perp(:,1), delta_B_perp(:,2), u0_par, grid)
  else
    call v_dot_grad_perp(grid,delta_B_perp,u0_par,dB_grad_u0)
  endif

  do i=1, grid%nnode
    ! Multiply by B_phi/B_tot
    dB_grad_u0(i)=dB_grad_u0(i)*grid%bfield(3,i)/grid%bfield(4,i)
  enddo

  return

end subroutine kink_drive

subroutine kink_drive2(grid,A_par,dB_grad_u0)
  use grid_class
  use sml_module, only: sml_e_charge, sml_old_grad_perp
  implicit none
  type(grid_type), intent(in) :: grid
  real (kind=8), intent(in) :: A_par(grid%nnode)
  real (kind=8), intent(out) :: dB_grad_u0(grid%nnode)
  integer :: i
  real (kind=8), dimension(grid%nnode) :: dadr, dadz

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!! Evaluate kink term dB_perp.grad(j_0/(e*B))
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !call grid_deriv(grid,a_par,dadr,dadz)
  !do i=1,grid%nnode
  !  dB_grad_u0(i)=grid%tearing_drive(i,1)*dadr(i) + grid%tearing_drive(i,2)*dadz(i)
  !enddo

  call mat_mult(grid%kink_mat,A_par,dB_grad_u0)

  return

end subroutine kink_drive2

subroutine grad_par_delta_pe(grid,tr,p,dx,bd,deltape,E_para_part)
  use grid_class
  use boundary_class
  use sml_module, only: sml_e_charge
  implicit none
  type(grid_type), intent(in) :: grid
  type(boundary2_type), intent(in) :: bd
#ifdef OLD_INIT_FF
  integer, intent(in) :: tr(2,grid%nnode)
  !rh note that with OLD_INIT_FF, dx=1/dl_parallel
  real (kind=8), intent(in) :: p(3,2,grid%nnode), dx(grid%nnode)
#else
  integer, intent(in) :: tr(grid%nnode,0:1)
  real (kind=8), intent(in) :: p(3,grid%nnode,0:1), dx(grid%nnode,0:1)
#endif
  real (kind=8), intent(in) :: deltape(grid%nnode)
  real (kind=8), intent(out) :: E_para_part(grid%nnode)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!! Evaluate b.grad(dP_e)/(e n_0)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Get the total E_para using the closure of E_para_tot=-grad_para Pe-

  ! get -grad_para deltape for psn%E_para_tot
  call GradParX(grid, tr, p, dx, bd, deltape, E_para_part)
  E_para_part = - E_para_part/(sml_e_charge*grid%dene) !-grad_para

  return

end subroutine grad_par_delta_pe


