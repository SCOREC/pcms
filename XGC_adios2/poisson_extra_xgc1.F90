#define XGC1 1
#ifdef XGC1
!!$#include <petscversion.h>
!!$#if PETSC_VERSION_LT(3,6,0)
!!$#include <finclude/petscsnesdef.h>
!!$#else
!!$#include <petsc-finclude/petscsnesdef.h>
!!$#endif

! send and receive other plane potential
subroutine send_recv_potential( dpot, recvr,sendl,nnode)
  use sml_module
  use perf_monitor
  implicit none
  include 'mpif.h'
  integer , intent(in) :: nnode
  real (kind=8) :: recvr(nnode), sendl(nnode), dpot(nnode,-1:2)
  integer :: icount, idest,isource,isendtag,irecvtag
  integer :: ierr
  integer :: istatus(mpi_status_size)

  ! receive 0
  sendl=dpot(:,1)
  recvr=0D0
  icount=nnode
  isource= modulo(sml_intpl_mype-1,sml_intpl_totalpe)
  idest  = modulo(sml_intpl_mype+1,sml_intpl_totalpe)
  isendtag=sml_intpl_mype
  irecvtag=isource

  call mpi_sendrecv(sendl,icount,MPI_DOUBLE_PRECISION,idest,isendtag,&
       recvr,icount,MPI_DOUBLE_PRECISION,isource,irecvtag,sml_intpl_comm,istatus,ierr)

  dpot(:,0) = recvr(:)

  ! receive -1
  sendl=dpot(:,0)
  recvr=0D0
  icount=nnode
  isource= modulo(sml_intpl_mype-1,sml_intpl_totalpe)
  idest  = modulo(sml_intpl_mype+1,sml_intpl_totalpe)

  isendtag=sml_intpl_mype
  irecvtag=isource

  call mpi_sendrecv(sendl,icount,MPI_DOUBLE_PRECISION,idest,isendtag,&
       recvr,icount,MPI_DOUBLE_PRECISION,isource,irecvtag,sml_intpl_comm,istatus,ierr)

  dpot(:,-1) = recvr(:)

  ! receive 1
  sendl=dpot(:,1)
  recvr=0D0
  icount=nnode
  isource= modulo(sml_intpl_mype+1,sml_intpl_totalpe)
  idest  = modulo(sml_intpl_mype-1,sml_intpl_totalpe)
  isendtag=sml_intpl_mype
  irecvtag=isource
  call mpi_sendrecv(sendl,icount,MPI_DOUBLE_PRECISION,idest,isendtag,&
       recvr,icount,MPI_DOUBLE_PRECISION,isource,irecvtag,sml_intpl_comm,istatus,ierr)

  dpot(:,2) = recvr(:)

end subroutine send_recv_potential


#ifndef EFIELD2
subroutine get_potential_grad(grid,psn)
  use psn_class
  use grid_class
  use sml_module
  use smooth_module
  use omp_module , only : split_indices
  use perf_monitor
  implicit none
  include 'mpif.h'
  type(grid_type) :: grid
  type(psn_type) :: psn
  integer :: i,iphi,itr,nodes(3),j
  real (kind=8) :: dp1,dp2
  integer :: iphi_p, iphi_m
  real (kind=8) :: dpot_p, dpot_m
  integer :: nodes_p(3),nodes_m(3),nd,sgn
  real (kind=8) :: area_sum, E0(2), E(2,0:1)
  integer :: irho, ipe
  ! mpi
  integer :: idest, isendtag, isource, irecvtag
  integer :: ierr
  integer, dimension(MPI_STATUS_SIZE) :: istatus
  ! omp
  integer :: itr_beg(sml_nthreads), itr_end(sml_nthreads)
  integer :: ith, i_beg(sml_nthreads), i_end(sml_nthreads)

  real (8), allocatable :: E_perp_tr(:,:,:), E_r_node(:,:),E_z_node(:,:), E_para(:,:)
  real (8), allocatable :: pot(:,:,:), E_rho(:,:,:,:)
  ! for electron simulation
  real (8), allocatable :: E_perp00_tr(:,:), E_r00_node(:,:), E_z00_node(:,:)
  character (len=30) :: filename
  real (8) :: turb_on, zero_on

  integer :: nphi, iphim1,iphip1,dir
  real(8) :: pot_mid, pot_i_0,pot_i_1

  integer, parameter :: idebug = 0

!pw  call t_startf("GET_POT_INIT")
  allocate(E_r_node(grid%nnode,0:1),E_z_node(grid%nnode,0:1), E_para(grid%nnode,0:1))
  allocate(E_rho(3,grid%nnode,0:1,0:grid%nrho))
  allocate(pot(grid%nnode,0:1,0:max(2,grid%nrho)))
  if(sml_electron_on) then
     allocate(E_r00_node(grid%nnode,0:1), E_z00_node(grid%nnode,0:1))
  endif

  if(sml_turb_efield) then
     turb_on=1D0
  else
     turb_on=0D0
  endif
  if(sml_00_efield) then
     zero_on=1D0
  else
     zero_on=0D0
  endif

  ! pot(:,:,1) as tmp variable
  do iphi=0,1
     pot(:,iphi,1)=zero_on*psn%pot0(:)+turb_on*psn%dpot(:,iphi)
  enddo
#ifdef  USE_CALC_GRADIENT
! ------------------------------------------------------
! copy and communicate psn%pot_phi_real(1:nnode,0:nphim1)
! which is used in computing the gradient
! E_phi_ff(1:3,0:1,1:nnode,0:nphim1)
! ------------------------------------------------------

  call t_startf("GATHER_POT_PHI_COPY")
  do iphi=0,sml_nphi_total-1
    if (sml_intpl_mype == iphi) then
      psn%pot_phi_real(:,iphi) = pot(:,0,1)
    endif
  enddo
  call t_stopf("GATHER_POT_PHI_COPY")

  call t_startf("GATHER_POT_PHI_GATHER")
  call mpi_allgather(MPI_IN_PLACE, grid%nnode,mpi_real8, &
                     psn%pot_phi_real(:,:),grid%nnode, mpi_real8,  &
                     sml_intpl_comm,ierr)
  if (ierr.ne.mpi_success) then
     write(*,*) 'mpi_allgather(psn%pot_phi_real) return ierr=',ierr
     stop '** error in  get_potential_grad'
  endif
  call t_stopf("GATHER_POT_PHI_GATHER")



    if ((idebug.ge.1).and.(sml_mype.eq.0)) then
      write(6,*) 'sml_totalpe ', sml_totalpe
      write(6,*) 'sml_intpl_totalpe ', sml_intpl_totalpe
      write(6,*) 'sml_plane_totalpe ', sml_plane_totalpe
      write(6,*) 'sml_nphi_total ', sml_nphi_total
      write(6,*) 'sml_grid_nrho ', sml_grid_nrho
      write(6,*) 'sml_intpl_mype ', sml_intpl_mype
      write(6,*) 'sml_plane_mype ', sml_plane_mype
      write(6,*) 'sml_mype ', sml_mype
      call flush(6)
    endif


#endif
!pw  call t_stopf("GET_POT_INIT")

  call t_startf("GET_POT_LOOPS")

  if (sml_old_grad_perp) then
    !rh Old FE perpendicular gradient operator

    allocate(E_perp_tr(2,grid%ntriangle,0:1))
    if(sml_electron_on) then
       allocate(E_perp00_tr(2,grid%ntriangle))
    endif


    ! Divide triangles among OpenMP threads
    call split_indices(grid%ntriangle, sml_nthreads, itr_beg, itr_end)

    !obtain perpendicular E-field -- total E-field
!$OMP PARALLEL DO &
!$OMP PRIVATE( ITH, ITR, NODES, DP1, DP2, IPHI)
    do ith=1,sml_nthreads

       ! store grad-2d pot0 + pot in E_perp
       do iphi=0, 1
          do itr=itr_beg(ith), itr_end(ith)
             nodes(:)=grid%nd(:,itr)

             dp1=pot(nodes(1),iphi,1)- pot(nodes(3),iphi,1)
             dp2=pot(nodes(2),iphi,1)- pot(nodes(3),iphi,1)

             E_perp_tr(:,itr,iphi)= &
                  -( dp1*grid%mapping(1,1:2,itr) + dp2*grid%mapping(2,1:2,itr) )
          enddo
       enddo

    enddo

    !obtain perpendicular 00 E-field -- required for weight evolution of electron
    if(sml_electron_on .and. sml_00_efield) then
!$OMP PARALLEL DO &
!$OMP PRIVATE( ITH, ITR, NODES, DP1, DP2)
       do ith=1,sml_nthreads

          ! store grad-2d pot0 + pot in E_perp
          do itr=itr_beg(ith), itr_end(ith)
             nodes(:)=grid%nd(:,itr)

             dp1=psn%pot0(nodes(1))- psn%pot0(nodes(3))
             dp2=psn%pot0(nodes(2))- psn%pot0(nodes(3))

             E_perp00_tr(:,itr)= &
                  -( dp1*grid%mapping(1,1:2,itr) + dp2*grid%mapping(2,1:2,itr) )
          enddo

       enddo
    endif



    ! Divide nodes among OpenMP threads
    call split_indices(grid%nnode, sml_nthreads, i_beg, i_end)

    ! Get area averge of perp E-field

!$OMP PARALLEL DO &
!$OMP PRIVATE( ITH, I, AREA_SUM, E, J, ITR )
    do ith=1,sml_nthreads

       do i=i_beg(ith),i_end(ith)
          area_sum=0D0
          E=0D0
          do j=1, grid%num_t_node(i)
             itr=grid%tr_node(j,i)
             area_sum=area_sum+grid%tr_area(itr)
             E(:,:)=E(:,:)+ E_perp_tr(:,itr,0:1) * grid%tr_area(itr)
          enddo
          E_r_node(i,:)=E(1,:)/area_sum
          E_z_node(i,:)=E(2,:)/area_sum
       enddo

    enddo

    if(sml_electron_on .and. sml_00_efield) then
!$OMP PARALLEL DO &
!$OMP PRIVATE( ITH, I, AREA_SUM, E, J, ITR )
       do ith=1,sml_nthreads

          do i=i_beg(ith),i_end(ith)
             area_sum=0D0
             E(:,1)=0D0
             do j=1, grid%num_t_node(i)
                itr=grid%tr_node(j,i)
                area_sum=area_sum+grid%tr_area(itr)
                E(:,1)=E(:,1)+ E_perp00_tr(:,itr) * grid%tr_area(itr)
             enddo
             E_r00_node(i,1)=E(1,1)/area_sum
             E_z00_node(i,1)=E(2,1)/area_sum
          enddo

       enddo
       E_r00_node(:,0)=E_r00_node(:,1)
       E_z00_node(:,0)=E_z00_node(:,1)

    endif

    deallocate(E_perp_tr)
    if (sml_electron_on) then
      deallocate(E_perp00_tr)
    endif

  else

    do iphi=0,1
      call grid_deriv(grid,pot(:,iphi,1),E_r_node(:,iphi),E_z_node(:,iphi))
      !rh we need -grad(pot) --->
      E_r_node(:,iphi)=-E_r_node(:,iphi)
      E_z_node(:,iphi)=-E_z_node(:,iphi)
    enddo
    if (sml_electron_on .and. sml_00_efield) then
      call grid_deriv(grid,psn%pot0(:),E_r00_node(:,1),E_z00_node(:,1))
      E_r00_node(:,1)=-E_r00_node(:,1)
      E_z00_node(:,1)=-E_z00_node(:,1)
      E_r00_node(:,0)=E_r00_node(:,1)
      E_z00_node(:,0)=E_z00_node(:,1)
    endif

  endif


#if defined(XGC1_EM) && defined(EM_NO_PHI00)
    !rh We can use b.grad(phi) in the electron fluid part
    !rh Here, E_para is still in real space
    psn%E_perp_node(1,:,:)=E_r_node(:,:)
    psn%E_perp_node(2,:,:)=E_z_node(:,:)
#endif


  !apply boundary condition of E-field
  !zero E-field (plus zero potential at boundary)

!  call set_boundary2_values(E_r_node(:,0),0D0,psn%pbd0_2)
!  call set_boundary2_values(E_r_node(:,1),0D0,psn%pbd0_2)
!  call set_boundary2_values(E_z_node(:,0),0D0,psn%pbd0_2)
!  call set_boundary2_values(E_z_node(:,1),0D0,psn%pbd0_2)
!  if(sml_electron_on)then
!  call set_boundary2_values(E_r00_node(:,0),0D0,psn%pbd0_2)
!  call set_boundary2_values(E_r00_node(:,1),0D0,psn%pbd0_2)
!  call set_boundary2_values(E_z00_node(:,0),0D0,psn%pbd0_2)
!  call set_boundary2_values(E_z00_node(:,1),0D0,psn%pbd0_2)
!  endif

  call t_stopf("GET_POT_LOOPS")

  ! convert field following coord
  call t_startf("GET_POT_CNVRT")
  call cnvt_grid_real2ff(grid,psn%ff_hdp_tr,psn%ff_hdp_p,pot(:,:,1),pot(:,:,0))
  call cnvt_grid_real2ff(grid,psn%ff_hdp_tr,psn%ff_hdp_p,E_r_node,pot(:,:,1))
  E_r_node=pot(:,:,1)
  call cnvt_grid_real2ff(grid,psn%ff_hdp_tr,psn%ff_hdp_p,E_z_node,pot(:,:,1))
  E_z_node=pot(:,:,1)
  if(sml_electron_on) then
     if(sml_00_efield) then
        call cnvt_grid_real2ff(grid,psn%ff_hdp_tr,psn%ff_hdp_p,E_r00_node,pot(:,:,1))
        E_r00_node=pot(:,:,1)
        call cnvt_grid_real2ff(grid,psn%ff_hdp_tr,psn%ff_hdp_p,E_z00_node,pot(:,:,1))
        E_z00_node=pot(:,:,1)
     else ! zero efield when sml_00_efield is .false.
        E_r00_node=0D0
        E_z00_node=0D0
     endif
  endif
  call t_stopf("GET_POT_CNVRT")

  ! E-parallel
  ! pot(:,0,1) is tmp variable
  ! obtain E_parallel at half position  --  1_dx is ~ 2 R dphi
#ifdef OLD_INIT_FF
  pot(:,0,1)= -sml_bt_sign*(pot(:,1,0)-pot(:,0,0))*psn%bfollow_1_dx(:)*2D0  ! ### replace bfollow_1_dx
#else
  pot(:,0,1)= -sml_bt_sign*(pot(:,1,0)-pot(:,0,0))/sum(psn%ff_hdp_dx,dim=2)
#endif


  call t_startf("GET_POT_SR")
  ! send left recv right
  idest=mod(sml_intpl_mype-1+sml_intpl_totalpe,sml_intpl_totalpe)
  isendtag=sml_intpl_mype
  isource=mod(sml_intpl_mype+1,sml_intpl_totalpe)
  irecvtag=isource
  call mpi_sendrecv(pot(:,0,1),grid%nnode,MPI_REAL8,  idest,isendtag, &
       E_para(:,1),grid%nnode,MPI_REAL8,isource,irecvtag, &
       sml_intpl_comm,istatus,ierr)

  ! send right recv left
  idest=mod(sml_intpl_mype+1,sml_intpl_totalpe)
  isendtag=sml_intpl_mype
  isource=mod(sml_intpl_mype-1+sml_intpl_totalpe,sml_intpl_totalpe)
  irecvtag=isource
  call mpi_sendrecv(pot(:,0,1),grid%nnode,MPI_REAL8,  idest,isendtag, &
       E_para(:,0),grid%nnode,MPI_REAL8,isource,irecvtag, &
       sml_intpl_comm,istatus,ierr)
  call t_stopf("GET_POT_SR")

  ! require converting of E_para(:,0) and E_perp(:,1) ????
  call t_startf("GET_POT_CNVRT")
  call cnvt_grid_real2ff(grid,psn%ff_1dp_tr,psn%ff_1dp_p,E_para(:,:),pot(:,:,2)) ! ???
  E_para(:,0)=0.5D0*(pot(:,0,2)+pot(:,0,1))
  E_para(:,1)=0.5D0*(pot(:,1,2)+pot(:,0,1))
  call t_stopf("GET_POT_CNVRT")

#if defined(XGC1_EM) && defined(EM_NO_PHI00)
  !rh We can use b.grad(phi) in the electron fluid part
  call cnvt_grid_ff2real(grid,psn%ff_hdp_tr,psn%ff_hdp_p,E_para(:,:),psn%E_para(:,:))
#endif

  ! obtain E_para_ff
  !E_para(:,0)=0.5D0*(E_para(:,0)+pot(:,0,1))
  !E_para(:,1)=0.5D0*(E_para(:,1)+pot(:,0,1))


  ! ----------- end of parallel E-field

  ! assign to E_rho
  E_rho(1,:,:,0)=E_r_node
  E_rho(2,:,:,0)=E_z_node
  E_rho(3,:,:,0)=E_para

  ! obtaion E_rho(:,:,n) and pot(:,:,n)
  if(sml_plane_mype<grid%nrho*2) then
     call t_startf("GET_POT_MAT_MULT")
     iphi=mod(sml_plane_mype,2)
     irho=sml_plane_mype/2  + 1

     call mat_mult(psn%gyro_avg_mat,pot(:,iphi,0),pot(:,iphi,irho))
     call mat_mult_tensor(psn%gyro_avg_mat,E_rho(:,:,iphi,0),3,E_rho(:,:,iphi,irho))
     call t_stopf("GET_POT_MAT_MULT")
  endif

  ! get all of E_rho and pot using bcast
  call t_startf("GET_POT_BCAST")
  do irho=1,grid%nrho
     do iphi=0,1
        ipe=(irho-1)*2 + iphi
        call mpi_bcast(pot(:,iphi,irho),grid%nnode,mpi_real8,ipe,sml_plane_comm,ierr)
        call mpi_bcast(E_rho(:,:,iphi,irho),grid%nnode*3,mpi_real8,ipe,sml_plane_comm,ierr)
     enddo
  enddo
  call t_stopf("GET_POT_BCAST")

#ifdef IDEN_DEBUG
  !for debuging
  if(sml_istep==200 .and. sml_ipc==1 .and. sml_mype==0) then
     do irho=1, grid%nrho
        do iphi=0,1
           write(filename,'("pdebug.",i2.2,".",i1.1,".txt")') irho,iphi
           open(unit=333,file=filename,status='replace')
           do i=1, grid%nnode
              write(333,1000) pot(i,iphi,irho), E_rho(:,i,iphi,irho)
           enddo
        enddo
     enddo
  endif
1000 format(18(e19.13,1x))
#endif


  ! change of indexing order
  call t_startf("GET_POT_IDX_REORD")
  do irho=0,grid%nrho
     do iphi=0,1
        psn%pot_rho_ff(iphi,irho,:)=pot(:,iphi,irho)
        psn%E_rho_ff(:,iphi,irho,:) = E_rho(:,:,iphi,irho)
     enddo
  enddo

  if(sml_electron_on) then
     do iphi=0,1
        psn%E00_ff(1,iphi,:)=E_r00_node(:,iphi)
        psn%E00_ff(2,iphi,:)=E_z00_node(:,iphi)
     enddo
  endif
  call t_stopf("GET_POT_IDX_REORD")

  ! for debugging ######
!  do iphi=0,1
!     psn%E_perp_node(1:2,:,iphi)=E_rho(1:2,:,iphi,0)
!     psn%E_para(:,iphi)     =E_rho(3  ,:,iphi,0)
!  enddo

  deallocate(E_r_node,E_z_node,E_para,E_rho,pot)
  if(sml_electron_on) deallocate(E_r00_node, E_z00_node)

end subroutine get_potential_grad

! on ifndef EFIELD2 -->
#else

subroutine get_potential_grad(grid,psn)
  use psn_class
  use grid_class
  use sml_module
  use smooth_module
  use omp_module , only : split_indices
  use perf_monitor
  implicit none
  include 'mpif.h'
  type(grid_type) :: grid
  type(psn_type) :: psn
  integer :: i,iphi,itr,nodes(3),j
  real (kind=8) :: dp1,dp2
  integer :: iphi_p, iphi_m
  real (kind=8) :: dpot_p, dpot_m
  integer :: nodes_p(3),nodes_m(3),nd,sgn
  real (kind=8) :: area_sum, E0(2), E(2,0:1)
  integer :: irho, ipe
  ! mpi
  integer :: idest, isendtag, isource, irecvtag
  integer :: ierr
  integer, dimension(MPI_STATUS_SIZE) :: istatus
  ! omp
  integer :: itr_beg(sml_nthreads), itr_end(sml_nthreads)
  integer :: ith, i_beg(sml_nthreads), i_end(sml_nthreads)

  real (8), allocatable :: E_perp_tr(:,:,:), E_r_node(:,:),E_z_node(:,:), E_para(:,:)
  real (8), allocatable :: pot(:,:,:), E_rho(:,:,:,:)
  ! for electron simulation
  real (8), allocatable :: E_perp00_tr(:,:), E_r00_node(:,:), E_z00_node(:,:)
  real (8) :: turb_on, zero_on


  allocate(E_r_node(grid%nnode,0:1),E_z_node(grid%nnode,0:1), E_para(grid%nnode,0:1))
  allocate(E_rho(3,grid%nnode,0:1,0:grid%nrho))
 allocate(pot(grid%nnode,0:1,0:grid%nrho))
  if(sml_electron_on) then
     allocate(E_r00_node(grid%nnode,0:1), E_z00_node(grid%nnode,0:1))
  endif


  if(sml_turb_efield) then
     turb_on=1D0
  else
     turb_on=0D0
  endif
  if(sml_00_efield) then
     zero_on=1D0
  else
     zero_on=0D0
  endif

  ! pot(:,:,1) as tmp variable
  do iphi=0,1
     pot(:,iphi,1)=zero_on*psn%pot0(:)+turb_on*psn%dpot(:,iphi)
  enddo

  if (sml_old_grad_perp) then

    allocate(E_perp_tr(2,grid%ntriangle,0:1))
    if(sml_electron_on) then
       allocate(E_perp00_tr(2,grid%ntriangle))
    endif

    ! Divide triangles among OpenMP threads
    call split_indices(grid%ntriangle, sml_nthreads, itr_beg, itr_end)

    !obtain perpendicular E-field -- total E-field
!$OMP PARALLEL DO &
!$OMP PRIVATE( ITH, ITR, NODES, DP1, DP2, IPHI)
    do ith=1,sml_nthreads

       ! store grad-2d pot0 + pot in E_perp
       do iphi=0, 1
          do itr=itr_beg(ith), itr_end(ith)
             nodes(:)=grid%nd(:,itr)

             dp1=pot(nodes(1),iphi,1)- pot(nodes(3),iphi,1)
             dp2=pot(nodes(2),iphi,1)- pot(nodes(3),iphi,1)

             E_perp_tr(:,itr,iphi)= &
                  -( dp1*grid%mapping(1,1:2,itr) + dp2*grid%mapping(2,1:2,itr) )
          enddo
       enddo

    enddo


    !obtain perpendicular 00 E-field -- required for weight evolution of electron
    if(sml_electron_on .and. sml_00_efield) then
!$OMP PARALLEL DO &
!$OMP PRIVATE( ITH, ITR, NODES, DP1, DP2)
       do ith=1,sml_nthreads

          ! store grad-2d pot0 + pot in E_perp
          do itr=itr_beg(ith), itr_end(ith)
             nodes(:)=grid%nd(:,itr)

             dp1=psn%pot0(nodes(1))- psn%pot0(nodes(3))
             dp2=psn%pot0(nodes(2))- psn%pot0(nodes(3))

             E_perp00_tr(:,itr)= &
                  -( dp1*grid%mapping(1,1:2,itr) + dp2*grid%mapping(2,1:2,itr) )
          enddo

       enddo
    endif




    ! Divide nodes among OpenMP threads
    call split_indices(grid%nnode, sml_nthreads, i_beg, i_end)

    ! Get area averge of perp E-field

!$OMP PARALLEL DO &
!$OMP PRIVATE( ITH, I, AREA_SUM, E, J, ITR )
    do ith=1,sml_nthreads

       do i=i_beg(ith),i_end(ith)
          area_sum=0D0
          E=0D0
          do j=1, grid%num_t_node(i)
             itr=grid%tr_node(j,i)
             area_sum=area_sum+grid%tr_area(itr)
             E(:,:)=E(:,:)+ E_perp_tr(:,itr,0:1) * grid%tr_area(itr)
          enddo
          E_r_node(i,:)=E(1,:)/area_sum
          E_z_node(i,:)=E(2,:)/area_sum
       enddo

    enddo

    if(sml_electron_on .and. sml_00_efield) then
!$OMP PARALLEL DO &
!$OMP PRIVATE( ITH, I, AREA_SUM, E, J, ITR )
       do ith=1,sml_nthreads

          do i=i_beg(ith),i_end(ith)
             area_sum=0D0
             E(:,1)=0D0
             do j=1, grid%num_t_node(i)
                itr=grid%tr_node(j,i)
                area_sum=area_sum+grid%tr_area(itr)
                E(:,1)=E(:,1)+ E_perp00_tr(:,itr) * grid%tr_area(itr)
             enddo
             E_r00_node(i,1)=E(1,1)/area_sum
             E_z00_node(i,1)=E(2,1)/area_sum
          enddo

       enddo
       E_r00_node(:,0)=E_r00_node(:,1)
       E_z00_node(:,0)=E_z00_node(:,1)

    endif

    deallocate(E_perp_tr)
    if (sml_electron_on) then
      deallocate(E_perp00_tr)
    endif

  else

    do iphi=0,1
      call grid_deriv(grid,pot(:,iphi,1),E_r_node(:,iphi),E_z_node(:,iphi))
      !rh we need -grad(pot) --->
      E_r_node(:,iphi)=-E_r_node(:,iphi)
      E_z_node(:,iphi)=-E_z_node(:,iphi)
    enddo
    if (sml_electron_on .and. sml_00_efield) then
      call grid_deriv(grid,psn%pot0(:),E_r00_node(:,1),E_z00_node(:,1))
      E_r00_node(:,1)=-E_r00_node(:,1)
      E_z00_node(:,1)=-E_z00_node(:,1)
      E_r00_node(:,0)=E_r00_node(:,1)
      E_z00_node(:,0)=E_z00_node(:,1)
    endif

  endif

#if defined(XGC1_EM) && defined(EM_NO_PHI00)
  !rh We can use b.grad(phi) in the electron fluid part
  !rh Here, E_para is still in real space
  psn%E_perp_node(1,:,:)=E_r_node(:,:)
  psn%E_perp_node(2,:,:)=E_z_node(:,:)
#endif


  ! convert field following coord
  call cnvt_grid_real2ff(grid,psn%ff_hdp_tr,psn%ff_hdp_p,pot(:,:,1),pot(:,:,0))
  call cnvt_grid_real2ff(grid,psn%ff_hdp_tr,psn%ff_hdp_p,E_r_node,pot(:,:,1))
  E_r_node=pot(:,:,1)
  call cnvt_grid_real2ff(grid,psn%ff_hdp_tr,psn%ff_hdp_p,E_z_node,pot(:,:,1))
  E_z_node=pot(:,:,1)
  if(sml_electron_on) then
     if(sml_00_efield) then
        call cnvt_grid_real2ff(grid,psn%ff_hdp_tr,psn%ff_hdp_p,E_r00_node,pot(:,:,1))
        E_r00_node=pot(:,:,1)
        call cnvt_grid_real2ff(grid,psn%ff_hdp_tr,psn%ff_hdp_p,E_z00_node,pot(:,:,1))
        E_z00_node=pot(:,:,1)
     else ! zero efield when sml_00_efield is .false.
        E_r00_node=0D0
        E_z00_node=0D0
     endif
  endif


  !------- Parallel field -- difference of EFIELD2 starts here


#ifdef OLD_INIT_FF
  if(sml_bt_sign>0D0) then
     sgn=1
  else
     sgn=-1
  endif
#else
  sgn=1
#endif

  ! store grad_para (phi0 + phi) in E_para

!$OMP PARALLEL DO &
!$OMP PRIVATE( ITH, I, IPHI, IPHI_P, IPHI_M, NODES_P, NODES_M, DPOT_P, DPOT_M, ND )
  do ith=1,sml_nthreads

     do i=i_beg(ith),i_end(ith)
        do iphi=0, 1
           iphi_p=iphi+sgn
           iphi_m=iphi-sgn
#ifdef OLD_INIT_FF
           nodes_p=grid%nd(:,psn%bfollow_tr(1,i))
           nodes_m=grid%nd(:,psn%bfollow_tr(2,i))
#else
           nodes_p=grid%nd(:,psn%ff_1dp_tr(i,1))
           nodes_m=grid%nd(:,psn%bfollow_tr(i,0))
#endif
           dpot_p=0D0
           dpot_m=0D0
           do nd=1,3
#ifdef OLD_INIT_FF
              dpot_p=dpot_p + psn%dpot(nodes_p(nd),iphi_p)*psn%bfollow_p(nd,1,i)
              dpot_m=dpot_m + psn%dpot(nodes_m(nd),iphi_m)*psn%bfollow_p(nd,2,i)
#else
              dpot_p=dpot_p + psn%dpot(nodes_p(nd),iphi_p)*psn%ff_1dp_p(nd,i,1)
              dpot_m=dpot_m + psn%dpot(nodes_m(nd),iphi_m)*psn%ff_1dp_p(nd,i,0)
#endif
           enddo
#ifdef OLD_INIT_FF
           E_para(i,iphi) = -(dpot_p-dpot_m)*psn%bfollow_1_dx(i)    !(psn%pot0(i)
#else
           E_para(i,iphi) = -sml_bt_sign*(dpot_p-dpot_m)/sum(psn%ff_1dp_dx(i,:),dim=2)
#endif
        enddo
     enddo

  enddo

#if defined(XGC1_EM) && defined(EM_NO_PHI00)
  !rh We can use b.grad(phi) in the electron fluid part
  !rh Here, E_para is still in real space
  psn%E_para(:,:)=E_para(:,:)
#endif

  call cnvt_grid_real2ff(grid,psn%ff_hdp_tr,psn%ff_hdp_p,E_para,pot(:,:,1))
  E_para=pot(:,:,1)


  !--------- end of parallel E-field


  ! assign to E_rho
  E_rho(1,:,:,0)=E_r_node
  E_rho(2,:,:,0)=E_z_node
  E_rho(3,:,:,0)=E_para

  ! obtaion E_rho(:,:,n) and pot(:,:,n)
  if(sml_plane_mype<grid%nrho*2) then
     iphi=mod(sml_plane_mype,2)
     irho=sml_plane_mype/2  + 1

     call mat_mult(psn%gyro_avg_mat,pot(:,iphi,0),pot(:,iphi,irho))
     call mat_mult_tensor(psn%gyro_avg_mat,E_rho(:,:,iphi,0),3,E_rho(:,:,iphi,irho))
  endif

  ! get all of E_rho and pot using bcast
  do irho=1,grid%nrho
     do iphi=0,1
        ipe=(irho-1)*2 + iphi
        call mpi_bcast(pot(:,iphi,irho),grid%nnode,mpi_real8,ipe,sml_plane_comm,ierr)
        call mpi_bcast(E_rho(:,:,iphi,irho),grid%nnode*3,mpi_real8,ipe,sml_plane_comm,ierr)
     enddo
  enddo

  ! change of indexing order
  do irho=0,grid%nrho
     do iphi=0,1
        psn%pot_rho_ff(iphi,irho,:)=pot(:,iphi,irho)
        psn%E_rho_ff(:,iphi,irho,:) = E_rho(:,:,iphi,irho)
     enddo
  enddo

  if(sml_electron_on) then
     do iphi=0,1
        psn%E00_ff(1,iphi,:)=E_r00_node(:,iphi)
        psn%E00_ff(2,iphi,:)=E_z00_node(:,iphi)
     enddo
  endif

  ! for debugging ######
!  do iphi=0,1
!     psn%E_perp_node(1:2,:,iphi)=E_rho(1:2,:,iphi,0)
!     psn%E_para(:,iphi)     =E_rho(3  ,:,iphi,0)
!  enddo

  deallocate(E_r_node,E_z_node,E_para,E_rho,pot)
  if(sml_electron_on) deallocate(E_r00_node, E_z00_node)

end subroutine get_potential_grad
#endif


subroutine get_dpot_ff(grid,psn,irk)
  use psn_class
  use grid_class
  use sml_module
  use smooth_module
  use omp_module , only : split_indices
  use perf_monitor
  implicit none
  include 'mpif.h'
  type(grid_type) :: grid
  type(psn_type) :: psn
  integer, intent(in) :: irk
  !
  real (8) :: tmp(grid%nnode,0:1)

  ! convert field following coord
  call cnvt_grid_real2ff(grid,psn%ff_hdp_tr,psn%ff_hdp_p,psn%dpot(:,0:1),psn%dpot_ff)

  !-------------------------------
  if(.true. .and. irk>0) then
     tmp=psn%dpot_ff

     psn%dpot_ff = psn%dpot_ff + (psn%save_dpot(:,:,irk)-psn%save_dpot0(:,:,irk))
  
     psn%save_dpot0(:,:,irk)=tmp
  endif

  
end subroutine get_dpot_ff
subroutine save_dpot(grid,psn,irk)
  use grid_class
  use psn_class
  use sml_module
  implicit none
  type(grid_type) :: grid
  type(psn_type) :: psn
  integer, intent(in) :: irk
  !

  psn%save_dpot(:,:,irk)=psn%dpot_ff

end subroutine save_dpot

! mode selection - works only the nodes between boundaries
! refine ngrid setting for generalization
subroutine mode_selection(nmode,grid,psn)
  use sml_module
  use grid_class
  use psn_class
  implicit none
  include 'mpif.h'
  integer, intent(in) :: nmode
  type(grid_type) :: grid
  type(psn_type) :: psn
  integer :: ngrid,ngrid2,grid_offset,ngrid_per_pe,nphi_total,nplane_per_pe,npe_per_plane 
  integer :: iphi_loc,iphi_offset, iseg,idest,isource,isendtag,irecvtag,istart
  integer :: i,inode,itmp1,itmp2,itmp3,istatus(mpi_status_size),ierr
  real (kind=8) ,allocatable :: buff1(:), buff2(:,:),buff3(:)
  complex (kind=8), allocatable :: ftc(:),iftc(:)
  complex (kind=8) :: cn

  ! 1. initialization
  nphi_total=sml_nphi_total
  nplane_per_pe=sml_plane_per_pe  ! one
  npe_per_plane=sml_pe_per_plane
  iphi_offset=grid%iphi_offset


  grid_offset=1
  ngrid=grid%nnode
 
  ! Parallelizaing - ngrid_per_pe : number of grid (inside boundary) per cpu 
  ! ngrid2 - number of grid total computed.  Node number which is  > ngrid will be ignored later.
  ngrid_per_pe=ngrid/sml_totalpe 
  if(ngrid==sml_totalpe*ngrid_per_pe) then   
     ngrid2=ngrid  
  else
     ngrid_per_pe=ngrid_per_pe+1
     ngrid2=sml_totalpe*ngrid_per_pe
  endif

  allocate( buff1(ngrid2),buff3(ngrid2),buff2(ngrid_per_pe,nphi_total) )
  allocate( ftc(nphi_total), iftc(nphi_total) )

  ! preparing (inverse) Fourier transform coefficient.
  do i=1, nphi_total
!     ftc(i)=cexp(cmplx(0D0,sml_2pi*real((i-1)*nmode,8)/real(nphi_total,8),8))
     ftc(i)=cexp(cmplx(0D0,sml_2pi*real((i-1)*nmode)/real(nphi_total)))
     iftc(i)=conjg(ftc(i))/real(nphi_total,8)
  enddo

  ! store delta potential value to the buff1 - with grid offset.
  buff1(1:ngrid)=psn%dpot(grid_offset:grid_offset-1+ngrid,1)
  buff1(ngrid+1:ngrid2)=0D0


  ! 2. send - receive field data
  ! send buff1 data (given plane field data) to buff2 data (poloidal angle (or grid number) localized data)
  do iseg=1, nphi_total
     iphi_loc=1+(iseg-1)/sml_totalpe   ! should be one
     idest  =mod(sml_mype+npe_per_plane*(iseg-1),sml_totalpe)
     isource=mod(sml_mype-npe_per_plane*(iseg-1),sml_totalpe)  
     isource=mod(isource+sml_totalpe,sml_totalpe) ! to preven minus
     isendtag=sml_mype
     irecvtag=isource
!     istart= 1 + ngrid_per_pe*mod( sml_mype+iseg-1 , sml_totalpe )
     istart= 1 + ngrid_per_pe*idest
     !sending buff1(istart:istart+ngrid_per_pe,iphi)
     !receiving buff2(1:ngrid_per_pe,1+modulo(iphi_offset+iseg-1,nphi_total))
     itmp1=istart
     itmp2=istart+ngrid_per_pe-1
!     itmp3=mod(  iphi_offset + nplane_per_pe*(1-iseg) + iphi_loc-1 , nphi_total   )
!     itmp3= 1 + mod(itmp3 + nphi_total, nphi_total)
     itmp3=iphi_loc + isource*nplane_per_pe/npe_per_plane
     ! send and receive data
     call MPI_SENDRECV( buff1(itmp1:itmp2), ngrid_per_pe, MPI_DOUBLE_PRECISION,idest,isendtag,&
          buff2(1:ngrid_per_pe,itmp3), ngrid_per_pe, MPI_DOUBLE_PRECISION,isource,irecvtag,&
          sml_comm,istatus,ierr)
  enddo


  do inode=1, ngrid_per_pe
     
     ! 3. Fourier transform. Single mode direct cal  
     cn=0D0
     do i=1, nphi_total     
        cn= cn + buff2(inode,i)*ftc(i)
     enddo

     ! 4. inverse fouirer transform. Single mode direct cal
     do i=1, nphi_total
        buff2(inode,i)=real(cn*iftc(i))*2D0 ! 2D0 for -nmode summation : nmode and (-nmode) are complex conjugate and should be summed.
     enddo

  enddo

  ! 5. send - receive filtered data
  buff1=0D0
  do iseg=1, nphi_total
     iphi_loc=1+(iseg-1)/sml_totalpe
     idest  =mod(sml_mype-npe_per_plane*(iseg-1),sml_totalpe)
     idest  =mod(idest+sml_totalpe,sml_totalpe) ! to preven minus
     isource=mod(sml_mype+npe_per_plane*(iseg-1),sml_totalpe)  
     isendtag=sml_mype
     irecvtag=isource
!     istart= 1 + ngrid_per_pe*mod( sml_mype+iseg-1 , sml_totalpe )
     istart= 1 + ngrid_per_pe*isource
     !sending buff1(istart:istart+ngrid_per_pe,iphi)
     !receiving buff2(1:ngrid_per_pe,1+modulo(iphi_offset+iseg-1,nphi_total))
     itmp1=istart
     itmp2=istart+ngrid_per_pe-1
!     itmp3=mod(  iphi_offset + nplane_per_pe*(1-iseg) + iphi_loc-1 , nphi_total   )
!     itmp3= 1 + mod(itmp3 + nphi_total, nphi_total)
     itmp3=iphi_loc + idest*nplane_per_pe/npe_per_plane

     ! send and receive data
     call MPI_SENDRECV( buff2(1:ngrid_per_pe,itmp3), ngrid_per_pe,MPI_DOUBLE_PRECISION,idest,isendtag,&
          buff1(itmp1:itmp2), ngrid_per_pe, MPI_DOUBLE_PRECISION,isource,irecvtag,&
          sml_comm,istatus,ierr)
  enddo

  ! 6. broad cast between group (for multiple cpu per plane)
  if(npe_per_plane > 1 ) then
     call MPI_ALLREDUCE(buff1,buff3,ngrid2,MPI_DOUBLE_PRECISION, MPI_SUM,sml_plane_comm,ierr)
  else
     buff3=buff1
  endif
  do iphi_loc=1, nplane_per_pe
     psn%dpot(grid_offset:grid_offset-1+ngrid,1)=buff3(1:ngrid)
  enddo
  ! 7. finalize
  deallocate(buff1,buff2,buff3,ftc,iftc)

end subroutine mode_selection

!rh subroutine q_evaluation(grid)
!rh   use grid_class
!rh   implicit none
!rh   type(grid_type) :: grid
!rh   stop 'calling dummmy subroutine q_evaluation'
!rh end subroutine q_evaluation

subroutine mode_selection_comb(nmode,grid,psn,dpot)
  use sml_module
  use grid_class
  use psn_class
  use eq_module, only : eq_rgrid, eq_zgrid,eq_mr,eq_mz
  implicit none
  include 'mpif.h'
  integer, intent(in) :: nmode
  type(grid_type) :: grid
  type(psn_type) :: psn
  integer :: nphi_total
  integer :: i,inode, mres, npsi, ip, k, ir, nq
  real (kind=8) :: epsi, psi0
  integer ierr
  integer , allocatable :: itheta0(:), ntheta0(:)
  real (kind=8), intent(inout) :: dpot(grid%nnode)
  real (kind=8) ,allocatable :: allpot(:), dpot_m(:)
  complex (kind=8), allocatable :: cn(:), fmtc(:,:),ifmtc(:,:), cm(:)
  complex (kind=8) ::  phiphase


  nphi_total=sml_nphi_total
  mres = 1 ! only keep 5 m modes near m=nq resonant surface
  phiphase = cexp(cmplx(0D0,-sml_bt_sign*sml_2pi*grid%iphi_offset*nmode/real(nphi_total)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   To filter out the resonant (n,m) modes  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !3. do Fourier transform for each flux surface and Select the resonant m modes.
  do ip = 1, npsi
     allocate(dpot_m(ntheta0(ip)), allpot(ntheta0(ip)))
     allocate(fmtc(ntheta0(ip),mres))
     allocate(ifmtc(ntheta0(ip),mres))
     allocate(cm(mres),cn(mres))
     cm=0D0
     cn=0D0
     dpot_m=0D0
     allpot=0D0
     fmtc=0D0
     ifmtc=0D0

     dpot_m = dpot(itheta0(ip):itheta0(ip)+ntheta0(ip)-1)
     !  SK-HELP grid%qsafety is in EM
#ifdef XGC1_EM
     ! theta integral
     do k =1, ntheta0(ip)
        do ir=1, mres
           nq = (nint(nmode*grid%qsafety(ip))-(ir-(mres+1)/2))
!           if(sml_mype==0 .and. sml_gstep==1) print *, 'nq=', nmode*grid%qsafety(ip), nq, ir
           if(nq>0)fmtc(k,ir)=cexp(cmplx(0D0,sml_2pi*real((k-1)*nq)/real(ntheta0(ip))))
           ifmtc(k,ir)=conjg(fmtc(k,ir))/real(ntheta0(ip),8)
           cm(ir) = cm(ir) + dpot_m(k)*fmtc(k,ir)*phiphase
        enddo
     enddo
#endif
     call MPI_AllREDUCE(cm, cn, mres, MPI_DOUBLE_COMPLEX, MPI_SUM, sml_intpl_comm,ierr)

     do k =1, ntheta0(ip)
        do ir=1, mres
           allpot(k)=allpot(k)+real(cn(ir)*ifmtc(k,ir)*conjg(phiphase)/real(nphi_total))*2D0
        enddo
     enddo

     do k=1, ntheta0(ip)
        dpot(itheta0(ip)+k-1) = allpot(k)
     enddo

     deallocate(dpot_m, allpot)
     deallocate(fmtc, ifmtc)
     deallocate(cm,cn)
  enddo

 end subroutine mode_selection_comb

subroutine extract_00mode(grid,phi)
  use grid_class
  use sml_module
  use smooth_module
  use omp_module , only : split_indices
  use perf_monitor
  implicit none
  type(grid_type) :: grid
  real (kind=8) :: phi(grid%nnode,-1:2)
  real (kind=8) :: tmp(grid%nnode), tmp2(grid%npsi00)
  integer :: i, j
  integer :: ith, i_beg(sml_nthreads), i_end(sml_nthreads)

  ! Divide nodes among OpenMP threads
  call split_indices(grid%nnode, sml_nthreads, i_beg, i_end)

!$OMP PARALLEL DO &
!$OMP PRIVATE( ITH, I )
  do ith=1,sml_nthreads

     do i=i_beg(ith),i_end(ith)
        tmp(i)=phi(i,1)
     enddo

  enddo
  ! sum-up
  call t_startf("EXTRACT_00MODE_RED")
  call my_mpi_allreduce(tmp,grid%rtmp1,grid%nnode) !rtmp1 is used for a temporory purpose. variable
  call t_stopf("EXTRACT_00MODE_RED")
  ! get density from charge summation
  !call smooth_pol0(grid, grid%rtmp1, smooth00)
  call convert_grid_2_001d(grid,grid%rtmp1,tmp2)
  call convert_001d_2_grid(grid,tmp2,grid%rtmp1)

!$OMP PARALLEL DO &
!$OMP PRIVATE( ITH, I )
  do ith=1,sml_nthreads

     do i=i_beg(ith),i_end(ith)
        phi(i,1)=phi(i,1) - grid%rtmp1(i)/(real(sml_totalpe))
     enddo

  enddo

end subroutine extract_00mode

subroutine smooth_tr_connect(grid,var)
  use sml_module, only : sml_nthreads
  use grid_class
  use omp_module, only : split_indices
  implicit none
  type(grid_type) :: grid
  real (8) :: var(grid%nnode)
  real (8), allocatable :: tmp(:)
  integer :: nodes(3)
  integer :: i,j, itr
  real (8) :: area_sum, nval
  integer :: ith, i_beg(sml_nthreads), i_end(sml_nthreads)

  integer :: ismooth
  integer, parameter :: nsmooth=1

  allocate(tmp(grid%ntriangle))

  do ismooth=1, nsmooth

     ! average value of each triangle

     call split_indices(grid%ntriangle, sml_nthreads, i_beg, i_end)

!$OMP PARALLEL DO &
!$OMP PRIVATE( ITH, I, NODES )
     do ith=1,sml_nthreads
        do itr=i_beg(ith), i_end(ith)
           nodes(:)=grid%nd(:,itr)
           ! get triangle average value
           tmp(itr)=(var(nodes(1))+var(nodes(2))+var(nodes(3)))/3D0
        enddo
     enddo


     ! angle average of triangle value at a node

     call split_indices(grid%nnode, sml_nthreads, i_beg, i_end)

!$OMP PARALLEL DO &
!$OMP PRIVATE( ITH, I, AREA_SUM, NVAL, J, ITR )
     do ith=1,sml_nthreads
        do i=i_beg(ith),i_end(ith)
           area_sum=0D0
           nval=0D0
           do j=1, grid%num_t_node(i)
              itr=grid%tr_node(j,i)
              area_sum=area_sum+         grid%tr_area(itr)
              nval    =nval    +tmp(itr)*grid%tr_area(itr)
           enddo
           var(i)=nval/area_sum
        enddo
     enddo
     !
  enddo

  ! finalize
  deallocate(tmp)

end subroutine smooth_tr_connect


#ifdef OLD_SMOOTH_POL
!rh This code is outdated and not used anywhere
!rh I removed this to prevent a conflict with the new
!rh subroutine smooth_pol (from XGCa)

!
! smooth_pol and related routines are out of date
!
subroutine smooth_pol(grid,var)
  use grid_class
  use smooth_module
  implicit none
  type(grid_type):: grid
  real (8) :: var(grid%nnode)
  real (8), allocatable :: tmp(:)
  logical :: first=.true.
  save first

  print *, 'smoth_pol is out of date. Do not use it'
  stop


  if(first) then
     call init_smooth_pol_mat(grid)
     first=.false.
  endif

  allocate(tmp(grid%nnode))

  call mat_mult(smooth_pol_mat,var,tmp)
  var=tmp

  deallocate(tmp)
end subroutine smooth_pol

subroutine init_smooth_pol_mat(grid)
  use sml_module
  use grid_class
  use smooth_module
  implicit none
  type(grid_type):: grid
  integer :: iphi, irho
  real (8) :: rho, phi0, phi
  integer :: i, larmor
  real (8) :: x(2), x_ring(2)
  integer, parameter :: n_gyro=8
  real (8) :: dx_unit(2,n_gyro), dx_unit2(2)
  real (8) :: angle, dpdr, dpdz, dp, cosa, sina
  integer :: itr, node, ip
  real (8) :: p(3), weight
  !
  real (8), external :: psi_interpol

  call new_mat(smooth_pol_mat,grid%nnode,n_gyro*3+1)

  rho=grid%drho*real(irho)
  phi0=0.5D0*grid%delta_phi

  do larmor=1, n_gyro
     angle=sml_2pi/real(n_gyro)*real(larmor-1)
     dx_unit(:,larmor)=(/cos(angle),sin(angle)/)
  enddo

  do i=1, grid%nnode
     x=grid%x(:,i)

     ! self : 0.5-- minimum diagonal
     weight=0.5D0
     call set_value(smooth_pol_mat,i,i,weight,1)


     do larmor=1, n_gyro
        !find 8-point with rhoi
        dx_unit2(1)= dx_unit(1,larmor)
        dx_unit2(2)= dx_unit(2,larmor)

        x_ring = x + rho* dx_unit2(:)

        ! find triangle
        call search_tr2(grid,x_ring,itr,p)

        if(itr>0) then
           do ip=1,3
              weight=p(ip)/real(n_gyro)*0.5D0
              node=grid%nd(ip,itr)
              call set_value(smooth_pol_mat,i,node,weight,1)
           enddo
        endif
     enddo

  enddo

end subroutine init_smooth_pol_mat
!rh on #ifdef OLD_SMOOTH_POL
#endif


subroutine get_potential_grad_epara(grid,psn)
  use psn_class
  use grid_class
  use sml_module
  use smooth_module
  use omp_module , only : split_indices
  use perf_monitor
  implicit none
  include 'mpif.h'
  type(grid_type) :: grid
  type(psn_type) :: psn
  integer :: i,iphi,itr,nodes(3),j
  real (kind=8) :: dp1,dp2
  integer :: iphi_p, iphi_m
  real (kind=8) :: dpot_p, dpot_m
  integer :: nodes_p(3),nodes_m(3),nd,sgn
  real (kind=8) :: area_sum, E0(2), E(2,0:1)
  integer :: irho, ipe
  ! mpi
  integer :: idest, isendtag, isource, irecvtag, icount
  integer ierr
  integer, dimension(MPI_STATUS_SIZE) :: istatus
  ! omp
  integer :: itr_beg(sml_nthreads), itr_end(sml_nthreads)
  integer :: ith, i_beg(sml_nthreads), i_end(sml_nthreads)

  real (8), allocatable :: E_para(:,:)
  real (8), allocatable :: E_para_rho(:,:,:)
  ! for electron simulation
  real (8) :: turb_on


  allocate(E_para(grid%nnode,0:1))
  allocate(E_para_rho(grid%nnode,0:1,0:grid%nrho))


  if(sml_turb_efield) then
     turb_on=1D0
  else
     turb_on=0D0
  endif



  ! Divide nodes among OpenMP threads
  call split_indices(grid%nnode, sml_nthreads, i_beg, i_end)


  ! send E_para(:,1) and recv E_para(:,:0)
  icount=grid%nnode
  isource= modulo(sml_intpl_mype-1,sml_intpl_totalpe)
  idest  = modulo(sml_intpl_mype+1,sml_intpl_totalpe)
  isendtag=sml_intpl_mype
  irecvtag=isource
  !  SK-HELP psn%E_para_tot is in EM
#ifdef XGC1_EM
  call mpi_sendrecv(psn%E_para_tot(:),icount,MPI_DOUBLE_PRECISION,idest,isendtag,&
      E_para(:,0),icount,MPI_DOUBLE_PRECISION,isource,irecvtag,sml_intpl_comm,istatus,ierr)
  E_para(:,1)=psn%E_para_tot
#endif
  call cnvt_grid_real2ff(grid,psn%ff_hdp_tr,psn%ff_hdp_p,E_para,E_para_rho(:,:,0))
  E_para = E_para_rho(:,:,0)

  !--------- end of parallel E-field


  ! assign to E_rho
  !E_para_rho(:,:,0)=E_para

  ! obtaion E_rho(:,:,n) and pot(:,:,n)
  if(sml_plane_mype<grid%nrho*2) then
     iphi=mod(sml_plane_mype,2)
     irho=sml_plane_mype/2  + 1

     call mat_mult(psn%gyro_avg_mat,E_para(:,iphi),E_para_rho(:,iphi,irho))
  endif

  ! get all of E_rho and pot using bcast
  do irho=1,grid%nrho
     do iphi=0,1
        ipe=(irho-1)*2 + iphi
        call mpi_bcast(E_para_rho(:,iphi,irho),grid%nnode,mpi_real8,ipe,sml_plane_comm,ierr)
     enddo
  enddo

  ! change of indexing order
  do irho=0,grid%nrho
     do iphi=0,1
        psn%E_rho_ff(3,iphi,irho,:) = E_para_rho(:,iphi,irho)
     enddo
  enddo

  ! for debugging ######
!  do iphi=0,1
!     psn%E_perp_node(1:2,:,iphi)=E_rho(1:2,:,iphi,0)
!     psn%E_para(:,iphi)     =E_rho(3  ,:,iphi,0)
!  enddo

  deallocate(E_para,E_para_rho)

end subroutine

!
! on ifdef XGC1 --->
#endif
