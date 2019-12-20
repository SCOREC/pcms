#include "switches.h"
#include "redef.h"
#include "intrinsic_sizes.h"
!> Module containing velocity space (vsp) diagnostics
!! All diagnostics with special focus on velocity space
!! analysis and only one direct space coordinate at max
MODULE diagnostics_vsp
  USE aux_fields
  USE axpy
  USE communications
  USE coordinates, only: mu, mu_weight, vp, xval, deli
  USE diagnostics_auxiliary, only: calc_dfields_dy, &
       &mat_20, mat_01, nmat_00, nmat_20, nmat_01, &
       &mat_30, mat_11, pTj_qjB0
  USE file_io
  USE fourier, only: initialize_fourier_x_1d,&
       &to_fourier_x_1d,to_real_x_1d,finalize_fourier_x_1d,&
       &initialize_fourier_y_1d,to_real_y_1d,finalize_fourier_y_1d
  USE geometry, only: geom, Bref, Lref, minor_r, rhostar
  USE par_mod
  USE profile_IO, only: Tref, mref, nref, omegatorref
  USE vel_space, only: mat_00, mat_10, fm, gyro_op_wrapper
  use compute_f
#ifdef WITHFUTILS
  USE futils
#endif

  IMPLICIT NONE

  PUBLIC :: initialize_all_diags_vsp, exec_all_diags_vsp, &
    finalize_all_diags_vsp, mem_est_all_diags_vsp, check_all_diags_vsp

  PUBLIC :: istep_vsp, istep_xvsp, xvsp_iz, xvsp_iy

  PRIVATE

  INTEGER :: istep_vsp = -1, istep_xvsp = -1, xvsp_iy = -1, xvsp_iz = -1
  CHARACTER(LEN=8) :: filestat='replace', filepos='rewind'
  INTEGER :: VSPFILE
  INTEGER, DIMENSION(:), ALLOCATABLE :: XVSPFILE
  !hdf5 output
#ifdef WITHFUTILS
  ! files, timestep identifiers, and coordinates variables
  integer, parameter :: bufsize = 20
  integer :: fidvsp_h5, isnap_vsp = 0, isnap_xvsp = 0
  integer, allocatable, dimension(:) :: fidxvsp_h5
  character (len=5), dimension(1:5) :: vsp_label =(/'G_es', 'G_em',&
       & 'Q_es', 'Q_em', '<f_>'/)
  character (len=7), dimension(1:1) :: xvsp_label=(/'delta_f'/)
  character (len=32), dimension(:), allocatable :: basexvsp_h5
  real, dimension(:), allocatable :: f_unit
#endif


  !administration
  REAL :: last_exec_diag_time = -3.14159265

CONTAINS

!!!******************************************************************!!!

  !>Give an estimate of the memory requirements of this module
  REAL FUNCTION mem_est_all_diags_vsp(mem_req_in)
    REAL :: mem_req_in
    REAL :: mem_loc = 0.0

    IF (istep_vsp .GT. 0) mem_loc = mem_est_diag_vsp(mem_loc)
    IF (istep_xvsp .GT. 0) mem_loc = mem_est_diag_xvsp(mem_loc)

    mem_est_all_diags_vsp = mem_req_in + mem_loc
  END FUNCTION mem_est_all_diags_vsp

!!!******************************************************************!!!

  SUBROUTINE check_all_diags_vsp

    IF (istep_vsp .LT. 0) istep_vsp = 0
    IF (istep_xvsp .LT. 0) istep_xvsp = 0

    IF ((istep_xvsp.GT.0).AND.(.not.nonlinear)) THEN
       if (mype==0) write(*,"(A)") &
            &"WARNING: xvsp deactivated - only possible for nonlinear runs"
       istep_xvsp = 0
    ENDIF

    IF ((istep_xvsp.GT.0).AND.(yx_order)) THEN
       if (mype==0) write(*,"(A)") &
            &"WARNING: xvsp deactivated - not yet available for yx_order"
       istep_xvsp = 0
    ENDIF

    IF ((istep_xvsp.GT.0).AND.(.not.xy_local)) THEN
       if (mype==0) write(*,"(A)") &
            &"WARNING: xvsp deactivated - not yet available for global runs"
       istep_xvsp = 0
    ENDIF

    IF (istep_xvsp.GT.0) THEN
       IF (xvsp_iy.LT.0) xvsp_iy = 0
       IF (xvsp_iz.LT.0) xvsp_iz = nz0/2
    ENDIF

    IF ((comp_type.eq.'EV').AND.(istep_vsp.GT.0)) istep_vsp = 1

  END SUBROUTINE check_all_diags_vsp


!!!******************************************************************!!!

  !!Each individual GENE diagnostic initialization should be called
  !!by this routine
  SUBROUTINE initialize_all_diags_vsp

    IF (istep_vsp .GT. 0) CALL initialize_diag_vsp
    IF (istep_xvsp .GT. 0) CALL initialize_diag_xvsp

  END SUBROUTINE initialize_all_diags_vsp

!!!******************************************************************!!!
  SUBROUTINE exec_all_diags_vsp(itime,time)
    INTEGER, INTENT(IN) :: itime
    REAL, INTENT(IN) :: time

    IF (istep_vsp.GT.0) THEN
         IF (MODULO(itime,istep_vsp).eq.0) CALL diag_vsp
      ENDIF
    IF (istep_xvsp.GT.0) THEN
         IF (MODULO(itime,istep_xvsp).eq.0) CALL diag_xvsp
      ENDIF

  END SUBROUTINE exec_all_diags_vsp

!!!******************************************************************!!!

  !>Finalizes all GENE internal diagnostics
  SUBROUTINE finalize_all_diags_vsp

    IF (istep_vsp .GT. 0) CALL finalize_diag_vsp
    IF (istep_xvsp .GT. 0) CALL finalize_diag_xvsp

  END SUBROUTINE finalize_all_diags_vsp

!!!******************************************************************!!!
!!!******************************************************************!!!

  !>Give an estimate of the memory requirements of diag_vsp
  Real Function mem_est_diag_vsp(mem_req_in)
    real:: mem_req_in
    real:: mem_loc=0

    !local variables
    !var
    mem_loc = lklmn0*5.*SIZE_OF_REAL_MB
    !fullarr
    mem_loc = mem_loc + nz0*nv0*nw0*n_spec*5*SIZE_OF_REAL_MB
    !dfields_dy
    mem_loc = mem_loc + lij0*n_fields*SIZE_OF_COMPLEX_MB
    !phi_bar, temporary, dens_with_corr
    mem_loc = mem_loc + (2.*li0+li0)*lj0*SIZE_OF_COMPLEX_MB

    mem_est_diag_vsp = mem_loc+mem_req_in

  End Function mem_est_diag_vsp

!!!******************************************************************!!!

  Subroutine initialize_diag_vsp
#ifdef WITHFUTILS
    integer, dimension(2) :: dims
    integer :: o, rank
#endif

    if (mype==0) then
       if(write_std) then
          call get_unit_nr(VSPFILE)
          OPEN(VSPFILE, file=trim(diagdir)//&
               &'/vsp'//trim(file_extension), form='unformatted', &
               status=filestat, position=filepos)
       END IF

#ifdef WITHFUTILS
       if (write_h5) then
          isnap_vsp = 0
          call creatf(trim(diagdir)//'/vsp'//trim(file_extension)//'.h5', &
               fidvsp_h5, "Velocity space diagnostics", 'd')
          call creatg(fidvsp_h5, '/vsp')
          do o = 1,5
             call creatg(fidvsp_h5, '/vsp/'//trim(vsp_label(o)))
          enddo
          rank = 0
          call creatd(fidvsp_h5, rank, dims, "/vsp/time", "time")
       end if
#endif
    endif

  End Subroutine Initialize_diag_vsp

  !>VSP (velocity space) diagnostic
  !!
  !!Calculates guiding-center-fluxes G_es/G_em and Q_es/em
  !!and sqrt(<|f_|^2>) averaged over kx and ky
  Subroutine diag_vsp
    Real:: var(lk1:lk2, ll1:ll2, lm1:lm2,ln1:ln2, 5)
    Real:: fullarr(0:nz0-1, 0:nv0-1, 0:nw0-1, 0:n_spec-1, 5)
    Integer :: k, l, m, n, o, i, ii, pni, ierr
    Real :: fnorm
    Complex, Dimension(li1:li2, lj1:lj2,1:n_fields):: dfields_dy
    COMPLEX, DIMENSION(li1:li2, lj1:lj2):: phi_bar,temporary
    ! dens_with_corr: density with FLR correction
    Complex, Dimension(li1:li2, lj1:lj2):: dens_with_corr
#ifdef WITHFUTILS
    character(len=FILENAME_MAX):: dset_name_vsp
#endif

    !PERFON('diag_vsp')
    pni=pn1
    do n=ln1,ln2
       if (pn0.gt.1) pni=n
       do m=lm1,lm2
          do l=ll1,ll2
             do k = lk1, lk2
                !!************************************************************!
                !!******************* calculate bx and vx ********************!
                !!************************************************************!
                do o=1,n_fields
                   call calc_dfields_dy(dfields_dy(:,:,o),k,o)
                enddo

                CALL gyro_op_wrapper(emfields(:,:,k,1),phi_bar,k,m,n,0)

                if (xy_local) then
                   temporary = f_(:,:,k,l,m,n)+spec(n)%charge/spec(n)%temp*&
                        & phi_bar*fm(pi1,pj1,k,l,m,pni)
                else
                   do i=li1,li2
                      temporary(i,:) = f_(i,:,k,l,m,n)+spec(n)%charge/spec(n)%temp*&
                           & phi_bar(i,:)*fm(i,pj1,k,l,m,n)
                   enddo
                endif

                Call gyro_op_wrapper(temporary,dens_with_corr,k,m,n,0)

                if (xy_local) then
                   call axpy_ij(lij0,-spec(n)%charge/spec(n)%temp*fm(pi1,pj1,k,l,m,pni),&
                        emfields(:,:,k,1),dens_with_corr)
                else
                   do i=li1,li2
                      dens_with_corr(i,:) = dens_with_corr(i,:) - &
                           & spec(n)%charge/spec(n)%temp*fm(i,pj1,k,l,m,n)*&
                           & emfields(i,:,k,1)
                   enddo
                endif

                !             the statements above calculate:
                !             dens_with_corr = jfac(:,:,k,m,n)*f_(li1:li2,:,k,l,m,n)-&
                !                fm(k,l,m)*(1-jfac(:,:,k,m,n)**2)*phi(li1:li2,:,k)*&
                !                spec(n)%charge/spec(n)%temp

                !!************************************************************!
                !!******************* calculate transport ********************!
                !!************************************************************!
                !... G_es = <vx ( n + FLR corr)>
                var(k,l,m,n,1) = -2.0*Sum(&
                     Real(Conjg(dfields_dy(:,:,1))*dens_with_corr))

                !... G_em = <bx u_parallel>
                if (n_fields.gt.1) then
                   var(k,l,m,n,2) = 2.0*Sum(&
                        Real(Conjg(dfields_dy(:,:,2))*dens_with_corr))
                else
                   var(k,l,m,n,2) = 0
                endif

                !... <f_>
                var(k,l,m,n,5) = 2.0*Sum(&
                     Real(Conjg(f_(li1:li2,lj1:lj2,k,l,m,n))*f_(li1:li2,lj1:lj2,k,l,m,n)))

                !substract kj=0 mode (it should not be counted twice)
                if (p_has_0_mode) then
                   if (xy_local.and.yx_order) then !here ky=0 is at i=li1
                      var(k,l,m,n,5) = var(k,l,m,n,5)&
                           -Real(Sum(Conjg(f_(li1,:,k,l,m,n)) * f_(li1,:,k,l,m,n)))
                   else
                      var(k,l,m,n,5) = var(k,l,m,n,5)&
                           -Real(Sum(Conjg(f_(li1:li2,lj1,k,l,m,n)) * f_(li1:li2,lj1,k,l,m,n)))
                   endif
                endif
             End Do
          End Do
       End Do
    End Do

    !... Q_es = <(3/2) p vx>
    var(:,:,:,:,3) = var(:,:,:,:,1)

    !... Q_em = <[(5/2) vpl + qpl + qpp] bx>
    var(:,:,:,:,4) = var(:,:,:,:,2)

    !   Calculate sum over n_procs_x, leave result in pex = 0.
    Call my_sum_to_0_real(var, Size(var), mpi_comm_x)
    !   Calculate sum over n_procs_y, leave result in pey = 0.
    Call my_sum_to_0_real(var, Size(var), mpi_comm_y)
!WARNING need to implement (x,y) dep. in the following part :
    do n=ln1,ln2
       var(:,:,:,n,1) = var(:,:,:,n,1) * mat_00(pi1,pj1,:,:,:) * spec(n)%dens
       var(:,:,:,n,2) = var(:,:,:,n,2) * mat_10(pi1,pj1,:,:,:,n) * spec(n)%dens
       var(:,:,:,n,3) = var(:,:,:,n,3) * (mat_20(pi1,pj1,:,:,:) + &
            &mat_01(pi1,pj1,:,:,:)) * spec(n)%dens * spec(n)%temp
       var(:,:,:,n,4) = var(:,:,:,n,4) * (mat_30(pi1,pj1,:,:,:,n) + &
            &mat_11(pi1,pj1,:,:,:,n)) * spec(n)%dens * spec(n)%temp
       var(:,:,:,n,5) = Sqrt( var(:,:,:,n,5) )
    enddo

    fullarr = 0.0

    Do ii = 1,5
       Do n = ln1,ln2
          Do m = lm1,lm2
             Do l = ll1,ll2
                CALL mpi_gather(var(lk1,l,m,n,ii),lk0, MPI_REAL_TYPE,&
                     fullarr(0,l-ll1,m-lm1,n-ln1,ii),lk0,MPI_REAL_TYPE,&
                     0, mpi_comm_z, ierr)
             Enddo
             CALL my_real_gather_to_0(fullarr(:,:,m-lm1,n-ln1,ii),&
                  &1,nz0*ll0,nz0*nv0,mpi_comm_v)
          Enddo
          CALL my_real_gather_to_0(fullarr(:,:,:,n-ln1,ii),&
               &1,nz0*nv0*lm0,nz0*nv0*nw0,mpi_comm_w)
       Enddo
       CALL my_real_gather_to_0(fullarr(:,:,:,:,ii),&
            &1,nz0*nv0*nw0*ln0,nz0*nv0*nw0*n_spec,mpi_comm_spec)
    Enddo

    if (xy_local) then
       fnorm = 1.0
    else
       fnorm = 1.0 / REAL(ni0)
    endif

    If (mype == 0) Then
       fullarr = fullarr * fnorm
       if(write_std) then
          Write(VSPFILE) time
          Write(VSPFILE) fullarr
         call flush(VSPFILE)
       end if
    End If

#ifdef WITHFUTILS
    if(write_h5) then
       if(mype.eq.0) then
          call append(fidvsp_h5, "/vsp/time", real(time,8))
          call attach(fidvsp_h5, "/vsp/time", "n_steps", isnap_vsp+1)
          do o = 1,5
             write(dset_name_vsp, "(A, '/', i10.10)") "/vsp/"//&
                  &trim(vsp_label(o)), isnap_vsp
             call putarr(fidvsp_h5, dset_name_vsp, fullarr(:,:,:,:,o)*fnorm)
             call attach(fidvsp_h5, dset_name_vsp, "time", time)
          enddo

          call flushh5(fidvsp_h5)
       end if
       isnap_vsp = isnap_vsp+1
    end if
#endif

    Call my_barrier()
    !PERFOFF
  End Subroutine diag_vsp

  Subroutine finalize_diag_vsp
    If (mype==0) then
       if (write_std) CLOSE(VSPFILE)
#ifdef WITHFUTILS
       if (write_h5) call closef(fidvsp_h5)
#endif
    endif
  End Subroutine finalize_diag_vsp

!!!******************************************************************!!!
!!!******************************************************************!!!

  !>Give an estimate of the memory requirements of diag_xvsp
  Real Function mem_est_diag_xvsp(mem_req_in)
    real:: mem_req_in
    real:: mem_loc=0

    !local variables (neglect arrays in initialization since they are
    !smaller
    !var
    mem_loc = li0*ll0*lm0*ln0*SIZE_OF_REAL_MB
    !fullarr
    mem_loc = mem_loc + nx0*nv0*nw0*n_spec*SIZE_OF_REAL_MB

    mem_est_diag_xvsp = mem_loc+mem_req_in

  End Function mem_est_diag_xvsp


!!!******************************************************************!!!

  Subroutine initialize_diag_xvsp
    integer :: n
#ifdef WITHFUTILS
    integer, dimension(2) :: dims
    integer :: rank, i, l, m, pni, ierr
    real, dimension(:,:,:), allocatable :: var
    real, dimension(:,:,:), allocatable :: fullarr
    real, dimension(:), allocatable :: Rpos, Zpos, x_axis
    real :: Blocal
    real :: Qref = 1.6022e-19, vp_unit, mu_unit, f0_unit
    isnap_xvsp = 0
#endif

    if ((write_std).and.(mype.eq.0)) then
       ALLOCATE(XVSPFILE(ln1:ln2))

       ! Moment files for each species in a different file
       DO n=0,n_spec-1
          call get_unit_nr(XVSPFILE(n))

          OPEN(XVSPFILE(n), file=trim(diagdir)//'/xvsp_'//&
               &trim(spec(n)%name)//''//trim(file_extension),&
               form='unformatted', status=filestat, position=filepos)
       END DO
    endif
#ifdef WITHFUTILS
    if (write_h5) then
       !collect F0 and geometry information
       allocate(var(ll1:ll2,lm1:lm2,ln1:ln2))
       allocate(x_axis(0:nx0-1),Rpos(0:nx0-1),Zpos(0:nx0-1))
       allocate(f_unit(0:n_spec-1))
       do i = 0, nx0-1
          x_axis(i) = -0.5*lx+i*deli
       enddo
       var = 0.0
       pni=pn1
       IF (my_pez == (xvsp_iz/lk0)) THEN
          Blocal = geom%Bfield(pi1,pj1,xvsp_iz)
          Rpos = geom%R(pi1,xvsp_iz)+geom%dxdR(pi1,xvsp_iz)/geom%gii(pi1,pj1,xvsp_iz)*&
               &x_axis*rhostar*minor_r*Lref
          Zpos = geom%Z(pi1,xvsp_iz)+geom%dxdZ(pi1,xvsp_iz)/geom%gii(pi1,pj1,xvsp_iz)*&
               &x_axis*rhostar*minor_r*Lref
          do n=ln1,ln2
             if (pn0.gt.1) pni=n
             do m=lm1,lm2
                do l=ll1,ll2
                   var(l,m,n) = fm(pi1,pj1,xvsp_iz,l,m,pni)
                enddo
             End Do
          End Do
       END IF
       !\TODO: reduce communication, send just to write_pe
       call mpi_bcast(var, ll0*lm0*ln0, MPI_REAL_TYPE, (xvsp_iz/lk0), mpi_comm_z, ierr)
       call mpi_bcast(Blocal, 1, MPI_REAL_TYPE, (xvsp_iz/lk0), mpi_comm_z, ierr)
       call mpi_bcast(Rpos, 1, MPI_REAL_TYPE, (xvsp_iz/lk0), mpi_comm_z, ierr)
       call mpi_bcast(Zpos, 1, MPI_REAL_TYPE, (xvsp_iz/lk0), mpi_comm_z, ierr)

       allocate(fullarr(0:nv0-1,0:nw0-1,0:n_spec-1))
       fullarr = 0.

       Do n = ln1,ln2
          Do m = lm1,lm2
             CALL mpi_gather(var(ll1,m,n),ll0, MPI_REAL_TYPE,&
                  fullarr(0,m-lm1,n-ln1),ll0,MPI_REAL_TYPE,&
                  0, mpi_comm_v, ierr)
          Enddo
          CALL my_real_gather_to_0(fullarr(:,:,n-ln1),&
               &1,nv0*lm0,nv0*nw0,mpi_comm_w)
       Enddo
       CALL my_real_gather_to_0(fullarr(:,:,:),&
            &1,nv0*nw0*ln0,nv0*nw0*n_spec,mpi_comm_spec)
       deallocate(var)

       if (mype.eq.0) then
          allocate(fidxvsp_h5(0:n_spec-1))
          allocate(basexvsp_h5(0:n_spec-1))
          do n=0, n_spec-1
             !compute some units
             vp_unit = SQRT(2.0*spec(n)%temp*Tref*1E3*Qref/(spec(n)%mass*mref*1.6726219E-27))
             mu_unit = spec(n)%temp*Tref*1E3*Qref/Bref
             f0_unit = spec(n)%dens*nref*1E19*vp_unit**(-3)
             f_unit(n) = rhostar*minor_r*f0_unit
             basexvsp_h5(n) = '/xvsp_'//trim(spec(n)%name)
             call creatf(trim(diagdir)//trim(basexvsp_h5(n))//trim(file_extension)//'.h5', &
                  fidxvsp_h5(n), "Radial velocity space diagnostics", 'd')
             call creatg(fidxvsp_h5(n), trim(basexvsp_h5(n)))
             call creatg(fidxvsp_h5(n), trim(basexvsp_h5(n))//'/'//trim(xvsp_label(1)))
             rank = 0
             call creatd(fidxvsp_h5(n), rank, dims, trim(basexvsp_h5(n))//"/time", "time")
             !---some additional axes information---
             call creatg(fidxvsp_h5(n), trim(basexvsp_h5(n))//'/axes')
             call putarr(fidxvsp_h5(n), trim(basexvsp_h5(n))//"/axes/Rpos_m", Rpos)
             call putarr(fidxvsp_h5(n), trim(basexvsp_h5(n))//"/axes/Zpos_m", Zpos)
             call putarr(fidxvsp_h5(n), trim(basexvsp_h5(n))//"/axes/mugene", mu)
             call putarr(fidxvsp_h5(n), trim(basexvsp_h5(n))//"/axes/muweights", mu_weight)
             call putarr(fidxvsp_h5(n), trim(basexvsp_h5(n))//"/axes/mu_Am2", mu*mu_unit)
             call putarr(fidxvsp_h5(n), trim(basexvsp_h5(n))//"/axes/muweights_SI", mu_weight*mu_unit)
             call putarr(fidxvsp_h5(n), trim(basexvsp_h5(n))//"/axes/vpargene", vp)
             call putarr(fidxvsp_h5(n), trim(basexvsp_h5(n))//"/axes/vpar_m_s", vp*vp_unit)
             call putarr(fidxvsp_h5(n), trim(basexvsp_h5(n))//"/axes/xgene_rhoref", x_axis)
             !---some additional general information---
             call creatg(fidxvsp_h5(n), trim(basexvsp_h5(n))//'/general information')
             call attach(fidxvsp_h5(n), trim(basexvsp_h5(n))//"/general information/","Bref,T", Bref)
             call attach(fidxvsp_h5(n), trim(basexvsp_h5(n))//"/general information/","Lref,m", Lref)
             call attach(fidxvsp_h5(n), trim(basexvsp_h5(n))//"/general information/","Tref,eV", Tref*1E3)
             call attach(fidxvsp_h5(n), trim(basexvsp_h5(n))//"/general information/","mref,kg", mref*1.6726219E-27)
             call attach(fidxvsp_h5(n), trim(basexvsp_h5(n))//"/general information/","nref,m^-3", nref*1E19)
             call attach(fidxvsp_h5(n), trim(basexvsp_h5(n))//"/general information/","omegatorref,(rad,s)", omegatorref)
             !---some additional misc information---
             call creatg(fidxvsp_h5(n), trim(basexvsp_h5(n))//'/misc')
             call putarr(fidxvsp_h5(n), trim(basexvsp_h5(n))//"/misc/F0", fullarr(:,:,n)*f0_unit)
             call attach(fidxvsp_h5(n), trim(basexvsp_h5(n))//"/misc/","Blocal,T",Blocal*Bref)
             call attach(fidxvsp_h5(n), trim(basexvsp_h5(n))//"/misc/","xvsp_iy",xvsp_iy)
             call attach(fidxvsp_h5(n), trim(basexvsp_h5(n))//"/misc/","xvsp_iz",xvsp_iz)
          enddo
       end if
       deallocate(fullarr, Rpos, Zpos)
    endif
#endif

    CALL initialize_fourier_x_1d
    CALL initialize_fourier_y_1d

  End Subroutine Initialize_diag_xvsp

!!!******************************************************************!!!

  !>XVSP (radially resolved velocity space) diagnostic
  !!
  !!writes real-valued f(x,vpar,mu) at user-defined (y,z) positions;
  !!current main purpose is interfacing with CECE forward modelling
  !!therefore, delta_f and F0 are written in unnormalized units(!)
  !!to h5 output - together with various additional information
  Subroutine diag_xvsp
    Real:: var(li1:li2,ll1:ll2,lm1:lm2,ln1:ln2)
    Real:: fullarr(0:nx0-1, 0:nv0-1, 0:nw0-1, 0:n_spec-1)
    Integer :: i, j, l, m, n, ierr
    COMPLEX, DIMENSION(li1:li2, lj1:lj2):: temp1,temp2
    COMPLEX, DIMENSION(lj1:lj2) :: tempc
    REAL, DIMENSION(0:2*nj0-3) :: tempr
#ifdef WITHFUTILS
    character(len=FILENAME_MAX):: dset_name_xvsp
#endif

    PERFON('diag_xvsp')
    var = 0.0

    IF (my_pez == (xvsp_iz/lk0)) THEN
       do n=ln1,ln2
          do m=lm1,lm2
             do l=ll1,ll2
                temp1 = f_(li1:li2,lj1:lj2,xvsp_iz,l,m,n)
                do j=lj1,lj2
                   call to_real_x_1d(temp1(li1:li2,j),temp2(li1:li2,j))
                enddo
                do i=li1,li2
                   tempc = temp2(i,lj1:lj2)
                   call to_real_y_1d(tempc,tempr)
                   var(i,l,m,n) = tempr(xvsp_iy)
                enddo
             End Do
          End Do
       End Do
    END IF

    !\TODO: reduce communication / send to write_pe only
    call mpi_bcast(var, li0*ll0*lm0*ln0,&
         MPI_REAL_TYPE, (xvsp_iz/lk0), mpi_comm_z, ierr)

    fullarr = 0.

    Do n = ln1,ln2
       Do m = lm1,lm2
          Do l = ll1,ll2
             CALL mpi_gather(var(li1,l,m,n),li0, MPI_REAL_TYPE,&
                  fullarr(0,l-ll1,m-lm1,n-ln1),li0,MPI_REAL_TYPE,&
                  0, mpi_comm_x, ierr)
          Enddo
          CALL my_real_gather_to_0(fullarr(:,:,m-lm1,n-ln1),&
               &1,nx0*ll0,nx0*nv0,mpi_comm_v)
       Enddo
       CALL my_real_gather_to_0(fullarr(:,:,:,n-ln1),&
            &1,nx0*nv0*lm0,nx0*nv0*nw0,mpi_comm_w)
    Enddo
    CALL my_real_gather_to_0(fullarr(:,:,:,:),&
         &1,nx0*nv0*nw0*ln0,nx0*nv0*nw0*n_spec,mpi_comm_spec)

    if (mype==0) then
       do n =0, n_spec-1
          if(write_std) then
             Write(XVSPFILE(n)) time
             Write(XVSPFILE(n)) fullarr(:,:,:,n)
             call flush(XVSPFILE(n))
          end if
#ifdef WITHFUTILS
          if(write_h5) then
             call append(fidxvsp_h5(n), trim(basexvsp_h5(n))//"/time", time)
             call attach(fidxvsp_h5(n), trim(basexvsp_h5(n))//"/time", "n_steps", isnap_xvsp+1)
             write(dset_name_xvsp, "(A, '/', i10.10)") trim(basexvsp_h5(n))//"/"//&
                  &trim(xvsp_label(1)), isnap_xvsp
             call putarr(fidxvsp_h5(n), dset_name_xvsp, fullarr(:,:,:,n)*f_unit(n))
             call flushh5(fidxvsp_h5(n))
          end if
#endif
       enddo
    endif


#ifdef WITHFUTILS
    isnap_xvsp = isnap_xvsp+1
#endif
    PERFOFF

  End Subroutine diag_xvsp

  Subroutine finalize_diag_xvsp
    integer :: n
    if (mype.eq.0) then
       if (write_std) then
          do n=0, n_spec-1
             CLOSE(XVSPFILE(n))
          enddo
          DEALLOCATE(XVSPFILE)
       endif
#ifdef WITHFUTILS
       if (write_h5) then
          do n=0, n_spec-1
             call closef(fidxvsp_h5(n))
          enddo
          deallocate(fidxvsp_h5,basexvsp_h5)
          deallocate(f_unit)
       endif
#endif
    endif

    CALL finalize_fourier_x_1d
    CALL finalize_fourier_y_1d

  End Subroutine finalize_diag_xvsp


END MODULE diagnostics_vsp
