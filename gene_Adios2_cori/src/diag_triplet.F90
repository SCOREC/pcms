#include "redef.h"
#include "intrinsic_sizes.h"
#include "switches.h"

! This modules contains the diagnostic for nonlinear energy transfer
! within a wavenumber triplet. This is analagous to the nlt_pod
! diagnostic, but the nonlinear energy transfers recorded are
! restricted to within the triplet

MODULE diagnostics_triplet
  USE par_mod
  USE file_io, ONLY: get_unit_nr
  USE vel_space, ONLY: fm, mat_00
  USE geometry
  USE aux_fields, ONLY: calc_aux_fields
  USE gyro_average_ff_mod, ONLY: gyro_average_ff
  USE diagnostics_extended, ONLY: get_k_indices
  USE mpi
  USE communications
  use axpy

  IMPLICIT NONE

  PUBLIC :: initialize_diag_triplet, diag_triplet, &
     &finalize_diag_triplet

  PRIVATE
  !IO units for electrostatic/magnetic output respectively, and mode loading from pod or ev
  INTEGER :: tplt_es_info,tplt_em_info,mode_io
  INTEGER :: ke_es_info,ke_em_info,uniq_tplts
  !array which stores loaded modes(dist. fn at kx/ky)
  !last index is positive/negative kx. needed because of complex conjugates
  COMPLEX, DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: df_modes,cfgamma_g_decomp
  COMPLEX, DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: g_1_temp,cfgamma_chi_temp,chi_temp,chi_em_temp
  INTEGER :: record_count_g,record_count_chi
  !index of unique mode per kx-ky mode
  INTEGER, DIMENSION(:), ALLOCATABLE :: tplt_k_indices
  INTEGER, DIMENSION(:), ALLOCATABLE :: DFS_PER_WAVENUMBER
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: kxkyset

CONTAINS

  SUBROUTINE finalize_diag_triplet
    IF(n_dfs_tplts.GT.0) THEN
      DEALLOCATE(df_modes)
    ENDIF
    IF((split_em_tplt).AND.(n_fields.GT.1)) THEN
      DEALLOCATE(chi_em_temp)
    ENDIF
    DEALLOCATE(cfgamma_g_decomp)
    DEALLOCATE(cfgamma_chi_temp)
    DEALLOCATE(chi_temp)
    DEALLOCATE(g_1_temp)
    DEALLOCATE(tplt_k_indices)
    IF(mype==0) THEN
     CLOSE(tplt_es_info)
     CLOSE(ke_es_info)
     IF((n_fields.GT.1).AND.(split_em_tplt)) CLOSE(tplt_em_info)
     IF((n_fields.GT.1).AND.(split_em_tplt)) CLOSE(ke_em_info)
    ENDIF
  END SUBROUTINE finalize_diag_triplet
!!!******************************************************************!!!

  SUBROUTINE initialize_diag_triplet
    INTEGER :: i,j,index,i_1,count,kx,ky,loaded_pods
    CHARACTER(len=4) :: kx_ind,ky_ind
    INTEGER :: ierr,dummy_i,pod_IND
    COMPLEX, DIMENSION(:,:,:,:), ALLOCATABLE :: mode_in_tplt
    LOGICAL :: found, I_EXIST

    record_count_g = 1
    record_count_chi = 1

    ALLOCATE(tplt_k_indices(3*n_tplts))
    ALLOCATE(kxkyset(2,3*n_tplts))

    loaded_pods=0
    count=0
    kxkyset=0
    tplt_k_indices=0

    i_exist = .false.

    !Find the unique modes by kx-ky so we don't have to double load the same vectors
    !saves on memory, nessesary with very large numbers of tplts
    !*****************************************
    !First find unique modes

    DO i=1,n_tplts
      DO j=1,3
        index=3*i-3+j
        kx=ABS(kx_triplets(index))
        ky=ky_triplets(index)
        found=.FALSE.
        DO i_1=1,count
          IF ((kxkyset(1,i_1) .EQ. kx) .AND. (kxkyset(2,i_1) .EQ. ky)) THEN
            tplt_k_indices(index)=i_1
            found=.TRUE.
          ENDIF
        ENDDO
        IF (.NOT. found) THEN
          count=count+1
          kxkyset(1,count)=ABS(kx)
          kxkyset(2,count)=ky
          tplt_k_indices(index)=count
        ENDIF
      ENDDO
    ENDDO
    uniq_tplts=count

    !df_modes: This is the pod mode at only kx=kx_center
    ALLOCATE(DFS_PER_WAVENUMBER(uniq_tplts))
    IF(n_dfs_tplts.GT.0) THEN
      ALLOCATE(df_modes(lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2,uniq_tplts,n_dfs_tplts,2))
      df_modes=CMPLX(0.0,0.0)
    ENDIF
    ALLOCATE(cfgamma_g_decomp(lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2,uniq_tplts,n_dfs_tplts+1,2))

    IF((n_fields.GT.1).AND.(split_em_tplt)) THEN
      ALLOCATE(chi_em_temp(lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2,uniq_tplts,2))
    ENDIF
    ALLOCATE(chi_temp(lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2,uniq_tplts,2))
    ALLOCATE(cfgamma_chi_temp(lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2,uniq_tplts,2))
    ALLOCATE(g_1_temp(lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2,uniq_tplts,2))

    IF(mype==0) WRITE(*,*) "istep_triplet = ", istep_triplet

    !CALL get_unit_nr(tplt_es_info)
    !IF((n_fields.GT.1).AND.(split_em_tplt)) CALL get_unit_nr(tplt_em_info)

    !CALL get_unit_nr(ke_es_info)
    !IF((n_fields.GT.1).AND.(split_em_tplt)) CALL get_unit_nr(ke_em_info)

    tplt_es_info = 8911
    tplt_em_info = 8912
    ke_es_info = 8913
    ke_em_info = 8914

    !some computers can want RECL=2
    OPEN(unit=tplt_es_info,action='write',file=TRIM(diagdir)//'/triplet_es.dat',status='unknown',&
       form='unformatted',access='direct',RECL=8)
    IF((n_fields.GT.1).AND.(split_em_tplt)) OPEN(unit=tplt_em_info,action='write',file=TRIM(diagdir)&
       //'/triplet_em.dat',status='unknown',form='unformatted',access='direct',RECL=8)

    OPEN(unit=ke_es_info,action='write',file=TRIM(diagdir)//'/triplet_ke_es.dat',status='unknown',&
       form='unformatted',access='direct',RECL=8)
    IF((n_fields.GT.1).AND.(split_em_tplt)) OPEN(unit=ke_em_info,action='write',file=TRIM(diagdir)&
       //'/triplet_ke_em.dat',status='unknown',form='unformatted',access='direct',RECL=8)

    ALLOCATE(mode_in_tplt(0:nz0-1,0:nv0-1,0:nw0-1,0:n_spec-1))

    DO i=1,n_tplts
      DO j=1,3
        index=3*i-3+j
        IF(tplt_k_indices(index).GT.loaded_pods) THEN
          loaded_pods=loaded_pods+1
          CALL get_unit_nr(mode_io)
          POD_IND=0
          IF(mype==0) THEN
            WRITE(kx_ind,'(i4.4)') ABS(kx_triplets(index))
            WRITE(ky_ind,'(i4.4)') ky_triplets(index)
            INQUIRE(FILE=TRIM(SVD_df_file_path)//'/EV_df_ky'//ky_ind//'kx'//kx_ind,EXIST=I_EXIST)
          ENDIF
          IF(.NOT.I_EXIST) GOTO 1
          OPEN(unit=mode_io,file=TRIM(SVD_df_file_path)//'/EV_df_ky'//ky_ind//'kx'//kx_ind,&
             form='unformatted',status='unknown')

          DO pod_ind = 0, n_dfs_tplts - 1
            IF(mype==0) THEN
              READ(mode_io,err=1,end=1) dummy_i
              READ(mode_io,err=1,end=1) mode_in_tplt
            ENDIF
            CALL MPI_BCAST(mode_in_tplt,nz0*nv0*nw0*n_spec,&
               MPI_COMPLEX_TYPE,0,my_mpi_comm_world,ierr)
            !Now give each local array the right data
            df_modes(:,:,:,:,tplt_k_indices(index),pod_ind+1,1)=&
               mode_in_tplt(lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
            mode_in_tplt(:,:,:,:)=mode_in_tplt(nz0-1:0:-1,nv0-1:0:-1,:,:)
            df_modes(:,:,:,:,tplt_k_indices(index),pod_ind+1,2)=&
               mode_in_tplt(lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
            !Normalize modes
            CALL s_norm(df_modes(:,:,:,:,tplt_k_indices(index),pod_ind+1,1))
            CALL s_norm(df_modes(:,:,:,:,tplt_k_indices(index),pod_ind+1,2))
          ENDDO
1         CALL MPI_BCAST(POD_IND,1,MPI_INT,0,my_mpi_comm_world,ierr)
          DFS_PER_WAVENUMBER(tplt_k_indices(index))=POD_IND
          IF(mype==0) CLOSE(mode_io)
        ENDIF !IF GREATER THAN LOADED PODS
      END DO!j subindex within tplt
    END DO!i triplets
    DEALLOCATE(mode_in_tplt)

    PRINT*, DFS_PER_WAVENUMBER

    IF(mype==0) PRINT*, "end initilize triplet diagnostic"

  END SUBROUTINE initialize_diag_triplet


  SUBROUTINE construct_g_nlt(i,j,p,nlt_g_tfr,em)
    !modes cfgamma
    INTEGER, INTENT(in) :: i,j,p
    LOGICAL, INTENT(in) :: em
    REAL, INTENT(out) :: nlt_g_tfr
    COMPLEX :: nlt_g_tfr_4d(lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
    REAL :: kx0,ky0,kxp,kyp,kxpp,kypp,ckkp2
    INTEGER :: index,ip,jp,ipp,jpp,indexp,indexpp
    LOGICAL :: conjg_p,conjg_pp
    INTEGER :: k_sign,kp_sign,kpp_sign

    index=3*i-3+j

    IF (kx_triplets(index).LT.0) THEN
      kx0=kx(kx_triplets(index)+nx0)
    ELSE
      kx0=kx(kx_triplets(index))
    ENDIF
    ky0=kymin*ky_triplets(index)

    IF(j.EQ.1) THEN
      IF(kx_triplets(index+1).LT.0) THEN
        kxp=kx(kx_triplets(index+1)+nx0)
      ELSE
        kxp=kx(kx_triplets(index+1))
      ENDIF
      kyp=kymin*ky_triplets(index+1)
      indexp=index+1
      indexpp=index+2
    ELSEIF(j.EQ.2) THEN
      IF(kx_triplets(index-1).LT.0) THEN
        kxp=kx(kx_triplets(index-1)+nx0)
      ELSE
        kxp=kx(kx_triplets(index-1))
      ENDIF
      kyp=kymin*ky_triplets(index-1)
      indexp=index-1
      indexpp=index+1
    ELSE
      IF(kx_triplets(index-2).LT.0) THEN
        kxp=kx(kx_triplets(index-2)+nx0)
      ELSE
        kxp=kx(kx_triplets(index-2))
      ENDIF
      kyp=kymin*ky_triplets(index-2)
      indexp=index-2
      indexpp=index-1
    ENDIF

    kxpp=kx0-kxp
    kypp=ky0-kyp

    ckkp2=kxpp*ky0-kx0*kypp

    CALL get_k_indices(kxp,kyp,ip,jp,conjg_p)        !get indices corresponding to kp
    !must take complex conjugate for ky>0 modes (indicated by conjg_p)
    CALL get_k_indices(kxpp,kypp,ipp,jpp,conjg_pp)   !get indices corresponding to kpp

    k_sign=MERGE(1,2,kx0.GE.0)
    kp_sign=MERGE(1,2,ip.LE.nx0/2)
    kpp_sign=MERGE(1,2,ipp.LE.nx0/2)

    IF (.NOT.em) THEN

      nlt_g_tfr_4d(:,:,:,:)=ckkp2*cfgamma_g_decomp(:,:,:,:,tplt_k_indices(index),p,k_sign)&

         *MERGE(chi_temp(:,:,:,:,tplt_k_indices(indexpp),kpp_sign),&
         CONJG(chi_temp(:,:,:,:,tplt_k_indices(indexpp),kpp_sign)),.NOT.conjg_pp)&

         *MERGE(g_1_temp(:,:,:,:,tplt_k_indices(indexp),kp_sign),&
         CONJG(g_1_temp(:,:,:,:,tplt_k_indices(indexp),kp_sign)),.NOT.conjg_p)

    ELSE

      nlt_g_tfr_4d(:,:,:,:)=ckkp2*cfgamma_g_decomp(:,:,:,:,tplt_k_indices(index),p,k_sign)&

         *MERGE(chi_em_temp(:,:,:,:,tplt_k_indices(indexpp),kpp_sign),&
         CONJG(chi_em_temp(:,:,:,:,tplt_k_indices(indexpp),kpp_sign)),.NOT.conjg_pp)&

         *MERGE(g_1_temp(:,:,:,:,tplt_k_indices(indexp),kp_sign),&
         CONJG(g_1_temp(:,:,:,:,tplt_k_indices(indexp),kp_sign)),.NOT.conjg_p)

    ENDIF

    IF(j.EQ.1) THEN
     IF(kx_triplets(index+2).LT.0) THEN
      kxp=kx(kx_triplets(index+2)+nx0)
     ELSE
      kxp=kx(kx_triplets(index+2))
     ENDIF
     kyp=kymin*ky_triplets(index+2)
     indexp=index+2
     indexpp=index+1
    ELSEIF(j.EQ.2) THEN
     IF(kx_triplets(index+1).LT.0) THEN
      kxp=-kx(kx_triplets(index+1)+nx0)
     ELSE
      kxp=-kx(kx_triplets(index+1))
     ENDIF
     kyp=-kymin*ky_triplets(index+1)
     indexp=index+1
     indexpp=index-1
    ELSE
     IF(kx_triplets(index-1).LT.0) THEN
      kxp=-kx(kx_triplets(index-1)+nx0)
     ELSE
      kxp=-kx(kx_triplets(index-1))
     ENDIF
     kyp=-kymin*ky_triplets(index-1)
     indexp=index-1
     indexpp=index-2
    ENDIF

    kxpp=kx0-kxp
    kypp=ky0-kyp
    ckkp2=kxpp*ky0-kx0*kypp

    CALL get_k_indices(kxp,kyp,ip,jp,conjg_p)        !get indices corresponding to kp
    !must take complex conjugate for ky>0 modes (indicated by conjg_p)
    CALL get_k_indices(kxpp,kypp,ipp,jpp,conjg_pp)   !get indices corresponding to kpp

    k_sign=MERGE(1,2,kx0.GE.0)
    kp_sign=MERGE(1,2,ip.LE.nx0/2)
    kpp_sign=MERGE(1,2,ipp.LE.nx0/2)

    IF (.NOT.em) THEN

      nlt_g_tfr_4d(:,:,:,:)=nlt_g_tfr_4d+ckkp2*cfgamma_g_decomp(:,:,:,:,tplt_k_indices(index),p,k_sign)&

         *MERGE(chi_temp(:,:,:,:,tplt_k_indices(indexpp),kpp_sign),&
         CONJG(chi_temp(:,:,:,:,tplt_k_indices(indexpp),kpp_sign)),.NOT.conjg_pp)&

         *MERGE(g_1_temp(:,:,:,:,tplt_k_indices(indexp),kp_sign),&
         CONJG(g_1_temp(:,:,:,:,tplt_k_indices(indexp),kp_sign)),.NOT.conjg_p)

    ELSE

      nlt_g_tfr_4d(:,:,:,:)=nlt_g_tfr_4d+ckkp2*cfgamma_g_decomp(:,:,:,:,tplt_k_indices(index),p,k_sign)&

         *MERGE(chi_em_temp(:,:,:,:,tplt_k_indices(indexpp),kpp_sign),&
         CONJG(chi_em_temp(:,:,:,:,tplt_k_indices(indexpp),kpp_sign)),.NOT.conjg_pp)&

         *MERGE(g_1_temp(:,:,:,:,tplt_k_indices(indexp),kp_sign),&
         CONJG(g_1_temp(:,:,:,:,tplt_k_indices(indexp),kp_sign)),.NOT.conjg_p)

    ENDIF

    CALL nlt_integral(nlt_g_tfr_4d,nlt_g_tfr)

  END SUBROUTINE construct_g_nlt

  !calculates nrg tfr to kx0 ky0 from kp,kpp
  !part from transfer of g, part from chi, chi usually smaller
  SUBROUTINE construct_chi_nlt(i,j,nlt_chi_tfr,em)
    !modes cfgamma
    INTEGER, INTENT(in) :: i,j
    LOGICAL, INTENT(in) :: em
    REAL,  INTENT(out) :: nlt_chi_tfr
    COMPLEX :: nlt_chi_tfr_4d(lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
    REAL :: kx0,ky0,kxp,kyp,kxpp,kypp,ckkp1
    INTEGER :: index,ip,jp,ipp,jpp,indexp,indexpp
    LOGICAL :: conjg_p,conjg_pp
    INTEGER :: k_sign,kp_sign,kpp_sign

    index=3*i-3+j

    IF (kx_triplets(index).LT.0) THEN
     kx0=kx(kx_triplets(index)+nx0)
    ELSE
     kx0=kx(kx_triplets(index))
    ENDIF
    ky0=kymin*ky_triplets(index)

    IF(j.EQ.1) THEN
     IF(kx_triplets(index+1).LT.0) THEN
      kxp=kx(kx_triplets(index+1)+nx0)
     ELSE
      kxp=kx(kx_triplets(index+1))
     ENDIF
     kyp=kymin*ky_triplets(index+1)
     indexp=index+1
     indexpp=index+2
    ELSEIF(j.EQ.2) THEN
     IF(kx_triplets(index-1).LT.0) THEN
      kxp=kx(kx_triplets(index-1)+nx0)
     ELSE
      kxp=kx(kx_triplets(index-1))
     ENDIF
     kyp=kymin*ky_triplets(index-1)
     indexp=index-1
     indexpp=index+1
    ELSE
     IF(kx_triplets(index-2).LT.0) THEN
      kxp=kx(kx_triplets(index-2)+nx0)
     ELSE
      kxp=kx(kx_triplets(index-2))
     ENDIF
     kyp=kymin*ky_triplets(index-2)
     indexp=index-2
     indexpp=index-1
    ENDIF

    kxpp=kx0-kxp
    kypp=ky0-kyp
    ckkp1=kxp*ky0-kx0*kyp

    CALL get_k_indices(kxp,kyp,ip,jp,conjg_p)        !get indices corresponding to kp
    !must take complex conjugate for ky>0 modes (indicated by conjg_p)
    CALL get_k_indices(kxpp,kypp,ipp,jpp,conjg_pp)   !get indices corresponding to kpp

    k_sign=MERGE(1,2,kx0.GE.0)
    kp_sign=MERGE(1,2,ip.LE.nx0/2)
    kpp_sign=MERGE(1,2,ipp.LE.nx0/2)

    IF (.NOT.em) THEN

      nlt_chi_tfr_4d=ckkp1*cfgamma_chi_temp(:,:,:,:,tplt_k_indices(index),k_sign)&

         *MERGE(chi_temp(:,:,:,:,tplt_k_indices(indexp),kp_sign),&
         CONJG(chi_temp(:,:,:,:,tplt_k_indices(indexp),kp_sign)),.NOT.conjg_p)&

         *MERGE(g_1_temp(:,:,:,:,tplt_k_indices(indexpp),kpp_sign),&
         CONJG(g_1_temp(:,:,:,:,tplt_k_indices(indexpp),kpp_sign)),.NOT.conjg_pp)

    ELSE

      nlt_chi_tfr_4d=ckkp1*cfgamma_chi_temp(:,:,:,:,tplt_k_indices(index),k_sign)&

         *MERGE(chi_em_temp(:,:,:,:,tplt_k_indices(indexp),kp_sign),&
         CONJG(chi_em_temp(:,:,:,:,tplt_k_indices(indexp),kp_sign)),.NOT.conjg_p)&

         *MERGE(g_1_temp(:,:,:,:,tplt_k_indices(indexpp),kpp_sign),&
         CONJG(g_1_temp(:,:,:,:,tplt_k_indices(indexpp),kpp_sign)),.NOT.conjg_pp)

    ENDIF

    IF(j.EQ.1) THEN
     IF(kx_triplets(index+2).LT.0) THEN
      kxp=kx(kx_triplets(index+2)+nx0)
     ELSE
      kxp=kx(kx_triplets(index+2))
     ENDIF
     kyp=kymin*ky_triplets(index+2)
     indexp=index+2
     indexpp=index+1
    ELSEIF(j.EQ.2) THEN
     IF(kx_triplets(index+1).LT.0) THEN
      kxp=-kx(kx_triplets(index+1)+nx0)
     ELSE
      kxp=-kx(kx_triplets(index+1))
     ENDIF
     kyp=-kymin*ky_triplets(index+1)
     indexp=index+1
     indexpp=index-1
    ELSE
     IF(kx_triplets(index-1).LT.0) THEN
      kxp=-kx(kx_triplets(index-1)+nx0)
     ELSE
      kxp=-kx(kx_triplets(index-1))
     ENDIF
     kyp=-kymin*ky_triplets(index-1)
     indexp=index-1
     indexpp=index-2
    ENDIF

    kxpp=kx0-kxp
    kypp=ky0-kyp
    ckkp1=kxp*ky0-kx0*kyp

    CALL get_k_indices(kxp,kyp,ip,jp,conjg_p)        !get indices corresponding to kp
    !must take complex conjugate for ky>0 modes (indicated by conjg_p)
    CALL get_k_indices(kxpp,kypp,ipp,jpp,conjg_pp)   !get indices corresponding to kpp

    k_sign=MERGE(1,2,kx0.GE.0)
    kp_sign=MERGE(1,2,ip.LE.nx0/2)
    kpp_sign=MERGE(1,2,ipp.LE.nx0/2)

    IF (.NOT.em) THEN

      nlt_chi_tfr_4d=nlt_chi_tfr_4d+ckkp1*cfgamma_chi_temp(:,:,:,:,tplt_k_indices(index),k_sign)&

         *MERGE(chi_temp(:,:,:,:,tplt_k_indices(indexp),kp_sign),&
         CONJG(chi_temp(:,:,:,:,tplt_k_indices(indexp),kp_sign)),.NOT.conjg_p)&

         *MERGE(g_1_temp(:,:,:,:,tplt_k_indices(indexpp),kpp_sign),&
         CONJG(g_1_temp(:,:,:,:,tplt_k_indices(indexpp),kpp_sign)),.NOT.conjg_pp)

    ELSE

      nlt_chi_tfr_4d=nlt_chi_tfr_4d+ckkp1*cfgamma_chi_temp(:,:,:,:,tplt_k_indices(index),k_sign)&

         *MERGE(chi_em_temp(:,:,:,:,tplt_k_indices(indexp),kp_sign),&
         CONJG(chi_em_temp(:,:,:,:,tplt_k_indices(indexp),kp_sign)),.NOT.conjg_p)&

         *MERGE(g_1_temp(:,:,:,:,tplt_k_indices(indexpp),kpp_sign),&
         CONJG(g_1_temp(:,:,:,:,tplt_k_indices(indexpp),kpp_sign)),.NOT.conjg_pp)

    ENDIF

    CALL nlt_integral(nlt_chi_tfr_4d,nlt_chi_tfr)

  END SUBROUTINE construct_chi_nlt

  !>This subroutine performs the mat_00 integral, z integral and sum over species in order to reduce from 4D to 0D
  SUBROUTINE nlt_integral(v4d,v0d)
    IMPLICIT NONE
    COMPLEX, DIMENSION(lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2), INTENT(in) :: v4d
    REAL, INTENT(out) :: v0d
    REAL :: v0d_temp
    COMPLEX, DIMENSION(lk1:lk2) :: v1d
    INTEGER :: k,l,m,n,ierr

    v1d=CMPLX(0.0,0.0)
    DO n=ln1,ln2
     DO m=lm1,lm2
      DO l=ll1,ll2
       DO k=lk1,lk2
        CALL axpy_ij(1,mat_00(pi1,pj1,k,l,m),v4d(k,l,m,n),v1d(k))
       ENDDO
      ENDDO
     ENDDO
    ENDDO

    CALL my_complex_sum_vwspec(v1d,SIZE(v1d))

    v0d=REAL(SUM(v1d(:)*geom%jacobian(pi1,pj1,:)/(REAL(nz0)*geom%avg_jaco),1))

    CALL mpi_allreduce(v0d,v0d_temp,1,MPI_REAL_TYPE,MPI_SUM,mpi_comm_z,ierr)
    v0d=v0d_temp

  END SUBROUTINE nlt_integral

  SUBROUTINE diag_triplet
    IMPLICIT NONE
    !i always of n_tplts, j of 1-3, p of 1-n_modes,pp second loop over modes
    INTEGER :: i,j,p,pp,index
    INTEGER :: k,l,m,n
    INTEGER :: i0,i0_dummy,j0,index_uniq_tplts
    REAL :: kx0,ky0
    COMPLEX :: cfgamma_g(lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
    COMPLEX :: cfgamma_chi(lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
    COMPLEX :: chi(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
    COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: chi_em
    COMPLEX :: fields_1(li1:li2,lj1:lj2,lbz:ubz,1:n_fields)
    COMPLEX :: f_1(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2)
    COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: h_1_tplt,apar_bar
    COMPLEX :: g_temp(lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
    COMPLEX :: hn_t_tplt,hn_all_tplt(n_dfs_tplts)
    REAL :: nlt_g_tfr,nlt_g_tfr2,nlt_chi_tfr,nlt_chi_tfr2  !ES/EM g/chi transfer

    IF((n_fields.GT.1).AND.(split_em_tplt)) ALLOCATE(chi_em(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2))

    IF(mype==0) THEN
      WRITE(tplt_es_info, REC=record_count_g) time
      WRITE(ke_es_info, REC=record_count_chi) time
      IF((n_fields.GT.1).AND.(split_em_tplt)) THEN
        WRITE(tplt_em_info, REC=record_count_g) time
        WRITE(ke_em_info, REC=record_count_chi) time
      ENDIF
      record_count_g=record_count_g+1
      record_count_chi=record_count_chi+1
    ENDIF

    ALLOCATE(h_1_tplt(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2))
    h_1_tplt=CMPLX(0.0,0.0)

    IF(n_fields.GT.1) ALLOCATE(apar_bar(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2))
    fields_1=0.0
    chi=0.0
    IF((n_fields.GT.1).AND.(split_em_tplt)) chi_em=0.0

    CALL calc_aux_fields(g_1,fields_1,f_1,.TRUE.,h_1_tplt)

    !construct Chi
    DO n=ln1,ln2
      DO m=lm1,lm2
        DO l=ll1,ll2
          DO k=lk1,lk2
            CALL gyro_average_ff(fields_1(:,:,k,1),chi(:,:,k,l,m,n),k,m,n)
            IF(n_fields.GT.1) CALL gyro_average_ff(fields_1(:,:,k,2),&
               apar_bar(:,:,k,l,m,n),k,m,n)
          END DO
        END DO
      END DO
    END DO
    IF(n_fields.GT.1) THEN
      DO n=ln1,ln2
        DO l=ll1,ll2
          IF(split_em_tplt) THEN
            chi_em(:,:,lk1:lk2,l,lm1:lm2,n)=0.0-(2.0*spec(n)%temp/spec(n)%mass)**0.5*&
               vp(l)*apar_bar(:,:,lk1:lk2,l,lm1:lm2,n)
          ELSE
            chi(:,:,lk1:lk2,l,lm1:lm2,n)=chi(:,:,lk1:lk2,l,lm1:lm2,n)-&
               (2.0*spec(n)%temp/spec(n)%mass)**0.5*vp(l)*apar_bar(:,:,lk1:lk2,l,lm1:lm2,n)
          ENDIF
        END DO
      END DO
    END IF

    IF(ALLOCATED(apar_bar)) DEALLOCATE(apar_bar)

    !build a array with all the triplet's g_1s, decomposition of g_1s, and chis
    !as in uniq_tplts, so we only compute them once instead of for every tplt
    g_1_temp=0.0
    chi_temp=0.0
    if((n_fields.gt.1).and.split_em_tplt) chi_em_temp=0.0
    cfgamma_g_decomp=0.0
    cfgamma_chi_temp=0.0

    DO index_uniq_tplts=1,uniq_tplts
      i0=kxkyset(1,index_uniq_tplts)
      j0=kxkyset(2,index_uniq_tplts)
      ky0=kymin*j0
      g_1_temp(:,:,:,:,index_uniq_tplts,1)=g_1(i0,j0,:,:,:,:)
      chi_temp(:,:,:,:,index_uniq_tplts,1)=chi(i0,j0,:,:,:,:)
      IF((n_fields.GT.1).AND.(split_em_tplt)) THEN
        chi_em_temp(:,:,:,:,index_uniq_tplts,1)=chi_em(i0,j0,:,:,:,:)
      ENDIF
      IF(i0.NE.0) THEN
        g_1_temp(:,:,:,:,index_uniq_tplts,2)=g_1(nx0-i0,j0,:,:,:,:)
        chi_temp(:,:,:,:,index_uniq_tplts,2)=chi(nx0-i0,j0,:,:,:,:)
        IF((n_fields.GT.1).AND.(split_em_tplt)) THEN
          chi_em_temp(:,:,:,:,index_uniq_tplts,2)=chi_em(nx0-i0,j0,:,:,:,:)
        ENDIF
        !not populating negative part for ky=0
      ENDIF
      CALL get_cfgamma_g_ij(g_1_temp(:,:,:,:,index_uniq_tplts,1),cfgamma_g)
      CALL get_cfgamma_h_ij(h_1_tplt(i0,j0,:,:,:,:),cfgamma_chi)

      cfgamma_chi_temp(:,:,:,:,index_uniq_tplts,1)=cfgamma_chi-cfgamma_g

      IF(i0.NE.0) THEN
        CALL get_cfgamma_g_ij(g_1_temp(:,:,:,:,index_uniq_tplts,2),cfgamma_g)
        CALL get_cfgamma_h_ij(h_1_tplt(nx0-i0,j0,:,:,:,:),cfgamma_chi)
        cfgamma_chi_temp(:,:,:,:,index_uniq_tplts,2)=cfgamma_chi-cfgamma_g
      ELSE
        CALL get_cfgamma_g_ij(g_1_temp(:,:,:,:,index_uniq_tplts,2),cfgamma_g)
        CALL get_cfgamma_h_ij(h_1_tplt(0,j0,:,:,:,:),cfgamma_chi)
        cfgamma_chi_temp(:,:,:,:,index_uniq_tplts,2)=cfgamma_chi-cfgamma_g
      ENDIF

      hn_all_tplt=0
      DO p = 1, n_dfs_tplts + 1
        IF(p.LE.DFS_PER_WAVENUMBER(index_uniq_tplts)) THEN
          g_temp=df_modes(:,:,:,:,index_uniq_tplts,p,1)
          CALL spro_tplt(g_temp,g_1_temp(:,:,:,:,index_uniq_tplts,1),hn_t_tplt)
          g_temp=hn_t_tplt*g_temp
          hn_all_tplt(p)=hn_t_tplt !saving how much overlap there was
          CALL get_cfgamma_g_ij(g_temp,cfgamma_g_decomp(:,:,:,:,index_uniq_tplts,p,1))
        ENDIF
        IF(p.EQ.n_dfs_tplts+1) THEN
          g_temp=g_1(i0,j0,:,:,:,:)!g_1_temp(:,:,:,:,index_uniq_tplts,1)
          DO pp=1,DFS_PER_WAVENUMBER(index_uniq_tplts)
            g_temp=g_temp-&
               hn_all_tplt(pp)*df_modes(:,:,:,:,index_uniq_tplts,pp,1)
          END DO
          CALL get_cfgamma_g_ij(g_temp,cfgamma_g_decomp(:,:,:,:,index_uniq_tplts,n_dfs_tplts+1,1))
        ENDIF
      ENDDO

      hn_all_tplt=0
      DO p = 1,n_dfs_tplts + 1
        IF(p.LE.DFS_PER_WAVENUMBER(index_uniq_tplts)) THEN
          g_temp=df_modes(:,:,:,:,index_uniq_tplts,p,2)
          CALL spro_tplt(g_temp,g_1_temp(:,:,:,:,index_uniq_tplts,2),hn_t_tplt)
          g_temp=hn_t_tplt*g_temp
          hn_all_tplt(p)=hn_t_tplt !saving how much overlap there was
          CALL get_cfgamma_g_ij(g_temp,cfgamma_g_decomp(:,:,:,:,index_uniq_tplts,p,2))
        ENDIF
        IF(p.EQ.n_dfs_tplts+1) THEN
          g_temp=g_1_temp(:,:,:,:,index_uniq_tplts,2)
          DO pp=1,DFS_PER_WAVENUMBER(index_uniq_tplts)
            g_temp=g_temp-&
               hn_all_tplt(pp)*df_modes(:,:,:,:,index_uniq_tplts,pp,2)
          END DO
          CALL get_cfgamma_g_ij(g_temp,cfgamma_g_decomp(:,:,:,:,index_uniq_tplts,n_dfs_tplts+1,2))
        ENDIF
      ENDDO
    ENDDO

    !now precomputation is done, can go through all triplets

   DO i=1,n_tplts
     DO j=1,3
      index=3*i-3+j
      i0=kx_triplets(index)
      j0=ky_triplets(index)
      IF((kx_triplets(index)).LT.0) THEN
       i0_dummy=kx_triplets(index)+nx0
      ELSE
       i0_dummy=kx_triplets(index)
      ENDIF
      kx0=kx(i0_dummy)
      ky0=kymin*ky_triplets(index)


      DO p = 1, n_dfs_tplts + 1
        nlt_g_tfr=0.0
        nlt_g_tfr2=0.0
        nlt_chi_tfr=0.0
        nlt_chi_tfr2=0.0
        IF(p.LE.DFS_PER_WAVENUMBER(tplt_k_indices(index))) THEN
          CALL construct_g_nlt(i,j,p,nlt_g_tfr,.FALSE.)
          IF((n_fields.GT.1).AND.(split_em_tplt)) THEN !seperately calculate electromagnetic energy transfer
            CALL construct_g_nlt(i,j,p,nlt_g_tfr2,.TRUE.)
          ENDIF

        ELSEIF(p.EQ.n_dfs_tplts+1) THEN!Now get nlt for g_1-sum(hn*gn) i.e. the residual distribution function
          CALL construct_g_nlt(i,j,p,nlt_g_tfr,.FALSE.)
          IF((n_fields.GT.1).AND.(split_em_tplt)) THEN
            CALL construct_g_nlt(i,j,p,nlt_g_tfr2,.TRUE.)
          ENDIF
          CALL construct_chi_nlt(i,j,nlt_chi_tfr,.FALSE.)
          IF((n_fields.GT.1).AND.(split_em_tplt)) THEN
            CALL construct_chi_nlt(i,j,nlt_chi_tfr2,.TRUE.)
          ENDIF
        END IF

        nlt_g_tfr=(nlt_g_tfr*2.0)/(spec(0)%temp*spec(0)%dens)
        nlt_chi_tfr=(nlt_chi_tfr*2.0)/(spec(0)%temp*spec(0)%dens)
        nlt_g_tfr2=(nlt_g_tfr2*2.0)/(spec(0)%temp*spec(0)%dens)
        nlt_chi_tfr2=(nlt_chi_tfr2*2.0)/(spec(0)%temp*spec(0)%dens)

        IF(mype==0) THEN
          WRITE(tplt_es_info,REC=record_count_g) nlt_g_tfr
          IF(p.EQ.n_dfs_tplts+1) WRITE(ke_es_info,REC=record_count_chi) nlt_chi_tfr
          IF((n_fields.GT.1).AND.(split_em_tplt)) THEN
            WRITE(tplt_em_info,REC=record_count_g) nlt_g_tfr2
            IF(p.EQ.n_dfs_tplts+1) WRITE(ke_em_info,REC=record_count_chi) nlt_chi_tfr2
          ENDIF
          IF(p.EQ.n_dfs_tplts+1) record_count_chi=record_count_chi+1
          record_count_g=record_count_g+1
        ENDIF
      END DO !n_modes_tplt
    END DO !1-3 w/i triplet
  END DO !n_tplts

  DEALLOCATE(h_1_tplt)
END SUBROUTINE diag_triplet

SUBROUTINE get_cfgamma_g_ij(g_ij_in,cfgamma_g_ij)
  IMPLICIT NONE
  COMPLEX,DIMENSION(lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2),INTENT(in) :: g_ij_in
  COMPLEX, DIMENSION(lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2), INTENT(out):: cfgamma_g_ij
  INTEGER :: k,l,m,n
  !just get_cfgamma_g cut down for an individual i,j as we don't need to recalculate everything
  DO n=ln1,ln2
    DO m=lm1,lm2
      DO l=ll1,ll2
        DO k=lk1,lk2
          cfgamma_g_ij(k,l,m,n)=spec(n)%temp*spec(n)%dens/fm(pi1,pj1,k,l,m,pn1)*&
             &(CONJG(g_ij_in(k,l,m,n)))
        END DO
      END DO
    END DO
  END DO
END SUBROUTINE get_cfgamma_g_ij

!>This subroutine returns Cj/F0*(hj*)=cfgamma_h, same as above but different boundaries
SUBROUTINE get_cfgamma_h_ij(h_ij_in,cfgamma_h_ij)
  IMPLICIT NONE
  COMPLEX,DIMENSION(lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2),INTENT(in) :: h_ij_in
  COMPLEX, DIMENSION(lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2), INTENT(out):: cfgamma_h_ij
  INTEGER :: k,l,m,n
  DO n=ln1,ln2
    DO m=lm1,lm2
      DO l=ll1,ll2
        DO k=lk1,lk2
          cfgamma_h_ij(k,l,m,n)=spec(n)%temp*spec(n)%dens/fm(pi1,pj1,k,l,m,pn1)*&
             CONJG(h_ij_in(k,l,m,n))
        END DO
      END DO
    END DO
  END DO
END SUBROUTINE get_cfgamma_h_ij

SUBROUTINE spro_tplt(vec1,vec2,sp)
  COMPLEX, INTENT(in), DIMENSION(lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)  :: vec1,vec2
  COMPLEX, INTENT(out) :: sp
  COMPLEX :: sptemp(1)
  INTEGER :: k,l,m,n,ierr

  sptemp(1)=CMPLX(0.0,0.0)

  DO n=ln1,ln2
    DO k=lk1,lk2
      DO m=lm1,lm2
        DO l=ll1,ll2
          sptemp(1)=sptemp(1)+CONJG(vec1(k,l,m,n))*vec2(k,l,m,n)*mat_00(pi1,pj1,k,l,m)*&
             geom%jacobian(pi1gl,pj1,k)/fm(pi1,pj1,k,l,m,pn1)
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  sptemp(1)=sptemp(1)/(REAL(nz0)*geom%avg_jaco)

  CALL my_complex_sum_vwspec(sptemp,1)
  CALL mpi_allreduce(sptemp(1),sp,1,MPI_COMPLEX_TYPE,MPI_SUM,mpi_comm_z,ierr)

END SUBROUTINE spro_tplt

SUBROUTINE s_norm(vec1)
  COMPLEX, INTENT(inout), DIMENSION(lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2) :: vec1
  COMPLEX :: norm
  CALL spro_tplt(vec1,vec1,norm)
  vec1=vec1/SQRT(norm)

END SUBROUTINE s_norm

END MODULE diagnostics_triplet
