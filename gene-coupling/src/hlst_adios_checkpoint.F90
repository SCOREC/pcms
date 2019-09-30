!
#ifndef HLST_HAC_USE_REAL
!#define HLST_HAC_USE_REAL 0
#error
#endif
#ifndef HLST_HAC_USE_CMSK
!#define HLST_HAC_USE_CMSK 1
#error
#endif
#ifndef HLST_HAC_USE_ADIM
!#define HLST_HAC_USE_ADIM 6
#error
#endif
!
! ##################### BEGIN USER DOCUMENTATION #######################
!> Fortran module for ADIOS based checkpoint of MPI applications;
!! developed by the EFDA High Level Support Team within the HLST-ITM-ADIOS project.
!!
!! Checkpoint / restart functionality is provided via the following routines, each documented individually.
!! \li #hac_init
!! \li #hac_write
!! \li #hac_info
!! \li #hac_read
!! \li #hac_exit
!!
!! This module is specialized in high throughput checkpointing of a single \e global array on a parallel filesystem.
!!
!!
!! The user is expected to specify the supported numerical type and number of
!! dimensions at module compile time via three preprocessor symbols, here described:
!! \li If \c HLST_HAC_USE_REAL is 1, will use REAL (otherwise COMPLEX) as type.
!! \li If \c HLST_HAC_USE_CMSK is 1, will not specify the type's KIND, leaving it compiler set;
!! otherwise a value of 8 will be used.
!! \li Value of \c HLST_HAC_USE_ADIM is the dimension of the supported array; either of 2,4,6.
!!
!! Before usage, the #hac_init initialization subroutine shall be called.
!!
!! The module supports a single global array per checkpoint file.
!! So to checkpoint two global arrays, two separate checkpoints are necessary.
!! To checkpoint a local array, this module is not appropriate.
!!
!! The checkpoint file format is specific to the ADIOS library.
!! It is possible to use tools distributed with ADIOS to manipulate that file; however this should not be necessary,
!! because this module contains all that is needed to read the data back, e.g. in a converter.
!!
!! Written metadata are kept at a minimum: program's metadata are expected to be saved by the user separately.
!! Only information about the global array is written to the checkpoint file: nothing about the originating
!! application instance.
!! A given checkpoint file can be read back from a number of tasks independent from the number at write time.
!!
!! A multidimensional global array shape is defined by its lower and upper bounds indices.
!! Size on each dimension is equal to the difference of upper and lower bounds plus one.
!! No dimension can be zero; there is no hardcoded upper size limit.
!! A global array is represented by \e slice arrays allocated on each participating MPI task.
!! Each slice array can have arbitrary lower and upper bounds and size (accessible via the LBOUND, UBOUND, SHAPE intrinsics).
!!
!! When writing (with #hac_write), lower bounds must specified explicitly by the user; the local slice
!! size will be assumed to be the array size unless explicitly provided.
!! It is responsibility of the user to provide non-overlapping array slices.
!!
!! The global array's upper and lower bounds will be computed from the array slices bounds.
!! The upper bounds will be assumed to be the maximal upper bounds across the array slices;
!! the lower bounds will be assumed to be the minimal lower bounds across the array slices.
!!
!! Two conventions are supported for the global array's lower bounds: 1-based or 0-based.
!! Either of the two conventions can be used, but shall be consistent on all of the dimensions.
!! When writing, the convention will be determined according to the minimal lower bound encountered.
!! When reading, 0 will be assumed to be the base unless otherwise specified.
!!
!! When writing or reading, it is not necessary for all of the array slices to be of the same size,
!! but in no case can a dimension be zero.
!! So, reads (with #hac_read) can be performed on any portion of the global array, as long as indices within
!! its lower and upper bounds are specified.
!! Before reading a checkpoint file, a user may inquire (with #hac_info) into the dimensions and type of the global array.
!!
!! After usage, the #hac_exit initialization subroutine shall be called.
!!
!! No environment variable is accessed directly by this module, exception made for the subroutine #hac_test.
!!
!! Official documentation consists of a man page (this document) and the following example files:
!!
!! \li \e hac_gene.F90
!! \li \e hac_gemza.F90
!! \li \e hac_nemorb.F90
!! \li \e hac_dump.F90
!! \li \e Makefile
!!
!!
!! \warning ADIOS-1.5 needs to allocate a buffer at initialization time or
!!     before the first I/O; as a consequence, one cannot allocate a
!!     larger buffer at the time of e.g. the second I/O.
!!     <em> So in case the second I/O is expected to be larger than the first, one should allocate enough at the
!!     beginning explicitly</em>.
!!     Otherwise, expect failures.
!!     ADIOS-1.5 and ADIOS-1.6 emit a warning when I/Os exceed 335 MiB of (local) size; the I/O should perform correctly.
!!
!! \remark Only \c hac_* prefixed Fortran identifiers are made \c PUBLIC by this module.
!! \remark This module uses ADIOS, so many \c adios_ prefixed C symbols may appear in your linked application.
!! \remark As a clean programming practice (to avoid any potential name clash), please don't name any of your
!! \remark identifiers as prefixed with \c 'hac_' or \c 'adios_'.
!! \remark ADIOS-1.5 emits a warning when the local to-save I/O array size exceeds the hard-coded amount of 704643072 bytes.
!!
!! \note
!! \li This module requires ADIOS-1.5.
!! \li In the current setup, a file (with name specified by the user) and a likely named (with ".dir" suffixed)
!! \li directory will be created. They shall always be used together, and addressed by the file path.
!! \li Currently using hardcoded Lustre filesystem related values in A_NO, A_SC, A_NA.
!! \li Local I/O times are measured without \c FSYNC() semantics.
!! \li Handling only of certain simple errors like e.g. 'no file found' is supported; failure to comply with the
!! \li specifications may go undetected and behaviour undefined.
!! \li One precision \c KIND and either \c REAL or \c COMPLEX is supported.
!! \li Assuming \c INTEGER is 4 bytes. Will break otherwise.
!! \li If the \c OPTIONAL error status variables are omitted, errors are fatal and handled via the Fortran \c STOP statement.
!! \li Reading a checkpoint file written with a different type or KIND may go undetected; you may end up reading garbage.
!! \li Description of ADIOS' internal format is in the "ADIOS Developer's Manual" available on the ADIOS web site.
!!
!! \version 20140608
!!
!! \author Michele Martone, EFDA High Level Support Team
!!
! ###################### END USER DOCUMENTATION ########################
! ############## WHAT FOLLOWS IS INTERNAL DOCUMENTATION ################
!
! Throughout this module, UPPERCASE text has been used for: PARAMETER
! identifiers, language tokens, intrinsic functions; an underscore is
! used trailing optional arguments identifiers.
! When USE'ing modules, relevant module identifiers follow in an ONLY
! list.
!
! A list of possible extensions/developments/TODO:
! - The "irregular" distribution is bugged; as an effect of the bug, the total size is too small.
! - In order to support multiple types, may declare one array for each type each type DGA,ZGA,SGA,CGA,
!   but keeping only one the 'big' one.
! - Arbitrary lower array bounds.
! - Draconian tests on the values of LB,UB.
!
MODULE hlst_adios_checkpoint
!
! What follows in this module is all internals, except for the PUBLIC interface mentioned before.
!
#if (HLST_HAC_USE_REAL==1)
#define HAC_FIELD REAL
#else
#define HAC_FIELD COMPLEX
#endif
#if (HLST_HAC_USE_CMSK==1)
#define HAC_KSPEC
#else
#define HAC_KSPEC (8)
#endif
#define HAC_NTYPE HAC_FIELD HAC_KSPEC
#define HLST_HAC_COMPILING_WITH_GNU_STANDARD 1
#define HLST_HAC_GOTO_END_ON_ERROR(LABEL) CALL hac_herr; IF (ierr.NE.OKV_.OR.aerr.NE.OKV_) GOTO LABEL
#define HLST_HAC_HANDLE_IERR(IERR_,IERR) IF (PRESENT(IERR_)) IERR_ = IERR
#define HAC_SET_IF_PRESENT(VAR_,VAR) IF(PRESENT(VAR_))VAR=VAR_
#define HAC_RET_IF_PRESENT(VAR_,VAR) IF(PRESENT(VAR_))VAR_=VAR
#define HAC_MINMAXSUM(VAR,MINVAR,MAXVAR,SUMVAR) MAXVAR=MAX(VAR,MAXVAR);MINVAR=MIN(VAR,MINVAR);SUMVAR=SUMVAR+ VAR
#define HAC_PCTD(F,T) 100.0*(((F)-(T))/(T))
!
#if   (HLST_HAC_USE_ADIM==6)
#define HAC_DSPEC DIMENSION(:,:,:,:,:,:)
#elif (HLST_HAC_USE_ADIM==4)
#define HAC_DSPEC DIMENSION(:,:,:,:)
#elif (HLST_HAC_USE_ADIM==2)
#define HAC_DSPEC DIMENSION(:,:)
#else
#error "Unsupported HLST_HAC_USE_ADIM case."
#endif
#define HLST_HAC_EXTRA_METADATA 1
!
  USE adios_write_mod, ONLY: adios_write, adios_close, adios_open,&
          & adios_write, adios_group_size, adios_define_var, &
          & adios_finalize , adios_select_method, adios_declare_group,&
          & adios_allocate_buffer
  ! Note: in ADIOS-1.5 adios_init_noxml is defined only in C sources.
  USE adios_defs_mod, ONLY: adios_double, adios_long, adios_complex, &
        & adios_unknown, adios_real
  USE adios_read_mod, ONLY: adios_read_init_method,&
        & adios_read_open_file, adios_get_scalar, ADIOS_READ_METHOD_BP,&
        & adios_get_scalar, adios_read_close, adios_read_open_file,&
        & adios_selection_boundingbox, adios_schedule_read, &
        & adios_perform_reads, adios_read_close, &
        & adios_selection_delete, adios_selection_delete,&
        & adios_read_finalize_method
  USE mpi, ONLY: MPI_INTEGER8, MPI_INTEGER, MPI_DOUBLE_PRECISION, MPI_COMM_NULL, MPI_SUM, &
          &MPI_WTIME, MPI_COMM_WORLD, MPI_MAX, MPI_MIN ! MPI_Allreduce is not here !?
  USE par_in, only: adios_mpi_agg, adios_mpi_ost
  USE par_other, only: itime
  USE ISO_FORTRAN_ENV, ONLY: ERROR_UNIT, OUTPUT_UNIT
  IMPLICIT NONE
  PRIVATE
  SAVE
  REAL HAC_KSPEC,PARAMETER,PRIVATE :: SAMPLEVAR = 0
! ############## BEGIN OF PUBLIC DECLARATIONS ##########################
  PUBLIC :: hac_init,hac_exit,hac_read,hac_write !< Init/exit/read/write
  PUBLIC :: hac_info, hac_test !
  INTEGER,PARAMETER,PUBLIC :: HAC_KIND = KIND(SAMPLEVAR) !< Supported KIND of the module's numerical type.
  INTEGER,PARAMETER,PUBLIC :: HAC_QUIET = 0   !< Verbosity level for hac_init: quiet value.
  INTEGER,PARAMETER,PUBLIC :: HAC_VERBOSE = 1 !< Verbosity level for hac_init: verbose value.
  INTEGER,PARAMETER,PUBLIC :: HAC_VERBOSE_ALL = 1000000000 !< Verbosity level for hac_init: verbose value,
!! output from all tasks.
  INTEGER,PARAMETER,PUBLIC :: HAC_DEBUG = -1 !< Verbosity level for hac_init: debug value.
  INTEGER,PARAMETER,PUBLIC :: HAC_DEBUG_ALL = -1000000000 !< Verbosity level for hac_init: debug value, output from all tasks.
! ############## END OF PUBLIC DECLARATIONS ############################
  PRIVATE hac_adios_lib_exit !< Here to force Doxygen-1.8.4 to document correctly.
  INTEGER,PARAMETER :: AMFNL_ = 128  !< Maximal file name length
  INTEGER,PARAMETER :: MEVVL_ = 128  !< Maximal environment variable value length
  INTEGER,PARAMETER :: MPSL_ = 128  !< Maximal Printable String Length
  ! Note: adios_double_complex(=HAC_ADIOS_DOUBLE_COMPLEX) is missing or not PUBLIC !
  INTEGER,PARAMETER :: HAC_ADIOS_DOUBLE_COMPLEX = 11 !< (16) ! type code for double complex
  INTEGER,PARAMETER :: OFF_K = 8 !< KIND for long integers
  INTEGER,PARAMETER :: HAC_A_INT_T = adios_long !< ADIOS type for integers
  INTEGER,PARAMETER :: HAC_M_INT_T = MPI_INTEGER8 !< MPI type for long integers
  INTEGER,PARAMETER :: HAC_4_INT_T = MPI_INTEGER  !< MPI type for normal integers
  INTEGER,PARAMETER :: TIME_K = 8 !< Time REAL kind
  INTEGER,PARAMETER :: AI_K=8 !< INTEGER kind for ADIOS interface
  INTEGER,PARAMETER :: k_ = 1000, MB_ = k_ **2 !< KB,  MB (see SI units, or man 7 units)
  INTEGER,PARAMETER :: Ki_= 1024, MiB_= Ki_**2 !< KiB, MiB(see SI units, or man 7 units)
  INTEGER :: A_SS = 4*MiB_, A_BS = 4*MiB_ !< ADIOS  Stripe Size, Block Size
  INTEGER,PARAMETER :: SUL_ = 5
#if 1
  INTEGER,PARAMETER :: XXKC = Ki_
  CHARACTER(LEN=*),PARAMETER :: XBS  =  "iB"
#else
  INTEGER,PARAMETER :: XXKC = K_
  CHARACTER(LEN=*),PARAMETER :: XBS  =  "B"
#endif
  INTEGER,PARAMETER :: OKV_ = 0 !< No error value
  INTEGER,PARAMETER :: ERRV_G = -1 !< General error value
  INTEGER,PARAMETER :: ERRV_M = -2 !< Dimensions mismatch value
  INTEGER,PARAMETER :: ERRV_F = -3 !< Invoking after post-finalization error value
  INTEGER,PARAMETER :: ERRV_I = -4 !< Multiple finalization error value
  INTEGER,PARAMETER :: ERRV_D = -5 !< Dimensions Limits related error (MAXLD trespassed)
  INTEGER,PARAMETER :: MAXLD = k_*k_*k_ !< Maximum Local Dimension
  INTEGER,PARAMETER :: RR = 0 !< Root Rank
  !INTEGER :: A_NA = 1024 !< ADIOS aggregators (output subfiles)
  !INTEGER :: A_NO = 128, A_SC=16 !< Lustre OST's and stripe count
  INTEGER :: A_NA = 1024 !< ADIOS aggregators (output subfiles)
  INTEGER :: A_NO = 128, A_SC=16 !< Lustre OST's and stripe count
  INTEGER :: eanf = 0 !< Errors Are Non Fatal (if 1)
  CHARACTER(LEN=*),PARAMETER :: A_GN = "all" !< ADIOS group name
  CHARACTER(LEN=*),PARAMETER :: DAC = "dac"  !< dimensions array count
  INTEGER,PARAMETER :: EU = ERROR_UNIT, OU=OUTPUT_UNIT !< output units
  INTEGER,PARAMETER :: MDIMS_ = 7, VIL = 3, ADT_ = 3 !< Maximum Dimensions, Variable Identifier Length, Array Dimension Types
  INTEGER,PARAMETER :: OFF_ = 1, GLO_ = 2, LOC_ = 3 !< Offset, Global, Local
  INTEGER,PARAMETER :: MDS_ = 100 !< Max Dumpable local Size (for debug)
  INTEGER :: DODEBUG_ = 0 !< Quiet mode if DODEBUG_.LT.1
  CHARACTER(LEN=VIL),PARAMETER :: AID = "GAS" !< Global Array Slice
  INTEGER(KIND=OFF_K),PARAMETER :: NDIMS = HLST_HAC_USE_ADIM !< dimensions of the I/O array
  INTEGER(KIND=OFF_K) :: mis !< Meta Information Size
  CHARACTER(LEN=*),PARAMETER :: BNR = "hlst_adios_checkpoint: " !< Banner
  CHARACTER(LEN=*),PARAMETER :: ARMP = "verbose=3"!read init method parms
  INTEGER :: verbose = 0 !< Quiet mode if .LT.1
  INTEGER :: atss, itss ! Array/Index type storage size
  REAL(KIND=TIME_K) :: wt, rt, yt !< Write time, read time, synchronization time
  INTEGER :: a_buf_mib = 0 ! ADIOS buffer size in MiB (in ADIOS sources adios.c it's multiplied by 1024**2
!! --- the PDF manual documents it wrong)
  INTEGER(KIND=AI_K) :: agid, avid !< ADIOS Group ID, ADIOS Variable ID
  CHARACTER(LEN=AMFNL_) :: fhn, aps, awm!< First host name, ADIOS write parameter string, ADIOS write method
  INTEGER :: acomm = MPI_COMM_NULL !< Communicator for ADIOS operations
  INTEGER :: asize = -1 !< ADIOS communicator size
  INTEGER :: arank !< Rank within the ADIOS communicator
  INTEGER :: aerr = 0 !< ADIOS routines error variable
  INTEGER :: ierr = 0 !< Generic error variable
  INTEGER :: fime = 0 !< Fill Memory (A Flush Technique)
  INTEGER :: hac_a_num_t = adios_unknown !< ADIOS array numerical type
  LOGICAL :: initialized = .FALSE., finalized = .FALSE. !< helper vars
  LOGICAL :: initialized_vars = .FALSE. !< helper vars
  LOGICAL :: abao = .FALSE. !< Allocate Buffer At Once
  LOGICAL :: amroot = .FALSE. !< Is This Rank Root ?
!
 CONTAINS
!
  SUBROUTINE hac_herr
!> Error handler. It performs collective operations.
   IMPLICIT NONE
   ! CHARACTER(LEN=AMFNL_) :: aes = ''
   CHARACTER(LEN=*),PARAMETER :: BBC = " *** hlst_adios_checkpoint *** "
   CHARACTER(LEN=*),PARAMETER :: ENR = " *** "
   INTEGER :: gierr = 0, ferr = 0
!
   CALL MPI_Allreduce(ierr, gierr, 1, MPI_INTEGER, MPI_SUM, acomm, ierr)
   ierr = gierr
   CALL MPI_Allreduce(aerr, gierr, 1, MPI_INTEGER, MPI_SUM, acomm, ierr)
   aerr = gierr
!
   IF ( amroot ) THEN
    SELECT CASE (ierr)
     CASE (ERRV_M)
     WRITE(EU,*)BBC,'Read Data Dimensions Mismatch !',ENR
     CASE (ERRV_D)
     WRITE(EU,*)BBC,'Read Data Dimensions Off-Limits !',ENR
     CASE (ERRV_F)
     WRITE(EU,*)BBC,'More than one finalization is not allowed !',ENR
     CASE (ERRV_I)
     WRITE(EU,*)BBC,'Initializing after finalization not allowed !',ENR
     !CASE (errv_a)
     !WRITE(EU,*)'Generic Error!'
     ! MPI Error ?!
    END SELECT
   END IF
   IF (ierr .NE. OKV_) THEN
    IF ( amroot ) THEN
     WRITE(EU,*)BBC,'Detected Error Code: ',ierr,ENR
    END IF
    ferr = ierr
   END IF
!
   IF (aerr .LT. OKV_) THEN
    IF ( amroot ) THEN
     ! CALL adios_errmsg(aes, aesl) ! This is not part of the interface.
     WRITE(EU,*) BBC,'Detected an ADIOS Error: ',aerr,ENR
     ! The following are not PUBLIC members of adios_defs_mod, with no adios_*
     ! prefixed symbol (as of ADIOS-1.5).
     SELECT CASE (aerr)
      CASE (-1)
       WRITE(EU,*) BBC,'Presumably a Memory Allocation Error!',ENR
      CASE (-2)
       WRITE(EU,*) BBC,'Presumably a File Open Error!',ENR
      CASE (-3)
       WRITE(EU,*) BBC,'Presumably a File Not Found Error!',ENR
      CASE (-24)
       WRITE(EU,*) BBC,'Presumably a Too Many Files Error!',ENR
      CASE DEFAULT
     END SELECT
    END IF
    ferr = aerr
   END IF
!
   IF ( eanf .EQ. 0 .AND. ferr .NE. 0 ) THEN
    IF ( amroot ) WRITE(EU,*)eanf
    IF ( amroot ) WRITE(EU,*)BBC,'Will abort job with code',ferr,'.',ENR
    CALL MPI_Abort(acomm, ferr, ierr)
   END IF
  END SUBROUTINE hac_herr
!
! Reset error
  SUBROUTINE hac_rerr
   aerr = OKV_
   ierr = OKV_
  END SUBROUTINE hac_rerr
!
  PURE SUBROUTINE hac_gadn ( vid, adi, adj)
!> ADIOS (array) Dimension(s) Name (string)
   IMPLICIT NONE
   CHARACTER(LEN=*),INTENT(INOUT) :: vid
   CHARACTER(LEN=VIL+1) :: vidt
   INTEGER,INTENT(IN) :: adi, adj
   INTEGER :: i,u,l
!
   IF (adi.LT.0) THEN
    l = 1
    u = - adi
   ELSE
    l = adi
    u = adi
   END IF
   WRITE (vid,'(a)') TRIM("")
   DO i = l, u
     IF (adj==OFF_) WRITE(vidt,'(a,i0)') TRIM("od"),i
     IF (adj==GLO_) WRITE(vidt,'(a,i0)') TRIM("gd"),i
     IF (adj==LOC_) WRITE(vidt,'(a,i0)') TRIM("ld"),i
     vid = TRIM(vid) // TRIM(vidt)
     IF (i.LT.u) vid = TRIM(vid) // TRIM(',')
   END DO
  END SUBROUTINE hac_gadn


  PURE SUBROUTINE dim_name ( vid, adi, adj)
   IMPLICIT NONE
   CHARACTER(LEN=*),INTENT(INOUT) :: vid
   INTEGER,INTENT(IN) :: adi, adj
   CHARACTER(len=10) :: tmp

     vid = ""

     if (adi==1) tmp = "x"
     if (adi==2) tmp = "ky"
     if (adi==3) tmp = "z"
     if (adi==4) tmp = "vz"
     if (adi==5) tmp = "vx"
     if (adi==6) tmp = "species"

     IF (adj==OFF_) WRITE(vid,'(a,a)') "global_offset_", TRIM(tmp)
     IF (adj==GLO_) WRITE(vid,'(a,a)') "global_num_", TRIM(tmp)
     IF (adj==LOC_) WRITE(vid,'(a,a)') "num_", TRIM(tmp)


  END SUBROUTINE dim_name

!
!< Define the necessary ADIOS variables.
!! Cannot be called before hac_init.
  SUBROUTINE hac_defv( A_NA_,A_NO_,A_SC_,A_SS_,A_BS_,ierr_ )
   IMPLICIT NONE
!
   INTEGER,INTENT(IN),OPTIONAL :: A_NA_,A_NO_,A_SC_,A_SS_,A_BS_,ierr_
!
   CHARACTER(LEN=VIL) :: vid
   INTEGER :: i,j
   CHARACTER(LEN=MDIMS_*(VIL+1)-1) :: lds,gds,ods
!
   IF ( initialized_vars .EQV. .FALSE. ) THEN
!
    HAC_SET_IF_PRESENT(A_NA_,A_NA)
    HAC_SET_IF_PRESENT(A_NO_,A_NO)
    HAC_SET_IF_PRESENT(A_SS_,A_SS)
    HAC_SET_IF_PRESENT(A_SC_,A_SC)
    HAC_SET_IF_PRESENT(A_BS_,A_BS)
!
    CALL adios_declare_group (agid, A_GN,"time0",1, aerr)
    aerr = 0 ! FIXME. This is here because of ADIOS non respecting its own specifications.
    CALL hac_herr
    IF (.FALSE.) THEN
     A_NA = 16
     A_NO = 64
     A_SC = 16
     A_SS = 4*MiB_
     A_BS = 4*MiB_
     WRITE(aps,'(a,i0,a,i0,a,i0,a,i0,a,i0)') &
     &TRIM("num_aggregators="),A_NA,&
     &TRIM(";num_ost="),A_NO,&
     &TRIM(",stripe_count="),A_SC,&
     &TRIM(",stripe_size="),A_SS,&
     &TRIM(",block_size="),A_BS
     WRITE(awm,'(a)') "MPI_AGGREGATE"
    ELSE
     WRITE(aps,'(a,i0,a,i0)') &
     &TRIM("num_aggregators="),A_NA,&
     &TRIM(";num_ost="),A_NO
     WRITE(awm,'(a)') "MPI_AMR"
    END IF
     !IF (DODEBUG_.NE.0.AND.amroot) THEN
     IF (amroot) THEN
      WRITE(OU,'(a,a,a)') BNR,'Using ADIOS method: ',TRIM(awm)
      WRITE(OU,'(a,a,a)') BNR,'Using ADIOS params: ',TRIM(aps)
     END IF
     CALL adios_select_method (agid,TRIM(awm),TRIM(aps),"",aerr)
     aerr = 0 ! FIXME. This is here because of ADIOS non respecting its own specifications.

     CALL adios_define_var (agid,DAC,"",HAC_A_INT_T,"","","",avid)

     DO j = 1, ADT_
      DO i = 1, NDIMS
       CALL hac_gadn(vid,i,j)
       CALL adios_define_var (agid, vid,"",HAC_A_INT_T,"","","",avid)

      END DO
     END DO
     CALL hac_herr
     CALL hac_gadn(ods,INT(-NDIMS),OFF_)
     CALL hac_gadn(gds,INT(-NDIMS),GLO_)
     CALL hac_gadn(lds,INT(-NDIMS),LOC_)
#if HLST_HAC_USE_REAL
     IF (KIND(SAMPLEVAR).EQ.4) THEN
      hac_a_num_t = adios_real
     ELSE
      hac_a_num_t = adios_double
     ENDIF
#else
     IF (KIND(SAMPLEVAR).EQ.4) THEN
      hac_a_num_t = adios_complex
     ELSE
      hac_a_num_t = HAC_ADIOS_DOUBLE_COMPLEX
     ENDIF
#endif
     !
     IF ( DODEBUG_ .GT. arank ) THEN
      WRITE (OU,'(a,a,i0,a,a,a,a,a,a,a,a,a)') BNR,('on rank '),&
      & arank," aid:'",aid,"' lds:'",TRIM(lds),"' gds:'",TRIM(gds)&
      &,"' ods:'",TRIM(ods),"'"
     ENDIF
     !

     CALL adios_define_var (agid, AID,"",hac_a_num_t,lds,gds,ods,avid)

#if (HLST_HAC_EXTRA_METADATA==1)
     CALL adios_define_var (agid, "st","",adios_double,"","","",avid)
     CALL adios_define_var (agid, "ts","",adios_double,"","","",avid)
#endif
     CALL hac_herr
     initialized_vars = .TRUE.
     HAC_SET_IF_PRESENT(ierr_,ierr)
    END IF
!
  END SUBROUTINE hac_defv
!
!> Initializes module hlst_adios_checkpoint.
!! It must be called once, after MPI_init.
!!
!! In verbose mode, supported arrays type/dimensions information will be printed out.
!!
  SUBROUTINE hac_init ( acomm_, verbose_, ierr_ )
   USE mpi, ONLY: MPI_Comm_dup
   IMPLICIT NONE
!
INTEGER,OPTIONAL,INTENT(IN) :: acomm_ !< Checkpoint specific MPI communicator. If not provided, MPI_COMM_WORLD will be used.
   INTEGER,OPTIONAL,INTENT(IN) :: verbose_ !< Requested verbosity level. Default is not verbose (0).
   INTEGER,OPTIONAL,INTENT(OUT) :: ierr_ !< Error status variable. If omitted, errors will be fatal.
!
   IF ( finalized .EQV. .TRUE. ) THEN
    ierr = ERRV_I
   END IF
!
   IF ( initialized .EQV. .FALSE. ) THEN
    IF ( PRESENT ( verbose_ ) ) THEN
     verbose = MAX(0, ABS(verbose_))
     IF ( verbose_ .LT. 0 ) DODEBUG_ = 1
    END IF
!
    IF ( PRESENT ( acomm_ ) ) THEN
     ! CALL MPI_Comm_dup (acomm_, acomm, ierr)
     acomm = acomm_
    ELSE
     CALL MPI_Comm_dup (MPI_COMM_WORLD, acomm, ierr)
     ! CALL hac_herr ! not allowed yet
    END IF
    CALL MPI_Comm_rank (acomm, arank, ierr)
    IF ( arank .EQ. RR ) amroot = .TRUE.
    CALL hac_herr
    CALL MPI_Comm_size (acomm, asize, ierr)
    CALL hac_herr

    !CALL adios_init_noxml (acomm, aerr)
    !CALL hac_herr
    !a_buf_mib = 0 ! Will be allocated at the next call of adios_allocate_buffer
    !call adios_set_max_buffer_size(a_buf_mib)
    !CALL hac_rerr ! FIXME: resetting aerr. This is tolerated by ADIOS-1.5.
    !CALL hac_herr

    ! CALL hac_defv
    itss = hac_isz()
    atss = hac_asz()
    mis = itss * (1 + NDIMS * ADT_)
    eanf = 0 ! errors are fatal
    fime = 0 !
#if HLST_HAC_COMPILING_WITH_GNU_STANDARD
!#ifdef __GFORTRAN__
     ! It seems like ifort stands this. Compatibility mode by default !?
     IF ( verbose .GT. 0 ) THEN
      CALL HOSTNM(fhn)
      fhn = TRIM (fhn)
#if (HLST_HAC_USE_REAL==1)
      IF (amroot) WRITE(OU,'(a,a,i0,a,i0,a)')BNR,'Using REAL*',atss ,&
#else
      IF (amroot) WRITE(OU,'(a,a,i0,a,i0,a)')BNR,'Using COMPLEX*',atss,&
#endif
              &' with ', HLST_HAC_USE_ADIM,' dimensions as array type.'
      IF (amroot) WRITE(OU,'(a,i0,a,a,a)') BNR//'First task of ',asize,&
              &' on host ',&
              &TRIM(fhn),'.'
      IF (arank==1) WRITE(OU,'(a,a,a,a)') BNR,'Second task on host ',&
              &TRIM(fhn),'.'
      IF (asize.GT.1.AND.arank.EQ.asize-1) WRITE(OU,'(a,a,i0,a,a,a)')&
              & BNR,'Last task of ',asize,' on host  ',TRIM(fhn),'.'
      IF ( amroot .AND. DODEBUG_ .NE. 0 ) WRITE(OU,'(a,a,i0,a,i0)')BNR,&
              & 'verbosity level: ',verbose_,', debug level:  ',dodebug_
     END IF
!#endif
#endif
    initialized = .TRUE.
    finalized = .FALSE.
    ierr = OKV_
    aerr = OKV_
   END IF
   HLST_HAC_HANDLE_IERR(ierr_,ierr)
  END SUBROUTINE hac_init
!
!!> \internal
  SUBROUTINE hac_adios_lib_exit
   IF ( finalized .EQV. .FALSE. ) THEN
    CALL adios_finalize (arank, ierr)
    CALL hac_herr
    finalized = .TRUE.
   ELSE
    ierr = ERRV_F
   END IF
   CALL hac_herr
  END SUBROUTINE hac_adios_lib_exit
!
!> Check the Dimensions Array offset
  INTEGER FUNCTION check_dda(LBA)
   IMPLICIT NONE
   INTEGER :: ierr = OKV_
   INTEGER(KIND=OFF_K),INTENT(IN),DIMENSION(:) :: LBA
!
   IF ( PRODUCT (INT8(LBA(1:NDIMS)) ) .EQ. 0) THEN
    ierr = ERRV_M
   END IF
!
   IF ( MINVAL(INT8(LBA(1:NDIMS)) ) .LT. 0) THEN
    ierr = ERRV_M
   END IF
!
   IF ( MINVAL (LBA(1:NDIMS) ) .GT. 1) THEN
    ierr = ERRV_M
   END IF
!
   IF ( MAXVAL (LBA(1:NDIMS) ) .GT. MAXLD ) THEN
    ierr = ERRV_D
   END IF
!
   check_dda = ierr
  END FUNCTION
!
!> Check the global Array's Dimensions Array
  SUBROUTINE check_ada(GASS, ada)
   IMPLICIT NONE
   INTEGER,INTENT(IN) :: GASS(:) !< Global Array Slice Shape
   INTEGER(KIND=OFF_K),INTENT(IN) :: ada(MDIMS_,ADT_) !< Array's Dimensions Array
!
   INTEGER :: lerr = OKV_
!
   lerr = check_dda(ada(1:NDIMS, OFF_))
   IF ( PRODUCT (ada(1:NDIMS, LOC_) ) .EQ. 0) THEN
    lerr = ERRV_M
   END IF
   IF ( PRODUCT (ada(1:NDIMS, GLO_) ) .EQ. 0) THEN
    lerr = ERRV_M
   END IF
   IF (PRODUCT (ada(1:NDIMS,LOC_)) .GT. PRODUCT(ada(1:NDIMS,GLO_))) THEN
    lerr = ERRV_M
   END IF
   IF ( SUM (GASS(1:NDIMS) - ada(1:NDIMS, LOC_) ) .NE. 0) THEN
    ierr = ERRV_M
   END IF
   ierr = lerr
  END SUBROUTINE check_ada
!
!> Integer type Size
INTEGER FUNCTION hac_isz()
  IMPLICIT NONE
  INTEGER :: sz
!
#if HLST_HAC_COMPILING_WITH_GNU_STANDARD
     sz = INT(SIZEOF(NDIMS))
#else
     ! This Requires Fortran-2008
     sz = STORAGE_SIZE(NDIMS)/ 8
#endif
   hac_isz = sz
  END FUNCTION hac_isz
!
!> Array numerical type Size
INTEGER FUNCTION hac_asz()
  IMPLICIT NONE
!
  HAC_NTYPE,PARAMETER :: XYE = 0.0 !< A sample value.
!
  INTEGER :: sz
!
#if HLST_HAC_COMPILING_WITH_GNU_STANDARD
     sz = INT(SIZEOF(XYE))
#else
     ! This Requires Fortran-2008
     sz = STORAGE_SIZE(XYE)  / 8
#endif
   hac_asz = sz
  END FUNCTION hac_asz
!
!> Numerical Type Character
CHARACTER FUNCTION hac_ntc()
  IMPLICIT NONE
  CHARACTER :: ntc = '?' !< Numerical Type Character.
!
   SELECT CASE (hac_asz())
#if (HLST_HAC_USE_REAL==1)
    CASE (4)
     ntc = 'S'
    CASE (8)
     ntc = 'D'
#else
    CASE (8)
     ntc = 'C'
    CASE (16)
     ntc = 'Z'
#endif
   END SELECT
   hac_ntc = ntc
  END FUNCTION hac_ntc
!
!> Reads header information from the checkpoint file.
!! Can be called repeatedly to get information out of a checkpoint file.
  SUBROUTINE hac_info(fname, ierr_, GADA_, ntc_, ts_, st_)
   IMPLICIT NONE
!
   INTEGER,INTENT(OUT),OPTIONAL :: GADA_(:) !< Global Array's Dimensions Array.
   !!The user can use a subset of these dimensions to ALLOCATE the right sized (or, SHAPE'd) local array.
   CHARACTER(LEN=*),INTENT(IN) :: fname !< Path to the checkpoint file. Any length is allowed.
   INTEGER,OPTIONAL,INTENT(OUT) :: ierr_ !< Error status variable. If omitted, errors will be fatal.
   CHARACTER,OPTIONAL,INTENT(OUT) :: ntc_ !< Numerical type character. Will be either S,D,C,Z.
#if (HLST_HAC_EXTRA_METADATA==1)
   REAL(KIND=8),OPTIONAL,INTENT(OUT) :: ts_ !< Meta-data scalar: Time Step. Can be useful to certain users.
   REAL(KIND=8),OPTIONAL,INTENT(OUT) :: st_ !< Meta-data scalar: Simulation time. Can be useful to certain users.
#endif
!
   CHARACTER :: ntc !< Numerical Type Character.
   INTEGER(KIND=OFF_K) :: ada(MDIMS_,ADT_) !< Array's Dimensions Array
   INTEGER :: i,j
   INTEGER(KIND=AI_K) :: fh ! File Handle
   CHARACTER(LEN=VIL) :: vid ! Variable Identifier
   INTEGER(KIND=OFF_K) :: ndims_v ! Number of dimensions (ADIOS) variable ID
#if (HLST_HAC_EXTRA_METADATA==1)
   REAL(KIND=8) :: ts = 0.0, st = 0.0 ! DOUBLE PRECISION is not supported by ADIOS
#endif
!
   IF (finalized) GOTO 9999
!
   CALL hac_defv
   rt = - MPI_Wtime()
   CALL adios_read_init_method (ADIOS_READ_METHOD_BP, acomm, ARMP, aerr)
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   CALL adios_read_open_file (fh,fname,ADIOS_READ_METHOD_BP,acomm, aerr)
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   CALL adios_get_scalar (fh, DAC, ndims_v, ierr)
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   IF ( ndims_v .GT. NDIMS ) THEN
    ierr = ERRV_M
   END IF
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   IF ( DODEBUG_ .GT. arank ) THEN
    WRITE (OU,*) "READ ",DAC,": ",NDIMS," on ", NDIMS
   END IF
   !IF ( PRODUCT (SHAPE (GAS) ) .EQ. 0 ) THEN
   ! ierr = ERRV_M
   !END IF
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   ada(1:NDIMS,OFF_) = 1 ! LBA(1:NDIMS) + 1 - ozb
   ierr = check_dda(ada(1:NDIMS,OFF_))
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   !CALL check_ada(SHAPE(GAS), ada)
   DO j = 1, ADT_
    DO i = 1, NDIMS
     IF (j.EQ.OFF_) CYCLE ! this variable should be locally set: read values won't match
     CALL hac_gadn (vid,i,j)
     CALL adios_get_scalar (fh, vid, ada(i,j), ierr)
     HLST_HAC_GOTO_END_ON_ERROR(9999)
     IF ( DODEBUG_ .GT. arank ) THEN
      WRITE (OU,*) "READ ",vid,": ",ada(i,j)," on ", arank
     END IF
    END DO
   END  DO
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   ada(1:NDIMS,OFF_) = 1 ! LBA(1:NDIMS) + 1 - ozb
#if (HLST_HAC_EXTRA_METADATA==1)
   CALL adios_get_scalar (fh, "ts", ts, ierr)
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   CALL adios_get_scalar (fh, "st", st, ierr)
   HLST_HAC_GOTO_END_ON_ERROR(9999)
#endif
   CALL adios_read_close (fh, ierr)
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   IF (.NOT.PRESENT(GADA_)) THEN
    ! WRITE (OU,*) "GADA ",ada(1:NDIMS,GLO_)
   ELSE
    GADA_(1:NDIMS) = INT(ada(1:NDIMS,GLO_)) ! FIXME: conversion overflow danger (no check on GADA_)
   END IF
!
   ntc = hac_ntc()
!
   IF (PRESENT(ntc_)) THEN
    ntc_ = ntc
   ELSE
    ! WRITE (OU,*) "TYPE: ",ntc, hac_asz()
   END IF
!
   HAC_RET_IF_PRESENT(ierr_,ierr)
#if (HLST_HAC_EXTRA_METADATA==1)
   HAC_RET_IF_PRESENT(ts_,ts)
   HAC_RET_IF_PRESENT(st_,st)
#endif
!
9999 HLST_HAC_HANDLE_IERR(ierr_,ierr)
  END SUBROUTINE hac_info
!
!> Get prefix
  SUBROUTINE hac_gp(p,ni_,nf_)
!
   CHARACTER(len=*),INTENT(INOUT) :: p
   INTEGER(KIND=OFF_K),OPTIONAL,INTENT(OUT) :: ni_
   REAL(KIND=8),       OPTIONAL,INTENT(OUT) :: nf_
!
   INTEGER(KIND=OFF_K),PARAMETER :: th = 1 ! threshold
   INTEGER(KIND=OFF_K),PARAMETER :: su = XXKC
!
   INTEGER(KIND=OFF_K) :: ni
   REAL(KIND=8) :: nf
!
   HAC_SET_IF_PRESENT(nf_,nf)
   HAC_SET_IF_PRESENT(ni_,ni)
   p = 'B'
   ni = INT8(nf)
   IF ( ni / su .LT. th ) THEN
    GOTO 9999
   END IF
!
   IF ( ni / su .GE. th ) THEN
    ni = ni / su
    nf= nf/ su
    p = 'K'
   END IF
   IF ( ni / su .GE. th ) THEN
    ni = ni / su
    nf= nf/ su
    p = 'M'
   END IF
   IF ( ni / su .GE. th ) THEN
    ni = ni / su
    nf= nf/ su
    p = 'G'
   END IF
   IF ( ni / su .GE. th ) THEN
    ni = ni / su
    nf= nf/ su
    p = 'T'
   END IF
   IF ( ni / su .GE. th ) THEN
    ni = ni / su
    nf= nf/ su
    p = 'P'
   END IF
!
   p = TRIM(p) // XBS
9999 HAC_RET_IF_PRESENT(nf_,nf)
   HAC_RET_IF_PRESENT(ni_,ni)
  END SUBROUTINE hac_gp
!
!> Prints Read/Write Statistics
  SUBROUTINE hac_prtsts(lcib,row,et,fname,tcib_)
   INTEGER(KIND=OFF_K),INTENT(IN) :: lcib ! local contribution in bytes
   LOGICAL, INTENT(IN) :: row ! read (.TRUE.) or write (.FALSE.)
   ! CHARACTER(LEN=*),INTENT(IN) :: ms ! message string
   CHARACTER(LEN=*),INTENT(IN) :: fname !< Path to the checkpoint file. Any length is allowed.
   REAL(KIND=TIME_K),INTENT(IN) :: et !< Elapsed time
!
   INTEGER(KIND=OFF_K) :: tcib ! local contribution in bytes/<x>bytes
   INTEGER(KIND=OFF_K),OPTIONAL,INTENT(INOUT) :: tcib_ !
   CHARACTER(LEN=SUL_) :: SS, XS
   REAL(KIND=8) :: tm, sm ! transfer and speed  measures
!
   tcib = 0
   CALL MPI_Allreduce( lcib, tcib, 1, HAC_M_INT_T, MPI_SUM, acomm, ierr)
   CALL hac_herr
!
   tm = tcib
   CALL hac_gp(XS,nf_=tm)
!
   IF ( verbose.GT.0 .AND. amroot ) THEN
    IF(row) THEN
     sm = (1.0/et)*tcib
     CALL hac_gp(SS,nf_=sm)
     WRITE(OU,'(a,a,f10.2,a,f10.3,a,f10.2,a,a)') BNR,"Read    ",&
             &tm," "//TRIM(XS)//" in ",et," s, at ",sm,&
             &" "//TRIM(SS)//"/s from ",TRIM(fname)
    ELSE
     sm = (1.0/et)*tcib
     CALL hac_gp(SS,nf_=sm)
     WRITE(OU,'(a,a,f10.2,a,f10.3,a,f10.2,a,a)') BNR,&
     &"Written ",tm," "//TRIM(XS)//" in ",et," s, at ",&
     &sm," "//TRIM(SS)//"/s  to  ",TRIM(fname)
    END IF
   END IF
!
  HAC_RET_IF_PRESENT(tcib_,tcib)
!
  END SUBROUTINE hac_prtsts
!
 SUBROUTINE hac_wrt_sda(SB,TDA,SA,wts_)
  INTEGER,INTENT(IN) :: TDA(:) !
  CHARACTER(LEN=*),INTENT(IN) :: SB,SA ! String  Before/After
  LOGICAL,INTENT(IN),OPTIONAL :: wts_ ! Want trailer string
!
  CHARACTER(LEN=MDS_) :: mts ! Modified trailer string
  INTEGER(KIND=OFF_K) :: nels
  REAL(KIND=8) :: nf
  CHARACTER(LEN=SUL_) :: XS
!
  mts = ''
  IF ( PRESENT(wts_) ) THEN
   nels = PRODUCT(TDA(1:))
   nf = REAL(nels*hac_asz())
   CALL hac_gp(XS,nf_=nf)
   WRITE (mts,'(a,i0,a,f0.1,a)')" [=",nels," -- ",nf," "//TRIM(XS)//" ]"
  END IF
!
   SELECT CASE (SIZE(TDA))
    CASE (2)
     WRITE(*,'(a,"[",&
      &i0,",",i0,"]",a)')SB,&
     &TDA(1),TDA(2),SA//TRIM(mts)
    CASE (4)
     WRITE(*,'(a,"[",&
      &i0,",",i0,",",i0,",",i0,"]",a)')SB,&
     &TDA(1),TDA(2),TDA(3),TDA(4),SA//TRIM(mts)
    CASE (6)
     WRITE(*,'(a,"[",&
      &i0,",",i0,",",i0,",",i0,",",i0,",",i0,"]",a)')SB,&
     &TDA(1),TDA(2),TDA(3),TDA(4),TDA(5),TDA(6),SA//TRIM(mts)
    CASE (8)
     WRITE(*,'(a,"[",&
      &i0,":",i0,",",i0,":",i0,",",i0,":",i0,",",i0,":",i0,"]",a)')SB,&
     &TDA(1),TDA(2),TDA(3),TDA(4),TDA(5),TDA(6),TDA(7),TDA(8),&
     &SA//TRIM(mts)
    CASE (12)
     WRITE(*,'(a,"[",&
      &i0,":",i0,",",i0,":",i0,",",i0,":",i0,",",i0,":",i0,","&
     &,i0,":",i0,",",i0,":",i0,"]",a)')SB,&
     &TDA(1),TDA(2),TDA(3),TDA(4),TDA(5),TDA(6),TDA(7),&
     &TDA(8),TDA(9),TDA(10),TDA(11),TDA(12),SA//TRIM(mts)
   END SELECT
 END SUBROUTINE hac_wrt_sda
!
!
!> Reads the distributed array from the checkpoint file.
!! An arbitrary array slice can be read from whatever rank.
!! For this reason unlike hac_write, here the index base cannot be unambiguously inferred from the input,
!! and it is advised to specify it.
!!
  SUBROUTINE hac_read(GAS, LBA, fname, LDA_, gaib_, ierr_)
   IMPLICIT NONE
!
   HAC_NTYPE,HAC_DSPEC,INTENT(INOUT) :: GAS !< Global Array Slice. \c SHAPE(GAS) will be assumed to be local
!! slice dimensions; LBOUND(GAS) its global lower bounds, unless \a LBA is provided. Expected to be already allocated.
   INTEGER,INTENT(IN),DIMENSION(:) :: LBA !< Lower Bounds Array. Either 1- or 0-based, on *all* the dimensions.
!! Offset of the local slice in GAS with respect to the global array.
   CHARACTER(LEN=*),INTENT(IN) :: fname !< Path to the checkpoint file. Any length is allowed.
   INTEGER,OPTIONAL,INTENT(IN),DIMENSION(:) :: LDA_ !< Local Dimensions Array. Dimensions of the local slice in GAS.
!! If omitted, SHAPE(GAS) will be used instead.
   INTEGER,OPTIONAL,INTENT(OUT) :: ierr_ !< Error status variable. If omitted, errors will be fatal.
   INTEGER,OPTIONAL,INTENT(IN) :: gaib_ !< Global Array Index Base; can be either 1 or 0; if not specified,
!! this will be taken to be the minimum value across LBA instances.
!
   CHARACTER(LEN=VIL) :: vid
   INTEGER :: i,j
   INTEGER(KIND=AI_K),DIMENSION(MDIMS_) :: off, rds !offset, read size
   INTEGER(KIND=AI_K) :: sel  ! ADIOS selection object
   INTEGER(KIND=AI_K) :: lcib  ! local contribution in bytes
   INTEGER(KIND=AI_K) :: fh !
   INTEGER(KIND=OFF_K) :: ndims_v
   INTEGER(KIND=OFF_K) :: ada(MDIMS_,ADT_) !< Array's Dimensions Array
   INTEGER(KIND=OFF_K) :: lmo, ozb ! Local Minimum Offset, One or Zero Base ?
!
   IF (finalized) GOTO 9999
!
   off = -1 ! initialization
   rds = -1 ! initialization
   ozb = 0
!
   IF (PRESENT(gaib_)) THEN
    lmo = gaib_ ! Either 1 or 0
     IF ( gaib_ .NE. 0 .AND. gaib_ .NE. 1 ) THEN
      ierr = ERRV_G
     ENDIF
    ELSE
    lmo = MINVAL(LBA(1:NDIMS) )
   END IF
   CALL hac_herr
   lmo = MIN(1_OFF_K,MAX(lmo,0_OFF_K))
   CALL MPI_Allreduce( lmo, ozb, 1, HAC_M_INT_T, MPI_MIN, acomm, ierr )
   CALL hac_defv
!
   IF (fime .NE. 0) CALL hac_flush_attempt()
!
   rt = - MPI_Wtime()
   CALL adios_read_init_method (ADIOS_READ_METHOD_BP, acomm, ARMP, aerr)
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   CALL adios_read_open_file (fh,fname,ADIOS_READ_METHOD_BP,acomm, aerr)
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   CALL adios_get_scalar (fh, DAC, ndims_v, ierr)
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   IF ( ndims_v .GT. NDIMS ) THEN
    ierr = ERRV_M
   END IF
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   IF ( DODEBUG_ .GT. arank ) THEN
    WRITE (OU,*) "READ ",DAC,": ",NDIMS," on ", NDIMS
   END IF
   IF ( PRODUCT (SHAPE (GAS) ) .EQ. 0 ) THEN
    ierr = ERRV_M
   END IF
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   ada(1:NDIMS,OFF_) = LBA(1:NDIMS) + 1 - ozb
   ierr = check_dda(ada(1:NDIMS,OFF_))
   CALL hac_herr
   CALL check_ada(SHAPE(GAS), ada)
   ada(1:NDIMS,LOC_) = SHAPE(GAS) ! The local SHAPE will determine how much we read
   DO j = 1, ADT_
    DO i = 1, NDIMS
     IF (j.EQ.OFF_) CYCLE ! this variable should be locally set: read values won't match
     IF (j.EQ.LOC_) CYCLE ! this variable should be locally set: depends on SHAPE(GAS)
     CALL hac_gadn (vid,i,j)
     CALL adios_get_scalar (fh, vid, ada(i,j), ierr)
     HLST_HAC_GOTO_END_ON_ERROR(9999)
     IF ( DODEBUG_ .GT. arank ) THEN
      WRITE (OU,*) "READ ",vid,": ",ada(i,j)," on ", arank
     END IF
    END DO
   END  DO
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   ada(1:NDIMS,OFF_) = LBA(1:NDIMS) + 1 - ozb
   CALL adios_read_close (fh, ierr)
   HLST_HAC_GOTO_END_ON_ERROR(9999)
!
   CALL check_ada(SHAPE(GAS), ada)
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   CALL adios_read_open_file (fh,fname,ADIOS_READ_METHOD_BP,acomm,ierr)
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   off(1:NDIMS) = ada(1:NDIMS, OFF_) - 1
   IF (PRESENT(LDA_)) THEN
    IF ( dodebug_.GT.0 .AND. amroot ) &
     &WRITE(OU,*)'Will read custom shape:', LDA_
    ada(1:NDIMS, LOC_) = LDA_(1:NDIMS) ! Overriding SHAPE(GAS)
   END IF
   rds(1:NDIMS) = ada(1:NDIMS, LOC_)
!
   CALL adios_selection_boundingbox (sel, INT(NDIMS), off, rds)
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   CALL adios_schedule_read (fh, sel, TRIM("/" // AID), 0, 1, GAS, ierr)
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   CALL adios_perform_reads (fh, ierr)
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   CALL adios_read_close (fh, ierr)
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   CALL adios_selection_delete (sel)
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   CALL adios_read_finalize_method (ADIOS_READ_METHOD_BP, ierr)
!
   IF (DODEBUG_.GT.arank.AND.PRODUCT(ada(:,LOC_))*asize.LT.MDS_ ) THEN
    WRITE (OU,*) "READ ",AID,": ",GAS," on ", arank
   END IF
!
   rt = rt + MPI_Wtime()
   lcib = mis + atss * PRODUCT(ada(1:NDIMS,LOC_))
!
   IF (DODEBUG_.GT.arank.AND.PRODUCT(ada(:,LOC_))*asize.LT.MDS_ ) THEN
    IF ( amroot ) WRITE (OU,*)"SIZEOF(NDIMS):", itss
    IF ( amroot ) WRITE (OU,*)"SIZEOF(TYPE):", atss
    ! IF ( amroot ) WRITE (OU,*)"SIZEOF(ada(1,1)):", itss
   END IF
!
   CALL hac_prtsts(lcib-mis,.TRUE.,rt,TRIM(fname))
!
9999 HLST_HAC_HANDLE_IERR(ierr_,ierr)
  END SUBROUTINE hac_read
!
!> Writes the global array to the checkpoint file.
  SUBROUTINE hac_write (GAS, LBA, fname, ts_, st_, a_buf_mib_,ws_,ierr_)
   IMPLICIT NONE
!
   HAC_NTYPE,HAC_DSPEC,INTENT(IN) :: GAS !< Global Array Slice. \c SHAPE(GAS) will be assumed to be the
   !!local slice dimensions.
   INTEGER,INTENT(IN),DIMENSION(:) :: LBA !< Lower Bounds Array. Either 1- or 0-based, on *all* the dimensions.
   !!Offset of the local slice in GAS with respect to the global array.
   CHARACTER(LEN=*),INTENT(IN) :: fname !< Path to the checkpoint file. Any length is allowed.
   INTEGER,OPTIONAL,INTENT(OUT) :: ierr_ !< Error status variable. If omitted, errors will be fatal.
   INTEGER,OPTIONAL,INTENT(IN) :: a_buf_mib_ !< ADIOS buffer size in MiB (1 MiB=1024*1024 bytes; see `man 7 units`).
   !!If omitted, it will be automatically set at the first write. If you foresee a subsequent write to be larger,
   !!submit an upper limit here at the first call.
#ifdef HLST_HAC_EXTRA_METADATA
   REAL(KIND=8),OPTIONAL,INTENT(IN) :: ts_ !< Meta-data scalar: Time Step. Can be useful to certain users.
   REAL(KIND=8),OPTIONAL,INTENT(IN) :: st_ !< Meta-data scalar: Simulation time. Can be useful to certain users.
!
   REAL(KIND=8) :: ts = 0.0, st = 0.0 ! DOUBLE PRECISION is not supported by ADIOS
#endif
   LOGICAL,OPTIONAL,INTENT(IN) :: ws_ !< If .TRUE., request synchronization. Off be default.
!
   INTEGER(KIND=AI_K) :: ah,ats,ags=0 !ADIOS handle, total size, group size
   INTEGER :: i,j

   !CHARACTER(LEN=VIL) :: vid
   CHARACTER(LEN=50) :: vid

   INTEGER(KIND=OFF_K) :: ada(MDIMS_,ADT_), adar(MDIMS_), adas(MDIMS_) !< Array's Dimensions Array (and recv/send buffer)
   INTEGER(KIND=OFF_K) :: lmo, ozb = 0 ! Local Minimum Offset, One or Zero Base
   INTEGER(KIND=OFF_K), PARAMETER :: hambs = 702545920  ! Hardcoded ADIOS maximal buffer size

   REAL(KIND=TIME_K) :: adios_close_time, max_adios_close_time
!
   IF(finalized) GOTO 9999
   wt = - MPI_Wtime()
   lmo = MINVAL(LBA(1:NDIMS) )
   lmo = MIN(1_OFF_K,MAX(lmo,0_OFF_K))
   CALL MPI_Allreduce( lmo, ozb, 1, HAC_M_INT_T, MPI_MIN, acomm, ierr )
!
!  TODO: here, draconian checks on LBA may be of use.
!
   ada(1:NDIMS,OFF_) = LBA(1:NDIMS) + 1 - ozb
   ierr = check_dda(ada(1:NDIMS,OFF_))
   HLST_HAC_GOTO_END_ON_ERROR(9999)

   !CALL hac_defv(A_NA_=adios_mpi_agg, A_NO_=adios_mpi_ost)
   !CALL adios_open (ah, A_GN, TRIM(fname), "w", acomm, aerr)

   CALL adios_open (ah, "restart", TRIM(fname), "w", acomm, aerr)
   initialized_vars = .TRUE.


   HLST_HAC_GOTO_END_ON_ERROR(9999)
   ada(:,:) = 1
   ada(1:NDIMS, OFF_) = LBA(1:NDIMS) + 1 - ozb
   ada(1:NDIMS, LOC_) = SHAPE(GAS)
   ada(1:NDIMS, GLO_) = ada(1:NDIMS,OFF_) + ada(1:NDIMS,LOC_) - 1
   adas(1:NDIMS) = ada(1:NDIMS,GLO_)
   CALL MPI_Allreduce(adas,adar,NDIMS,HAC_M_INT_T,MPI_MAX,acomm,ierr)
   ada(1:NDIMS,GLO_) = adar(1:NDIMS)
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   ags = mis + atss * PRODUCT( ada(1:NDIMS, LOC_))
!
! FIXME
!
! Need a statement in the documentation about the
! WARN : MPI_AMR method (BG): The max allowed aggregation buffer is 704643072 bytes.
! But this ADIOS method needs 1814080644 bytes for aggregation warning.
!
!  Need to document max possible I/O size. Maybe hard-code this. In both ADIOS-1.5 and 1.6,in adios_mpi_amr.c,
!  #define MAX_AGG_BUF 704643072  so max array is half of that: 2*335*1024*1024 = 702545920.
!   IF( ags .LT. hambs ) STOP
!
   IF ( DODEBUG_ .GT. arank ) THEN
    IF ( amroot )WRITE (OU,*)"WRITE NDIMS: ",NDIMS," ags: ",&
            & ags, " PROD(dims):", PRODUCT(ada(1:NDIMS,LOC_)),&
            & " offsets:", (ada(1:NDIMS,OFF_)),&
            & " ozb:", ozb,&
            & " offsets:", (ada(1:NDIMS,OFF_))
   END IF
   IF (abao.EQV..FALSE.) THEN
    IF (PRESENT(a_buf_mib_)) THEN
      a_buf_mib = INT( 1 + ( a_buf_mib_ ) ) ! The extra 1 is to round up.
    ELSE
      a_buf_mib = INT( 1 + ( ags / MiB_ ) )
    END IF
    call adios_set_max_buffer_size(2*a_buf_mib)

    ! ADIOS-1.5 allows adios_allocate_buffer only once.
    abao = .TRUE.
   END IF
   CALL adios_group_size (ah, ags, ats, aerr)
   HLST_HAC_GOTO_END_ON_ERROR(9999)

   !CALL adios_write (ah, DAC, NDIMS, aerr)
   CALL adios_write (ah, "num_dims", NDIMS, aerr)

   HLST_HAC_GOTO_END_ON_ERROR(9999)
   CALL check_ada(SHAPE(GAS), ada)
   ada(:,OFF_) = ada(:,OFF_) - 1 ! for ADIOS
   DO j = 1, ADT_
    DO i = 1, NDIMS

     !CALL hac_gadn(vid,i,j)
     !CALL adios_write (ah, vid, ada(i,j), aerr)

     CALL dim_name(vid,i,j)
     CALL adios_write (ah, vid, ada(i,j), aerr)

     HLST_HAC_GOTO_END_ON_ERROR(9999)
    END DO
   END DO
#if (HLST_HAC_EXTRA_METADATA==1)
   IF (PRESENT(ts_)) ts = ts_
   IF (PRESENT(st_)) st = st_
   !CALL adios_write (ah, "ts", ts, aerr)
   !CALL adios_write (ah, "st", st, aerr)

   CALL adios_write (ah, "timestep", ts, aerr)
   CALL adios_write (ah, "time", st, aerr)
#endif
   ada(:,OFF_) = ada(:,OFF_) + 1 ! for ADIOS
!
   !CALL adios_write (ah, AID, GAS, aerr)

   CALL adios_write (ah, "g", GAS, aerr)

   HLST_HAC_GOTO_END_ON_ERROR(9999)

   adios_close_time = - MPI_Wtime()
   CALL adios_close (ah, aerr)
   adios_close_time = adios_close_time + MPI_Wtime()

   HLST_HAC_GOTO_END_ON_ERROR(9999)
   IF (DODEBUG_.GT.arank.AND.PRODUCT(ada(:,LOC_))*asize.LT.MDS_) THEN
    WRITE (OU,*) "WRITE ",AID,": ",GAS," on ", arank
   END IF
!
   wt = wt + MPI_Wtime()
   yt = 0
   IF ( PRESENT(ws_) ) THEN
    IF ( ws_ ) THEN
     yt = - MPI_Wtime()
     CALL hac_sync(TRIM(fname))
     yt = yt + MPI_Wtime()
    END IF
   END IF
   CALL hac_prtsts(ags-mis,.FALSE.,wt,TRIM(fname))

   if (verbose.gt.0) then
       call MPI_Reduce (adios_close_time, max_adios_close_time, 1, MPI_DOUBLE_PRECISION, MPI_MAX, RR, acomm, aerr)
       if (amroot) then
           write (OU, *) new_line('A')//&
               "ADIOS restart write, time step: ", itime, new_line('A')//&
               "groupsize: ", ags, new_line('A')//&
               "close time: ", max_adios_close_time, new_line('A')
       endif
   endif

!
   IF ( verbose.GT.RR .AND. amroot ) THEN
   END IF
9999 HLST_HAC_HANDLE_IERR(ierr_,ierr)
  END SUBROUTINE hac_write
!
!> Finalizes module hlst_adios_checkpoint.
!! It must be called once, before MPI_Finalize.
!! \note #hac_exit reports successful operation as default.
  SUBROUTINE hac_exit(ierr_)
    IMPLICIT NONE
!
   INTEGER,OPTIONAL,INTENT(OUT) :: ierr_ !< Error status variable. If omitted, errors will be fatal.
!
    CALL hac_adios_lib_exit
    IF ( PRESENT(ierr_) ) THEN
     ierr_ = OKV_
    END IF
    HLST_HAC_HANDLE_IERR(ierr_,ierr)
  END SUBROUTINE hac_exit
!
!> Gets an INTEGER environment variable value.
INTEGER FUNCTION get_env_opt_i(name,dval_)
  IMPLICIT NONE
  CHARACTER(LEN=*),INTENT(IN) :: name
  CHARACTER(LEN=MEVVL_) :: value
  INTEGER :: i, status = 0
  INTEGER,INTENT(IN),OPTIONAL :: dval_
!
  ! CALL GETENV(name,value)
  CALL GET_ENVIRONMENT_VARIABLE(name,value,status=status)
  READ (value,'(i10)') i
  IF ( PRESENT(dval_) ) THEN
    !IF (status.EQ.1 .AND. i .EQ. 0 ) i = dval_
    IF (status.EQ.1 ) i = dval_ ! if var unset, override
    ! IF (i==0) i = dval_
  END IF
  get_env_opt_i = i
END FUNCTION get_env_opt_i
!
!> Gets a CHARACTER string.
SUBROUTINE get_env_opt_s(name,val,sval_)
  IMPLICIT NONE
  CHARACTER(LEN=*),INTENT(IN) :: name
  INTEGER :: status = 0
  CHARACTER(LEN=*),INTENT(INOUT) :: val
!
  CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: sval_
!
  CALL GET_ENVIRONMENT_VARIABLE(name,val,status=status)
  IF ( PRESENT(sval_) ) THEN
    IF (status.EQ.1 ) val = sval_ ! if var unset, override
  END IF
END SUBROUTINE get_env_opt_s
!
!> Local Array Allocate (and initialize).
SUBROUTINE hac_laa(AS,LB,UB,erank,diodz)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: LB(NDIMS) !< Lower Bounds
  INTEGER,INTENT(IN) :: UB(NDIMS) !< Upper Bounds
  INTEGER,INTENT(IN) :: erank !
  LOGICAL,INTENT(IN) :: diodz ! Do init or do zero ?
  HAC_NTYPE,HAC_DSPEC,ALLOCATABLE,INTENT(INOUT) :: AS !< Array Slice
!
#if   (HLST_HAC_USE_ADIM==6)
  INTEGER,PARAMETER :: MF = 10
#elif   (HLST_HAC_USE_ADIM==4)
  INTEGER,PARAMETER :: MF = 100
#elif   (HLST_HAC_USE_ADIM==2)
  INTEGER,PARAMETER :: MF = 1000
#else
#error !
#endif
  INTEGER :: i1, i2, i3, i4, i5, i6
!
#if   (HLST_HAC_USE_ADIM==6)
 ALLOCATE(AS(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3)&
           &,LB(4):UB(4),LB(5):UB(5),LB(6):UB(6)))
#elif   (HLST_HAC_USE_ADIM==4)
 ALLOCATE(AS(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3)&
           &,LB(4):UB(4)))
#elif   (HLST_HAC_USE_ADIM==2)
 ALLOCATE(AS(LB(1):UB(1),LB(2):UB(2)))
#else
#error !
#endif
   AS =        10**8 * (erank + 1)
#if (HLST_HAC_USE_REAL!=1)
   AS = AS + CMPLX( 0.0, erank + 1 )
#endif
   IF ( diodz ) THEN
#if   (HLST_HAC_USE_ADIM==6)
   FORALL(i1=LB(1):UB(1),i2=LB(2):UB(2),i3=LB(3):UB(3),i4=LB(4):UB(4)&
   &,i5=LB(5):UB(5),i6=LB(6):UB(6)) &
   &AS(i1,i2,i3,i4,i5,i6) = &
   & MF**0*i1+MF**1*i2+MF**2*i3+MF**3*i4+MF**4*i5+MF**6*i6
#elif   (HLST_HAC_USE_ADIM==4)
   FORALL(i1=LB(1):UB(1),i2=LB(2):UB(2),i3=LB(3):UB(3),i4=LB(4):UB(4))&
                   &AS(i1,i2,i3,i4) = i1+MF*i2+i3*MF**2+MF**3*i4
#elif   (HLST_HAC_USE_ADIM==2)
   FORALL(i1=LB(1):UB(1),i2=LB(2):UB(2)) AS(i1,i2) = i1+MF*i2
#else
#error !
#endif
!
   ELSE
    AS = 0
   ENDIF
!
END SUBROUTINE hac_laa
!
!> Performs a basic correctness test run.
!! Single write and read with the specified local lower/upper bounds.
 SUBROUTINE hac_swar(ecomm,LB,UB,fname,reps_)
  IMPLICIT NONE
!
  INTEGER,INTENT(INOUT) :: ecomm
  INTEGER,INTENT(IN) :: LB(HLST_HAC_USE_ADIM)
  INTEGER,INTENT(IN) :: UB(HLST_HAC_USE_ADIM)
  INTEGER,OPTIONAL :: reps_  ! samples
!
  CHARACTER, PARAMETER :: s = ' '
  REAL(KIND=TIME_K), PARAMETER :: IBT = 1000.0*1000.0*1000.0 ! Impossibly Big Time
  CHARACTER(LEN=AMFNL_),INTENT(IN) :: fname
!
  INTEGER :: LADA(NDIMS),GADA(NDIMS) !< Local/Global Array's Dimensions Array
  CHARACTER :: ntc !< Numerical Type Character.
  CHARACTER(LEN=MPSL_) :: ds
  INTEGER :: erank, esize, ii, ti
  HAC_NTYPE,HAC_DSPEC,ALLOCATABLE :: GAS,GCP !< Global Array Slice
  REAL(KIND=TIME_K) :: twt, trt, tyt !< Total Write time, total read time, total synchronized write time
  REAL(KIND=TIME_K) :: lwt, lrt, lyt, uwt, urt, uyt, awt, art, ayt !<
  REAL(KIND=TIME_K), PARAMETER :: zero = 0.0
  INTEGER :: vbackup
  INTEGER(KIND=OFF_K) :: tcib !
!
  tcib = 0 !
  ti = 1
!
  tyt = 0.0
  twt = 0.0
  trt = 0.0
  ayt = 0.0
  awt = 0.0
  art = 0.0
  lyt = IBT
  lwt = IBT
  lrt = IBT
  uyt = 0.0
  uwt = 0.0
  urt = 0.0
!
  HAC_SET_IF_PRESENT(reps_,ti)
!
  CALL MPI_Comm_rank (ecomm, erank, ierr)
  CALL MPI_Comm_size (ecomm, esize, ierr)
!
  CALL hac_laa(GAS,LB,UB,erank,.TRUE.)
  CALL hac_laa(GCP,LB,UB,erank,.TRUE.)
!
  GCP = GAS ! A copy (for recheck after read).
  IF ( verbose.GT.erank ) THEN
   WRITE(OU,'(a,a)') BNR//''
   CALL hac_wrt_sda(BNR//'Local Array Dimensions    :',SHAPE(GAS),&
           &'',.TRUE.)
  END IF
!
  LADA = LBOUND(GAS)+SHAPE(GAS)
  CALL MPI_Allreduce(LADA,GADA,NDIMS,HAC_4_INT_T,MPI_MAX,acomm,ierr)
!
  IF ( verbose.GT.erank.AND.amroot ) THEN
   CALL hac_wrt_sda(BNR//'Global Array Dimensions: ',GADA,'')
  END IF
!
  vbackup = verbose
  verbose = HAC_QUIET
  DO ii = 1, ti
   CALL hac_write (GAS, LBOUND(GAS), TRIM(fname), &
           ts_=zero,st_=zero,ws_=.TRUE.&
           &)
   HAC_MINMAXSUM(wt,lwt,uwt,twt)
   HAC_MINMAXSUM(yt,lyt,uyt,tyt)
   GAS = 0
   CALL hac_info (fname, ierr)
   CALL hac_read (GAS, LBOUND(GAS), TRIM(fname))
   HAC_MINMAXSUM(rt,lrt,urt,trt)
   IF (SUM(GAS-GCP).NE.0.0) STOP "Critical error in restart --- &
           &read values do not match written ones !"
  END DO
  verbose = vbackup
  awt = twt/ti
  ayt = tyt/ti
  art = trt/ti
!
  CALL hac_info (fname, ierr, GADA_=GADA )
  IF ( verbose.GT.erank.AND.amroot ) THEN
  ! CALL hac_wrt_sda(BNR//'Global Array Dimensions: ',GADA,'')
  ! CALL hac_wrt_sda(BNR//'Global Array Lower Bounds: ',LBOUND(GAS),'')
  END IF
  CALL hac_info (fname, ierr, ntc_=ntc)
  IF ( verbose.GT.0 ) THEN
   CALL hac_gds(ecomm,LB,UB,ds)
   IF ( verbose.GT.erank .AND. amroot ) THEN
    WRITE(OU,'(a,a)') BNR//'Descriptor String: ',TRIM(ds)
   ! WRITE(OU,'(a,a)') BNR//'Numerical Type Character: ',NTC
   END IF
   IF ( ti .GT. 1 .AND. verbose.GT.erank .AND. amroot ) THEN
    WRITE(OU,'(a,i0,a)') BNR//"Average of the ",ti," measured samples:"
    WRITE(OU,'(a,f10.2,a,f10.2,a,f10.2,a,f10.2,a)') BNR//&
    &  "Diff in read time:",&
    &HAC_PCTD(lrt,art),&
    &"%/",HAC_PCTD(urt,art),&
    &"%, diff in write time:",&
    &HAC_PCTD(lwt,awt),&
    &"%/",HAC_PCTD(uwt,awt),"%"
    yt = 0.0
    wt = 0.0
    rt = 0.0
   END IF
  END IF
!
  CALL hac_prtsts(INT8(atss*PRODUCT(UB-LB+1)),.FALSE.,awt,(fname),tcib)
  CALL hac_prtsts(INT8(atss*PRODUCT(UB-LB+1)),.TRUE. ,art,(fname))
!
  IF ( verbose.GT.erank .AND. amroot ) THEN
  IF (.TRUE.) THEN
  WRITE(OU,"(a)") BNR//"# marker: desc samples nt&
          & minwt avgwt maxwt&
          & minrt avgrt maxrt&
          & minyt avgyt maxyt&
          & totMiB"
  WRITE(OU,&
   &"(a,i0,a,i0,&
   &a,f10.2,a,f10.2,a,f10.2,&
   &a,f10.2,a,f10.2,a,f10.2,&
   &a,f10.2,a,f10.2,a,f10.2,&
   &a,f10.2)")&
   &BNR//"HAC:"//s//TRIM(ds)//' '&
   &,ti,s,esize,&
   &s,lwt,s,awt,s,uwt,&
   &s,lrt,s,art,s,urt,&
   &s,lyt,s,ayt,s,uyt,&
   &s,DBLE(tcib/MiB_)
  END IF
  END IF
!
  DEALLOCATE(GAS)
  DEALLOCATE(GCP)
!
 END SUBROUTINE hac_swar
!
! Generate Descriptor String
! Descriptor chars are:
! 'r': replicated
! 'd': distributed, regular
! 'i': distributed, irregular
 SUBROUTINE hac_gds(ecomm, LBA, UBA, ds)
  USE mpi
  IMPLICIT NONE
  INTEGER, INTENT(INOUT) :: ecomm
  INTEGER,INTENT(IN),DIMENSION(:) :: LBA, UBA ! Lower/Upper Bounds Array
  CHARACTER(LEN=MPSL_),INTENT(INOUT) :: ds ! description string
!
  INTEGER,DIMENSION(NDIMS) :: GSA ! Global Shape Array
  INTEGER,ALLOCATABLE,DIMENSION(:,:) :: ALBA, AUBA ! All Lower/Upper Bounds Arrays
  INTEGER :: asize, ierr, erank, MIT, i
  CHARACTER(LEN=MPSL_) :: hs
!
  MIT = MPI_INTEGER
!
  CALL MPI_Comm_size (acomm, asize, ierr)
  CALL MPI_Comm_rank (ecomm, erank, ierr)
!
  WRITE (ds,'(a)') TRIM('')
  ALLOCATE(ALBA(NDIMS,asize),AUBA(NDIMS,asize))
  CALL MPI_Gather(LBA,NDIMS,MIT,ALBA,NDIMS,MIT,RR,ecomm,ierr)
  CALL MPI_Gather(UBA,NDIMS,MIT,AUBA,NDIMS,MIT,RR,ecomm,ierr)
  IF (erank.EQ.0) THEN
   GSA = (/(( MAXVAL(AUBA(i,:)) - MINVAL(ALBA(i,:)) + 1 ), i=1,NDIMS)/)
   WRITE (ds,'(a,i0,a)') hac_ntc(),NDIMS,'_'
   DO i = 1, NDIMS
    WRITE(hs,'(i0)') GSA(i)
    ds = TRIM(ds) // TRIM(hs)
    IF ( MINVAL(ALBA(i,:)) .NE. MAXVAL(ALBA(i,:)) ) THEN ! distributed
     IF (MINVAL(AUBA(i,:)-ALBA(i,:)).NE.MAXVAL(AUBA(i,:)-ALBA(i,:)))THEN ! irregular local dimensions
      WRITE(hs,'(a)') TRIM("i")
     ELSE
      WRITE(hs,'(a)') TRIM("d")
     END IF
    ELSE
     WRITE(hs,'(a)') TRIM("r")
    END IF
    ds = TRIM(ds) // TRIM(hs)
    IF (i .NE. NDIMS) THEN
     WRITE(hs,'(a)') TRIM("_")
     ds = TRIM(ds) // TRIM(hs)
    END IF
   END DO
   DEALLOCATE(ALBA,AUBA)
  END IF
!
 END SUBROUTINE hac_gds
!
! General test driver
 SUBROUTINE hac_gtd(ecomm, fname, miub_, maub_, reps_,&
        & loib_, upib_, fddi_, lddi_, ti_, lszl_, uszl_, mibpt_, sdims_)
  IMPLICIT NONE
  INTEGER, INTENT(INOUT) :: ecomm
  INTEGER, INTENT(IN), OPTIONAL :: miub_ , maub_ ! minimal/maximal unbalance
  INTEGER, INTENT(IN), OPTIONAL :: reps_ ! repetitions
  INTEGER, INTENT(IN), OPTIONAL :: mibpt_ !
  INTEGER, INTENT(IN), OPTIONAL :: lszl_ , uszl_ ! lower/upper size limit
  INTEGER, INTENT(IN), OPTIONAL :: loib_ , upib_ ! lower/upper index base
  INTEGER, INTENT(IN), OPTIONAL :: fddi_ , lddi_ ! first/last distributed dimension index
  INTEGER, INTENT(IN), OPTIONAL :: sdims_ ! sub-dimensions
  INTEGER, INTENT(OUT), OPTIONAL :: ti_ ! test iterations
  CHARACTER(LEN=AMFNL_),INTENT(IN) :: fname
!
  INTEGER :: miub, maub ! minimal/maximal unbalance
  INTEGER :: reps ! repetitions
  INTEGER :: mibpt !
  INTEGER :: lszl, uszl ! lower/upper size limit
  INTEGER :: loib, upib ! lower/upper index base
  INTEGER :: fddi, lddi ! first/last distributed dimension index
!
  INTEGER :: acub ! actual unbalance
  INTEGER :: i,tdim
  INTEGER :: erank
  INTEGER :: ti ! test iterations
  INTEGER :: losz, esize ! , rep
  INTEGER :: LB(NDIMS), UB(NDIMS)
  INTEGER :: sdims ! sub-dimensions
  INTEGER :: ib ! Index Base
!
  INTEGER, DIMENSION(NDIMS+1) :: MDIMS, GDIMS ! , LB, UB
  INTEGER, DIMENSION(NDIMS) :: LCRDS
  INTEGER, DIMENSION(NDIMS+1) :: DDIMS, LDIMS
!
  miub = 0
  mibpt = 0
  maub = 0
  reps = 1
  lszl = 1
  uszl = 1
  loib = 0
  upib = 1
  fddi = 1
  lddi = HLST_HAC_USE_ADIM ! first/last distributed dimension index
  sdims = NDIMS
!
  HAC_SET_IF_PRESENT(miub_,miub)
  HAC_SET_IF_PRESENT(maub_,maub)
  HAC_SET_IF_PRESENT(mibpt_,mibpt)
  HAC_SET_IF_PRESENT(reps_,reps)
  HAC_SET_IF_PRESENT(lszl_,lszl)
  HAC_SET_IF_PRESENT(uszl_,uszl)
  HAC_SET_IF_PRESENT(loib_,loib)
  HAC_SET_IF_PRESENT(upib_,upib)
  HAC_SET_IF_PRESENT(fddi_,fddi)
  HAC_SET_IF_PRESENT(lddi_,lddi)
  HAC_SET_IF_PRESENT(sdims_,sdims)
  sdims = MAX(MIN(sdims,INT(NDIMS)),0)
!
  !sdims = MIN(sdims,5) ! FIXME: this seems to prevent a bug.
!
  IF( sdims .EQ. 0 ) sdims = INT(NDIMS)
!
  IF ( mibpt .NE. 0 ) THEN
    lszl = 0
    uszl = 0
  END IF
!
  CALL MPI_Comm_rank (ecomm, erank, ierr)
  CALL MPI_Comm_size (ecomm, esize, ierr)
!
  ti = 0
!
  IF ( mibpt .NE. 0 ) THEN
  DO ib = loib, upib
  DO acub = miub, maub !
!
   IF ( acub .EQ. 0 ) THEN
    LDIMS = 0
    DDIMS = 0
    MDIMS = 0
!
    LDIMS(sdims+1:) = 1
    DDIMS(sdims+1:) = 1
    MDIMS(sdims+1:) = 1
    ! On OpenMPI-1.6.5 MPI_Dims_create is horribly slow above say, 2**24, hence the following trick of triple call and multiplying.
    CALL MPI_Dims_create((MiB_/hac_asz()), sdims, MDIMS, ierr)
    CALL MPI_Dims_create(mibpt,            sdims, LDIMS, ierr)
    CALL MPI_Dims_create(esize,            sdims, DDIMS, ierr)
    LDIMS(1:sdims) = LDIMS(sdims:1:-1) ! Reversal
    LDIMS(1:sdims) = LDIMS(1:sdims)*MDIMS(1:sdims)
    GDIMS(1:sdims) = LDIMS(1:sdims)*DDIMS(1:sdims)
    MDIMS(1:sdims) = 0
!
    DO i = 1, sdims
     LCRDS(i) = MOD(erank, PRODUCT(DDIMS(i:sdims)))
     LCRDS(i) = LCRDS(i) / PRODUCT(DDIMS(i+1:sdims))
    END DO
!
    LB = ib
    UB = ib
    LB(1:sdims) = ib +    LCRDS(1:sdims)  * LDIMS(1:sdims)
    UB(1:sdims) = ib + (1+LCRDS(1:sdims)) * LDIMS(1:sdims) - 1
!
   ELSE ! ( acub .NE. 0 )
    LDIMS = 0
    DDIMS = 0
    MDIMS = 0
    ! On OpenMPI-1.6.5 MPI_Dims_create is horribly slow above say, 2**24, hence the following trick of triple call and multiplying.
    MDIMS(NDIMS) = 1
    CALL MPI_Dims_create((Ki_/hac_asz()), NDIMS-1, MDIMS, ierr)
    LDIMS(1:NDIMS) = 1
    LDIMS(NDIMS) = Ki_*mibpt
    CALL MPI_Dims_create(esize,            NDIMS, DDIMS, ierr)
    LDIMS(1:NDIMS) = LDIMS(NDIMS:1:-1) ! Reversal
    LDIMS(1:NDIMS) = LDIMS(1:NDIMS)*MDIMS(1:sdims)
    LDIMS(1:NDIMS) = LDIMS(1:NDIMS)*DDIMS(1:NDIMS)
    GDIMS = 0
    CALL MPI_Allreduce(LDIMS,GDIMS,NDIMS,HAC_4_INT_T,MPI_MAX,acomm,ierr)
    write(*,*) erank, LDIMS(1:NDIMS)
    MDIMS(1:sdims) = 0
    STOP
    ! FIXME: need to set LB, UP correctly here.
    ! Either doing a prefix sum or using another trick.
!
!    DO i = 1, NDIMS
!     LCRDS(i) = MOD(erank, PRODUCT(DDIMS(i:NDIMS)))
!     LCRDS(i) = LCRDS(i) / PRODUCT(DDIMS(i+1:NDIMS))
!    END DO
!
    LB = ib +    LCRDS  * LDIMS(1:NDIMS)
    UB = ib + (1+LCRDS) * LDIMS(1:NDIMS) - 1
!
!     LB(tdim) = ib + hac_liscab(losz,esize,erank-1,1,1+acub,.TRUE.) - 0
!     UB(tdim) = ib + hac_liscab(losz,esize,erank+0,1,1+acub,.TRUE.) - 1
   END IF
   CALL hac_swar(ecomm,LB,UB,fname,reps_=reps)
  END DO
  END DO
!
  ELSE ! ( mibpt .EQ. 0 )
  ! DO rep = 1, reps
  DO tdim = fddi, lddi
  DO losz = lszl, uszl !
  DO acub = miub, maub !
  DO ib = loib, upib
    ! All dimensions are replicated and fixed.
    LB = (/(ib,       i=1,NDIMS)/)
    UB = (/(ib+losz-1,i=1,NDIMS)/)
    ! ... except this one which is distributed and if acub.GT.0, skewed.
    IF ( acub .EQ. 0 ) THEN
     LB(tdim) = ib + (erank+0) * losz - 0
     UB(tdim) = ib + (erank+1) * losz - 1
    ELSE
     ! LB(tdim) = ib + (erank+0) * losz * (1+acub * (erank+0)) - 0
     ! UB(tdim) = ib + (erank+1) * losz * (1+acub * (erank+1)) - 1
     LB(tdim) = ib + hac_liscab(losz,esize,erank-1,1,1+acub,.TRUE.) - 0
     UB(tdim) = ib + hac_liscab(losz,esize,erank+0,1,1+acub,.TRUE.) - 1
    END IF
    CALL hac_swar(ecomm,LB,UB,fname,reps_=reps)
    IF (verbose.GT.0) THEN
     ! IF (erank.EQ.0) WRITE(OU,*)' ',losz,'/',100,' ib=',ib
    END IF
    ti = ti + 1
  END DO
  END DO
  END DO
!  IF(erank.EQ.0.AND.MOD(rep,reps/10).EQ.0)WRITE(OU,'(a,i0,a,i0,a,i0)')&
!  'at ',rep,' of ',reps,', sz=',losz
  END DO
  END IF
  !END DO
!
  HAC_RET_IF_PRESENT(ti_,ti)
!
 END SUBROUTINE hac_gtd
!
 SUBROUTINE hac_stress(ecomm,fname)
  IMPLICIT NONE
  INTEGER, INTENT(INOUT) :: ecomm
  CHARACTER(LEN=AMFNL_),INTENT(IN) :: fname
!
  INTEGER, PARAMETER :: miub = 0, maub = 2 ! minimal/maximal unbalance
  INTEGER, PARAMETER :: reps = 1 ! repetitions
  INTEGER, PARAMETER :: lszl = 1, uszl = 4 ! lower/upper size limit
  INTEGER, PARAMETER :: loib = 0, upib = 1 ! lower/upper index base
  INTEGER, PARAMETER :: fddi = 1, lddi = HLST_HAC_USE_ADIM ! first/last distributed dimension index
!
  INTEGER :: erank, ti
!
  CALL MPI_Comm_rank (ecomm, erank, ierr)
!
  IF (erank .EQ. 0) THEN
   WRITE(OU,'(a)')BNR//'Stress testing...'
   WRITE(OU,'(a)')BNR//'(may be damaging for a mechanical hard drive).'
  END IF
  ! verbose = HAC_QUIET
  ! verbose = HAC_VERBOSE
  ! verbose = 3
  CALL hac_gtd(ecomm,fname,miub, maub, reps, &
          & loib, upib, fddi, lddi, ti, lszl, uszl)
!
  IF (erank .EQ. 0) THEN
   WRITE(OU,'(a,i0,a)')&
    &BNR//'Stress testing done: ',ti,' tries successful.'
   WRITE(OU,'(a)')BNR//'END Benchmarking'
  END IF
!
  END SUBROUTINE hac_stress
!
 INTEGER PURE &
        &RECURSIVE FUNCTION hac_liscab(sz,esize,erank,miu,mau,wantcum_)&
        & RESULT(res)
  INTEGER,INTENT(IN) :: sz,miu,mau,erank,esize
  LOGICAL,OPTIONAL,INTENT(IN) :: wantcum_
!
  DOUBLE PRECISION :: p ! proportion
!
  IF ( erank .LT. 0 ) THEN
   res = 0
   RETURN
  END IF
  IF ( PRESENT(wantcum_) ) THEN
   res =        hac_liscab(sz,esize,erank+0,miu,mau) +&
           &    hac_liscab(sz,esize,erank-1,miu,mau,wantcum_)
   RETURN
  END IF
  IF (esize .EQ. 1 .OR. mau .LE. miu) THEN
   res = sz
  ELSE
   p = (DBLE(miu)/mau) + (((DBLE(mau-miu))/mau)*erank)/(esize-1)
   res = FLOOR( DBLE(sz) * p )
  END IF
 END FUNCTION hac_liscab
!
 SUBROUTINE hac_bench(ecomm,fname,mibpt_,reps_)
  ! USE MPI ! , ONLY: MPI_Dims_create, MPI_Comm_rank, MPI_Comm_size
  IMPLICIT NONE
  INTEGER, INTENT(INOUT) :: ecomm
  INTEGER, OPTIONAL, INTENT(IN) :: mibpt_,reps_  ! MiB's per task
  CHARACTER(LEN=AMFNL_),INTENT(IN) :: fname
!
  INTEGER, PARAMETER :: loib = 0, upib = 0 ! lower/upper index base
  INTEGER, PARAMETER :: fddi = 1, lddi = HLST_HAC_USE_ADIM ! first/last distributed dimension index
  INTEGER, PARAMETER :: ib = 0
!
  INTEGER :: miub, maub ! minimal/maximal unbalance
  INTEGER :: reps ! repetitions
  INTEGER :: lszl, uszl  ! lower/upper size limit
  INTEGER :: erank, esize, ti
  INTEGER :: mibpt ! MiB's of data per task, default value.
  INTEGER :: sdims ! subdimensions [based test]
  INTEGER :: odpt ! one dimension per time [based test]
!
  mibpt = 1
  miub = 0
  maub = 0
  reps = 1
  lszl = 0
  uszl = 0
!
  HAC_SET_IF_PRESENT(mibpt_,mibpt)
  HAC_SET_IF_PRESENT(reps_, reps )
!
  CALL MPI_Comm_rank (ecomm, erank, ierr)
  CALL MPI_Comm_size (ecomm, esize, ierr)
!
  IF (erank .EQ. 0) WRITE(OU,'(a)')BNR//'Begin benchmark...'
!
! First test: reasonably balanced multidimensional distribution.
  sdims = get_env_opt_i("HAC_MDIMSTDC",INT(NDIMS))
  !IF (.TRUE.) THEN
  IF ( sdims .NE. 0 ) THEN
!  IF (.FALSE.) THEN
!
   CALL hac_gtd(ecomm,fname, miub, maub, reps,&
          & loib, upib, fddi, lddi, ti, mibpt_=mibpt, sdims_=sdims )
!
  END IF
!
! Second test: only one dimension distributed; the others replicated, no skew.
  ! IF (.TRUE.) THEN
  odpt = get_env_opt_i("HAC_ODPT",1)
  IF ( odpt .NE. 0 ) THEN
!  IF (.FALSE.) THEN
   uszl = ((mibpt *MiB_)/hac_asz())
   uszl = INT( DBLE(uszl)**(1.0/DBLE(NDIMS)) ) ! NDIMS'th root
   lszl = uszl
   CALL hac_gtd(ecomm,fname, miub, maub, reps,&
          & loib, upib, fddi, lddi, ti, lszl, uszl)
  END IF
!
! Third test: only one dimension distributed; the others replicated, with skew.
  IF (.FALSE.) THEN ! FIXME: this case is broken, a dead line.
!  IF (.TRUE.) THEN
   uszl = ((mibpt *MiB_)/hac_asz())
   uszl = uszl / esize
   uszl = INT( DBLE(uszl)**(1.0/DBLE(NDIMS)) ) ! NDIMS'th root
   lszl = uszl
   miub = 1
   maub = 1
   WRITE(*,*) esize,mibpt,lszl,uszl
   CALL hac_gtd(ecomm,fname, miub, maub, reps, loib, upib, fddi, lddi,&
           & ti, lszl, uszl)
  END IF
!
!  IF (.TRUE.) THEN
  IF (.FALSE.) THEN ! FIXME: this case is broken, under rework.
   miub = 1
   maub = 1
   CALL hac_gtd(ecomm,fname, miub, maub, reps, &
          & loib, upib, fddi, lddi, ti, mibpt_=mibpt)
  END IF
!
  IF (erank .EQ. 0) WRITE(OU,'(a)')BNR//'End benchmark...'
!
  END SUBROUTINE hac_bench
!
!> Performs basic error testing.
 SUBROUTINE hac_ex_ckp_err_test
  IMPLICIT NONE
!
  eanf = 1 !
!
! CALL hac_info ('/dev/zero', ierr) ! This crashes ADIOS. No workaround.
! CALL hac_info ('/dev/non-existent-file', ierr) ! This does not crash ADIOS but leads it to an inconsistent state.
  ierr = OKV_
  CALL hac_rerr
  eanf = 1 - eanf
!
 END SUBROUTINE hac_ex_ckp_err_test
!
!> This is a complete MPI program and can only be called once, e.g.: in an empty program body.
 SUBROUTINE hac_test
  IMPLICIT NONE
  INTEGER :: loff, losz, hac_mb, reps ! local offset, size, MiB of local size for stress tester, repetitions
  INTEGER :: ib ! Index Base
  INTEGER :: erank, ecomm ! rank and communicator
  HAC_NTYPE,HAC_DSPEC,ALLOCATABLE :: GAS,GCP !< Global Array Slice and copy
#if HLST_HAC_USE_REAL
  INTEGER,PARAMETER :: HAC_1MB = 1024*1024/(1*HAC_KIND)
#else
  INTEGER,PARAMETER :: HAC_1MB = 1024*1024/(2*HAC_KIND)
#endif
!
  CHARACTER(AMFNL_) :: fname
  INTEGER :: tdim, ld(HLST_HAC_USE_ADIM)
  CHARACTER :: ntc !< Numerical Type Character.
  INTEGER :: GADA(NDIMS) !< Global Array's Dimensions Array
  INTEGER :: i
  INTEGER :: LB(NDIMS)
  REAL(KIND=TIME_K) :: ttt
!
  loff = 0
  losz = 0
  ecomm = MPI_COMM_WORLD
  CALL MPI_Init (ierr)
  CALL MPI_Comm_rank (ecomm, erank, ierr)
  ttt = - MPI_Wtime()
!
  IF (erank .EQ. 0) THEN
   WRITE(OU,'(a)')'+--------------------------------------------------+'
   WRITE(OU,'(a)')'| MODULE hlst_adios_checkpoint test program        |'
   WRITE(OU,'(a)')'+--------------------------------------------------+'
   WRITE(OU,'(a)')'| # You can use the HAC_MB env. var.:              |'
   WRITE(OU,'(a)')'| export HAC_MB=1  # if>0,  size of I/O tests      |'
   WRITE(OU,'(a)')'| export HAC_MB=-1 # if<=1, auto choice            |'
   WRITE(OU,'(a)')'| export HAC_ST=0 # stress test; default is 1      |'
   WRITE(OU,'(a)')'| export HAC_ET=0 # error test; default is 1       |'
   WRITE(OU,'(a)')'| export HAC_VL=0 # verbosity level; default is 1  |'
   WRITE(OU,'(a)')'| export HAC_CHECKPOINT_FILENAME="checkpoint.dat"  |'
   WRITE(OU,'(a)')'| # you can increase it, negative values will      |'
   WRITE(OU,'(a)')'| # request extra (debug) output.                  |'
   WRITE(OU,'(a)')'| # Or also:                                       |'
   WRITE(OU,'(a)')'| export HAC_BENCH=1.. # Default benchmark         |'
   WRITE(OU,'(a)')'| export HAC_BENCH=3.. # Benchmark with 3 MiB..    |'
   WRITE(OU,'(a)')'| # And additionally (ADIOS parameters--INTERNALS):|'
   WRITE(OU,'(a)')'| export HAC_SAMPLES=... #  (for HAC_BENCH>1)      |'
   WRITE(OU,'(a)')'| export HAC_MEFIBFL=... #  (for HAC_BENCH>1)      |'
   WRITE(OU,'(a)')'| export HAC_MDIMSTDC=... #                        |'
   WRITE(OU,'(a)')'| export HAC_ODPT=... #                            |'
   WRITE(OU,'(a)')'| export HAC_A_NA=... #                            |'
   WRITE(OU,'(a)')'| export HAC_A_NO=... #                            |'
   WRITE(OU,'(a)')'| export HAC_A_SC=... #                            |'
   WRITE(OU,'(a)')'| export HAC_A_SS=... #                            |'
   WRITE(OU,'(a)')'| export HAC_A_BS=... #                            |'
   WRITE(OU,'(a)')'+--------------------------------------------------+'
  END IF
!
  CALL get_env_opt_s("HAC_CHECKPOINT_FILENAME",fname,"checkpoint.dat")
  CALL hac_init (ecomm, get_env_opt_i("HAC_VL",HAC_VERBOSE), ierr)
!
  CALL hac_defv( &
         & A_NA_ = get_env_opt_i("HAC_A_NA",A_NA),&
         & A_NO_ = get_env_opt_i("HAC_A_NO",A_NO),&
         & A_SC_ = get_env_opt_i("HAC_A_SC",A_SC),&
         & A_SS_ = get_env_opt_i("HAC_A_SS",A_SS),&
         & A_BS_ = get_env_opt_i("HAC_A_BS",A_BS),&
         & ierr_=ierr)
!
  CALL hac_herr
  ! HLST_HAC_GOTO_END_ON_ERROR(9999)
!
  IF ( get_env_opt_i("HAC_BENCH",0) .NE. 0 ) THEN
   reps = get_env_opt_i("HAC_SAMPLES",1)
   fime = get_env_opt_i("HAC_MEFIBFL",0)
   IF ( get_env_opt_i("HAC_BENCH",0) .EQ. 1 ) THEN
    CALL hac_bench(ecomm,fname,reps_=reps)
   ELSE
    CALL hac_bench(ecomm,fname,&
            &mibpt_=get_env_opt_i("HAC_BENCH",0),reps_=reps)
   END IF
   GOTO 9999
  END IF
!
  hac_mb = get_env_opt_i("HAC_MB",1)
  IF (hac_mb .LE. 0) THEN
   hac_mb = 1 ! Start here to improve the 'auto choice' option.
  END IF
  losz = hac_mb * HAC_1MB
  IF ( erank .EQ. 0 ) &
   & WRITE(OU,'(a,i0)') BNR//'BEGIN test with HAC_MB = ',hac_mb
  DO ib = 0, 1
  DO tdim = 2, HLST_HAC_USE_ADIM
   ld = 1 + ib - 1
   loff = erank + ib
   ld(tdim) = losz + ib - 1 !
   ld(1) = loff             !
!
   LB = (/ld(1),(ib,i=2,NDIMS)/)
   CALL hac_laa(GAS,LB,ld,erank,.TRUE.)
   CALL hac_laa(GCP,LB,ld,erank,.TRUE.)
!
   GCP = GAS
   CALL hac_write (GAS, LBOUND(GAS), TRIM(fname), st_=1.0_8, ts_=1.0_8) !
   GAS = 0
   CALL hac_info (fname, ierr, GADA_=GADA )
   ! IF (erank .EQ. 0) WRITE(OU,*)'GADA ',GADA
   IF(erank.EQ.0) THEN
     CALL hac_wrt_sda(BNR//'Global Array Dimensions:',GADA,'')
     CALL hac_wrt_sda(BNR//'Global Array Local Bounds: ',LBOUND(GAS),'')
   END IF
   CALL hac_info (fname, ierr, ntc_=ntc)
   IF (erank.EQ.0)WRITE(OU,'(a,a)')BNR//'Numerical Type Character: ',NTC
   CALL hac_info (fname, ierr)
   CALL hac_read (GAS, LBOUND(GAS), TRIM(fname))
   CALL hac_info (fname, ierr)
   CALL hac_read (GAS, LBOUND(GAS), TRIM(fname), gaib_=ib)
   IF (SUM(GAS-GCP).NE.0.0) STOP "Critical error in restart !"
   DEALLOCATE(GAS)
   DEALLOCATE(GCP)
  END DO
  END DO
!
  IF ( get_env_opt_i("HAC_ST",1) .NE. 0 ) THEN
   CALL hac_stress(ecomm,fname)
  END IF
!
  IF ( get_env_opt_i("HAC_ET",1) .NE. 0 ) THEN
   IF(erank.EQ.0) WRITE(OU,'(a)')BNR//'BEGIN Error test.'
   CALL hac_ex_ckp_err_test
   IF(erank.EQ.0) WRITE(OU,'(a)')BNR//'END Error test.'
  END IF
!
9999 ttt = ttt + MPI_Wtime()
  IF(erank.EQ.0) WRITE(OU,'(a,f10.2,a)')BNR//'Done (took ',ttt,' s).'
  CALL hac_exit(ierr)
  CALL MPI_Finalize (ierr)
 END SUBROUTINE hac_test
!
#undef HAC_FIELD
#undef HAC_KSPEC
#undef HAC_NTYPE
#undef HLST_HAC_USE_ADIM
#undef HLST_HAC_USE_CMSK
!
  SUBROUTINE hac_sync(fname, ierr_)
   USE ISO_C_BINDING
   IMPLICIT NONE
   CHARACTER(LEN=*), TARGET,INTENT(IN) :: fname !< Path to the checkpoint file. Any length is allowed.
   INTEGER, OPTIONAL, INTENT(OUT) :: ierr_ !< Error status variable. If omitted, errors will be fatal.
!
   INTEGER, PARAMETER :: O_RDONLY=0 ! Do not count on this; check fcntl.h
   INTEGER, TARGET :: FD, cerr
!
  INTERFACE
   INTEGER(C_INT) FUNCTION C_FDATASYNC (FD) BIND(c,NAME = 'fdatasync')
   USE ISO_C_BINDING
   INTEGER(C_INT), VALUE :: FD
   END FUNCTION C_FDATASYNC
  END INTERFACE
!
  INTERFACE
   SUBROUTINE C_SYNC () BIND(c,NAME = 'sync')
   ! "sync() causes all buffered modifications to file metadata and data to be written to the underlying file systems."
   USE ISO_C_BINDING
   END SUBROUTINE C_SYNC
  END INTERFACE
!
  INTERFACE
   INTEGER(C_INT) FUNCTION C_SYNCFS (FD) BIND(c,NAME = 'syncfs') ! Linux Specific, _GNU_SOURCE
   ! "syncfs() is like sync(), but synchronizes just the file system containing file referred to by the open file descriptor fd."
   USE ISO_C_BINDING
   INTEGER(C_INT), VALUE :: FD
   END FUNCTION C_SYNCFS
  END INTERFACE
!
  INTERFACE
   INTEGER(C_INT) FUNCTION C_CLOSE (FD) BIND(c,NAME = 'close')
   USE ISO_C_BINDING
   INTEGER(C_INT), VALUE :: FD
   END FUNCTION C_CLOSE
  END INTERFACE
!
  INTERFACE
   INTEGER(C_INT) FUNCTION C_OPEN (PATHNAME,FLAGS) BIND(c,NAME = 'open')
   USE ISO_C_BINDING
   !CHARACTER(C_CHAR), DIMENSION(*) :: PATHNAME
   CHARACTER(C_CHAR), INTENT(IN) :: PATHNAME(*)
   !TYPE(C_PTR),VALUE :: PATHNAME
   !CHARACTER(C_CHAR), TARGET :: PATHNAME
   INTEGER(C_INT), VALUE :: FLAGS
   END FUNCTION C_OPEN
  END INTERFACE
!
  INTERFACE
   SUBROUTINE C_PRINTF (STR) BIND(c,NAME = 'printf')
   USE ISO_C_BINDING
   CHARACTER(C_CHAR), INTENT(IN) :: STR(*)
   END SUBROUTINE C_PRINTF
  END INTERFACE
!
   FD = 0
   cerr = 0
!
!   CALL MPI_Barrier (acomm, cerr)
!
   IF (.TRUE.) THEN
    ! The following is not guaranteed to do the job.
    CALL C_SYNC ()
    GOTO 99
   END IF
!
   IF (arank .EQ. RR) THEN
     FD = C_OPEN(TRIM(fname)//C_NULL_CHAR,O_RDONLY)
    IF ( FD .eq. -1 ) THEN
     WRITE(EU,"(a)")BNR//"C_OPEN error"
     cerr = FD
    ELSE
     ! The following is of little use: this file is only meta data.
     cerr = C_FDATASYNC(FD)
     IF ( cerr /= 0 ) WRITE(EU,"(a)")BNR//"C_FDATASYNC error"
     !
     ! The following is not ifort friendly.
     ! cerr = C_SYNCFS (FD)
     ! IF ( cerr /= 0 ) WRITE(EU,"(a)")BNR//"C_SYNCFS error"
     !
     cerr = C_CLOSE(FD)
     IF ( cerr /= 0 ) WRITE(EU,"(a)")BNR//"C_CLOSE error"
    END IF
    IF ( cerr .NE. 0 ) THEN
     WRITE(EU,"(a,a)")BNR//"synchronization failed on file ",fname
    ELSE
!     WRITE(EU,"(a,a)")BNR//"synchronization success on file ",fname
    END IF
   END IF
99 HAC_RET_IF_PRESENT(ierr_,cerr)
!
  END SUBROUTINE hac_sync
!
  SUBROUTINE hac_flush_attempt()
   ! Note: this subroutine is experimental.
   ! To improve it, one such test per node only should be performed.
   ! But writing code for this is quite tricky.
   ! And in any case, this test is very slow.
   CALL MPI_Barrier (acomm,ierr)
   CALL hac_max_alloc_test()
   ierr = OKV_
  END SUBROUTINE hac_flush_attempt
!
  SUBROUTINE hac_max_alloc_test()
   INTEGER(KIND=8) :: se,sa
   INTEGER :: ret
!
   DO se = 0,16
    sa = 2**se
    ! WRITE(OU,"(a,i0,a)")BNR//"allocating ",sa," MiB"
    IF ( hac_alloc_test(sa) .NE. 0 ) THEN
     !IF (DODEBUG_.NE.0.AND.amroot) THEN
     IF (amroot) THEN
      WRITE(OU,"(a,i0,a)")BNR//&
      &"Stopping flush-triggering allocation attempts after ",sa," MiB"
     END IF
     !END IF
     GOTO 99
    END IF
   END DO
99 ret = hac_alloc_test((sa/2)*3) ! up to 3/4 of the max
   RETURN
  END SUBROUTINE hac_max_alloc_test
!
  INTEGER FUNCTION hac_alloc_test(mta)
   USE ISO_C_BINDING
   IMPLICIT NONE
   INTEGER(KIND=8),INTENT(IN) :: mta ! MiB to alloc
!
   INTEGER(KIND=8) :: i
   INTEGER :: res
   REAL :: r
   INTEGER(C_INT),DIMENSION(:),ALLOCATABLE :: IA ! integer array
!
   ALLOCATE(IA(MAX(mta*(Mib_/4),INT8(1))),stat=res)
   IF ( res .EQ. 0 ) THEN
    DO i = 1, SIZE(IA)
     ! We only care for the following loop to exist in object code.
     ! If not, the optimizer will swallow the whole ALLOCATE ...
     ! CALL RANDOM_NUMBER(r)
     r = SQRT(REAL(i))
     IA(i) = INT(i)-1+INT(r)
    END DO
    DEALLOCATE(IA)
   END IF
!
   hac_alloc_test = res
  END FUNCTION hac_alloc_test
!
END MODULE hlst_adios_checkpoint
!
