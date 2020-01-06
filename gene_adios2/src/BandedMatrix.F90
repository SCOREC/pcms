#include "intrinsic_sizes.h"
#include "redef.h"
#undef DEBUG
#ifdef DEBUG
#define DEBSTART(string) print*,"== START ",string," =="
#define DEBEND(string) print*,"== END   ",string," =="
#else
#define DEBSTART(string)
#define DEBEND(string)
#endif

!> This is a frontend module for the usage of a BandedMatrix. It more or less has the interfaces to the
!! underlying module, which does the real storage. At the moment the underlying module is
!! StoreBandedMatrixModule, which handles the storage.
MODULE BandedMatrixModule
  USE ProcessGridModule
  USE StoreBandedMatrixModule
  USE MatrixModule
  use VectorModule
  !USE mpi
  IMPLICIT NONE
  PRIVATE

  public :: my_sum,static_size_of_BandedMatrix
  public :: initialize_BandedMatrix_module,finalize_BandedMatrix_module

!> Just knows the mathematical extent of the matrix, that is the
!! number of rows and columns. It also has a data object embedded
!! which contains the real data in form of a StoreBandedMatrixObject.
  TYPE,public ::  BandedMatrix
     INTEGER :: NCols
     INTEGER :: NRows
     type(StoreBandedMatrixObject) :: data
   contains

     PROCEDURE,private :: bm_initialize_matrix
     PROCEDURE,private :: bm_initialize_vector
     generic,public :: initialize => bm_initialize_matrix, bm_initialize_vector

     final :: bm_finalize_matrix

     PROCEDURE :: allocate => bm_allocate
     PROCEDURE,private :: bm_set_complex_value
     PROCEDURE,private :: bm_set_real_value
     generic,public :: set_value => bm_set_complex_value, bm_set_real_value
     PROCEDURE,private ::  bm_add_complex_value
     PROCEDURE,private ::  bm_add_real_value
     generic,public :: add_value => bm_add_complex_value, bm_add_real_value
     procedure :: commit_values => bm_commit_values
     PROCEDURE :: set_zero => bm_set_zero
     procedure :: get_global_matrix_locally => bm_get_global_matrix_locally
     !procedure :: convert_banded_to_full => bm_convert_Banded_to_Full
     !procedure :: convert_Full_to_banded => bm_convert_Full_to_Banded
     procedure :: get_local_abs_square_sum => bm_get_local_abs_square_sum
     PROCEDURE :: mat_get_value => bm_get_value
     PROCEDURE :: mat_get_row_pointer => bm_get_row_pointer
     procedure :: autotune => bm_autotune
     PROCEDURE,private :: bm_dot_multiply
     PROCEDURE,private :: bm_matmat_dot_multiply
     procedure,private :: bm_dot_multiply_Banded_with_Full
     generic,public :: dot_multiply => bm_dot_multiply, bm_matmat_dot_multiply,&
         &bm_dot_multiply_Banded_with_Full
     procedure :: square_matrix => bm_square_matrix
     PROCEDURE,private ::  bm_show_on_screen
     PROCEDURE,private ::  bm_show_in_file
     generic,public :: show => bm_show_on_screen, bm_show_in_file
     PROCEDURE,private ::  bm_add_to_matrix
     PROCEDURE,private ::  bm_add_matrix
     generic,public :: add_matrix => bm_add_to_matrix, bm_add_matrix
     PROCEDURE,private ::  bm_subtract_matrix
     PROCEDURE,private ::  bm_subtract_from_matrix
     generic,public :: subtract_matrix =>bm_subtract_matrix, bm_subtract_from_matrix
     procedure,private :: bm_multiply_matrix_with_real
     procedure,private :: bm_scale_matrix_by_real
     generic,public :: multiply_matrix_with_scalar => bm_multiply_matrix_with_real,&
          &bm_scale_matrix_by_real
     procedure,private :: bm_assign_matrix
     generic,public :: assignment(=) => bm_assign_matrix
     procedure :: output_data => bm_output_data_matrix
     procedure :: isInitialized => bm_isInitialized
     PROCEDURE :: print_storage_details => bm_print_storage_details
     procedure :: get_number_of_bands =>bm_get_number_of_bands

     procedure :: row_axpy => bm_row_axpy
     procedure :: transpose_and_conjugate => bm_transpose_and_conjugate
     procedure :: transpose_storage => bm_transpose_storage
     procedure :: LU_factor => bm_LU_factor
     procedure :: LU_factor_ok => bm_LU_factor_ok
     procedure :: LU_solve => bm_LU_solve
  END TYPE BandedMatrix

  INTERFACE initialize_BandedMatrix_module
    MODULE PROCEDURE bm_initialize_BandedMatrix_module
  END INTERFACE

  INTERFACE finalize_BandedMatrix_module
     MODULE PROCEDURE bm_finalize_BandedMatrix_module
  END INTERFACE

  !interface BandedMatrix
  !   module procedure bm_initialize_matrix,bm_initialize_vector
  !end interface BandedMatrix

  interface my_sum
     module procedure bm_sum_0D, bm_sum_2D, bm_sum_3D
  end interface my_sum



  ! Use the same process grid for all Matrices instantiated whithin this module
  TYPE(ProcessGrid),SAVE :: PG

CONTAINS
  FUNCTION static_size_of_BandedMatrix() RESULT(memory_need)
    integer :: memory_need

    memory_need = 2*SIZE_OF_INTEGER &
         & + static_size_of_storeBandedmatrixobject()

  END FUNCTION static_size_of_BandedMatrix

  !function bm_initialize_matrix(n_rows,n_cols,transp) result(mat)
  !type(BandedMatrix) :: mat
  subroutine bm_initialize_matrix(mat,n_rows,n_cols,transp)
    Class(BandedMatrix), INTENT(INOUT) :: mat
    INTEGER :: n_rows, n_cols
    LOGICAL, OPTIONAL :: transp

    DEBSTART("bm_initialize_matrix")
    IF (PG%isInitialized) THEN
       mat%NRows = n_rows
       mat%NCols = n_cols
       !if (.not.associated(mat%Data)) then
       !   print*,"initialize: mat%Data not allocated."
       !   allocate(mat%Data)
       !end if
       IF (PRESENT(transp)) THEN
          CALL initialize(mat%Data,n_rows,n_cols, PG,transp)
       ELSE
          CALL initialize(mat%Data,n_rows,n_cols, PG)
       END IF
       call allocate(mat%Data)
    ELSE
       PRINT*,"Before initializing a matrix, the process grid has to be initialized."
       STOP
    END IF
    !if (.not.associated(mat%Data)) then
    !   print*,"initialize: mat%Data immer noch nicht alloziert."
    !end if
    DEBEND("bm_initialize_matrix")
  end subroutine bm_initialize_matrix

  subroutine bm_initialize_vector(this,n_rows)
    class(BandedMatrix) :: this
    integer :: n_rows
    !type(BandedMatrix) :: vec

    call this%initialize(n_rows,1)
    !vec = BandedMatrix(n_rows,1)
  end subroutine bm_initialize_vector

  SUBROUTINE bm_initialize_BandedMatrix_module
    CALL get_processgrid(PG)
  END SUBROUTINE bm_initialize_BandedMatrix_module

  LOGICAL FUNCTION bm_isInitialized(mat)
    class(BandedMatrix) :: mat

    bm_isInitialized = isInitialized(mat%Data)
  END FUNCTION bm_isInitialized

  SUBROUTINE bm_allocate(mat)
    Class(BandedMatrix), intent(INOUT) :: mat

    CALL allocate(mat%Data)
  END SUBROUTINE bm_allocate

  SUBROUTINE bm_finalize_matrix(mat)
    type(BandedMatrix) :: mat
    call finalize(mat%Data)
    !deallocate(mat%Data)
  END SUBROUTINE bm_finalize_matrix

  SUBROUTINE bm_finalize_BandedMatrix_module()
    !CALL finalize(PG)
  END SUBROUTINE bm_finalize_BandedMatrix_module

  SUBROUTINE bm_set_complex_value(mat,irow,icol,value)
    Class(BandedMatrix),INTENT(INOUT) :: mat
    INTEGER,INTENT(IN) :: irow,icol
    COMPLEX, INTENT(IN) :: value

    !PRINT*,"In routine set_value of MatrixModule."
    CALL set_value(mat%Data,irow,icol,value)
  END SUBROUTINE bm_set_complex_value

  SUBROUTINE bm_set_real_value(mat,irow,icol,value)
    Class(BandedMatrix),INTENT(INOUT) :: mat
    INTEGER,INTENT(IN) :: irow,icol
    REAL, INTENT(IN) :: value

    !PRINT*,"In routine set_value of MatrixModule."
    CALL set_value(mat%Data,irow,icol,CMPLX(value,0.0,KIND(value)))
  END SUBROUTINE bm_set_real_value

  SUBROUTINE bm_add_complex_value(mat,irow,icol,value)
    Class(BandedMatrix),INTENT(INOUT) :: mat
    INTEGER,INTENT(IN) :: irow,icol
    COMPLEX, INTENT(IN) :: value

    CALL add_value(mat%Data,irow,icol,value)
  END SUBROUTINE bm_add_complex_value

  SUBROUTINE bm_add_real_value(mat,irow,icol,value)
    Class(BandedMatrix),INTENT(INOUT) :: mat
    INTEGER,INTENT(IN) :: irow,icol
    REAL, INTENT(IN) :: value

    CALL add_value(mat%Data,irow,icol,CMPLX(value,0.0,KIND(value)))
  END SUBROUTINE bm_add_real_value


  SUBROUTINE bm_commit_values(mat)
    Class(BandedMatrix) :: mat

    CALL commit_values(mat%Data)
  END SUBROUTINE bm_commit_values

  subroutine bm_autotune(mat,vec,res)
    Class(BandedMatrix),intent(IN) :: mat
    Type(Vector),intent(IN) :: vec
    Type(Vector),intent(INOUT) :: res

    IF ((mat%NCols.EQ.vec%NRows).AND. &
         &(mat%NRows.EQ.res%NRows)) THEN
       call autotune(mat%Data,vec%Data,res%Data)
    ELSE
       PRINT*,"BandedMatrix shapes for multiplication with Vector do not match, ",&
            &mat%NCols,mat%NRows,vec%NRows,res%NRows
       stop
    END IF

  end subroutine bm_autotune

  SUBROUTINE bm_dot_multiply(mat,vec, res)
    Class(BandedMatrix),intent(IN) :: mat
    Type(Vector),intent(IN) :: vec
    Type(Vector),intent(INOUT) :: res

    !PRINT*,vec%NCols,vec%NRows,res%NCols,res%NRows
    IF ((mat%NCols.EQ.vec%NRows).AND. &
         &(mat%NRows.EQ.res%NRows)) THEN
       CALL dot_multiply(mat%Data,vec%Data,res%Data)
    ELSE
       PRINT*,"BandedMatrix shapes for multiplication with Vector do not match, ",&
            &mat%NCols,mat%NRows,vec%NRows,res%NRows
       stop
    END IF
  END SUBROUTINE bm_dot_multiply

  subroutine bm_matmat_dot_multiply(mat,mat2,res,transA_in,transB_in)
    class(BandedMatrix),intent(IN) :: mat, mat2
    class(BandedMatrix),intent(INOUT) :: res
    CHARACTER, intent(IN), optional :: transA_in, transB_in

    IF ((mat%NCols.EQ.mat2%NRows).AND. &
         &(mat%NRows.EQ.res%NRows)) THEN
       if (.not.(present(transA_in).or.present(transB_in))) then
          CALL dot_multiply(mat%Data,mat2%Data,res%Data)
       elseif (present(transA_in).and.present(transB_in)) then
          CALL dot_multiply(mat%Data,mat2%Data,res%Data,transA_in,transB_in)
       elseif (present(transA_in)) then
          CALL dot_multiply(mat%Data,mat2%Data,res%Data,transA_in)
       endif
    ELSE
       PRINT*,"BandedMatrix shapes for multiplication with Matrix do not match, ",&
            &mat%NCols,mat%NRows,mat2%NRows,res%NRows
       stop
    END IF
  end subroutine bm_matmat_dot_multiply

  subroutine bm_dot_multiply_Banded_with_Full(bmat,mat,res,transA_in,transB_in)
    class(BandedMatrix),intent(IN) :: bmat
    type(Matrix), intent(IN) :: mat
    type(Matrix),intent(INOUT) :: res
    CHARACTER, intent(IN), optional :: transA_in, transB_in

    IF ((bmat%NCols.EQ.mat%NRows).AND. &
         &(bmat%NRows.EQ.res%NRows)) THEN
       if (.not.(present(transA_in).or.present(transB_in))) then
          CALL dot_multiply(bmat%Data,mat%Data,res%Data)
       elseif (present(transA_in).and.present(transB_in)) then
          CALL dot_multiply(bmat%Data,mat%Data,res%Data,transA_in,transB_in)
       elseif (present(transA_in)) then
          CALL dot_multiply(bmat%Data,mat%Data,res%Data,transA_in)
       endif
    ELSE
       PRINT*,"BandedMatrix shapes for multiplication with Matrix do not match, ",&
            &bmat%NCols,bmat%NRows,mat%NRows,res%NRows
       stop
    END IF
  end subroutine bm_dot_multiply_Banded_with_Full

  SUBROUTINE bm_square_matrix(bmat,res_bmat)
    Class(BandedMatrix),intent(IN) :: bmat
    Class(BandedMatrix),intent(INOUT) :: res_bmat

    CALL square_matrix(bmat%Data,res_bmat%Data)
  END SUBROUTINE bm_square_matrix

  SUBROUTINE bm_set_zero(mat)
    class(BandedMatrix) :: mat

    DEBSTART("bm_set_zero")
    call set_zero(mat%Data)
    DEBEND("bm_set_zero")
  END SUBROUTINE bm_set_zero

  SUBROUTINE bm_get_global_matrix_locally(mat,localfullmat)
    Class(BandedMatrix),INTENT(IN) :: mat
    COMPLEX, DIMENSION(:,:),INTENT(OUT) :: localfullmat

    CALL get_global_matrix_locally(mat%Data,localfullmat)
  END SUBROUTINE bm_get_global_matrix_locally

#ifdef USE_OLD_CONVERT_ROUTINES
  SUBROUTINE bm_convert_Banded_to_Full(bmat, fmat)
    Class(BandedMatrix), intent(IN) :: bmat
    TYPE(Matrix), intent(INOUT) :: fmat

    INTEGER :: irow, icol, RowsPerRank

    IF ((fmat%NRows.NE.bmat%NRows).OR.(fmat%NCols.NE.bmat%NCols)) THEN
       PRINT*,"Number of Rows and Number of Cols must be equal for full and banded matrix."
       stop
    END IF

    RowsPerRank = bmat%NRows/PG%NProcs

    IF (bmat%data%isTransposed) THEN
       DO irow=PG%Rank*RowsPerRank+1,(PG%Rank+1)*RowsPerRank
          !DO irow=1,bmat%NRows
          DO icol=1,bmat%NCols
             CALL set_value(fmat,irow,icol,mat_get_value(bmat,irow,icol))
          END DO
       END DO
    ELSE
       PRINT*, "Can only convert banded to full matrix if banded matrix is in transposed storage mode."
       stop
    END IF
  END SUBROUTINE bm_convert_Banded_to_Full

  SUBROUTINE bm_convert_Full_to_Banded(fmat, bmat)
    TYPE(Matrix), intent(IN) :: fmat
    Class(BandedMatrix), intent(INOUT) :: bmat

    INTEGER :: irow, icol, RowsPerRank

    IF ((fmat%NRows.NE.bmat%NRows).OR.(fmat%NCols.NE.bmat%NCols)) THEN
       PRINT*,"Number of Rows and Number of Cols must be equal for full and banded matrix."
       stop
    END IF

    call set_zero(bmat)
    RowsPerRank = bmat%NRows/PG%NProcs

    IF (bmat%data%isTransposed) THEN
       DO irow=PG%Rank*RowsPerRank+1,(PG%Rank+1)*RowsPerRank
          DO icol=1,bmat%NCols
             CALL set_value(bmat,irow,icol,mat_get_value(fmat,irow,icol))
          END DO
       END DO
       CALL commit_values(bmat)
    ELSE
       PRINT*, "Can only convert full to banded matrix if banded matrix is in transposed storage mode."
       stop
    END IF
  END SUBROUTINE bm_convert_Full_to_Banded
#endif

  FUNCTION bm_get_value(mat,irow,icol,trans_in) result(value)
    class(BandedMatrix) :: mat
    INTEGER :: irow, icol
    complex :: value
    character(len=1),optional,intent(IN) :: trans_in

    if (present(trans_in)) then
       value = mat_get_value(mat%Data,irow,icol,trans_in)
    else
       value = mat_get_value(mat%Data,irow,icol,'N')
    end if
  END FUNCTION bm_get_value

  subroutine bm_get_row_pointer(mat,global_row,ptr_row,s_global_col,e_global_col)
    Class(BandedMatrix) :: mat
    INTEGER,intent(IN) :: global_row
    complex,dimension(:),pointer,intent(OUT) :: ptr_row
    integer, INTENT(OUT) :: s_global_col, e_global_col

    call mat_get_row_pointer(mat%Data,global_row,ptr_row,s_global_col,e_global_col)
  end subroutine bm_get_row_pointer

  FUNCTION bm_get_local_abs_square_sum(mat) RESULT(value)
    class(BandedMatrix) :: mat
    real :: value

    value = get_local_abs_square_sum(mat%Data)
  END FUNCTION bm_get_local_abs_square_sum

  SUBROUTINE bm_show_on_screen(mat)
    class(BandedMatrix) :: mat

    CALL show(mat%Data)
  END SUBROUTINE bm_show_on_screen

  SUBROUTINE bm_show_in_file(mat,filename)
    class(BandedMatrix) :: mat
    CHARACTER(len=*) :: filename

    CALL show(mat%Data,filename)
  END SUBROUTINE bm_show_in_file

  SUBROUTINE bm_output_data_matrix(this,basename)
    class(BandedMatrix) :: this
    CHARACTER(len=*) :: basename

    character(len=FILENAME_MAX) :: filename

    WRITE(filename,"(3A)") "./",TRIM(basename),".dat"
    PRINT*,"Writing to file ",TRIM(filename)
    CALL this%show(filename)
  END SUBROUTINE bm_output_data_matrix

  ! =======================================
  ! == Define some operators on matrices ==
  ! =======================================
  SUBROUTINE bm_add_to_matrix(this,mat)
    Class(BandedMatrix), intent(INOUT) :: this
    Class(BandedMatrix),INTENT(IN) :: mat
    CALL add_matrix(this%Data,mat%Data)
  END SUBROUTINE bm_add_to_matrix

  SUBROUTINE bm_add_matrix(this,mat,res)
    Class(BandedMatrix), INTENT(IN) :: this,mat
    Class(BandedMatrix), INTENT(INOUT) :: res
    !PRINT*,"in bm_add_matrix of MatrixModule"
    CALL add_matrix(this%Data,mat%Data,res%Data)
  END SUBROUTINE bm_add_matrix

  SUBROUTINE bm_subtract_matrix(this,mat,res)
    Class(BandedMatrix),INTENT(IN) :: this,mat
    Class(BandedMatrix),intent(INOUT) :: res
    ! calculate res = this - mat

    CALL subtract_matrix(this%Data,mat%Data, res%Data)
  END SUBROUTINE bm_subtract_matrix

  SUBROUTINE bm_subtract_from_matrix(this,mat)
    Class(BandedMatrix),INTENT(INOUT) :: this
    Class(BandedMatrix), intent(IN) :: mat
    ! calculate this = this - mat

    CALL subtract_matrix(this%Data,mat%Data)
  END SUBROUTINE bm_subtract_from_matrix

  SUBROUTINE bm_multiply_matrix_with_real(this, scalar, res)
    Class(BandedMatrix), intent(IN) :: this
    REAL, intent(IN) :: scalar
    Class(BandedMatrix),intent(INOUT) :: res

    CALL multiply_matrix_with_scalar(this%Data,scalar, res%Data)
  END SUBROUTINE bm_multiply_matrix_with_real

  SUBROUTINE bm_scale_matrix_by_real(this, scalar)
    Class(BandedMatrix), intent(INOUT) :: this
    REAL, intent(IN) :: scalar

    CALL multiply_matrix_with_scalar(this%Data,scalar)
  END SUBROUTINE bm_scale_matrix_by_real

  SUBROUTINE bm_assign_matrix(lmat,rmat)
    Class(BandedMatrix), intent(INOUT) :: lmat
    Type(BandedMatrix), intent(IN) :: rmat

    lmat%Nrows = rmat%Nrows
    lmat%Ncols = rmat%Ncols

    lmat%Data = rmat%Data
       !allocate(lmat%data,source=rmat%Data)
    !   PRINT*, "Left and right matrices in an assignment must have the same dimensions."
    !   write(*,"(4(A,I4))") "lmat,rmat: ",lmat%Nrows,"x",lmat%Ncols," .neq. ",rmat%NRows,"x",rmat%NCols
  END SUBROUTINE bm_assign_matrix

  SUBROUTINE bm_print_storage_details(bmat)
    Class(BandedMatrix), INTENT(IN) :: bmat

    CALL print_storage_details(bmat%Data)
  END SUBROUTINE bm_print_storage_details

  FUNCTION bm_get_number_of_bands(mat)
    Class(BandedMatrix), intent(IN) :: mat
    integer :: bm_get_number_of_bands

    bm_get_number_of_bands = get_number_of_bands(mat%Data)
  END FUNCTION bm_get_number_of_bands

  ! ===================================
  ! == Define communication routines ==
  ! ===================================
  SUBROUTINE bm_sum_0D(localsum,communicator)
!    INTEGER, intent(IN) :: len
    Class(BandedMatrix), INTENT(INOUT) :: localsum
    integer :: communicator

!    IF (len.NE.1) THEN
!       PRINT*,"BandedMatrix:bm_sum_0D: A zero-dimensional array of banded matrices must have a len of 1."
!       stop
!    END IF

    CALL my_sum(localsum%Data,communicator)
  END SUBROUTINE bm_sum_0D

  SUBROUTINE bm_sum_2D(localsum,communicator)
    Class(BandedMatrix), DIMENSION(:,:), INTENT(INOUT) :: localsum
    integer :: communicator

    integer :: i,j

    do j=1,size(localsum,2)
       do i=1,size(localsum,1)
          call my_sum(localsum(i,j)%Data,communicator)
       end do
    end do
    !CALL bm_sum_generic(localsum,len,communicator)
  END SUBROUTINE bm_sum_2D

  SUBROUTINE bm_sum_3D(localsum,communicator)
    Class(BandedMatrix), DIMENSION(:,:,:), INTENT(INOUT) :: localsum
    integer :: communicator

    integer :: k

    do k=1,size(localsum,3)
       call my_sum(localsum(:,:,k),communicator)
    end do
    !CALL bm_sum_generic(localsum,len,communicator)
  END SUBROUTINE bm_sum_3D

#if 0
  SUBROUTINE bm_sum_generic(localsum,len,communicator)
    INTEGER, intent(IN) :: len
    Class(BandedMatrix), dimension(*),intent(INOUT) :: localsum
    integer :: communicator

    integer :: iMatrix
    DO iMatrix=1,len
       CALL my_sum(localsum(iMatrix)%Data,communicator)
    END DO
  END SUBROUTINE bm_sum_generic
#endif

  SUBROUTINE bm_row_axpy(this,this_row,mat,scalar)
    Class(BandedMatrix),intent(INOUT) :: this
    INTEGER, INTENT(IN) :: this_row
    Class(BandedMatrix), intent(IN) :: mat
    COMPLEX, intent(IN) :: scalar

    IF ((mat%Nrows.NE.this%NRows).OR.(mat%NCols.NE.this%NCols)) THEN
       PRINT*,"row_axpy: mat and this must have the same shape."
       STOP
    END IF

    IF ((this_row.GE.1).AND.(this_row.LE.this%NRows)) THEN
       CALL row_axpy(this%data,this_row,mat%data,scalar)
    ELSE
       PRINT*,"row_axpy: this_row must be between 1 and ",this%NRows," but is ",this_row
       STOP
    END IF
  END SUBROUTINE bm_row_axpy

  subroutine bm_transpose_and_conjugate(mat,dagger_mat)
    class(BandedMatrix) :: mat,dagger_mat

    call transpose_and_conjugate(mat%data,dagger_mat%data)
  end subroutine bm_transpose_and_conjugate

  subroutine bm_transpose_storage(mat,mat_t)
    class(BandedMatrix) :: mat,mat_t

    call transpose_storage(mat%data,mat_t%data)
  end subroutine bm_transpose_storage

  SUBROUTINE bm_LU_factor(mat)
    Class(BandedMatrix),intent(INOUT) :: mat

    CALL LU_factor(mat%DATA)

  END SUBROUTINE bm_LU_factor

  LOGICAL FUNCTION bm_LU_factor_ok(mat)
    Class(BandedMatrix),intent(IN) :: mat

    bm_LU_factor_ok = LU_factor_ok(mat%DATA)

  END FUNCTION bm_LU_factor_ok

  SUBROUTINE bm_LU_solve(mat,vec,res)
    Class(BandedMatrix),intent(IN) :: mat
    Type(Vector),intent(IN) :: vec
    Type(Vector),intent(INOUT) :: res

    CALL LU_solve(mat%DATA,vec%data,res%data)

  END SUBROUTINE bm_LU_solve

#ifdef ALL_ROUTINES
  ! ====================================================
  ! == Define some mathematical routines for matrices ==
  ! ====================================================
  SUBROUTINE bm_invert_matrix(mat)
    Class(BandedMatrix) :: mat

    CALL invert(mat%Data)
  END SUBROUTINE bm_invert_matrix
#endif

END MODULE BandedMatrixModule
