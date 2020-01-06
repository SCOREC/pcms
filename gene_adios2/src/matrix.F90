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

MODULE MatrixModule
  USE ProcessGridModule
  USE StoreFullMatrixModule
  !USE mpi
  IMPLICIT NONE
  ! setting default accessibility to private
  private

  public :: get_processgrid,static_size_of_matrix,my_sum
  public :: initialize_matrix_module, finalize_matrix_module

  TYPE,public :: Matrix
     INTEGER :: NCols
     INTEGER :: NRows

     class(StoreFullMatrixObject),allocatable :: data
   contains
     procedure,private :: initialize_matrix
     procedure,private :: initialize_vector
     generic,public :: initialize => initialize_matrix,initialize_vector
     procedure :: finalize => mp_finalize_matrix
     procedure :: allocate => mp_allocate

     procedure,private :: attach_vector
     procedure,private :: attach_matrix
     generic,public :: attach => attach_vector,attach_matrix

     procedure,private :: set_real_value
     procedure,private :: set_complex_value
     generic :: set_value => set_complex_value,set_real_value

     procedure,private :: add_real_value
     procedure,private :: add_complex_value
     generic,public :: add_value => add_real_value, add_complex_value

     procedure :: add_to_row => mat_add_to_row
     procedure :: commit_values => mp_commit_values
     procedure :: set_zero => mp_set_zero
     procedure :: get_global_matrix_locally => mp_get_global_matrix_locally
     procedure :: get_local_abs_square_sum => mp_get_local_abs_square_sum
     procedure :: mat_get_value => mp_get_value
     procedure :: dot_multiply => mp_dot_multiply_MatMat
     procedure,private :: show_on_screen => mp_show_on_screen
     procedure,private :: show_in_file => mp_show_in_file
     generic,public :: show => show_on_screen,show_in_file
     procedure,private :: add_to_matrix => mp_add_to_matrix
     procedure,private :: mp_add_matrix
     generic,public :: add_matrix => add_to_matrix,mp_add_matrix
     procedure,private :: mp_subtract_matrix
     procedure,private :: subtract_from_matrix
     generic,public :: subtract_matrix => mp_subtract_matrix,subtract_from_matrix
     procedure,private :: scale_matrix_by_real => mp_scale_matrix_by_real
     procedure,private :: multiply_matrix_with_real => mp_multiply_matrix_with_real
     generic,public :: multiply_matrix_with_scalar => multiply_matrix_with_real,scale_matrix_by_real
     procedure :: invert => mp_invert_matrix
     procedure :: isInitialized => mp_isInitialized
     procedure,private :: assign_matrix => mp_assign_matrix
     generic,public :: assignment(=) => assign_matrix
     procedure :: output_data => mp_output_data_matrix
     procedure :: LU_factor => mp_LU_factor
     procedure :: LU_solve => mp_LU_solve
  END TYPE matrix


  INTERFACE initialize_matrix_module
     MODULE PROCEDURE mp_initialize_matrix_module
  END INTERFACE

  INTERFACE finalize_matrix_module
     MODULE PROCEDURE mp_finalize_module
  END INTERFACE

  interface my_sum
     module PROCEDURE my_matrix_sum_2D, my_matrix_sum_3D, my_matrix_sum_0D
  end interface my_sum

  ! Use the same process grid for all Matrices instantiated whithin this module
  TYPE(ProcessGrid),SAVE :: PG

CONTAINS
  FUNCTION static_size_of_matrix() RESULT(memory_need)
    integer :: memory_need

    memory_need = 2*SIZE_OF_INTEGER &
         & + static_size_of_storefullmatrixobject()

  END FUNCTION static_size_of_matrix

  SUBROUTINE initialize_matrix(this,n_rows,n_cols)
    Class(Matrix), INTENT(INOUT) :: this
    INTEGER :: n_rows, n_cols

    DEBSTART("mp_initialize_matrix")
    IF (PG%isInitialized) THEN
       this%NRows = n_rows
       this%NCols = n_cols
       if (.not.allocated(this%Data)) allocate(this%Data)
       CALL initialize(this%Data,n_rows,n_cols, PG)
    ELSE
       PRINT*,"Before initializing a matrix, the process grid has to be initialized."
       STOP
    END IF
    DEBEND("mp_initialize_matrix")
  END SUBROUTINE initialize_matrix

  SUBROUTINE get_processgrid(ao_PG)
    type(ProcessGrid) :: ao_PG
    ao_PG = PG
  END SUBROUTINE get_processgrid

  SUBROUTINE initialize_vector(this,n_rows)
    class(matrix) :: this
    integer :: n_rows

    CALL this%initialize(n_rows,1)
  END SUBROUTINE initialize_vector

  SUBROUTINE mp_initialize_matrix_module(comm)
    INTEGER :: comm

    DEBSTART("mp_initialize_matrix_module")
    CALL initialize(PG,MAT_BLOCKS_OF_ROWS, comm)
    !CALL initialize(PG_row,MAT_BLOCKS_OF_ROWS, comm)
    DEBEND("mp_initialize_matrix_module")
  END SUBROUTINE mp_initialize_matrix_module

  LOGICAL FUNCTION mp_isInitialized(this)
    class(matrix) :: this

    mp_isInitialized = isInitialized(this%Data)
  END FUNCTION mp_isInitialized

  SUBROUTINE mp_allocate(this)
    Class(Matrix), intent(INOUT) :: this

    CALL allocate(this%Data)
  END SUBROUTINE mp_allocate

  SUBROUTINE attach_vector(this,localarray)
    Class(Matrix),INTENT(INOUT) :: this
    COMPLEX, DIMENSION(:),INTENT(IN),TARGET :: localarray

    DEBSTART("mp_attach_vector")
    CALL attach(this%Data,localarray)
    DEBEND("mp_attach_vector")
  END SUBROUTINE attach_vector

  SUBROUTINE attach_matrix(this,localarray)
    Class(Matrix),INTENT(INOUT) :: this
    COMPLEX, DIMENSION(:,:),INTENT(IN),TARGET :: localarray

    CALL attach(this%Data,localarray)
  END SUBROUTINE attach_matrix

  SUBROUTINE mp_finalize_matrix(this)
    class(matrix) :: this
    call finalize(this%Data)
    deallocate(this%Data)
  END SUBROUTINE mp_finalize_matrix

  SUBROUTINE mp_finalize_module()
    call finalize(PG)
  END SUBROUTINE mp_finalize_module

  SUBROUTINE set_complex_value(this,irow,icol,value)
    Class(Matrix),INTENT(INOUT) :: this
    INTEGER,INTENT(IN) :: irow,icol
    COMPLEX, INTENT(IN) :: value

    !PRINT*,"In routine set_value of MatrixModule."
    CALL set_value(this%Data,irow,icol,value)
  END SUBROUTINE set_complex_value

  SUBROUTINE set_real_value(this,irow,icol,value)
    Class(Matrix),INTENT(INOUT) :: this
    INTEGER,INTENT(IN) :: irow,icol
    REAL, INTENT(IN) :: value

    !PRINT*,"In routine set_value of MatrixModule."
    CALL set_value(this%Data,irow,icol,CMPLX(value,0.0))
  END SUBROUTINE set_real_value

  SUBROUTINE add_complex_value(this,irow,icol,value)
    Class(Matrix),INTENT(INOUT) :: this
    INTEGER,INTENT(IN) :: irow,icol
    COMPLEX, INTENT(IN) :: value

    CALL add_value(this%Data,irow,icol,value)
  END SUBROUTINE add_complex_value

  SUBROUTINE add_real_value(this,irow,icol,value)
    Class(Matrix),INTENT(INOUT) :: this
    INTEGER,INTENT(IN) :: irow,icol
    REAL, INTENT(IN) :: value

    CALL add_value(this%Data,irow,icol,CMPLX(value,0.0))
  END SUBROUTINE add_real_value

  SUBROUTINE mat_add_to_row(this,irow,whole_row)
    Class(Matrix),INTENT(INOUT) :: this
    INTEGER,INTENT(IN) :: irow
    COMPLEX, dimension(:),INTENT(IN) :: whole_row

    CALL add_to_row(this%Data,irow,whole_row)
  END SUBROUTINE mat_add_to_row


  SUBROUTINE mp_commit_values(this)
    Class(Matrix) :: this

    CALL commit_values(this%Data)
  END SUBROUTINE mp_commit_values

  SUBROUTINE mp_dot_multiply_MatMat(this,mat, resmat, transA_in, transB_in)
    Class(Matrix),intent(IN) :: this
    Class(Matrix),intent(IN) :: mat
    Class(Matrix),intent(INOUT) :: resmat
    CHARACTER, intent(IN), optional :: transA_in, transB_in

    !PRINT*,mat%NCols,mat%NRows,resmat%NCols,resmat%NRows
    IF ((this%NCols.EQ.mat%NRows).AND. &
         &(this%NRows.EQ.resmat%NRows).AND.&
         &(mat%NCols.EQ.resmat%NCols)) THEN
       if (.not.(present(transA_in).or.present(transB_in))) then
          CALL dot_multiply(this%Data,mat%Data,resmat%Data)
       elseif (present(transA_in).and.present(transB_in)) then
          CALL dot_multiply(this%Data,mat%Data,resmat%Data,transA_in,transB_in)
       elseif (present(transA_in)) then
          CALL dot_multiply(this%Data,mat%Data,resmat%Data,transA_in)
       endif
    ELSE
       PRINT*,"Matrix shapes for multiplication with Matrix do not match, ",&
            &this%NCols,this%NRows,mat%NCols,mat%NRows,resmat%NCols,resmat%NRows
       stop
    END IF
  END SUBROUTINE mp_dot_multiply_MatMat

  SUBROUTINE mp_set_zero(this)
    class(matrix) :: this

    ! Local variables
    INTEGER :: icol, irow

    DEBSTART("mp_set_zero")
    DO irow=1,this%NRows
       DO icol=1,this%NCols
          CALL this%set_value(irow,icol,CMPLX(0.0,0.0))
       END DO
    END DO
    CALL this%commit_values()
    DEBEND("mp_set_zero")
  END SUBROUTINE mp_set_zero

  SUBROUTINE mp_get_global_matrix_locally(mat,localfullmat)
    Class(Matrix),INTENT(IN) :: mat
    COMPLEX, DIMENSION(:,:),INTENT(OUT) :: localfullmat

    CALL get_global_matrix_locally(mat%Data,localfullmat)
  END SUBROUTINE mp_get_global_matrix_locally

  FUNCTION mp_get_value(mat,irow,icol) result(value)
    class(matrix) :: mat
    INTEGER :: irow, icol
    complex :: value

    value = mat_get_value(mat%Data,irow,icol)
  END FUNCTION mp_get_value

  FUNCTION mp_get_local_abs_square_sum(mat) RESULT(value)
    class(matrix) :: mat
    real :: value

    value = get_local_abs_square_sum(mat%Data)
  END FUNCTION mp_get_local_abs_square_sum

  SUBROUTINE mp_show_on_screen(mat)
    class(matrix) :: mat

    CALL show(mat%Data)
  END SUBROUTINE mp_show_on_screen

  SUBROUTINE mp_show_in_file(mat,filename)
    class(matrix) :: mat
    CHARACTER(len=*) :: filename

    CALL show(mat%Data,filename)
  END SUBROUTINE mp_show_in_file

  SUBROUTINE mp_output_data_matrix(this,basename)
    class(matrix) :: this
    CHARACTER(len=*) :: basename

    character(len=100) :: filename

    WRITE(filename,"(3A)") "./",TRIM(basename),".dat"
    PRINT*,"Writing to file ",TRIM(filename)
    CALL this%show(filename)
  END SUBROUTINE mp_output_data_matrix

  ! =======================================
  ! == Define some operators on matrices ==
  ! =======================================
  SUBROUTINE mp_add_to_matrix(this,mat)
    Class(Matrix), intent(INOUT) :: this
    Class(Matrix),INTENT(IN) :: mat
    CALL add_matrix(this%Data,mat%Data)
  END SUBROUTINE mp_add_to_matrix

  SUBROUTINE mp_add_matrix(this,mat,res)
    Class(Matrix), INTENT(IN) :: this,mat
    Class(Matrix), INTENT(INOUT) :: res
    !PRINT*,"in mp_add_matrix of MatrixModule"
    CALL add_matrix(this%Data,mat%Data,res%Data)
  END SUBROUTINE mp_add_matrix

  SUBROUTINE mp_subtract_matrix(this,mat,res)
    Class(Matrix),INTENT(IN) :: this,mat
    Class(Matrix),intent(INOUT) :: res
    ! calculate res = this - mat

    CALL subtract_matrix(this%Data,mat%Data, res%Data)
  END SUBROUTINE mp_subtract_matrix

  SUBROUTINE subtract_from_matrix(this,mat)
    Class(Matrix),INTENT(INOUT) :: this
    Class(Matrix), intent(IN) :: mat
    ! calculate this = this - mat

    CALL subtract_matrix(this%Data,mat%Data)
  END SUBROUTINE subtract_from_matrix

  SUBROUTINE mp_multiply_matrix_with_real(this, scalar, res)
    Class(Matrix), intent(IN) :: this
    REAL, intent(IN) :: scalar
    Class(Matrix),intent(INOUT) :: res

    CALL multiply_matrix_with_scalar(this%Data,scalar, res%Data)
  END SUBROUTINE mp_multiply_matrix_with_real

  SUBROUTINE mp_scale_matrix_by_real(this, scalar)
    Class(Matrix), intent(INOUT) :: this
    REAL, intent(IN) :: scalar

    CALL multiply_matrix_with_scalar(this%Data,scalar)
  END SUBROUTINE mp_scale_matrix_by_real

  SUBROUTINE mp_assign_matrix(lmat,rmat)
    Class(Matrix), intent(INOUT) :: lmat
    Class(Matrix), intent(IN) :: rmat

    IF ((lmat%Nrows.EQ.rmat%Nrows).AND.(lmat%Ncols.EQ.rmat%Ncols)) THEN
       lmat%Data = rmat%Data
    ELSE
       PRINT*, "Left and right matrices in an assignment must have the same dimensions."
    END IF
  END SUBROUTINE mp_assign_matrix

  ! ====================================================
  ! == Define some mathematical routines for matrices ==
  ! ====================================================
  SUBROUTINE mp_invert_matrix(mat)
    Class(Matrix) :: mat

    CALL invert(mat%Data)
  END SUBROUTINE mp_invert_matrix

  SUBROUTINE mp_LU_factor(mat)
    Class(Matrix) :: mat

    IF (mat%NRows.EQ.mat%NCOLS) THEN
      CALL LU_factor(mat%Data)
    ELSE
      PRINT*,"Matrix for _LU_factor is not square"
      STOP
    ENDIF
  END SUBROUTINE mp_LU_factor

  SUBROUTINE mp_LU_solve(this, mat, resmat, phase)
    Class(Matrix),intent(INOUT) :: this
    Class(Matrix),intent(IN) :: mat
    Class(Matrix),intent(INOUT) :: resmat
    INTEGER,intent(IN) :: phase

    IF ( this%NRows.EQ.this%NCols .AND. &
         this%NRows.EQ.mat%NRows .AND.&
         this%NRows.EQ.resmat%NRows .AND.&
         mat%NCols.EQ.resmat%NCols ) THEN
       CALL LU_solve(this%Data,mat%Data,resmat%Data,phase)
    ELSE
       PRINT*,"Matrix shapes for LU_solve do not match, ",&
            &this%NCols,this%NRows,mat%NCols,mat%NRows,resmat%NCols,resmat%NRows
       STOP
    END IF
  END SUBROUTINE mp_LU_solve

  ! ===================================
  ! == Define communication routines ==
  ! ===================================
  SUBROUTINE my_matrix_sum_0D(localsum,communicator)
    !INTEGER, intent(IN) :: len
    Class(Matrix), INTENT(INOUT) :: localsum
    integer :: communicator

    CALL my_sum(localsum%Data,communicator)
  END SUBROUTINE my_matrix_sum_0D

  SUBROUTINE my_matrix_sum_2D(localsum,communicator)
    Class(Matrix), DIMENSION(:,:), INTENT(INOUT) :: localsum
    integer :: communicator

    integer :: i,j

    do j=1,size(localsum,2)
       do i=1,size(localsum,1)
          call my_sum(localsum(i,j)%data,communicator)
       end do
    end do
    !CALL my_matrix_sum_generic(localsum,len,communicator)
  END SUBROUTINE my_matrix_sum_2D

  SUBROUTINE my_matrix_sum_3D(localsum,communicator)
    Class(Matrix), DIMENSION(:,:,:), INTENT(INOUT) :: localsum
    integer :: communicator

    integer :: k

    do k=1,size(localsum,3)
       call my_sum(localsum(:,:,k),communicator)
    end do
    !CALL my_matrix_sum_generic(localsum,len,communicator)
  END SUBROUTINE my_matrix_sum_3D
END MODULE MatrixModule
