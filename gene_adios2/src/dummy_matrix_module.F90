#include "switches.h"
#include "intrinsic_sizes.h"
#include "redef.h"
!****h* /MatrixModule
! DESCRIPTION
!****
MODULE MatrixModule
  !USE StoreFullMatrixModule
  !USE mpi
  IMPLICIT NONE

  private

  INTEGER, PARAMETER,PUBLIC :: MAT_BLOCKS_OF_ROWS=1000, MAT_BLOCKS_OF_COLS=1001, MAT_BOTH=1002

  TYPE,public :: Matrix
     INTEGER :: NCols
     INTEGER :: NRows

     !class(StoreFullMatrixObject),allocatable :: data
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

  public :: initialize_matrix_module,finalize_matrix_module,my_sum
  public :: static_size_of_matrix
CONTAINS
  FUNCTION static_size_of_matrix() RESULT(memory_need)
    integer :: memory_need

    memory_need = 2*SIZE_OF_INTEGER

  END FUNCTION static_size_of_matrix

  SUBROUTINE initialize_matrix(this,n_rows,n_cols)
    Class(Matrix), INTENT(INOUT) :: this
    INTEGER :: n_rows, n_cols

  END SUBROUTINE initialize_matrix

  SUBROUTINE initialize_vector(this,n_rows)
    class(matrix) :: this
    integer :: n_rows

  END SUBROUTINE initialize_vector

  SUBROUTINE mp_initialize_matrix_module(comm)
    INTEGER :: comm

  END SUBROUTINE mp_initialize_matrix_module

  LOGICAL FUNCTION mp_isInitialized(this)
    class(matrix) :: this

    mp_isInitialized = .true.
  END FUNCTION mp_isInitialized

  SUBROUTINE mp_allocate(this)
    Class(Matrix), intent(INOUT) :: this

  END SUBROUTINE mp_allocate

  SUBROUTINE attach_vector(this,localarray)
    Class(Matrix),INTENT(INOUT) :: this
    COMPLEX, DIMENSION(:),INTENT(IN),TARGET :: localarray

  END SUBROUTINE attach_vector

  SUBROUTINE attach_matrix(this,localarray)
    Class(Matrix),INTENT(INOUT) :: this
    COMPLEX, DIMENSION(:,:),INTENT(IN),TARGET :: localarray

  END SUBROUTINE attach_matrix

  SUBROUTINE mp_finalize_matrix(this)
    class(matrix) :: this
  END SUBROUTINE mp_finalize_matrix

  SUBROUTINE mp_finalize_module()
  END SUBROUTINE mp_finalize_module

  SUBROUTINE set_complex_value(this,irow,icol,value)
    Class(Matrix),INTENT(INOUT) :: this
    INTEGER,INTENT(IN) :: irow,icol
    COMPLEX, INTENT(IN) :: value

  END SUBROUTINE set_complex_value

  SUBROUTINE set_real_value(this,irow,icol,value)
    Class(Matrix),INTENT(INOUT) :: this
    INTEGER,INTENT(IN) :: irow,icol
    REAL, INTENT(IN) :: value

  END SUBROUTINE set_real_value

  SUBROUTINE add_complex_value(this,irow,icol,value)
    Class(Matrix),INTENT(INOUT) :: this
    INTEGER,INTENT(IN) :: irow,icol
    COMPLEX, INTENT(IN) :: value
  END SUBROUTINE add_complex_value

  SUBROUTINE add_real_value(this,irow,icol,value)
    Class(Matrix),INTENT(INOUT) :: this
    INTEGER,INTENT(IN) :: irow,icol
    REAL, INTENT(IN) :: value
  END SUBROUTINE add_real_value

  SUBROUTINE mat_add_to_row(this,irow,whole_row)
    Class(Matrix),INTENT(INOUT) :: this
    INTEGER,INTENT(IN) :: irow
    COMPLEX, dimension(:),INTENT(IN) :: whole_row
  END SUBROUTINE mat_add_to_row

  SUBROUTINE mp_commit_values(this)
    Class(Matrix) :: this
  END SUBROUTINE mp_commit_values

  SUBROUTINE mp_dot_multiply_MatMat(this,mat, resmat, transA_in, transB_in)
    Class(Matrix),intent(IN) :: this
    Class(Matrix),intent(IN) :: mat
    Class(Matrix),intent(INOUT) :: resmat
    CHARACTER, intent(IN), optional :: transA_in, transB_in
  END SUBROUTINE mp_dot_multiply_MatMat

  SUBROUTINE mp_set_zero(this)
    class(matrix) :: this
  END SUBROUTINE mp_set_zero

  SUBROUTINE mp_get_global_matrix_locally(mat,localfullmat)
    Class(Matrix),INTENT(IN) :: mat
    COMPLEX, DIMENSION(:,:),INTENT(OUT) :: localfullmat

    localfullmat = cmplx(0.0,0.0)
  END SUBROUTINE mp_get_global_matrix_locally

  FUNCTION mp_get_value(mat,irow,icol) result(value)
    class(matrix) :: mat
    INTEGER :: irow, icol
    complex :: value
  END FUNCTION mp_get_value

  FUNCTION mp_get_local_abs_square_sum(mat) RESULT(value)
    class(matrix) :: mat
    real :: value
  END FUNCTION mp_get_local_abs_square_sum

  SUBROUTINE mp_show_on_screen(mat)
    class(matrix) :: mat
  END SUBROUTINE mp_show_on_screen

  SUBROUTINE mp_show_in_file(mat,filename)
    class(matrix) :: mat
    CHARACTER(len=*) :: filename
  END SUBROUTINE mp_show_in_file

  SUBROUTINE mp_output_data_matrix(this,basename)
    class(matrix) :: this
    CHARACTER(len=*) :: basename
  END SUBROUTINE mp_output_data_matrix

  ! =======================================
  ! == Define some operators on matrices ==
  ! =======================================
  SUBROUTINE mp_add_to_matrix(this,mat)
    Class(Matrix), intent(INOUT) :: this
    Class(Matrix),INTENT(IN) :: mat
  END SUBROUTINE mp_add_to_matrix

  SUBROUTINE mp_add_matrix(this,mat,res)
    Class(Matrix), INTENT(IN) :: this,mat
    Class(Matrix), INTENT(INOUT) :: res
  END SUBROUTINE mp_add_matrix

  SUBROUTINE mp_subtract_matrix(this,mat,res)
    Class(Matrix),INTENT(IN) :: this,mat
    Class(Matrix),intent(INOUT) :: res
  END SUBROUTINE mp_subtract_matrix

  SUBROUTINE subtract_from_matrix(this,mat)
    Class(Matrix),INTENT(INOUT) :: this
    Class(Matrix), intent(IN) :: mat
  END SUBROUTINE subtract_from_matrix

  SUBROUTINE mp_multiply_matrix_with_real(this, scalar, res)
    Class(Matrix), intent(IN) :: this
    REAL, intent(IN) :: scalar
    Class(Matrix),intent(INOUT) :: res
  END SUBROUTINE mp_multiply_matrix_with_real

  SUBROUTINE mp_scale_matrix_by_real(this, scalar)
    Class(Matrix), intent(INOUT) :: this
    REAL, intent(IN) :: scalar
  END SUBROUTINE mp_scale_matrix_by_real

  SUBROUTINE mp_assign_matrix(lmat,rmat)
    Class(Matrix), intent(INOUT) :: lmat
    Class(Matrix), intent(IN) :: rmat
  END SUBROUTINE mp_assign_matrix

  ! ====================================================
  ! == Define some mathematical routines for matrices ==
  ! ====================================================
  SUBROUTINE mp_invert_matrix(mat)
    Class(Matrix) :: mat
  END SUBROUTINE mp_invert_matrix

  SUBROUTINE mp_LU_factor(mat)
    Class(Matrix) :: mat
  END SUBROUTINE mp_LU_factor

  SUBROUTINE mp_LU_solve(this, mat, resmat, phase)
    Class(Matrix),intent(INOUT) :: this
    Class(Matrix),intent(IN) :: mat
    Class(Matrix),intent(INOUT) :: resmat
    INTEGER,intent(IN) :: phase
  END SUBROUTINE mp_LU_solve

  ! ===================================
  ! == Define communication routines ==
  ! ===================================
  SUBROUTINE my_matrix_sum_0D(localsum,communicator)
    !INTEGER, intent(IN) :: len
    Class(Matrix), INTENT(INOUT) :: localsum
    integer :: communicator
  END SUBROUTINE my_matrix_sum_0D

  SUBROUTINE my_matrix_sum_2D(localsum,communicator)
    Class(Matrix), DIMENSION(:,:), INTENT(INOUT) :: localsum
    integer :: communicator
  END SUBROUTINE my_matrix_sum_2D

  SUBROUTINE my_matrix_sum_3D(localsum,communicator)
    Class(Matrix), DIMENSION(:,:,:), INTENT(INOUT) :: localsum
    integer :: communicator
  END SUBROUTINE my_matrix_sum_3D

END MODULE MatrixModule

!=========================== VectorModule dummy ===============================

MODULE VectorModule
  implicit none

  private
  PUBLIC :: Vector,initialize,attach, finalize, isInitialized, &
       & initialize_vector_module, finalize_vector_module, &
       & static_size_of_vector

  TYPE Vector
     INTEGER :: NRows
  END TYPE Vector

  INTERFACE initialize
     MODULE PROCEDURE mp_initialize_vector
  END INTERFACE

  INTERFACE attach
     MODULE PROCEDURE mp_attach_vector
  END INTERFACE

  INTERFACE finalize
     MODULE PROCEDURE mp_finalize_vector
  END INTERFACE

  INTERFACE isInitialized
     module procedure mp_isInitialized
  END INTERFACE

contains
  FUNCTION static_size_of_vector() RESULT(memory_need)
    integer :: memory_need

    memory_need = SIZE_OF_INTEGER

  END FUNCTION static_size_of_vector

  SUBROUTINE initialize_vector_module
  END SUBROUTINE initialize_vector_module

  SUBROUTINE finalize_vector_module
  END SUBROUTINE finalize_vector_module

  SUBROUTINE mp_initialize_vector(vec,n_rows)
    TYPE(Vector), INTENT(INOUT) :: vec
    INTEGER :: n_rows

  END SUBROUTINE mp_initialize_vector

  SUBROUTINE mp_finalize_vector(vec)
    type(Vector) :: vec
  END SUBROUTINE mp_finalize_vector

  SUBROUTINE mp_attach_vector(vec,localarray)
    TYPE(Vector),INTENT(INOUT) :: vec
    COMPLEX, DIMENSION(:),INTENT(IN),TARGET :: localarray
  END SUBROUTINE mp_attach_vector

  LOGICAL FUNCTION mp_isInitialized(vec)
    type(Vector) :: vec

    mp_isInitialized = .false.
  END FUNCTION mp_isInitialized

END MODULE VectorModule

!=========================== BandedMatrixModule dummy ========================

MODULE BandedMatrixModule
  use MatrixModule
  use VectorModule
  implicit none

  PRIVATE

  public :: my_sum,static_size_of_BandedMatrix
  public :: initialize_BandedMatrix_module,finalize_BandedMatrix_module

!> Just knows the mathematical extent of the matrix, that is the
!! number of rows and columns. It also has a data object embedded
!! which contains the real data in form of a StoreBandedMatrixObject.
  TYPE,public ::  BandedMatrix
     INTEGER :: NCols
     INTEGER :: NRows
     !type(StoreBandedMatrixObject) :: data
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

  interface my_sum
     module procedure bm_sum_0D, bm_sum_2D, bm_sum_3D
  end interface my_sum

CONTAINS
  FUNCTION static_size_of_BandedMatrix() RESULT(memory_need)
    integer :: memory_need

    memory_need = 2*SIZE_OF_INTEGER
  END FUNCTION static_size_of_BandedMatrix

  subroutine bm_initialize_matrix(mat,n_rows,n_cols,transp)
    Class(BandedMatrix), INTENT(INOUT) :: mat
    INTEGER :: n_rows, n_cols
    LOGICAL, OPTIONAL :: transp
  end subroutine bm_initialize_matrix

  subroutine bm_initialize_vector(this,n_rows)
    class(BandedMatrix) :: this
    integer :: n_rows
  end subroutine bm_initialize_vector

  SUBROUTINE bm_initialize_BandedMatrix_module
  !  CALL get_processgrid(PG)
  END SUBROUTINE bm_initialize_BandedMatrix_module

  LOGICAL FUNCTION bm_isInitialized(mat)
    class(BandedMatrix) :: mat
    bm_isInitialized = .true.
  END FUNCTION bm_isInitialized

  SUBROUTINE bm_allocate(mat)
    Class(BandedMatrix), intent(INOUT) :: mat
  END SUBROUTINE bm_allocate

  SUBROUTINE bm_finalize_matrix(mat)
    type(BandedMatrix) :: mat
  END SUBROUTINE bm_finalize_matrix

  SUBROUTINE bm_finalize_BandedMatrix_module()
  END SUBROUTINE bm_finalize_BandedMatrix_module

  SUBROUTINE bm_set_complex_value(mat,irow,icol,value)
    Class(BandedMatrix),INTENT(INOUT) :: mat
    INTEGER,INTENT(IN) :: irow,icol
    COMPLEX, INTENT(IN) :: value
  END SUBROUTINE bm_set_complex_value

  SUBROUTINE bm_set_real_value(mat,irow,icol,value)
    Class(BandedMatrix),INTENT(INOUT) :: mat
    INTEGER,INTENT(IN) :: irow,icol
    REAL, INTENT(IN) :: value
  END SUBROUTINE bm_set_real_value

  SUBROUTINE bm_add_complex_value(mat,irow,icol,value)
    Class(BandedMatrix),INTENT(INOUT) :: mat
    INTEGER,INTENT(IN) :: irow,icol
    COMPLEX, INTENT(IN) :: value
  END SUBROUTINE bm_add_complex_value

  SUBROUTINE bm_add_real_value(mat,irow,icol,value)
    Class(BandedMatrix),INTENT(INOUT) :: mat
    INTEGER,INTENT(IN) :: irow,icol
    REAL, INTENT(IN) :: value
  END SUBROUTINE bm_add_real_value


  SUBROUTINE bm_commit_values(mat)
    Class(BandedMatrix) :: mat
  END SUBROUTINE bm_commit_values

  subroutine bm_autotune(mat,vec,res)
    Class(BandedMatrix),intent(IN) :: mat
    Type(Vector),intent(IN) :: vec
    Type(Vector),intent(INOUT) :: res
  end subroutine bm_autotune

  SUBROUTINE bm_dot_multiply(mat,vec, res)
    Class(BandedMatrix),intent(IN) :: mat
    Type(Vector),intent(IN) :: vec
    Type(Vector),intent(INOUT) :: res
  END SUBROUTINE bm_dot_multiply

  subroutine bm_matmat_dot_multiply(mat,mat2,res,transA_in,transB_in)
    class(BandedMatrix),intent(IN) :: mat, mat2
    class(BandedMatrix),intent(INOUT) :: res
    CHARACTER, intent(IN), optional :: transA_in, transB_in
  end subroutine bm_matmat_dot_multiply

  subroutine bm_dot_multiply_Banded_with_Full(bmat,mat,res,transA_in,transB_in)
    class(BandedMatrix),intent(IN) :: bmat
    type(Matrix), intent(IN) :: mat
    type(Matrix),intent(INOUT) :: res
    CHARACTER, intent(IN), optional :: transA_in, transB_in
  end subroutine bm_dot_multiply_Banded_with_Full

  SUBROUTINE bm_square_matrix(bmat,res_bmat)
    Class(BandedMatrix),intent(IN) :: bmat
    Class(BandedMatrix),intent(INOUT) :: res_bmat
  END SUBROUTINE bm_square_matrix

  SUBROUTINE bm_set_zero(mat)
    class(BandedMatrix) :: mat
  END SUBROUTINE bm_set_zero

  SUBROUTINE bm_get_global_matrix_locally(mat,localfullmat)
    Class(BandedMatrix),INTENT(IN) :: mat
    COMPLEX, DIMENSION(:,:),INTENT(OUT) :: localfullmat
    localfullmat = cmplx(0.0,0.0)
  END SUBROUTINE bm_get_global_matrix_locally

  FUNCTION bm_get_value(mat,irow,icol,trans_in) result(value)
    class(BandedMatrix) :: mat
    INTEGER :: irow, icol
    complex :: value
    character(len=1),optional,intent(IN) :: trans_in
  END FUNCTION bm_get_value

  subroutine bm_get_row_pointer(mat,global_row,ptr_row,s_global_col,e_global_col)
    Class(BandedMatrix) :: mat
    INTEGER,intent(IN) :: global_row
    complex,dimension(:),pointer,intent(OUT) :: ptr_row
    integer, INTENT(OUT) :: s_global_col, e_global_col
    ptr_row =>NULL()
    s_global_col=0
    e_global_col=0
  end subroutine bm_get_row_pointer

  FUNCTION bm_get_local_abs_square_sum(mat) RESULT(value)
    class(BandedMatrix) :: mat
    real :: value
  END FUNCTION bm_get_local_abs_square_sum

  SUBROUTINE bm_show_on_screen(mat)
    class(BandedMatrix) :: mat
  END SUBROUTINE bm_show_on_screen

  SUBROUTINE bm_show_in_file(mat,filename)
    class(BandedMatrix) :: mat
    CHARACTER(len=*) :: filename
  END SUBROUTINE bm_show_in_file

  SUBROUTINE bm_output_data_matrix(this,basename)
    class(BandedMatrix) :: this
    CHARACTER(len=*) :: basename
  END SUBROUTINE bm_output_data_matrix

  ! =======================================
  ! == Define some operators on matrices ==
  ! =======================================
  SUBROUTINE bm_add_to_matrix(this,mat)
    Class(BandedMatrix), intent(INOUT) :: this
    Class(BandedMatrix),INTENT(IN) :: mat
  END SUBROUTINE bm_add_to_matrix

  SUBROUTINE bm_add_matrix(this,mat,res)
    Class(BandedMatrix), INTENT(IN) :: this,mat
    Class(BandedMatrix), INTENT(INOUT) :: res
  END SUBROUTINE bm_add_matrix

  SUBROUTINE bm_subtract_matrix(this,mat,res)
    Class(BandedMatrix),INTENT(IN) :: this,mat
    Class(BandedMatrix),intent(INOUT) :: res
  END SUBROUTINE bm_subtract_matrix

  SUBROUTINE bm_subtract_from_matrix(this,mat)
    Class(BandedMatrix),INTENT(INOUT) :: this
    Class(BandedMatrix), intent(IN) :: mat
  END SUBROUTINE bm_subtract_from_matrix

  SUBROUTINE bm_multiply_matrix_with_real(this, scalar, res)
    Class(BandedMatrix), intent(IN) :: this
    REAL, intent(IN) :: scalar
    Class(BandedMatrix),intent(INOUT) :: res
  END SUBROUTINE bm_multiply_matrix_with_real

  SUBROUTINE bm_scale_matrix_by_real(this, scalar)
    Class(BandedMatrix), intent(INOUT) :: this
    REAL, intent(IN) :: scalar
  END SUBROUTINE bm_scale_matrix_by_real

  SUBROUTINE bm_assign_matrix(lmat,rmat)
    Class(BandedMatrix), intent(INOUT) :: lmat
    Type(BandedMatrix), intent(IN) :: rmat
  END SUBROUTINE bm_assign_matrix

  SUBROUTINE bm_print_storage_details(bmat)
    Class(BandedMatrix), INTENT(IN) :: bmat
  END SUBROUTINE bm_print_storage_details

  FUNCTION bm_get_number_of_bands(mat)
    Class(BandedMatrix), intent(IN) :: mat
    integer :: bm_get_number_of_bands
    bm_get_number_of_bands=1
  END FUNCTION bm_get_number_of_bands

  ! ===================================
  ! == Define communication routines ==
  ! ===================================
  SUBROUTINE bm_sum_0D(localsum,communicator)
    Class(BandedMatrix), INTENT(INOUT) :: localsum
    integer :: communicator
  END SUBROUTINE bm_sum_0D

  SUBROUTINE bm_sum_2D(localsum,communicator)
    Class(BandedMatrix), DIMENSION(:,:), INTENT(INOUT) :: localsum
    integer :: communicator
  END SUBROUTINE bm_sum_2D

  SUBROUTINE bm_sum_3D(localsum,communicator)
    Class(BandedMatrix), DIMENSION(:,:,:), INTENT(INOUT) :: localsum
    integer :: communicator
  END SUBROUTINE bm_sum_3D

  SUBROUTINE bm_row_axpy(this,this_row,mat,scalar)
    Class(BandedMatrix),intent(INOUT) :: this
    INTEGER, INTENT(IN) :: this_row
    Class(BandedMatrix), intent(IN) :: mat
    COMPLEX, intent(IN) :: scalar
  END SUBROUTINE bm_row_axpy

  subroutine bm_transpose_and_conjugate(mat,dagger_mat)
    class(BandedMatrix) :: mat,dagger_mat
  end subroutine bm_transpose_and_conjugate

  subroutine bm_transpose_storage(mat,mat_t)
    class(BandedMatrix) :: mat,mat_t
  end subroutine bm_transpose_storage

  SUBROUTINE bm_LU_factor(mat)
    Class(BandedMatrix),intent(INOUT) :: mat
  END SUBROUTINE bm_LU_factor

  LOGICAL FUNCTION bm_LU_factor_ok(mat)
    Class(BandedMatrix),intent(IN) :: mat
    bm_LU_factor_ok = .true.
  END FUNCTION bm_LU_factor_ok

  SUBROUTINE bm_LU_solve(mat,vec,res)
    Class(BandedMatrix),intent(IN) :: mat
    Type(Vector),intent(IN) :: vec
    Type(Vector),intent(INOUT) :: res
  END SUBROUTINE bm_LU_solve

#ifdef ALL_ROUTINES
  ! ====================================================
  ! == Define some mathematical routines for matrices ==
  ! ====================================================
  SUBROUTINE bm_invert_matrix(mat)
    Class(BandedMatrix) :: mat
  END SUBROUTINE bm_invert_matrix
#endif

END MODULE BandedMatrixModule

! ======================= Derivative Matrix Module Dummy =====================
MODULE DerivativeMatrixModule
  USE Grid1DModule
  USE MatrixModule
  use BandedMatrixModule

  implicit none

  public :: static_size_of_DerivativeMatrix

  type,extends(BandedMatrix),public ::  DerivativeMatrix
     integer :: derivative_order !the order of the finite differences formulas for the derivatives
     logical :: periodic_boundary
     !TYPE(Grid1D), pointer :: grid
   contains
     procedure,private ::  dm_initialize_matrix
     generic,public :: initialize => dm_initialize_matrix
     final :: dm_finalize_matrix

     procedure,private :: dm_assign_matrix
     generic,public :: assignment(=) => dm_assign_matrix

     procedure :: calculate
  END TYPE DerivativeMatrix


CONTAINS

  FUNCTION static_size_of_DerivativeMatrix() RESULT(memory_need)
    integer :: memory_need

    memory_need = static_size_of_BandedMatrix()&
         & + SIZE_OF_INTEGER&
         & + SIZE_OF_LOGICAL

  END FUNCTION static_size_of_DerivativeMatrix

  SUBROUTINE mp_initialize_module
  END SUBROUTINE mp_initialize_module

  subroutine dm_initialize_matrix(this,grid,p_derivative_order,transposed)
    class(DerivativeMatrix),intent(inout) :: this !> derivative matrix
    TYPE(Grid1D),TARGET,intent(in) :: grid     !> grid definition
    INTEGER,intent(in) :: p_derivative_order   !> derivative order
    LOGICAL, OPTIONAL :: transposed
  end subroutine dm_initialize_matrix

  subroutine dm_finalize_matrix(this)
    type(DerivativeMatrix) :: this
  end subroutine dm_finalize_matrix

  SUBROUTINE dm_assign_matrix(lmat,rmat)
    Class(DerivativeMatrix), intent(INOUT) :: lmat
    type(DerivativeMatrix), intent(IN) :: rmat
   end SUBROUTINE dm_assign_matrix

  SUBROUTINE calculate(this, which_derivative, rad_bc_type)
    Class(DerivativeMatrix) :: this             !> matrix to contain the derivative matrix
    INTEGER, intent(IN) :: which_derivative    !> number of the derivative to calculate the matrix for
    !                   can be 1 or 2 for the first or second derivative
    INTEGER,intent(in) :: rad_bc_type          !> radial boundary condition
  END SUBROUTINE calculate
END MODULE DerivativeMatrixModule
