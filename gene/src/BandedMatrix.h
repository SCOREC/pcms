!> Just knows the mathematical extend of the matrix, that is the 
!! number of rows and columns. It also has a data object embedded
!! which contains the real data in form of a StoreBandedMatrixObject.
  TYPE BandedMatrix
     INTEGER :: NCols
     INTEGER :: NRows
     TYPE(StoreBandedMatrixObject) :: DATA
  END TYPE BandedMatrix

  INTERFACE initialize
    MODULE PROCEDURE mp_initialize_matrix, mp_initialize_vector
  END INTERFACE

  INTERFACE initialize_BandedMatrix_module
    MODULE PROCEDURE mp_initialize_BandedMatrix_module
  END INTERFACE

  INTERFACE allocate
    MODULE PROCEDURE mp_allocate
  END INTERFACE

  INTERFACE finalize
    MODULE PROCEDURE mp_finalize_matrix
  END INTERFACE

  INTERFACE finalize_BandedMatrix_module
     MODULE PROCEDURE mp_finalize_BandedMatrix_module
  END INTERFACE

  INTERFACE set_value
    MODULE PROCEDURE mp_set_complex_value, mp_set_real_value
  END INTERFACE

  INTERFACE add_value
    MODULE PROCEDURE mp_add_complex_value, mp_add_real_value
  END INTERFACE

  INTERFACE commit_values
    module procedure mp_commit_values
  END INTERFACE

  INTERFACE set_zero
    MODULE PROCEDURE mp_set_zero
  END INTERFACE

  INTERFACE get_global_matrix_locally
    module procedure mp_get_global_matrix_locally
  END INTERFACE

  INTERFACE convert_Banded_to_Full
    module procedure mp_convert_Banded_to_Full
  END INTERFACE

  INTERFACE convert_Full_to_Banded
    module procedure mp_convert_Full_to_Banded
  END INTERFACE

  INTERFACE get_local_abs_square_sum
    module procedure mp_get_local_abs_square_sum
  END INTERFACE

  INTERFACE mat_get_value
    MODULE PROCEDURE mp_get_value
  END INTERFACE

  INTERFACE mat_get_row_pointer
     MODULE PROCEDURE bm_get_row_pointer
  END INTERFACE

  interface autotune
     module procedure bm_autotune
  end interface

  INTERFACE dot_multiply
    MODULE PROCEDURE mp_dot_multiply, bm_matmat_dot_multiply
    module procedure bm_dot_multiply_Banded_with_Full
  END INTERFACE

  INTERFACE square_matrix
     module procedure bm_square_matrix
  END INTERFACE

  INTERFACE show
    MODULE PROCEDURE mp_show_on_screen, mp_show_in_file
  END INTERFACE

  INTERFACE add_matrix
    MODULE PROCEDURE mp_add_to_matrix, mp_add_matrix!, mp_add_Matrix_to_BandedMatrix
  END INTERFACE

  interface subtract_matrix
  !INTERFACE OPERATOR(-)
     MODULE PROCEDURE mp_subtract_matrix, mp_subtract_from_matrix
  END INTERFACE

  interface multiply_matrix_with_scalar
  !INTERFACE OPERATOR(*)
    module procedure mp_multiply_matrix_with_real
    module procedure mp_scale_matrix_by_real
  END INTERFACE

  INTERFACE ASSIGNMENT(=)
    module procedure mp_assign_matrix
  END INTERFACE

  INTERFACE output_data
    module procedure mp_output_data_matrix
  END INTERFACE

  INTERFACE isInitialized
    module procedure mp_isInitialized
  END INTERFACE

  INTERFACE print_storage_details
     MODULE PROCEDURE bm_print_storage_details
  END INTERFACE

  INTERFACE get_number_of_bands
     module procedure bm_get_number_of_bands
  END INTERFACE

  INTERFACE my_sum
     MODULE PROCEDURE bm_sum_generic, bm_sum_2D, bm_sum_3D, bm_sum_0D
  END INTERFACE

  INTERFACE row_axpy
     module procedure bm_row_axpy
  END INTERFACE

  interface transpose_and_conjugate
     module procedure bm_transpose_and_conjugate
  end interface

  interface transpose_storage
     module procedure bm_transpose_storage
  end interface

  INTERFACE LU_factor
     module procedure bm_LU_factor
  END INTERFACE

  INTERFACE LU_factor_ok
     module procedure bm_LU_factor_ok
  END INTERFACE

  INTERFACE LU_solve
     module procedure bm_LU_solve
  END INTERFACE

#ifdef ALL_ROUTINES
  INTERFACE invert
    module procedure mp_invert_matrix
  END INTERFACE
#endif
