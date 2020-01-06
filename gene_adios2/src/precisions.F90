MODULE precisions
  USE,intrinsic :: iso_c_binding
  implicit none

! Module precisions
 integer, parameter :: sp = selected_real_kind( 6, 35) ! single
 integer, parameter :: dp = selected_real_kind(12,307) ! double
 !INTEGER, PARAMETER :: qp = SELECTED_REAL_KIND(30,307) ! quad

 integer, parameter :: ip = kind(1) ! default integer

#ifdef DOUBLE_PREC
 integer, parameter :: rp = dp
 integer, parameter :: cp = dp
 INTEGER, parameter :: c_rp = C_DOUBLE
 INTEGER, parameter :: c_cp = C_DOUBLE_COMPLEX
#else
 ! single precision
 INTEGER, parameter :: rp = sp
 INTEGER, parameter :: cp = sp
 INTEGER, parameter :: c_rp = C_FLOAT
 INTEGER, parameter :: c_cp = C_FLOAT_COMPLEX
#endif

END MODULE precisions
