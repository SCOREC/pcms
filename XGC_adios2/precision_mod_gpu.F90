module precision_mod_gpu

! single precision
integer, parameter :: single_p = kind(1.0)

! double precision
integer, parameter :: double_p = kind(1.0d0)

! working precision, set wp = dp or  wp = sp

integer, parameter :: work_p = double_p
integer, parameter :: eq_mpsi = 151
integer, parameter :: eq_mr = 151
integer, parameter :: eq_mz = 151
end module precision_mod_gpu

