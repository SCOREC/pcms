! This interfaces to external C bicub evaluation routines,
! which see some speedup because they use 2D cudaArrays for
! the acoeff array.
module ex_bicub_mod
  use cudafor
  interface
    subroutine init_external_bicub(rc_cub, zc_cub, acoeff, ncoeff, nr, nz, &
        rmin, zmin, dr_inv, dz_inv) &
      bind(C,name="init_cuda_bicub")
      use iso_c_binding
      implicit none
      real (c_double), dimension(*) :: rc_cub, zc_cub, acoeff
      integer (c_int), value :: ncoeff, nr, nz
      real (c_double), value :: rmin, zmin, dr_inv, dz_inv
    end subroutine init_external_bicub

    subroutine destroy_external_bicub() &
      bind(C,name="destroy_cuda_bicub")
      use iso_c_binding
    end subroutine destroy_external_bicub

    attributes(device) &
    subroutine fetch_external_bicub(rc_cub, zc_cub, acoeff, ir, iz) &
      bind(C,name="fetch_cuda_bicub")
      use iso_c_binding
      use cudafor
      implicit none
      real (c_double), device ::  rc_cub, zc_cub
      real (c_double), dimension(*), device :: acoeff
      integer (c_int), value :: ir, iz
    end subroutine fetch_external_bicub

    attributes(device) &
    subroutine eval_0_external(x, y, f00) &
      bind(C,name="eval_0_cuda")
      use iso_c_binding
      use cudafor
      implicit none
      real (c_double), value ::  x, y
      real (c_double), device :: f00
    end subroutine eval_0_external

    attributes(device) &
    subroutine eval_1_external(x, y, f00, f10, f01) &
      bind(C,name="eval_1_cuda")
      use iso_c_binding
      use cudafor
      implicit none
      real (c_double), value ::  x, y
      real (c_double), device :: f00, f10, f01
    end subroutine eval_1_external

    attributes(device) &
    subroutine eval_2_external(x, y, f00, f10, f01, f11, f20, f02) &
      bind(C,name="eval_2_cuda")
      use iso_c_binding
      use cudafor
      implicit none
      real (c_double), value ::  x, y
      real (c_double), device :: f00, f10, f01, f11, f20, f02
    end subroutine eval_2_external
  end interface

end module ex_bicub_mod
