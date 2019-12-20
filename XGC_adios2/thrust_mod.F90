      module thrust_mod

      interface thrustsort
        subroutine sort_int(input,N)                                     &
     &             bind(C,name="sort_int_wrapper")
        use iso_c_binding
        integer(c_int),device:: input(*)
        integer(c_int),value:: N
        end subroutine

        subroutine sort_float(input,N)                                   &
     &             bind(C,name="sort_float_wrapper")
        use iso_c_binding
        real(c_float),device:: input(*)
        integer(c_int),value:: N
        end subroutine

        subroutine sort_double( input,N)                                 &
     &             bind(C,name="sort_double_wrapper")
        use iso_c_binding
        real(c_double),device:: input(*)
        integer(c_int),value:: N
        end subroutine

      end interface


      interface thrustscan
        subroutine scan_int(input,N)                                     &
     &             bind(C,name="scan_int_wrapper")
        use iso_c_binding
        integer(c_int),device:: input(*)
        integer(c_int),value:: N
        end subroutine

        subroutine scan_float(input,N)                                   &
     &             bind(C,name="scan_float_wrapper")
        use iso_c_binding
        real(c_float),device:: input(*)
        integer(c_int),value:: N
        end subroutine

        subroutine scan_double( input,N)                                 &
     &             bind(C,name="scan_double_wrapper")
        use iso_c_binding
        real(c_double),device:: input(*)
        integer(c_int),value:: N
        end subroutine
      end interface thrustscan

      end module thrust_mod


