        subroutine geomhash_rz(n,rz, xstart, list )
        use eq_module
#if defined(CAM_TIMERS)
        use perf_mod
#endif
        implicit none

        integer, intent(in) :: n
        real(kind=8), dimension(1:2,n), intent(in) :: rz
        integer, dimension((eq_mr-1)*(eq_mz-1)+1), intent(inout) :: xstart
        integer, dimension(n), intent(inout) :: list



        real(kind=8) :: dr,dz,dr_inv,dz_inv,rmin,zmin
        integer :: nr, nz
#if defined(CAM_TIMERS)
        call t_startf('geomhash_rz')
#endif
        dr = (eq_rgrid(2) - eq_rgrid(1))
        dz = (eq_zgrid(2) - eq_zgrid(1))
        dr_inv = 1.0d0/dr
        dz_inv = 1.0d0/dz
        nr = eq_mr - 1
        nz = eq_mz - 1
        rmin = eq_rgrid(1)
        zmin = eq_zgrid(1)

        call geomhash(n,rz, nr,nz, rmin,dr_inv,zmin,dz_inv,xstart,list)
#if defined(CAM_TIMERS)
        call t_stopf('geomhash_rz')
#endif
        return
        end subroutine 
