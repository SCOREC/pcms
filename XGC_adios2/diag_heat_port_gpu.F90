attributes(device) &
subroutine diag_heat_port_gpu(w, pot, epara, eperp, ct, old_ph, new_ph, dphi, stype, ith)
    use sml_module_gpu
    use ptl_module_gpu
    use diag_module_gpu !, only : diag_heat_nvar
    use precision_mod_gpu

    implicit none
    real (kind=work_p), intent(in) :: w, pot, epara, eperp
    real (kind=work_p), intent(in) :: ct(ptl_nconst), old_ph(ptl_nphase), new_ph(ptl_nphase), dphi
    integer, intent(in) :: stype, ith
    real (kind=work_p) :: wp, r, z, psi, rn, zn, pn, wr,wz, ws, v(5)
    integer :: ir, iz, ip, itmp

    real (kind=work_p) :: dval0, dummy0,dval1,dummy1, dval2, dummy2, dval3, dummy3
    integer :: lb,ub,i,i0,ith0
    real (kind=work_p) :: x(2), phi, xff(2), phi_mid

    ! chracteristic r and z - use mean value for simplicity
!    r=(old_ph(pir)+new_ph(pir))*0.5_work_p
!    z=(old_ph(piz)+new_ph(piz))*0.5_work_p

    !get field following
    x=new_ph(1:2)
    phi=new_ph(3)
    phi_mid=(floor(phi/dphi) + 0.5_work_p) * dphi
    call field_following_pos2_gpu(x,phi,phi_mid,xff)
    r=xff(1)
    z=xff(2)

    wp=w*ct(piw0)

    v(1)=1_work_p
    v(2)=epara
    v(3)=eperp
    v(4)=pot
    v(5)=sqrt((old_ph(pir)-new_ph(pir))**2+(old_ph(piz)-old_ph(piz))**2)

    ! for all sections
    do i=1, diag_heat_nsection
        ! check range
        if(diag_heat_rmin(i) < r .and. r < diag_heat_rmax(i) &
            .and. diag_heat_zmin(i) < z .and. z < diag_heat_zmax(i)) then

        !r index
        rn=(r-diag_heat_rmin(i))/diag_heat_dr(i)
        ir=floor(rn)+1
        if(ir<1 .or. diag_heat_nr<= ir) cycle
        wr=1.0_work_p - rn + real(ir-1,kind=work_p)

        !z index
        zn=(z-diag_heat_zmin(i))/diag_heat_dz(i)
        iz=floor(zn)+1
        if(iz<1 .or. diag_heat_nz<= iz) cycle
        wz=1.0_work_p - zn + real(iz-1,kind=work_p)

#ifdef ORIGINAL
        diag_heat_pv(:,ir  ,iz  ,i,stype,ith)=diag_heat_pv(:,ir  ,iz  ,i,stype,ith) + v*wp*wr      *wz
        diag_heat_pv(:,ir+1,iz  ,i,stype,ith)=diag_heat_pv(:,ir+1,iz  ,i,stype,ith) + v*wp*(1.0_work_p-wr)*wz
        diag_heat_pv(:,ir  ,iz+1,i,stype,ith)=diag_heat_pv(:,ir  ,iz+1,i,stype,ith) + v*wp*wr      *(1.0_work_p-wz)
        diag_heat_pv(:,ir+1,iz+1,i,stype,ith)=diag_heat_pv(:,ir+1,iz+1,i,stype,ith) + v*wp*(1.0_work_p-wr)*(1.0_work_p-wz)
#else
        ith0 = 1 + mod(ith,size(diag_heat_pv,6))
        lb = lbound(diag_heat_pv,1)
        ub = ubound(diag_heat_pv,1)

        do itmp=lb,ub
           i0 = (itmp-lb) + lbound(v,1)
           dval0 = v(i0)*wp*wr      *wz
           dummy0 = atomicAdd( diag_heat_pv(itmp,ir,iz,i,stype,ith0), dval0 )

           dval1 = v(i0)*wp*(1.0_work_p-wr)*wz
           dummy1 = atomicAdd( diag_heat_pv(itmp,ir+1,iz,i,stype,ith0), dval1 )

           dval2 = v(i0)*wp*wr      *(1.0_work_p-wz)
           dummy2 = atomicAdd( diag_heat_pv(itmp,ir,iz+1,i,stype,ith0), dval2 )

           dval3 = v(i0)*wp*(1.0_work_p-wr)*(1.0_work_p-wz)
           dummy3 = atomicAdd( diag_heat_pv(itmp,ir+1,iz+1,i,stype,ith0), dval3 )
        enddo
#endif
        !psi 
        psi=psi_interpol_gpu(r,z,0,0)
        pn=(psi-diag_heat_pmin(i))/diag_heat_dp(i)
        ip=floor(pn)+1
        if(ip<1 .or. diag_heat_npsi<= ip) cycle
        ws=1.0_work_p - pn + real(ip-1,kind=work_p)

#ifdef ORIGINAL
        diag_heat_pv_psi(:,ip  ,i,stype,ith)=diag_heat_pv_psi(:,ip  ,i,stype,ith) + v*wp*ws
        diag_heat_pv_psi(:,ip+1,i,stype,ith)=diag_heat_pv_psi(:,ip+1,i,stype,ith) + v*wp*(1.0_work_p-ws)
#else
        ith0 = 1 + mod(ith,size(diag_heat_pv_psi,5))
        lb = lbound(diag_heat_pv_psi,1)
        ub = ubound(diag_heat_pv_psi,1)

        do itmp=lb,ub
           i0 = (itmp-lb) + lbound(v,1)
           dval0 = v(i0)*wp*ws
           dummy0 = atomicAdd( diag_heat_pv_psi(itmp,ip,i,stype,ith0), dval0 )

           dval1 = v(i0)*wp*(1.0_work_p-ws)
           dummy1 = atomicAdd( diag_heat_pv_psi(itmp,ip+1,i,stype,ith0), dval1 )
        enddo
#endif 
        endif
    enddo 
return
end subroutine diag_heat_port_gpu
