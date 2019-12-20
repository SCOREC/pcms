!!****************************************************************************
!!> calculate field of a given point using interpolation funtions
!!  adopted from xorbit
!!
!!  first created : 2000/10/19
!!  last modified : 2006/02/01
!!  B->-B routine added (commented out)
!!  time dependance of field is added (2002/6/19)
!!  2006/02/01 fld module update
!!  2002/09/10 code optimization for speed
!!  2002/11/18 code modification for gxc
!!****************************************************************************
attributes(device) &
subroutine field_gpu(fld,t,rz_outside)
    use sml_module_gpu, only :sml_bp_sign, sml_time
    use eq_module_gpu, only : eq_min_r,eq_max_r,eq_min_z,eq_max_z, eq_x_psi, eq_x_z
    use fld_module
    use precision_mod_gpu
    implicit none
    type(fld_type),intent(inout) :: fld  !! Field information
    real (kind=work_p), intent(in) :: t  !! time
    logical , intent(out) :: rz_outside

    real (kind=work_p) :: psi, dpsi_dr, dpsi_dz, d2psi_d2r, d2psi_drdz, d2psi_d2z
    real (kind=work_p) :: r,z, ret(6)
    !real (kind=8) , external :: I_interpol
    real (kind=work_p) :: r2, over_r,over_r2 !! variables for opimization
    real (kind=work_p) :: cos_rip, sin_rip, ripp, dripp_dr, dripp_dz
    real (kind=work_p) :: rippbphi, drippbphi_dr, drippbphi_dz, drippbphi_dphi

    r=fld%r
    z=fld%z
    r2=r**2
    over_r=1/r
    over_r2=over_r**2
    if(r<eq_min_r) then
        r=eq_min_r
        rz_outside=.true.
    else if (r>eq_max_r)then
        r=eq_max_r
        rz_outside=.true.
    endif
    if(z<eq_min_z) then
        z=eq_min_z
        rz_outside=.true.
    else if (z>eq_max_z)then
       z=eq_max_z
       rz_outside=.true.
    else
       rz_outside=.false.
    endif

    call psi_der_all_gpu(r,z,ret)
    psi        =ret(1)
    dpsi_dr    =ret(2)
    dpsi_dz    =ret(3)
    d2psi_d2r  =ret(4)
    d2psi_drdz =ret(5)
    d2psi_d2z  =ret(6)

    fld%psi=psi
    fld%dpsidr=dpsi_dr
    fld%dpsidz=dpsi_dz
    ! added 2001/06/01  - lower bound of psi
    if(psi<0D0) then
        psi=0D0
    endif

    !fld_q=q_interpol(psi,0)  --> no need
    
!    if(psi<eq_x_psi .AND. z<eq_x_z) then
    if(.not. is_rgn12_gpu(r,z,psi) ) then
        fld%I=I_interpol_gpu(psi,0,3)
        fld%dIdpsi = I_interpol_gpu(psi,1,3)
    else
        fld%I=I_interpol_gpu(psi,0,1)
        fld%dIdpsi = I_interpol_gpu(psi,1,1)
    endif

    
    fld%br=- dpsi_dz *over_r * sml_bp_sign
    fld%bz= dpsi_dr *over_r  * sml_bp_sign
    fld%bphi=fld%I *over_r   
        
    !derivativs
    fld%dbrdr= (dpsi_dz *over_r2 - d2psi_drdz *over_r) * sml_bp_sign
    fld%dbrdz=- d2psi_d2z *over_r                   * sml_bp_sign
    fld%dbrdp=0D0                                   * sml_bp_sign
    
    fld%dbzdr= (- dpsi_dr * over_r2 + d2psi_d2r *over_r) * sml_bp_sign
    fld%dbzdz= d2psi_drdz *over_r                      * sml_bp_sign
    fld%dbzdp=0D0                                      * sml_bp_sign
           
    fld%dbpdr= dpsi_dr * fld%dIdpsi *over_r - fld%I *over_r2
    fld%dbpdz= fld%dIdpsi * dpsi_dz *over_r
    fld%dbpdp=0D0
!    call efield(r,z,fld_psi,t)  ! set fld_er, ez, ephi value for a give r,z value  -> move to derivs

end subroutine field_gpu
