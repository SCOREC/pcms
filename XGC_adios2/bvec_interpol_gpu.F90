attributes(device) &
subroutine bvec_interpol_gpu(r,z,phi,br,bz,bphi)
  use eq_module_gpu
  use sml_module_gpu, only : sml_bp_sign
  use precision_mod_gpu
  use bicub_mod_gpu, only : bicub_interpol1
  implicit none
  real (kind=work_p) , intent(in) :: r,z,phi
  real (kind=work_p)              :: psi,dpsi_dr,dpsi_dz,br,bz,bphi,fi
  real (kind=work_p)              :: ripp, dripp_dr, dripp_dz

 call bicub_interpol1(r,z,psi,dpsi_dr,dpsi_dz)

 if(psi<eq_x_psi .AND. z<eq_x_z) then
    fi=I_interpol_gpu(psi,0,3)
 else
    fi=I_interpol_gpu(psi,0,1)
 endif

 br=- dpsi_dz / r * sml_bp_sign
 bz= dpsi_dr / r  * sml_bp_sign
 bphi=fi / r

end subroutine bvec_interpol_gpu


