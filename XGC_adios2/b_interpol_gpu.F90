attributes(device) &
  function b_interpol_gpu(r,z,phi)
  use eq_module_gpu, only : eq_x_psi, eq_x_z
  use bicub_mod_gpu, only : bicub_interpol1
  use precision_mod_gpu

  implicit none
  real (kind=work_p)              :: b_interpol_gpu
  real (kind=work_p) , intent(in) :: r,z,phi
  real (kind=work_p)              :: psi,dpsi_dr,dpsi_dz,br,bz,bphi,fi
  real (kind=work_p)              :: ripp,dripp_dr,dripp_dz


  
  call bicub_interpol1(r,z,psi,dpsi_dr,dpsi_dz)
  
  if(psi<eq_x_psi .AND. z<eq_x_z) then
     fi=I_interpol_gpu(psi,0,3)
  else
     fi=I_interpol_gpu(psi,0,1)
  endif
  
  br=- dpsi_dz / r   ! sign ignored   sml_bp_sign
  bz= dpsi_dr / r  
  bphi=fi / r
  
  b_interpol_gpu= sqrt(br**2+bz**2+bphi**2)

end function b_interpol_gpu
