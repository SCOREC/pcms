attributes(device) &
  subroutine derivs_gpu(x,phi,dx)
    use sml_module_gpu
    use precision_mod_gpu

    implicit none
    real (kind=work_p), intent(in)  :: x(2), phi
    real (kind=work_p), intent(out) :: dx(2)
    real (kind=work_p) :: b(2), bphi, r,z

    r=min(max(x(1),sml_bd_min_r),sml_bd_max_r)
    z=min(max(x(2),sml_bd_min_z),sml_bd_max_z)


    call bvec_interpol_gpu(r,z,phi,b(1),b(2),bphi)
    dx = b/bphi*x(1)

  end subroutine derivs_gpu

