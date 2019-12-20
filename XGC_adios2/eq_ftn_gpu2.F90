  ! function evaluation
  attributes(device) &
    function eq_ftn_gpu2(p,r,z,ftn)
    use eq_module_gpu, only : eq_x_psi, eq_x_z,eq_ftn_type
    use precision_mod_gpu
    !use sml_module
    implicit none
    real(kind=work_p) :: eq_ftn_gpu2
    type(eq_ftn_type) :: ftn
    real(kind=work_p) :: p,r, z 


    real(kind=work_p) :: tmp, tmp2, tmp3, tmp4

   if(is_rgn12_gpu(r,z,p)) then
       tmp=p
    else
       tmp=eq_x_psi
    endif

    !tmp=min(max(p,sml_pmin), sml_pmax)

    select case(ftn%shape)
    case(0) ! constant
       eq_ftn_gpu2 = ftn%iny(1)
    case(1) ! hyperbolic tanh
       eq_ftn_gpu2 = ftn%sv(2)*tanh((ftn%inx(1)-tmp)*ftn%sv(3))+ftn%sv(1)
    case(2) ! linear
       eq_ftn_gpu2 = ftn%sv(1)*tmp + ftn%sv(2)
    case(3) ! a exp(-B w tanh( (r-r0)/w ))
       eq_ftn_gpu2 = ftn%iny(1)*exp(ftn%sv(1)*tanh((ftn%inx(1)-tmp)*ftn%sv(2)))
    case(4) 
       eq_ftn_gpu2 = ftn%sv(2)*tanh((ftn%inx(1)-tmp)*ftn%sv(3))+ftn%sv(1)
       if(tmp < ftn%inx(1)-0.5_work_p*ftn%inx(2) ) then
          eq_ftn_gpu2 = eq_ftn_gpu2 + ftn%sv(4)*sqrt(tmp) + ftn%sv(5)
       endif
    case(5)
       tmp2=ftn%sv(3)*( ftn%inx(1)-sqrt(ftn%inx(1)*tmp) )  ! z
       tmp3=exp(tmp2)                           ! expz
       tmp4=exp(-tmp2)                          ! expmz
       ! A * ( (1+z*slope)*expz - expmz )/(expz + expmz) + B
       eq_ftn_gpu2 = ftn%sv(2)*( (1+tmp2*ftn%sv(4))*tmp3 - tmp4 )/(tmp3+tmp4) + ftn%sv(1)
    case(-1)
       ! print*,'shape number -1 not implemented in eq_ftn_gpu2'
       ! tmp=min(max(tmp,ftn_min),ftn_max)
       ! call ftn_evaluation(ftn_spl,tmp,eq_ftn_gpu2)
       tmp=min(max(tmp,ftn%min),ftn%max)
       eq_ftn_gpu2 = huge(0.0_work_p)
!       print *, "GPU code is not properly implemented for this case"
    case default
       ! print *, 'Invalid shape number in eq_ftn_gpu2', ftn%shape
       ! stop
    end select
    
  end function eq_ftn_gpu2
