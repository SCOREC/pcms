  ! function evaluation
  attributes(device) &
    function eq_dftn_gpu2(p,r,z,  ftn)
         
    use eq_module_gpu, only :  eq_x_z, eq_x_psi,eq_ftn_type
    use precision_mod_gpu
    ! use sml_module
    implicit none
    real(kind=work_p) :: eq_dftn_gpu2
    type(eq_ftn_type) :: ftn
    real(kind=work_p) :: p,r, z 
    real(kind=work_p) :: tmp


    ! region searching and enforcing region 3 value to x-point value
    if(is_rgn12_gpu(r,z,p)) then
       tmp=p
    else
       tmp=eq_x_psi
       eq_dftn_gpu2=0.0_work_p
       return
    endif


    select case(ftn%shape)
    case(0) ! constant
       eq_dftn_gpu2=0
    case(1) ! hyperbolic tanh
       eq_dftn_gpu2 = -ftn%sv(3)*ftn%sv(2)*(1.0_work_p/cosh((ftn%inx(1)-tmp)*ftn%sv(3)))**2
    case(2) ! linear
       eq_dftn_gpu2 = ftn%sv(1)
    case(3)
       eq_dftn_gpu2 = - ftn%iny(1)*exp(ftn%sv(1)*tanh((ftn%inx(1)-tmp)*ftn%sv(2))) &
            /((cosh((ftn%inx(1)-tmp)*ftn%sv(2)))**2*ftn%inx(3))
       
       
    case(-1)
       if(ftn%min < tmp .and. tmp < ftn%max) then
          eq_dftn_gpu2=-huge(0.0_work_p)
!          print *, "GPU code is not properly implemented for this case"
       else
          ! EFD not sure how to work around this
          ! call ftn_derivative(ftn_spl,tmp,1,eq_dftn_gpu2)
          eq_dftn_gpu2 = -0.0_work_p
          ! print*,'eqn_dftn_gpu2: ftn_derivative not implemented'
       endif
    case default
       ! print *, 'Invalid shape number in eq_dftn_gpu2', ftn%shape
       ! stop
    end select
    
    ! case  4 and case 5 does not have derivative function, yet. --> case 5 is required

  end function eq_dftn_gpu2
