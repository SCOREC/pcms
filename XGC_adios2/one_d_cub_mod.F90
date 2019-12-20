
!module for 1D I interpolation without pspline
module one_d_cub_mod
   implicit none

   !real (8), allocatable, dimension(:) :: one_d_cub_psi
   real (8), allocatable, dimension(:,:) :: one_d_cub_psi_acoef
   real (8) :: one_d_cub_dpsi_inv, one_d_cub_psimin



contains
   subroutine setup_psi_one_d_cub(eq_mpsi, eq_psi_grid)
      use sml_module, only : sml_mype
      use EZspline_obj
      use EZspline
      use itp_module, only : spl_psi

      implicit none
      integer, intent(in) :: eq_mpsi
      real (8), dimension(eq_mpsi), intent(in) :: eq_psi_grid

      integer, parameter :: ndeg = 3
      real (8), dimension(0:ndeg) :: acoef
      real (8), dimension(0:ndeg) :: fval

      real (8) :: psi, dpsi
      integer :: npsi, i, j, ideriv

      integer :: ier
      real (kind=8) :: r8value

      !test
      integer, parameter :: test_n=320
      real (kind=8) :: psp_value0(test_n), oned_value0(test_n)
      real (kind=8) :: psp_value1(test_n), oned_value1(test_n)
      real (kind=8) :: ratio_value0(test_n)
      real (kind=8) :: ratio_value1(test_n), inputpsi(test_n)



      npsi = eq_mpsi-1
      dpsi = eq_psi_grid(2)-eq_psi_grid(1)

      one_d_cub_dpsi_inv=dble(1)/dpsi
      one_d_cub_psimin = eq_psi_grid(1)

      allocate(one_d_cub_psi_acoef(0:ndeg, npsi))

      do i=1, npsi
         do j=0, 1 !j=0 : left point, j=1 : right point
            psi=eq_psi_grid(i+j)
            do ideriv=0, 1 !k=0 : I value, k=1 : dI/dpsi
               call EZspline_derivative(spl_psi,ideriv,psi,r8value,ier)
               call EZspline_error(ier)
               fval(2*j+ideriv)=r8value
            enddo
         enddo
         call gencoef_one_d_cub(dpsi,fval,acoef)
         one_d_cub_psi_acoef(:,i)=acoef
      enddo

      !test routine
      do i=1,test_n
         inputpsi(i)=real(i)/real(test_n)*eq_psi_grid(eq_mpsi)
         ideriv=0
         call EZspline_derivative(spl_psi,ideriv,inputpsi(i),r8value,ier)
         call EZspline_error(ier)
         psp_value0(i)=r8value

         call I_interpol_wo_pspline(inputpsi(i), ideriv, r8value)
         oned_value0(i)=r8value

         ideriv=1
         call EZspline_derivative(spl_psi,ideriv,inputpsi(i),r8value,ier)
         call EZspline_error(ier)
         psp_value1(i)=r8value

         call I_interpol_wo_pspline(inputpsi(i), ideriv, r8value)
         oned_value1(i)=r8value

         ratio_value0(i)=psp_value0(i)/oned_value0(i)
         ratio_value1(i)=psp_value1(i)/oned_value1(i)
      enddo
   !   if (sml_mype==0) then
   !      open(598, file='ondd_cub_test.txt', status='replace')
   !      do i=1,300
   !         write(598,2100) psp_value0




    !     2100 format(30(e19.13,' '))
    !  endif


   end subroutine setup_psi_one_d_cub


   subroutine gencoef_one_d_cub(dpsi,fval,acoef)
      implicit none
      !generate coefficients of 1D cubic polynomials
      !cf) http://mathworld.wolfram.com/CubicSpline.html

      integer, parameter :: ndeg = 3
      real (8),intent(in) :: dpsi
      real (8),dimension(0:ndeg),intent(in) :: fval
      real (8),dimension(0:ndeg),intent(inout)::acoef

      acoef(0)=fval(0)
      acoef(1)=fval(1)*dpsi
      acoef(2)=3D0*(fval(2)-fval(0))-(2D0*fval(1)+fval(3))*dpsi
      acoef(3)=2D0*(fval(0)-fval(2))+(fval(1)+fval(3))*dpsi

   end subroutine gencoef_one_d_cub

   subroutine I_interpol_wo_pspline(psi, ideriv, ivalue)
      implicit none

      integer, parameter :: ndeg = 3
      real (8), intent(in) :: psi
      integer, intent(in) :: ideriv
      real (8), intent(inout) :: ivalue

      real (8) :: pn, wp, acoef(0:ndeg)
      integer :: ip

      pn=psi*one_d_cub_dpsi_inv
      ip=floor(pn)+1
      ip=min(max(ip,1),ubound(one_d_cub_psi_acoef,2))
      wp=pn-real(ip-1,8)

      acoef=one_d_cub_psi_acoef(:,ip)

      if (ideriv==0) then
         ivalue=acoef(0)+(acoef(1)+(acoef(2)+acoef(3)*wp)*wp)*wp
      elseif (ideriv==1) then
         ivalue=(acoef(1)+(2D0*acoef(2)+3D0*acoef(3)*wp)*wp)*one_d_cub_dpsi_inv
      else
         print *, 'ideriv in I_interpol_wo_pspline should be 0 or 1'
         stop
      endif
   end subroutine I_interpol_wo_pspline


end module one_d_cub_mod
