#define USE_COPY_EQ_FTN_TYPE 1

module eq_module_gpu
use precision_mod_gpu
real (kind=work_p) :: eq_min_r,eq_max_r,eq_min_z,eq_max_z 
real (kind=work_p) :: eq_x_psi, eq_x_z, eq_axis_r, eq_axis_z, eq_x_slope, eq_x_r
real (kind=work_p) :: eq_x2_r, eq_x2_z, eq_x2_slope
attributes(constant) :: eq_min_r,eq_max_r,eq_min_z,eq_max_z 
attributes(constant) :: eq_x_psi, eq_x_z, eq_axis_r, eq_axis_z, eq_x_slope, eq_x_r
attributes(constant) :: eq_x2_r, eq_x2_z, eq_x2_slope


  ! data structure for equilbirum profile
  type eq_ftn_type
     integer :: shape
     real (kind=work_p):: inx(3), iny(3)
     real (kind=work_p):: sv(6)
     ! for arbitrary profile function - use pspline
     ! character (len=100) :: filename
     ! type (EZspline1_r8) :: spl
     real (kind=work_p) :: min, max
  end type eq_ftn_type

type(eq_ftn_type) :: eq_tempi, eq_tempe, eq_den, eq_flowi, eq_flowe

attributes(constant) :: eq_tempi, eq_tempe, eq_den, eq_flowi, eq_flowe

contains

#if USE_COPY_EQ_FTN_TYPE
  attributes(host) &
  subroutine copy_eq_ftn_type(eq_tempi_gpu, eq_tempi_host)
  ! --------------------------------------------------
  ! copy structure type(eq_ftn_type)  from host to gpu
  ! --------------------------------------------------
  use eq_module, only: eq_ftn_type_host => eq_ftn_type
  implicit none

  type(eq_ftn_type), intent(inout) :: eq_tempi_gpu
  type(eq_ftn_type_host),intent(in) :: eq_tempi_host

  attributes(constant) :: eq_tempi_gpu

  type(eq_ftn_type) :: p_eq_tempi


  p_eq_tempi%inx = eq_tempi_host%inx
  p_eq_tempi%iny = eq_tempi_host%iny
  p_eq_tempi%sv  = eq_tempi_host%sv
  p_eq_tempi%shape = eq_tempi_host%shape
  p_eq_tempi%min = eq_tempi_host%min
  p_eq_tempi%max = eq_tempi_host%max

  eq_tempi_gpu = p_eq_tempi


  return
  end subroutine copy_eq_ftn_type
#endif

  attributes(host) &
  subroutine update_device_eq()
  use sml_module, only : sml_mype

  use eq_module, only : &
  eq_ftn_type_host => eq_ftn_type, &
  eq_axis_r_host => eq_axis_r, &
  eq_axis_z_host => eq_axis_z, &
  eq_min_r_host => eq_min_r, &
  eq_max_r_host => eq_max_r, &
  eq_min_z_host => eq_min_z, &
  eq_max_z_host => eq_max_z, &
  eq_x_psi_host => eq_x_psi, &
  eq_x_r_host => eq_x_r,     &
  eq_x_z_host => eq_x_z,     &
  eq_tempi_host => eq_tempi, &
  eq_tempe_host => eq_tempe, &
  eq_den_host => eq_den,     &
  eq_flowi_host => eq_flowi, &
  eq_flowe_host => eq_flowe, &
  eq_x2_z_host => eq_x2_z,   &
  eq_x2_r_host => eq_x2_r,   &
  eq_x2_slope_host => eq_x2_slope,  &
  eq_x_slope_host => eq_x_slope

  implicit none

  integer, parameter :: idebug = 0

  eq_axis_r = eq_axis_r_host
  eq_min_r = eq_min_r_host
  eq_max_r = eq_max_r_host
  eq_min_z = eq_min_z_host
  eq_max_z = eq_max_z_host
  eq_x_psi = eq_x_psi_host
  eq_x_r = eq_x_r_host
  eq_x_z = eq_x_z_host
  eq_x2_r = eq_x2_r_host
  eq_x2_z = eq_x2_z_host
  eq_x2_slope = eq_x2_slope_host
  eq_x_slope = eq_x_slope_host

  if (idebug >= 1) then
!    if (sml_mype == 0) print*,'before copy_eq_ftn_type()'
  endif

#if USE_COPY_EQ_FTN_TYPE
  call copy_eq_ftn_type(eq_tempi, eq_tempi_host)
  call copy_eq_ftn_type(eq_tempe, eq_tempe_host)
  call copy_eq_ftn_type(eq_den, eq_den_host)
  call copy_eq_ftn_type(eq_flowi, eq_flowi_host)
  call copy_eq_ftn_type(eq_flowe, eq_flowe_host)
#else
  eq_tempi%inx = eq_tempi_host%inx
  eq_tempi%iny = eq_tempi_host%iny
  eq_tempi%sv  = eq_tempi_host%sv
  eq_tempi%shape = eq_tempi_host%shape
  eq_tempi%min = eq_tempi_host%min
  eq_tempi%max = eq_tempi_host%max


  eq_tempi%inx = eq_tempi_host%inx
  eq_tempi%iny = eq_tempi_host%iny
  eq_tempi%sv  = eq_tempi_host%sv
  eq_tempi%shape = eq_tempi_host%shape
  eq_tempi%min = eq_tempi_host%min
  eq_tempi%max = eq_tempi_host%max

  
  eq_den%inx = eq_den_host%inx
  eq_den%iny = eq_den_host%iny
  eq_den%sv  = eq_den_host%sv
  eq_den%shape = eq_den_host%shape
  eq_den%min = eq_den_host%min
  eq_den%max = eq_den_host%max

  eq_den%inx = eq_den_host%inx
  eq_den%iny = eq_den_host%iny
  eq_den%sv  = eq_den_host%sv
  eq_den%shape = eq_den_host%shape
  eq_den%min = eq_den_host%min
  eq_den%max = eq_den_host%max

  eq_flowi%inx = eq_flowi_host%inx
  eq_flowi%iny = eq_flowi_host%iny
  eq_flowi%sv  = eq_flowi_host%sv
  eq_flowi%shape = eq_flowi_host%shape
  eq_flowi%min = eq_flowi_host%min
  eq_flowi%max = eq_flowi_host%max


  eq_flowe%inx = eq_flowe_host%inx
  eq_flowe%iny = eq_flowe_host%iny
  eq_flowe%sv  = eq_flowe_host%sv
  eq_flowe%shape = eq_flowe_host%shape
  eq_flowe%min = eq_flowe_host%min
  eq_flowe%max = eq_flowe_host%max

#endif

  end subroutine update_device_eq


end module eq_module_gpu
