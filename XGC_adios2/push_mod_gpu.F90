! collection of codes related to GPU
! try bottom up approach
!
#include "precision_mod_gpu.F90"


module util_mod_gpu
contains
#include "get_gpu_streamid.F90"
#include "setval_gpu.F90"
#include "sum2_gpu.F90"
#include "sum4_gpu.F90"
end module util_mod_gpu

! --------------------------------
! code related to particle sorting 
! --------------------------------
#include "thrust_mod.F90"
#include "reorder_gpu_mod.F90"
#include "gen_perm_gpu_mod.F90"



#include "dimensions_mod_gpu.F90"
#include "sml_module_gpu.F90"
#include "neu_module_gpu.F90"
#include "eq_module_gpu.F90"
#include "itp_module_gpu.F90"
#include "one_d_cub_mod_gpu.F90"
#include "external_bicub.F90"
#include "bicub_mod_gpu.F90"
#include "ptl_module_gpu.F90"
#include "grid_class_gpu.F90"
#include "psn_class_gpu.F90"
#include "boundary_class_gpu.F90"
#include "diag_module_gpu.F90"
#include "bnc_module_gpu.F90"

module push_mod_gpu
use cudafor
use util_mod_gpu

use dimensions_mod_gpu
use sml_module_gpu
use neu_module_gpu
use fld_module
use eq_module_gpu
use itp_module_gpu
use one_d_cub_mod_gpu
use bicub_mod_gpu
use ptl_module_gpu
use grid_class_gpu
use boundary_class_gpu
use psn_class_gpu
use diag_module_gpu
use bnc_module_gpu
use precision_mod_gpu
implicit none

type tbuf
   real (kind=work_p) :: ph(ptl_nphase)
   real (kind=work_p) :: ct(ptl_nconst)
#ifdef USE_TR_CHECK
   integer :: tr
#endif
end type tbuf

#ifdef USE_GPU_EMU      
      interface atomicAdd
        module procedure atomicAdd_d, atomicAdd_i
      end interface
#endif

contains

#ifdef USE_GPU_EMU
      attributes(device)                                                 &
     &real*8 function atomicAdd_d( x, value ) 
      real*8, intent(in) :: value
      real*8, intent(inout) :: x

      real*8 :: xold

      xold = x
      x = x + value
      atomicAdd_d = xold
      return
      end function atomicAdd_d


      attributes(device)                                                 &
     &integer function atomicAdd_i( x, value ) 
      integer, intent(in) :: value
      integer, intent(inout) :: x

      integer :: xold

      xold = x
      x = x + value
      atomicAdd_i = xold
      return
      end function atomicAdd_i

#endif


! -------------------------------------------------
! all code must be in the same module, a limitation
! of PGI cuda compiler
! -------------------------------------------------

#include "psi_interpol_gpu.F90"
#include "I_interpol_gpu.F90"
#include "psi_der_all_gpu.F90"
#include "field_gpu.F90"
#include "b_interpol_gpu.F90"
#include "bvec_interpol_gpu.F90"
#include "rho_mu_to_ev_pitch2_gpu.F90"
#include "remove_particle_gpu.F90"

#include "z_psi_min_gpu.F90"
#include "b_interpol_sym_gpu.F90"

#include "eq_ftn_gpu2.F90"
#include "eq_dftn_gpu2.F90"
#include "is_rgn1_gpu.F90"
#include "is_rgn12_gpu.F90"

#include "derivs_sp_elec_gpu.F90"
#include "efield_gpu.F90"
#include "diag_1d_port1_gpu.F90"
#include "diag_heat_port_gpu.F90"
#include "derivs_single_gpu.F90"
#include "restrict_weight_gpu.F90"
#include "derivs_single_with_e_gpu.F90"
#include "push_single_gpu.F90"
#include "pushe_single_gpu.F90"
#include "bounce_gpu.F90"

#include "derivs_gpu.F90"

#include "field_following_pos2_gpu.F90"


#include "init_push_mod_gpu.F90"
#include "guess_gpu.F90"
#include "push_update_host_gpu.F90"
#include "push_update_device_gpu.F90"
#include "pushe_gpu.F90"
#include "pushe_kernel_gpu.F90"
#include "gen_perm.F90"
#include "gen_perm_tri.F90"
#include "gen_perm_tri_gpu.F90"
#include "sheath_calculation_gpu.F90"
#include "ignore_factor_gpu.F90"

#ifdef USE_CALC_GRADIENT
#include "calc_gradxy_gpu.F90"
#include "calc_Er_Ez_gpu.F90"
#include "calc_E_para_gpu.F90"
#include "calc_E_phi_ff_gpu.F90"
#endif


end module push_mod_gpu



