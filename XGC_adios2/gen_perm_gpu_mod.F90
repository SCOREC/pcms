      module gen_perm_gpu_mod
      use cudafor
      use precision_mod_gpu
      implicit none

      contains

#include "isetval_gpu.F90"
#include "isum_col_gpu.F90"
#include "iprefix_sum_gpu.F90"
#include "setup_lcount_gpu.F90"

#include "gen_perm_gpu_pass1.F90"
#include "gen_perm_gpu_pass2.F90"
#include "gen_perm_gpu.F90"

       end module gen_perm_gpu_mod
