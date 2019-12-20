       module reorder_gpu_mod
       use cudafor
       implicit none


       interface reorder_gpu
         module procedure ireorder_gpu, sreorder_gpu,                    &
     &          dreorder_gpu, lreorder_gpu
       end interface 

       interface reorder1d_gpu
         module procedure ireorder1d_gpu, sreorder1d_gpu,                    &
     &          dreorder1d_gpu, lreorder1d_gpu
       end interface 

       interface reorder2d_gpu
         module procedure ireorder2d_gpu, sreorder2d_gpu,                    &
     &          dreorder2d_gpu, lreorder2d_gpu
       end interface 

       contains

#include "lreorder_gpu.F90"
#include "ireorder_gpu.F90"
#include "sreorder_gpu.F90"
#include "dreorder_gpu.F90"

#include "lreorder1d_gpu.F90"
#include "ireorder1d_gpu.F90"
#include "sreorder1d_gpu.F90"
#include "dreorder1d_gpu.F90"

#include "lreorder2d_gpu.F90"
#include "ireorder2d_gpu.F90"
#include "sreorder2d_gpu.F90"
#include "dreorder2d_gpu.F90"

       end module reorder_gpu_mod
