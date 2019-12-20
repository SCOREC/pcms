module dimensions_mod_gpu
! integer, parameter :: ptl_maxnum_dim = 22*100*100 + 5001
! integer, parameter :: ptl_nummax_dim = ptl_maxnum_dim
! integer, parameter :: maxnum_dim = ptl_maxnum_dim
! integer, parameter :: nummax_dim = maxnum_dim
! integer, parameter :: lost_maxnum_dim = 100 + (ptl_maxnum_dim / 10)
! integer, parameter :: lost_nummax_dim = lost_maxnum_dim

! integer, parameter :: nnode_dim = 20694
! integer, parameter :: ntriangle_dim = 41087

! integer, parameter :: guess_table_n1_dim = 512
! integer, parameter :: guess_table_n2_dim = 512
! integer, parameter :: guess_list_size = 7005393

 ! -----------------------------------
 ! nseg_dim used in boundary_class_gpu
 ! -----------------------------------
 integer, parameter :: iseg_dim = 100

 ! ------------------------------
 ! nrho_dim used in psn_class_gpu
 ! ------------------------------
! integer, parameter :: nrho_dim = 8
! integer, parameter :: nphi_dim = 32

 ! -------------------------------------
 ! nthreads_dim used in diag_module_gpu
 ! -------------------------------------
#ifdef ORIGINAL
 integer, parameter :: grid_dim = 64
 integer, parameter :: block_dim = 64
#else
 integer, parameter :: grid_dim = 384
 integer, parameter :: block_dim = 64
#endif

 integer, parameter :: nthreads_dim = grid_dim * block_dim
 integer, parameter :: texture_max = 2**27

end module dimensions_mod_gpu

