.SUFFIXES: .o .F90

%.o: %.mod

#include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/variables

#USER_LIB_DIR=/project/projectdirs/m499/rhager/software/cori_haswell/intel
#PETSC_DIR=$(USER_LIB_DIR)/petsc-cori-haswell-opt64-intel
#include ${PETSC_DIR}/lib/petsc/conf/variables
#PETSC_INC=$(PETSC_FC_INCLUDES)
#PETSC_LIB=$(PETSC_KSP_LIB)

OBJ= module.o search.o module_psn.o pol_decomp.o f0module.o adios_comm_mod.o adios2_comm_mod.o \
	elliptics.o collisionf.o lbal_mod.o \
	$(EXTRA_OBJ) psmooth.o bicub_mod.o one_d_cub_mod.o initial_perturbation_GM.o \
        interfaces.o charge.o diagnosis.o poisson_extra_xgc1.o qevaluateandtrapped.o setup.o   \
        read.o gen_perm.o sort_particles.o push.o pushe.o load.o \
	 efield.o interpolation.o $(MPI_OBJ) \
	limiter.o bounce.o diagnosis2.o collision.o collision2.o  \
        collisionf2.o  heat.o \
	turbulence.o neutral.o  neutral2.o  neutral3.o linearsolver.o fem2d.o \
	fem_ops.o poisson_extra_common.o petsc_solve.o
#coupling_core_edge.o
#OBJ_EM=module-em.o bicub_mod.o one_d_cub_mod.o  search-em.o pol_decomp.o f0module.o em_advance_petsc.o  interfaces-em.o  em_main.o em_hyb.o em_poisson.o charge-em.o diagnosis-em.o poisson_extra_xgc1-em.o collisionf.o setup-em.o 
OBJ_EM=em_main.o em_hyb.o em_poisson.o em_advance_petsc.o 

#OBJ_ES=module.o bicub_mod.o one_d_cub_mod.o search.o pol_decomp.o f0module.o interfaces.o es_main.o es_poisson.o charge.o diagnosis.o poisson_extra_xgc1.o collisionf.o setup.o 
OBJ_ES=es_main.o es_poisson.o  

OBJ_GPU=push_mod_gpu.o cscan.o csort.o cbicub.o 

#
# These are added into EXTRA_OBJ & MPI_OBJ in makefiles, funny dependancy
#
IMSL_OBJ=my_imsl.o
PORT_OBJ=bspline90_22.o taus88.o derf.o datanh.o pppack.o fmin.o 
SER_OBJ=mpisingle.o
PAR_OBJ=mpi.o

GPU_SRC =   push_mod_gpu.F90 \
        dimensions_mod_gpu.F90  \
        precision_mod_gpu.F90  \
        sml_module_gpu.F90  \
        eq_module_gpu.F90  \
        itp_module_gpu.F90  \
        one_d_cub_mod_gpu.F90  \
        bicub_mod_gpu.F90  \
        ptl_module_gpu.F90  \
        grid_class_gpu.F90  \
        boundary_class_gpu.F90  \
        psn_class_gpu.F90  \
        diag_module_gpu.F90  \
        bnc_module_gpu.F90  \
        push_update_device_gpu.F90  \
        psi_interpol_gpu.F90  \
        I_interpol_wo_pspline_gpu.F90  \
        I_interpol_gpu.F90  \
        psi_der_all_gpu.F90  \
        field_gpu.F90  \
        b_interpol_gpu.F90  \
        rho_mu_to_ev_pitch2_gpu.F90  \
        remove_particle_gpu.F90  \
        efield_gk_gpu.F90  \
        derivs_sp_gpu.F90  \
        efield_gpu.F90  \
        diag_1d_port1_gpu.F90  \
        derivs_single_gpu.F90  \
        derivs_single_with_e_gpu.F90  \
        push_single_gpu.F90  \
        z_psi_min_gpu.F90  \
        b_interpol_sym_gpu.F90  \
        bounce_gpu.F90  \
        push_kernel_gpu.F90  \
        push_update_host_gpu.F90  \
        guess_gpu.F90 \
        search_tr2_gpu.F90 \
        bvec_interpol_gpu.F90 \
        derivs_gpu.F90 \
        field_following_pos2_gpu.F90 \
        pushe_gpu.F90 \
        pushe_kernel_gpu.F90 \
        t_coeff_gpu.F90 \
        dreorder1d_gpu.F90 \
        dreorder2d_gpu.F90 \
        dreorder_gpu.F90 \
        ireorder1d_gpu.F90 \
        ireorder2d_gpu.F90 \
        ireorder_gpu.F90 \
        lreorder1d_gpu.F90 \
        lreorder2d_gpu.F90 \
        lreorder_gpu.F90 \
        reorder_gpu_mod.F90 \
        sreorder1d_gpu.F90 \
        sreorder2d_gpu.F90 \
        sreorder_gpu.F90 \
        gen_perm_gpu.F90 \
        gen_perm_gpu_mod.F90 \
        gen_perm_gpu_pass1.F90 \
        gen_perm_gpu_pass2.F90 \
        isetval_gpu.F90 \
        isum_col_gpu.F90 \
        iprefix_sum_gpu.F90 \
        setup_lcount_gpu.F90

