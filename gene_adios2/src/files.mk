###########################################################################
### The files and directories                                           ###
###########################################################################

###########################################################################
### The source files	                                                ###
###########################################################################

F77SRC  =

F90PSRC = \
	adiabatic_response.F90\
	all_rhs_terms.F90\
	antenna.F90\
	arrays.F90\
	aux_fields.F90\
	aux_func.F90\
	axpy.F90\
	blockindex.F90\
	boundary.F90\
	boundary_exchange_general.F90\
	boundary_exchange_vw.F90\
	boundary_exchange_x.F90\
	boundary_exchange_z.F90\
	BoundaryDescription.F90\
	box_data_module.F90\
	calc_rhs.F90\
	charge_curr_dens.F90\
	chease.F90\
	check_parameters.F90\
	checkpoint.F90\
	circular.F90\
	codemods.F90\
	collisions.F90\
	comm.F90\
	compute_f.F90\
	convergence_monitoring.F90\
	coordinates.F90\
	coordinates_adptv.F90\
	dchidxy_terms.F90\
	dchidz_term.F90\
	debug_output.F90\
	dfdxy_terms.F90\
	dfdzv_terms.F90\
	dzv_terms.F90\
	dgdxy_terms.F90\
	discretization.F90\
	discretization_adptv.F90\
	diag.F90\
	diag_svd.F90 \
        eigen_direct.F90\
	diag_aux.F90\
	diag_df.F90\
	diag_finit.F90\
	diag_fmsvd.F90\
	diag_vsp.F90\
	eigenvalue_comp.F90\
	eigen_parameters.F90\
	equilibrium_fields.F90\
	external_contr.F90\
	exchange_z_ff.F90\
	f0_term.F90\
	field_solver.F90\
	field_solve_ff.F90\
	field_solve_df.F90 \
	file_io.F90 \
	flr_corr.F90 \
	fourier_$(FFTLIB).F90\
	full_rhs_mat.F90\
	fullmatrix_aux.F90\
	gauss.quadrature.F90\
	gene_scan.F90\
	gene_subroutine.F90\
	geometry.F90\
	grid1d.F90\
	gyro_average_ff.F90\
	gyro_average_df.F90 \
        hybrid.F90\
	Gyro_LES.F90 \
	initial_value_comp.F90\
	init_cond.F90\
	init_physical.F90\
	lag_interp.F90 \
	localpolynombase.F90 \
	localpolynombase2d.F90 \
	miller_geometry.F90 \
	mtrandom.F90\
	neo_equilibrium.F90\
	numerical_damping.F90\
	parallel_nonlin.F90\
	parameters_IO.F90\
	par_geom.F90\
	par_in.F90\
	par_mod.F90\
	par_other.F90\
	perf_opt.F90\
	petsc_precond.F90\
	prefactors.F90\
        profile_IO.F90\
	profile_smoothing.F90\
	profiles.F90\
	poly_roots.F90\
	quasilin_model.F90\
	reset_mod.F90\
	rhs_computations.F90\
	RK_coefficients.F90\
	sources_mod.F90\
	spatial_averages.F90\
	spectype.F90\
	spline_interp.F90\
	compute_dt.F90\
	island.F90\
	time_averages.F90\
	time_scheme.F90\
	tracer.F90\
	tracer_aux.F90\
	tracer_IO.F90\
	tracer_rk5.F90\
	tracer_rk5_util.F90\
	vel_space.F90\
	x_derivatives.F90 \
	rhs_term.F90 \
	nonlinear_term.F90 \
	df_nonlinear_term.F90 \
	ff_nonlinear_term.F90 \
	ff_nonlinear_term_h.F90 \
	df_arakawa_nonlinear_term.F90 \
	df_nonlinear_prefactor.F90 \
	df_epc_nl_prefactor.F90\
	pinning.F90\
	par_poloidal_planes.F90\
	poloidal_planes_aux.F90\
	poloidal_planes_nonuni.F90\
	poloidal_planes.F90\
	create_nonuni_grid.F90\
	import_XGC_grid.F90\
	boundary_exchange_XGC.F90


ifneq ($(COUPLE_XGC),yes)

F90PSRC += \
	diag_neoclass.F90\
	diag_Gyro_LES.F90\
	diag_Gyro_LES_common.F90\
	diag_Gyro_LES_fetime.F90\
	diag_Gyro_LES_spec2d.F90\
        diag_Gyro_LES_real_3d.F90\
	diag_Gyro_LES_transfer.F90\
	diag_Gyro_LES_transfer3d.F90\
	diag_energy.F90\
	diag_extended.F90\
	diag_triplet.F90\
	diag_fsa_moments.F90\
	stellopt_interface.F90

else

    F90PSRC += diag_energy.F90

endif


CSRC = vectorized_kernels.c

ifeq ($(SLEPC),yes)
F90PSRC += eigen_iterative.F90
F90PSRC += impl_petsc.F90
F90PSRC += petsc_aux.F90
F90PSRC += slepc_aux.F90
endif

ifeq ($(SCALAPACK),yes)
F90PSRC += impl_scalapack.F90
F90PSRC += scalapack_aux.F90
#F90PSRC += matrix_module.F90
F90PSRC += matrix.F90
F90PSRC += storefullmatrix.F90
F90PSRC += processgrid.F90
F90PSRC += derivative_matrix.F90
F90PSRC += BandedMatrix.F90
F90PSRC += Vector.F90
F90PSRC += storebandedmatrix.F90
F90PSRC += storevector.F90
F90PSRC += ListObject.F90
else
F90PSRC += dummy_matrix_module.F90
endif

ifeq ($(WITH_CUTILS),no)
F90PSRC += dummy_c_utils.F90
else
CSRC += c_utils.c
endif

ifeq ($(USE_C_NONLIN),yes)
F90PSRC += df_nonlinear_c.F90
CSRC += c_nonlinearity.c
CSRC += c_fourier_fftw.c
endif

ifeq ($(USE_SPLIT_NONLIN),yes)
F90PSRC += df_split_nonlinear_term.F90
endif

ifeq ($(USE_CUDA_NONLIN),yes)
F90PSRC += df_nonlinear_cuda.F90
F90PSRC += page_locked_memory.F90
CUDASRC += cuda_nonlinearity.cu
CUDASRC += cuda_fourier_cufft.cu
CUDASRC += cuda_overlap.cu
CUDASRC += cuda_general.cu
endif

ifeq ($(USE_MIC_NONLIN),yes)
F90PSRC += df_mic_nonlinear_term.F90
F90PSRC += fourier_mkl_mic.F90
endif

ifeq ($(OPENMP),yes)
	ifeq ($(USE_OMPXY_NONLIN),yes)
		F90PSRC += df_ompxy_nonlinear_term.F90
	else
		F90PSRC += df_omp_nonlinear_term.F90
	endif
endif

ifeq ($(USE_PERFLIB),vtune)
F90PSRC += itt_api.F90
CSRC += itt_from_fortran.c
endif


MPIDEP += BandedMatrix.F90\
	BoundaryDescription.F90\
	boundary_exchange_general.F90\
	boundary.F90\
	checkpoint.F90\
	comm.F90\
	df_arakawa_nonlinear_term.F90\
	df_mic_nonlinear_term.F90\
	df_nonlinear_term.F90\
	df_nonlinear_term.F90\
	df_omp_nonlinear_term.F90\
	diag.F90\
	diag_svd.F90\
	diag_fmsvd.F90\
	dummy_matrix_module.F90\
	ff_nonlinear_term.F90\
	ff_nonlinear_term_h.F90\
	field_solve_df.F90\
	FR_perf.F90\
	gene.F90\
	gene_scan.F90\
	gyro_average_df.F90\
	libITMgene.F90\
	libtringene.F90\
	matrix.F90\
	parallel_nonlin.F90\
	storebandedmatrix.F90\
	storefullmatrix.F90\
	storevector.F90\
	test_libtringene.F90\
	test_libtringene.F90\
	Vector.F90\
	boundary_exchange_XGC.F90



ifeq ($(ADIOS),yes)
    F90PSRC += adios_io.F90 adios_comm_mod.F90
    MPIDEP += adios_io.F90 poloidal_plane.F90
    XMLFILE = $(SRCDIR)/gene.xml
endif


ifeq ($(ADIOS2),yes)
    F90PSRC += adios2_comm_mod.F90
    MPIDEP  += poloidal_plane.F90
endif

ifeq ($(COUPLE_XGC),yes)
    F90_TMP = coupling_core_gene.F90 coupling_base_setup.F90
    F90PSRC += $(F90_TMP)
    MPIDEP  += $(F90_TMP)
endif

# If GPTL isn't available perf_monitor just doesn't do anything
F90PSRC += perf_monitor.F90
MPIDEP += perf_monitor.F90


F90FULLSRC = $(addprefix $(SRCDIR)/,$(F90PSRC))
F90PRENAMES = $(addprefix $(PPDIR)/,$(subst .F90,.f90,$(F90PSRC)))
F90POBJ = $(addprefix $(OBJDIR)/,$(subst .F90,.o,$(F90PSRC)))
F77OBJ  = $(addprefix $(OBJDIR)/,$(subst .f,.o,$(F77SRC)))
COBJ = $(addprefix $(OBJDIR)/,$(subst .c,.o,$(CSRC)))
CUDAOBJ = $(addprefix $(OBJDIR)/,$(subst .cu,.o,$(CUDASRC)))
MPIOBJ = $(addprefix $(OBJDIR)/,$(subst .F90,.o,$(MPIDEP)))


#VPATH = $(SRCDIR)
OBJLIST = $(F77OBJ) $(F90POBJ)
COBJLIST = $(COBJ) $(CUDAOBJ)

RPL_STR =-WF,
CPREPROC += $(subst $(RPL_STR),,$(PREPROC))

#For libraries
LIBELEMENTS = $(OBJLIST) $(COBJLIST) $(OBJDIR)/$(LIBN).o
LIBN2=$(subst lib,,$(LIBN))

# to make some conditionals in the rules simpler, we map the case
# with USE_PERFLIB=none on the case when USE_PERFLIB is undefined or empty
USE_PERFLIB := $(if $(findstring none,$(USE_PERFLIB)),,$(strip $(USE_PERFLIB)))


###########################################################################
### The dependencies  (might be checked/updated using ../tools/makedep) ###
###########################################################################

$(OBJLIST) $(COBJLIST):		$(SRCDIR)/redef.h \
			$(SRCDIR)/intrinsic_sizes.h \
			$(BASEDIR)/makefile $(MKFILE) \
			$(MACHMKFILE)


$(OBJDIR)/adiabatic_response.o:	$(SRCDIR)/switches.h \
			$(OBJDIR)/discretization.o\
			$(OBJDIR)/par_in.o\
			$(OBJDIR)/par_other.o\
			$(OBJDIR)/coordinates.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/geometry.o\
			$(OBJDIR)/aux_func.o\
			$(OBJDIR)/gyro_average_ff.o\
			$(OBJDIR)/spatial_averages.o\
			$(OBJDIR)/hybrid.o

$(OBJDIR)/all_rhs_terms.o: $(OBJDIR)/nonlinear_term.o

$(OBJDIR)/antenna.o:	$(OBJDIR)/par_mod.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/coordinates.o\
			$(OBJDIR)/blockindex.o\
			$(OBJDIR)/vel_space.o\
			$(OBJDIR)/gyro_average_ff.o\
			$(OBJDIR)/mtrandom.o

$(OBJDIR)/arrays.o:	$(OBJDIR)/par_mod.o\
			$(OBJDIR)/discretization.o

$(OBJDIR)/aux_fields.o:	$(SRCDIR)/switches.h \
			$(OBJDIR)/par_mod.o\
			$(OBJDIR)/field_solve_ff.o\
			$(OBJDIR)/field_solve_df.o\
			$(OBJDIR)/field_solver.o\
			$(OBJDIR)/gyro_average_df.o\
			$(OBJDIR)/gyro_average_ff.o\
			$(OBJDIR)/prefactors.o\
			$(OBJDIR)/boundary.o\
			$(OBJDIR)/charge_curr_dens.o\
			$(OBJDIR)/compute_f.o\
			$(OBJDIR)/geometry.o\
			$(OBJDIR)/fourier_$(FFTLIB).o

$(OBJDIR)/aux_func.o:

$(OBJDIR)/axpy.o:	$(OBJDIR)/comm.o\
			$(OBJDIR)/discretization.o

$(OBJDIR)/blockindex.o:	$(OBJDIR)/discretization.o

$(OBJDIR)/boundary.o:	$(OBJDIR)/par_mod.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/boundary_exchange_general.o\
			$(OBJDIR)/boundary_exchange_x.o\
			$(OBJDIR)/boundary_exchange_z.o\
			$(OBJDIR)/boundary_exchange_vw.o\
			$(OBJDIR)/BoundaryDescription.o\
			$(OBJDIR)/box_data_module.o\
			$(OBJDIR)/exchange_z_ff.o\
			$(OBJDIR)/blockindex.o

$(OBJDIR)/BoundaryDescription.o:

$(OBJDIR)/boundary_exchange_general.o:	$(OBJDIR)/BoundaryDescription.o

$(OBJDIR)/boundary_exchange_vw.o:	$(OBJDIR)/BoundaryDescription.o\
					$(OBJDIR)/boundary_exchange_general.o

$(OBJDIR)/boundary_exchange_x.o:	$(OBJDIR)/discretization.o\
			$(OBJDIR)/BoundaryDescription.o\
			$(OBJDIR)/boundary_exchange_general.o\
			$(OBJDIR)/par_other.o

$(OBJDIR)/boundary_exchange_z.o:	$(OBJDIR)/BoundaryDescription.o\
			$(OBJDIR)/boundary_exchange_general.o\
			$(OBJDIR)/box_data_module.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/fourier_$(FFTLIB).o\
			$(OBJDIR)/discretization.o\
			$(OBJDIR)/par_in.o\
			$(OBJDIR)/par_other.o

$(OBJDIR)/box_data_module.o:	$(OBJDIR)/grid1d.o

$(OBJDIR)/calc_rhs.o:	$(SRCDIR)/switches.h \
			$(OBJDIR)/aux_fields.o\
			$(OBJDIR)/par_mod.o\
			$(OBJDIR)/prefactors.o\
			$(OBJDIR)/boundary.o\
			$(OBJDIR)/parallel_nonlin.o\
			$(OBJDIR)/f0_term.o\
			$(OBJDIR)/dchidxy_terms.o\
			$(OBJDIR)/dgdxy_terms.o\
			$(OBJDIR)/dfdxy_terms.o\
			$(OBJDIR)/dchidz_term.o\
			$(OBJDIR)/dfdzv_terms.o\
			$(OBJDIR)/dzv_terms.o\
			$(OBJDIR)/collisions.o\
			$(OBJDIR)/aux_fields.o\
			$(OBJDIR)/gyro_average_df.o\
			$(OBJDIR)/blockindex.o\
			$(OBJDIR)/sources_mod.o\
			$(OBJDIR)/RK_coefficients.o\
			$(OBJDIR)/axpy.o\
			$(OBJDIR)/numerical_damping.o\
			$(OBJDIR)/comm.o \
			$(OBJDIR)/nonlinear_term.o \
			$(OBJDIR)/df_nonlinear_term.o \
			$(OBJDIR)/ff_nonlinear_term.o \
			$(OBJDIR)/ff_nonlinear_term_h.o \
			$(OBJDIR)/df_arakawa_nonlinear_term.o \
			$(OBJDIR)/x_derivatives.o\
			$(OBJDIR)/all_rhs_terms.o

$(OBJDIR)/codemods.o:	$(OBJDIR)/discretization.o\
			$(OBJDIR)/par_in.o\
			$(OBJDIR)/par_other.o\
			$(OBJDIR)/file_io.o

$(OBJDIR)/nonlinear_term.o:	$(OBJDIR)/discretization.o \
				$(OBJDIR)/par_other.o \
				$(OBJDIR)/rhs_term.o

$(OBJDIR)/df_nonlinear_term.o:	$(SRCDIR)/switches.h\
				$(OBJDIR)/nonlinear_term.o \
				$(OBJDIR)/df_nonlinear_prefactor.o \
				$(OBJDIR)/df_epc_nl_prefactor.o \
				$(OBJDIR)/par_mod.o \
				$(OBJDIR)/par_other.o \
				$(OBJDIR)/comm.o \
				$(OBJDIR)/discretization.o\
				$(OBJDIR)/fourier_$(FFTLIB).o


ifeq ($(OPENMP),yes)
ifeq ($(USE_OMPXY_NONLIN),yes)
$(OBJDIR)/calc_rhs.o:	$(OBJDIR)/df_ompxy_nonlinear_term.o
else
$(OBJDIR)/calc_rhs.o:	$(OBJDIR)/df_omp_nonlinear_term.o
endif

$(OBJDIR)/df_omp_nonlinear_term.o:	$(OBJDIR)/df_nonlinear_term.o\
			$(OBJDIR)/par_mod.o\
			$(OBJDIR)/par_other.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/discretization.o\
			$(OBJDIR)/fourier_$(FFTLIB).o

$(OBJDIR)/df_ompxy_nonlinear_term.o:	$(OBJDIR)/df_nonlinear_term.o\
					$(OBJDIR)/par_mod.o\
					$(OBJDIR)/par_other.o\
					$(OBJDIR)/comm.o\
					$(OBJDIR)/discretization.o\
					$(OBJDIR)/fourier_$(FFTLIB).o
endif

$(OBJDIR)/df_nonlinear_prefactor.o:	$(OBJDIR)/discretization.o \
					$(OBJDIR)/par_mod.o \
					$(OBJDIR)/par_in.o \
					$(OBJDIR)/geometry.o\
					$(OBJDIR)/vectorized_kernels.o

$(OBJDIR)/df_nonlinear_prefactor_dummy.o:	$(OBJDIR)/discretization.o \
						$(OBJDIR)/df_nonlinear_prefactor.o

$(OBJDIR)/df_epc_nl_prefactor.o:	$(OBJDIR)/df_nonlinear_prefactor.o\
			$(OBJDIR)/par_in.o\
			$(OBJDIR)/blockindex.o\
			$(OBJDIR)/par_other.o\
			$(OBJDIR)/geometry.o\
			$(OBJDIR)/par_mod.o\
			$(OBJDIR)/coordinates.o

$(OBJDIR)/df_arakawa_nonlinear_term.o: $(SRCDIR)/switches.h\
				$(SRCDIR)/redef.h\
				$(SRCDIR)/intrinsic_sizes.h\
				$(OBJDIR)/df_nonlinear_term.o\
				$(OBJDIR)/par_mod.o \
				$(OBJDIR)/coordinates.o \
				$(OBJDIR)/comm.o \
				$(OBJDIR)/discretization.o\
				$(OBJDIR)/blockindex.o\
				$(OBJDIR)/fourier_$(FFTLIB).o \
				$(OBJDIR)/prefactors.o\
				$(OBJDIR)/x_derivatives.o\
				$(OBJDIR)/geometry.o


$(OBJDIR)/ff_nonlinear_term.o:	$(SRCDIR)/switches.h\
				$(OBJDIR)/nonlinear_term.o \
				$(OBJDIR)/par_mod.o \
				$(OBJDIR)/par_other.o \
				$(OBJDIR)/par_in.o \
				$(OBJDIR)/coordinates.o \
				$(OBJDIR)/comm.o \
				$(OBJDIR)/discretization.o\
				$(OBJDIR)/blockindex.o\
				$(OBJDIR)/fourier_$(FFTLIB).o \
				$(OBJDIR)/prefactors.o\
				$(OBJDIR)/x_derivatives.o\
				$(OBJDIR)/geometry.o

$(OBJDIR)/ff_nonlinear_term_h.o:	$(SRCDIR)/switches.h\
				$(OBJDIR)/ff_nonlinear_term.o \
				$(OBJDIR)/par_mod.o \
				$(OBJDIR)/par_other.o \
				$(OBJDIR)/par_in.o \
				$(OBJDIR)/coordinates.o \
				$(OBJDIR)/comm.o \
				$(OBJDIR)/discretization.o\
				$(OBJDIR)/blockindex.o\
				$(OBJDIR)/fourier_$(FFTLIB).o \
				$(OBJDIR)/prefactors.o\
				$(OBJDIR)/x_derivatives.o\
				$(OBJDIR)/geometry.o


$(OBJDIR)/charge_curr_dens.o:	$(SRCDIR)/switches.h \
			$(OBJDIR)/comm.o\
			$(OBJDIR)/discretization.o\
			$(OBJDIR)/gyro_average_df.o\
			$(OBJDIR)/gyro_average_ff.o\
			$(OBJDIR)/hybrid.o \
			$(OBJDIR)/par_in.o \
			$(OBJDIR)/par_other.o \
			$(OBJDIR)/vel_space.o

$(OBJDIR)/chease.o:	$(OBJDIR)/discretization.o\
			$(OBJDIR)/coordinates.o\
			$(OBJDIR)/par_geom.o\
			$(OBJDIR)/lag_interp.o

$(OBJDIR)/check_parameters.o:	$(SRCDIR)/switches.h \
			$(OBJDIR)/par_in.o\
			$(OBJDIR)/par_other.o\
			$(OBJDIR)/coordinates.o\
			$(OBJDIR)/discretization.o\
			$(OBJDIR)/equilibrium_fields.o\
			$(OBJDIR)/geometry.o\
			$(OBJDIR)/compute_dt.o\
			$(OBJDIR)/external_contr.o\
			$(OBJDIR)/eigen_parameters.o\
			$(OBJDIR)/petsc_precond.o\
			$(OBJDIR)/prefactors.o\
			$(OBJDIR)/fourier_$(FFTLIB).o\
			$(OBJDIR)/collisions.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/file_io.o\
			$(OBJDIR)/numerical_damping.o\
			$(OBJDIR)/field_solver.o\
			$(OBJDIR)/time_averages.o

$(OBJDIR)/checkpoint.o:	$(OBJDIR)/par_mod.o\
			$(OBJDIR)/file_io.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/lag_interp.o\
			$(OBJDIR)/gauss.quadrature.o\
			$(OBJDIR)/boundary.o\
			$(OBJDIR)/par_other.o

$(OBJDIR)/circular.o:	$(OBJDIR)/discretization.o\
			$(OBJDIR)/coordinates.o\
			$(OBJDIR)/par_geom.o

$(OBJDIR)/collisions.o:	$(SRCDIR)/switches.h\
			$(OBJDIR)/par_mod.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/vel_space.o\
			$(OBJDIR)/boundary.o\
			$(OBJDIR)/geometry.o\
			$(OBJDIR)/aux_func.o\
			$(OBJDIR)/profiles.o\
			$(OBJDIR)/gauss.quadrature.o\
			$(OBJDIR)/gyro_average_ff.o\
			$(OBJDIR)/gyro_average_df.o

$(OBJDIR)/comm.o:	$(OBJDIR)/par_mod.o\
			$(OBJDIR)/discretization.o\
			$(OBJDIR)/par_other.o\
			$(OBJDIR)/par_in.o

$(OBJDIR)/compute_dt.o:	$(OBJDIR)/par_mod.o\
			$(OBJDIR)/discretization.o\
			$(OBJDIR)/coordinates.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/collisions.o\
			$(OBJDIR)/init_physical.o\
			$(OBJDIR)/aux_fields.o\
			$(OBJDIR)/eigen_parameters.o\
			$(OBJDIR)/RK_coefficients.o\
			$(OBJDIR)/calc_rhs.o\
			$(OBJDIR)/external_contr.o\
			$(OBJDIR)/numerical_damping.o\
			$(OBJDIR)/x_derivatives.o

$(OBJDIR)/compute_f.o:	$(OBJDIR)/par_mod.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/vel_space.o\
			$(OBJDIR)/gyro_average_df.o\
			$(OBJDIR)/gyro_average_ff.o\
			$(OBJDIR)/antenna.o\
			$(OBJDIR)/boundary.o\
			$(OBJDIR)/island.o

$(OBJDIR)/convergence_monitoring.o:	$(OBJDIR)/discretization.o\
			$(OBJDIR)/coordinates.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/par_in.o\
			$(OBJDIR)/par_other.o\
			$(OBJDIR)/geometry.o\
			$(OBJDIR)/aux_fields.o\
			$(OBJDIR)/RK_coefficients.o\
			$(OBJDIR)/file_io.o\
			$(OBJDIR)/quasilin_model.o\
			$(OBJDIR)/vel_space.o\
			$(OBJDIR)/flr_corr.o

$(OBJDIR)/coordinates.o:	$(OBJDIR)/box_data_module.o \
			$(OBJDIR)/discretization.o\
			$(OBJDIR)/gauss.quadrature.o\
			$(OBJDIR)/par_in.o

$(OBJDIR)/coordinates_adptv.o:	$(SRCDIR)/redef.h \
			$(OBJDIR)/discretization_adptv.o\
			$(OBJDIR)/par_in.o

$(OBJDIR)/dchidxy_terms.o:	$(SRCDIR)/switches.h \
			$(OBJDIR)/par_mod.o\
			$(OBJDIR)/discretization.o\
			$(OBJDIR)/blockindex.o\
			$(OBJDIR)/coordinates.o\
			$(OBJDIR)/prefactors.o\
			$(OBJDIR)/external_contr.o\
			$(OBJDIR)/gyro_average_ff.o\
			$(OBJDIR)/gyro_average_df.o\
			$(OBJDIR)/x_derivatives.o

$(OBJDIR)/dchidz_term.o:	$(OBJDIR)/par_mod.o\
			$(OBJDIR)/discretization.o\
			$(OBJDIR)/blockindex.o\
			$(OBJDIR)/prefactors.o\
			$(OBJDIR)/gyro_average_ff.o\
			$(OBJDIR)/gyro_average_df.o\
			$(OBJDIR)/axpy.o\
			$(OBJDIR)/geometry.o\
			$(OBJDIR)/vel_space.o

$(OBJDIR)/dfdxy_terms.o:	$(SRCDIR)/switches.h \
			$(OBJDIR)/par_mod.o\
			$(OBJDIR)/discretization.o\
			$(OBJDIR)/numerical_damping.o\
			$(OBJDIR)/blockindex.o\
			$(OBJDIR)/boundary.o

$(OBJDIR)/dfdzv_terms.o:	$(OBJDIR)/par_mod.o\
			$(OBJDIR)/discretization.o\
			$(OBJDIR)/blockindex.o\
			$(OBJDIR)/prefactors.o\
			$(OBJDIR)/boundary.o\
			$(OBJDIR)/axpy.o\
			$(OBJDIR)/geometry.o\
			$(OBJDIR)/numerical_damping.o\
			$(OBJDIR)/par_geom.o\
			$(OBJDIR)/equilibrium_fields.o

$(OBJDIR)/dgdxy_terms.o:	$(SRCDIR)/switches.h \
			$(OBJDIR)/par_mod.o\
			$(OBJDIR)/discretization.o\
			$(OBJDIR)/blockindex.o\
			$(OBJDIR)/coordinates.o\
			$(OBJDIR)/prefactors.o\
			$(OBJDIR)/external_contr.o\
			$(OBJDIR)/x_derivatives.o

$(OBJDIR)/diag.o:	$(SRCDIR)/switches.h \
			$(OBJDIR)/par_mod.o\
			$(OBJDIR)/file_io.o\
			$(OBJDIR)/vel_space.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/aux_fields.o\
			$(OBJDIR)/antenna.o\
			$(OBJDIR)/flr_corr.o\
			$(OBJDIR)/aux_func.o\
			$(OBJDIR)/boundary.o\
			$(OBJDIR)/diag_aux.o\
			$(OBJDIR)/diag_df.o\
			$(OBJDIR)/diag_finit.o\
			$(OBJDIR)/diag_vsp.o\
			$(OBJDIR)/convergence_monitoring.o\
			$(OBJDIR)/spatial_averages.o\
			$(OBJDIR)/time_averages.o\
			$(OBJDIR)/checkpoint.o\
			$(OBJDIR)/eigen_parameters.o\
			$(OBJDIR)/geometry.o\
			$(OBJDIR)/axpy.o

ifneq ($(COUPLE_XGC),yes)
  $(OBJDIR)/paramters_IO.o:	$(OBJDIR)/diag_fsa_moments.o

  $(OBJDIR)/diag.o:	$(OBJDIR)/diag_energy.o\
			$(OBJDIR)/diag_extended.o\
			$(OBJDIR)/diag_Gyro_LES.o\
			$(OBJDIR)/diag_fsa_moments.o\
			$(OBJDIR)/diag_triplet.o


  $(OBJDIR)/diag_Gyro_LES.o:	$(SRCDIR)/switches.h\
			$(OBJDIR)/par_mod.o\
			$(OBJDIR)/diag_Gyro_LES_common.o\
			$(OBJDIR)/diag_Gyro_LES_fetime.o\
			$(OBJDIR)/diag_Gyro_LES_spec2d.o\
                        $(OBJDIR)/diag_Gyro_LES_real_3d.o\
			$(OBJDIR)/diag_Gyro_LES_transfer.o\
			$(OBJDIR)/diag_Gyro_LES_transfer3d.o

  $(OBJDIR)/diag_Gyro_LES_common.o:	$(SRCDIR)/switches.h\
			$(OBJDIR)/par_mod.o\
			$(OBJDIR)/antenna.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/vel_space.o\
			$(OBJDIR)/geometry.o\
			$(OBJDIR)/aux_fields.o\
			$(OBJDIR)/gyro_average_ff.o\
			$(OBJDIR)/aux_func.o

  $(OBJDIR)/diag_Gyro_LES_fetime.o:	$(SRCDIR)/switches.h\
			$(OBJDIR)/par_mod.o\
			$(OBJDIR)/file_io.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/antenna.o\
			$(OBJDIR)/collisions.o\
			$(OBJDIR)/diag_energy.o\
			$(OBJDIR)/dfdzv_terms.o\
			$(OBJDIR)/dzv_terms.o\
			$(OBJDIR)/dfdxy_terms.o\
			$(OBJDIR)/dgdxy_terms.o\
			$(OBJDIR)/dchidz_term.o\
			$(OBJDIR)/dchidxy_terms.o\
			$(OBJDIR)/spatial_averages.o\
			$(OBJDIR)/prefactors.o\
			$(OBJDIR)/numerical_damping.o\
			$(OBJDIR)/diag_Gyro_LES_common.o\
			$(OBJDIR)/gyro_average_ff.o\
			$(OBJDIR)/calc_rhs.o

  $(OBJDIR)/diag_Gyro_LES_spec2d.o:	$(SRCDIR)/switches.h\
			$(OBJDIR)/par_mod.o\
			$(OBJDIR)/file_io.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/antenna.o\
			$(OBJDIR)/collisions.o\
			$(OBJDIR)/diag_energy.o\
			$(OBJDIR)/dfdzv_terms.o\
			$(OBJDIR)/dzv_terms.o\
			$(OBJDIR)/dfdxy_terms.o\
			$(OBJDIR)/dgdxy_terms.o\
			$(OBJDIR)/dchidz_term.o\
			$(OBJDIR)/dchidxy_terms.o\
			$(OBJDIR)/spatial_averages.o\
			$(OBJDIR)/prefactors.o\
			$(OBJDIR)/numerical_damping.o\
			$(OBJDIR)/diag_Gyro_LES_common.o\
			$(OBJDIR)/calc_rhs.o

  $(OBJDIR)/diag_Gyro_LES_real_3d.o:      $(SRCDIR)/switches.h\
                        $(OBJDIR)/par_mod.o\
                        $(OBJDIR)/file_io.o\
                        $(OBJDIR)/fourier_$(FFTLIB).o\
                        $(OBJDIR)/comm.o\
                        $(OBJDIR)/antenna.o\
                        $(OBJDIR)/collisions.o\
                        $(OBJDIR)/diag_energy.o\
                        $(OBJDIR)/dfdzv_terms.o\
                        $(OBJDIR)/dzv_terms.o\
                        $(OBJDIR)/dfdxy_terms.o\
                        $(OBJDIR)/dgdxy_terms.o\
                        $(OBJDIR)/dchidz_term.o\
                        $(OBJDIR)/dchidxy_terms.o\
                        $(OBJDIR)/spatial_averages.o\
                        $(OBJDIR)/prefactors.o\
                        $(OBJDIR)/numerical_damping.o\
                        $(OBJDIR)/diag_Gyro_LES_common.o\
                        $(OBJDIR)/calc_rhs.o

  $(OBJDIR)/diag_Gyro_LES_transfer.o:	$(SRCDIR)/switches.h\
			$(OBJDIR)/par_mod.o\
			$(OBJDIR)/file_io.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/collisions.o\
			$(OBJDIR)/diag_energy.o\
			$(OBJDIR)/dfdzv_terms.o\
			$(OBJDIR)/dzv_terms.o\
			$(OBJDIR)/dfdxy_terms.o\
			$(OBJDIR)/dgdxy_terms.o\
			$(OBJDIR)/dchidz_term.o\
			$(OBJDIR)/dchidxy_terms.o\
			$(OBJDIR)/spatial_averages.o\
			$(OBJDIR)/prefactors.o\
			$(OBJDIR)/numerical_damping.o\
			$(OBJDIR)/diag_Gyro_LES_common.o\
			$(OBJDIR)/calc_rhs.o

  $(OBJDIR)/diag_Gyro_LES_transfer3d.o:	$(SRCDIR)/switches.h\
			$(OBJDIR)/par_mod.o\
			$(OBJDIR)/file_io.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/diag_energy.o\
			$(OBJDIR)/spatial_averages.o\
			$(OBJDIR)/prefactors.o\
			$(OBJDIR)/diag_Gyro_LES_common.o\
			$(OBJDIR)/calc_rhs.o

  $(OBJDIR)/diag_extended.o:	$(SRCDIR)/switches.h \
			$(OBJDIR)/par_mod.o\
			$(OBJDIR)/file_io.o\
			$(OBJDIR)/vel_space.o\
			$(OBJDIR)/calc_rhs.o\
			$(OBJDIR)/geometry.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/aux_fields.o\
			$(OBJDIR)/gyro_average_ff.o\
			$(OBJDIR)/gyro_average_df.o\
			$(OBJDIR)/collisions.o\
			$(OBJDIR)/dfdzv_terms.o\
			$(OBJDIR)/dzv_terms.o\
			$(OBJDIR)/dfdxy_terms.o\
			$(OBJDIR)/dgdxy_terms.o\
			$(OBJDIR)/dchidz_term.o\
			$(OBJDIR)/dchidxy_terms.o\
			$(OBJDIR)/diag_energy.o\
			$(OBJDIR)/sources_mod.o\
			$(OBJDIR)/blockindex.o\
			$(OBJDIR)/axpy.o\
			$(OBJDIR)/prefactors.o\
			$(OBJDIR)/numerical_damping.o\
			$(OBJDIR)/boundary.o

  $(OBJDIR)/diag_fsa_moments.o:	$(SRCDIR)/switches.h\
			$(OBJDIR)/par_mod.o\
			$(OBJDIR)/file_io.o\
			$(OBJDIR)/vel_space.o\
			$(OBJDIR)/calc_rhs.o\
			$(OBJDIR)/geometry.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/aux_fields.o\
			$(OBJDIR)/collisions.o\
			$(OBJDIR)/dfdzv_terms.o\
			$(OBJDIR)/dgdxy_terms.o\
			$(OBJDIR)/dzv_terms.o\
			$(OBJDIR)/dfdxy_terms.o\
			$(OBJDIR)/dchidz_term.o\
			$(OBJDIR)/dchidxy_terms.o\
			$(OBJDIR)/f0_term.o\
			$(OBJDIR)/sources_mod.o\
			$(OBJDIR)/spatial_averages.o\
			$(OBJDIR)/blockindex.o\
			$(OBJDIR)/axpy.o\
			$(OBJDIR)/prefactors.o\
			$(OBJDIR)/numerical_damping.o\
			$(OBJDIR)/boundary.o


endif



$(OBJDIR)/diag_aux.o:	\
			$(OBJDIR)/geometry.o\
			$(OBJDIR)/par_mod.o\
			$(OBJDIR)/vel_space.o\
			$(OBJDIR)/x_derivatives.o\
			$(OBJDIR)/aux_fields.o

$(OBJDIR)/diag_df.o:	$(OBJDIR)/par_mod.o\
			$(OBJDIR)/file_io.o\
			$(OBJDIR)/vel_space.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/boundary.o\
			$(OBJDIR)/aux_fields.o\
			$(OBJDIR)/geometry.o\
			$(OBJDIR)/sources_mod.o\
			$(OBJDIR)/f0_term.o\
			$(OBJDIR)/coordinates.o\
			$(OBJDIR)/profile_smoothing.o\
			$(OBJDIR)/lag_interp.o

$(OBJDIR)/diag_energy.o:	$(SRCDIR)/switches.h \
			$(OBJDIR)/par_mod.o\
			$(OBJDIR)/file_io.o\
			$(OBJDIR)/vel_space.o\
			$(OBJDIR)/calc_rhs.o\
			$(OBJDIR)/geometry.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/aux_fields.o\
			$(OBJDIR)/antenna.o\
			$(OBJDIR)/field_solve_ff.o\
			$(OBJDIR)/gyro_average_ff.o\
			$(OBJDIR)/gyro_average_df.o\
			$(OBJDIR)/collisions.o\
			$(OBJDIR)/dfdzv_terms.o\
			$(OBJDIR)/dzv_terms.o\
			$(OBJDIR)/dfdxy_terms.o\
			$(OBJDIR)/dgdxy_terms.o\
			$(OBJDIR)/dchidz_term.o\
			$(OBJDIR)/dchidxy_terms.o\
			$(OBJDIR)/sources_mod.o\
			$(OBJDIR)/spatial_averages.o\
			$(OBJDIR)/blockindex.o\
			$(OBJDIR)/axpy.o\
			$(OBJDIR)/prefactors.o\
			$(OBJDIR)/numerical_damping.o\
			$(OBJDIR)/boundary.o\
			$(OBJDIR)/f0_term.o\
			$(OBJDIR)/spatial_averages.o



ifneq ($(COUPLE_XGC),yes)
$(OBJDIR)/diag_triplet.o:	$(OBJDIR)/par_mod.o\
			$(OBJDIR)/file_io.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/vel_space.o\
			$(OBJDIR)/geometry.o\
			$(OBJDIR)/aux_fields.o\
			$(OBJDIR)/gyro_average_ff.o\
			$(OBJDIR)/diag_extended.o
endif

$(OBJDIR)/diag_finit.o:	$(OBJDIR)/par_mod.o\
			$(OBJDIR)/file_io.o\
			$(OBJDIR)/vel_space.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/aux_fields.o

$(OBJDIR)/diag_fmsvd.o:	$(OBJDIR)/par_mod.o\
			$(OBJDIR)/file_io.o\
			$(OBJDIR)/geometry.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/init_physical.o\
			$(OBJDIR)/arrays.o

ifneq ($(COUPLE_XGC),yes)
$(OBJDIR)/diag_neoclass.o:	$(OBJDIR)/par_mod.o\
			$(OBJDIR)/aux_fields.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/geometry.o\
			$(OBJDIR)/file_io.o\
			$(OBJDIR)/spatial_averages.o\
			$(OBJDIR)/vel_space.o
endif
$(OBJDIR)/diag_svd.o:	$(OBJDIR)/par_mod.o\
			$(OBJDIR)/init_physical.o\
			$(OBJDIR)/file_io.o\
			$(OBJDIR)/geometry.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/vel_space.o\
			$(OBJDIR)/diag.o\
			$(OBJDIR)/arrays.o

$(OBJDIR)/diag_vsp.o:	$(OBJDIR)/aux_fields.o\
			$(OBJDIR)/axpy.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/coordinates.o\
			$(OBJDIR)/diag_aux.o\
			$(OBJDIR)/file_io.o\
			$(OBJDIR)/geometry.o\
			$(OBJDIR)/par_mod.o\
			$(OBJDIR)/profile_IO.o\
			$(OBJDIR)/vel_space.o

$(OBJDIR)/discretization.o:	$(SRCDIR)/switches.h \
				$(OBJDIR)/par_in.o\
				$(OBJDIR)/par_other.o

$(OBJDIR)/discretization_adptv.o:	$(OBJDIR)/discretization.o\
					$(OBJDIR)/coordinates.o\
					$(OBJDIR)/geometry.o

$(OBJDIR)/dzv_terms.o:	$(SRCDIR)/switches.h \
			$(OBJDIR)/par_mod.o\
			$(OBJDIR)/discretization.o\
			$(OBJDIR)/discretization_adptv.o\
			$(OBJDIR)/blockindex.o\
			$(OBJDIR)/prefactors.o\
			$(OBJDIR)/boundary.o\
			$(OBJDIR)/geometry.o\
			$(OBJDIR)/vel_space.o\
			$(OBJDIR)/boundary_exchange_z.o\
			$(OBJDIR)/axpy.o\
			$(OBJDIR)/lag_interp.o\
			$(OBJDIR)/gyro_average_ff.o\
			$(OBJDIR)/fourier_$(FFTLIB).o\
			$(OBJDIR)/equilibrium_fields.o\
			$(OBJDIR)/numerical_damping.o

$(OBJDIR)/eigen_parameters.o:	$(OBJDIR)/discretization.o\
			$(OBJDIR)/par_in.o\
			$(OBJDIR)/coordinates.o

$(OBJDIR)/eigenvalue_comp.o:	$(OBJDIR)/parameters_IO.o\
				$(OBJDIR)/eigen_parameters.o\
				$(OBJDIR)/aux_fields.o\
				$(OBJDIR)/init_physical.o

$(OBJDIR)/equilibrium_fields.o:	$(OBJDIR)/par_in.o\
			$(OBJDIR)/discretization.o\
			$(OBJDIR)/external_contr.o\
			$(OBJDIR)/geometry.o\
			$(OBJDIR)/comm.o

$(OBJDIR)/exchange_z_ff.o:	$(OBJDIR)/par_in.o\
				$(OBJDIR)/par_other.o\
				$(OBJDIR)/discretization.o\
				$(OBJDIR)/comm.o\
				$(OBJDIR)/fourier_$(FFTLIB).o

$(OBJDIR)/external_contr.o:	$(OBJDIR)/par_mod.o\
			$(OBJDIR)/coordinates.o\
			$(OBJDIR)/aux_func.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/profiles.o\
			$(OBJDIR)/fourier_$(FFTLIB).o\
			$(OBJDIR)/geometry.o\
			$(OBJDIR)/file_io.o\
			$(OBJDIR)/lag_interp.o

$(OBJDIR)/f0_term.o:	$(OBJDIR)/par_mod.o\
			$(OBJDIR)/discretization.o\
			$(OBJDIR)/prefactors.o\
			$(OBJDIR)/blockindex.o\
			$(OBJDIR)/geometry.o\
			$(OBJDIR)/sources_mod.o\
			$(OBJDIR)/vel_space.o

$(OBJDIR)/field_solve_df.o:	$(OBJDIR)/spectype.o\
			$(OBJDIR)/discretization.o\
			$(OBJDIR)/coordinates.o\
			$(OBJDIR)/par_in.o\
			$(OBJDIR)/par_other.o\
			$(OBJDIR)/geometry.o\
			$(OBJDIR)/vel_space.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/gyro_average_df.o\
			$(OBJDIR)/debug_output.o\
			$(OBJDIR)/grid1d.o\
			$(OBJDIR)/flr_corr.o\
			$(OBJDIR)/adiabatic_response.o\
			$(OBJDIR)/spatial_averages.o\
			$(OBJDIR)/hybrid.o\
			$(OBJDIR)/island.o


$(OBJDIR)/field_solve_ff.o:	$(OBJDIR)/par_mod.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/coordinates.o\
			$(OBJDIR)/aux_func.o\
			$(OBJDIR)/gyro_average_ff.o\
			$(OBJDIR)/vel_space.o\
			$(OBJDIR)/geometry.o\
			$(OBJDIR)/charge_curr_dens.o\
			$(OBJDIR)/adiabatic_response.o\
			$(OBJDIR)/hybrid.o\
			$(OBJDIR)/antenna.o\
			$(OBJDIR)/equilibrium_fields.o\
			$(OBJDIR)/spatial_averages.o\
			$(OBJDIR)/island.o

$(OBJDIR)/field_solver.o:	$(SRCDIR)/switches.h \
			$(OBJDIR)/discretization.o\
			$(OBJDIR)/par_in.o\
			$(OBJDIR)/par_other.o\
			$(OBJDIR)/charge_curr_dens.o\
			$(OBJDIR)/adiabatic_response.o\
			$(OBJDIR)/antenna.o\
			$(OBJDIR)/gyro_average_df.o\
			$(OBJDIR)/gyro_average_ff.o\
			$(OBJDIR)/field_solve_ff.o\
			$(OBJDIR)/field_solve_df.o\
			$(OBJDIR)/island.o\
			$(OBJDIR)/poloidal_planes.o


$(OBJDIR)/file_io.o:	$(OBJDIR)/discretization.o\
			$(OBJDIR)/comm.o

$(OBJDIR)/flr_corr.o:	$(OBJDIR)/par_mod.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/gyro_average_df.o\
			$(OBJDIR)/geometry.o\
			$(OBJDIR)/vel_space.o\
			$(OBJDIR)/aux_func.o\
			$(OBJDIR)/hybrid.o\
			$(OBJDIR)/equilibrium_fields.o


$(OBJDIR)/fourier_$(FFTLIB).o:	$(OBJDIR)/discretization.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/par_in.o\
			$(OBJDIR)/par_poloidal_planes.o

$(OBJDIR)/fourier_mkl.o:	$(OBJDIR)/mkl_dfti.o

$(OBJDIR)/fullmatrix_aux.o:	$(OBJDIR)/par_mod.o\
			$(OBJDIR)/eigen_parameters.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/calc_rhs.o\
			$(OBJDIR)/file_io.o\
			$(OBJDIR)/init_physical.o\
			$(OBJDIR)/aux_fields.o\
			$(OBJDIR)/par_other.o

$(OBJDIR)/full_rhs_mat.o:	$(OBJDIR)/par_mod.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/fullmatrix_aux.o\
			$(OBJDIR)/calc_rhs.o\
			$(OBJDIR)/RK_coefficients.o

$(OBJDIR)/gauss.quadrature.o:	$(OBJDIR)/discretization.o

$(OBJDIR)/gene.o:	$(SRCDIR)/switches.h \
			$(SRCDIR)/intrinsic_sizes.h \
			$(OBJDIR)/gene_scan.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/par_mod.o\
			$(OBJDIR)/parameters_IO.o\
			$(OBJDIR)/check_parameters.o\
			$(OBJDIR)/gene_subroutine.o\
			$(OBJDIR)/file_io.o\
			$(OBJDIR)/diag.o\
			$(OBJDIR)/pinning.o

ifeq ($(USE_MEMORY_CHECKER),own)
$(OBJDIR)/gene.o:	$(OBJDIR)/memory_check_module.o
endif

$(OBJDIR)/gene_scan.o:	$(SRCDIR)/switches.h \
			$(SRCDIR)/intrinsic_sizes.h \
			$(OBJDIR)/comm.o\
			$(OBJDIR)/parameters_IO.o\
			$(OBJDIR)/file_io.o\
			$(OBJDIR)/gene_subroutine.o\
			$(OBJDIR)/diag.o\
			$(OBJDIR)/checkpoint.o\
			$(OBJDIR)/discretization.o\
			$(OBJDIR)/par_in.o\
			$(OBJDIR)/par_other.o

$(OBJDIR)/gene_subroutine.o:	$(OBJDIR)/par_mod.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/perf_opt.o\
			$(OBJDIR)/eigenvalue_comp.o\
			$(OBJDIR)/initial_value_comp.o\
			$(OBJDIR)/compute_dt.o\
			$(OBJDIR)/codemods.o\
			$(OBJDIR)/parameters_IO.o\
			$(OBJDIR)/fullmatrix_aux.o\
			$(OBJDIR)/neo_equilibrium.o\
			$(OBJDIR)/diag.o\
			$(OBJDIR)/time_averages.o\
			$(OBJDIR)/diag_df.o\
			$(OBJDIR)/checkpoint.o\
			$(OBJDIR)/geometry.o\
			$(OBJDIR)/profiles.o\
			$(OBJDIR)/collisions.o\
			$(OBJDIR)/diag_fmsvd.o

$(OBJDIR)/geometry.o:	$(OBJDIR)/par_mod.o\
			$(OBJDIR)/par_geom.o\
			$(OBJDIR)/file_io.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/chease.o\
			$(OBJDIR)/circular.o\
			$(OBJDIR)/tracer.o\
			$(OBJDIR)/lag_interp.o\
			$(OBJDIR)/profile_IO.o\
			$(OBJDIR)/miller_geometry.o

$(OBJDIR)/grid1d.o:


$(OBJDIR)/gyro_average_df.o:    $(OBJDIR)/localpolynombase.o\
			$(OBJDIR)/comm.o\
                        $(OBJDIR)/discretization.o\
                        $(OBJDIR)/coordinates.o\
                        $(OBJDIR)/par_other.o\
                        $(OBJDIR)/par_in.o\
                        $(OBJDIR)/geometry.o\
                        $(OBJDIR)/BoundaryDescription.o\
			$(OBJDIR)/grid1d.o\
			$(OBJDIR)/lag_interp.o\
                        $(OBJDIR)/file_io.o

$(OBJDIR)/gyro_average_ff.o:	$(OBJDIR)/discretization.o\
			$(OBJDIR)/coordinates.o\
			$(OBJDIR)/par_in.o\
			$(OBJDIR)/par_other.o\
			$(OBJDIR)/aux_func.o\
			$(OBJDIR)/boundary.o\
			$(OBJDIR)/geometry.o

$(OBJDIR)/Gyro_LES.o:   $(SRCDIR)/switches.h \
			$(OBJDIR)/par_mod.o \
			$(OBJDIR)/comm.o\
			$(OBJDIR)/aux_fields.o\
			$(OBJDIR)/vel_space.o\
			$(OBJDIR)/geometry.o\
			$(OBJDIR)/discretization.o\
			$(OBJDIR)/file_io.o\
			$(OBJDIR)/prefactors.o\
			$(OBJDIR)/blockindex.o\
			$(OBJDIR)/numerical_damping.o\
			$(OBJDIR)/axpy.o\
			$(OBJDIR)/calc_rhs.o

$(OBJDIR)/hybrid.o:	$(OBJDIR)/comm.o\
			$(OBJDIR)/geometry.o\
			$(OBJDIR)/discretization.o\
			$(OBJDIR)/coordinates.o

$(OBJDIR)/init_cond.o:	$(OBJDIR)/par_mod.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/aux_func.o\
			$(OBJDIR)/geometry.o\
			$(OBJDIR)/gyro_average_ff.o\
			$(OBJDIR)/aux_fields.o\
			$(OBJDIR)/fourier_$(FFTLIB).o\
			$(OBJDIR)/parallel_nonlin.o\
			$(OBJDIR)/discretization.o\
			$(OBJDIR)/vel_space.o\
			$(OBJDIR)/neo_equilibrium.o\
			$(OBJDIR)/antenna.o\
			$(OBJDIR)/mtrandom.o\
			$(OBJDIR)/island.o


$(OBJDIR)/initial_value_comp.o:	$(OBJDIR)/par_mod.o\
			$(OBJDIR)/parameters_IO.o\
			$(OBJDIR)/diag.o\
			$(OBJDIR)/diag_df.o\
			$(OBJDIR)/diag_svd.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/time_scheme.o\
			$(OBJDIR)/aux_fields.o\
			$(OBJDIR)/antenna.o\
			$(OBJDIR)/checkpoint.o\
			$(OBJDIR)/reset_mod.o\
			$(OBJDIR)/init_cond.o\
			$(OBJDIR)/init_physical.o\
			$(OBJDIR)/external_contr.o\
			$(OBJDIR)/file_io.o\
			$(OBJDIR)/sources_mod.o\
			$(OBJDIR)/convergence_monitoring.o\
			$(OBJDIR)/time_averages.o \
			$(OBJDIR)/itt_api.o


$(OBJDIR)/init_physical.o: \
			$(OBJDIR)/arrays.o\
			$(OBJDIR)/blockindex.o\
			$(OBJDIR)/boundary.o\
			$(OBJDIR)/box_data_module.o\
			$(OBJDIR)/checkpoint.o\
			$(OBJDIR)/collisions.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/coordinates.o\
			$(OBJDIR)/discretization.o\
			$(OBJDIR)/dzv_terms.o\
			$(OBJDIR)/equilibrium_fields.o\
			$(OBJDIR)/external_contr.o\
			$(OBJDIR)/geometry.o\
			$(OBJDIR)/grid1d.o\
			$(OBJDIR)/lag_interp.o\
			$(OBJDIR)/numerical_damping.o\
			$(OBJDIR)/par_geom.o\
			$(OBJDIR)/par_in.o\
			$(OBJDIR)/par_mod.o\
			$(OBJDIR)/par_other.o\
			$(OBJDIR)/prefactors.o\
			$(OBJDIR)/profiles.o\
			$(OBJDIR)/vel_space.o\
			$(OBJDIR)/x_derivatives.o

$(OBJDIR)/lag_interp.o:

$(OBJDIR)/localpolynombase.o:	$(OBJDIR)/grid1d.o

$(OBJDIR)/localpolynombase2d.o:   $(OBJDIR)/grid1d.o\

$(OBJDIR)/circular.o:	$(OBJDIR)/discretization.o\
			$(OBJDIR)/coordinates.o\
			$(OBJDIR)/par_geom.o\
			$(OBJDIR)/fourier_$(FFTLIB).o

$(OBJDIR)/miller_geometry.o:	$(OBJDIR)/discretization.o\
			$(OBJDIR)/coordinates.o\
			$(OBJDIR)/par_geom.o\
			$(OBJDIR)/fourier_$(FFTLIB).o\
			$(OBJDIR)/par_in.o\
			$(OBJDIR)/par_other.o\
			$(OBJDIR)/lag_interp.o

$(OBJDIR)/mtrandom.o:

$(OBJDIR)/neo_equilibrium.o:	$(OBJDIR)/par_mod.o\
			$(OBJDIR)/discretization.o\
			$(OBJDIR)/f0_term.o\
			$(OBJDIR)/init_physical.o\
			$(OBJDIR)/aux_fields.o\
			$(OBJDIR)/diag.o\
			$(OBJDIR)/calc_rhs.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/checkpoint.o

ifeq ($(USE_C_NONLIN),yes)
$(OBJDIR)/calc_rhs.o:		$(OBJDIR)/df_nonlinear_c.o
$(OBJDIR)/df_nonlinear_c.o:	$(OBJDIR)/df_nonlinear_term.o \
				$(OBJDIR)/c_nonlinearity.o \
				$(OBJDIR)/c_fourier_fftw.o
endif

ifeq ($(USE_CUDA_NONLIN),yes)
$(OBJDIR)/calc_rhs.o:		$(OBJDIR)/df_nonlinear_cuda.o
$(OBJDIR)/df_nonlinear_cuda.o:	$(OBJDIR)/df_nonlinear_term.o \
				$(OBJDIR)/cuda_nonlinearity.o
$(OBJDIR)/cuda_nonlinearity.o:	$(SRCDIR)/cuda_overlap.h \
				$(SRCDIR)/cuda_kernels.cu \
				$(SRCDIR)/reduction_kernel.cu \
				$(OBJDIR)/cuda_overlap.o \
				$(OBJDIR)/cuda_fourier_cufft.o
$(OBJDIR)/cuda_fourier_cufft.o:	$(SRCDIR)/cuda_overlap.h
$(OBJDIR)/cuda_overlap.o:	$(SRCDIR)/cuda_overlap.h
$(OBJDIR)/time_scheme.o:	$(OBJDIR)/page_locked_memory.o
$(OBJDIR)/page_locked_memory.o: $(OBJDIR)/cuda_nonlinearity.o
$(OBJDIR)/initial_value_comp.o: $(OBJDIR)/cuda_general.o
endif

ifeq ($(USE_SPLIT_NONLIN),yes)
$(OBJDIR)/calc_rhs.o:	$(OBJDIR)/df_split_nonlinear_term.o
$(OBJDIR)/df_split_nonlinear_term.o:	$(OBJDIR)/df_nonlinear_term.o \
					$(OBJDIR)/df_nonlinear_cuda.o
endif

ifeq ($(USE_MIC_NONLIN),yes)
$(OBJDIR)/calc_rhs.o:	$(OBJDIR)/df_mic_nonlinear_term.o
$(OBJDIR)/df_mic_nonlinear_term.o: $(OBJDIR)/fourier_mkl_mic.o
endif

$(OBJDIR)/numerical_damping.o:	$(OBJDIR)/discretization.o\
			$(OBJDIR)/coordinates.o\
			$(OBJDIR)/geometry.o\
			$(OBJDIR)/par_other.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/par_in.o

$(OBJDIR)/parallel_nonlin.o:	$(SRCDIR)/switches.h \
			$(OBJDIR)/par_mod.o\
			$(OBJDIR)/par_other.o\
			$(OBJDIR)/par_in.o\
			$(OBJDIR)/coordinates.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/discretization.o\
			$(OBJDIR)/blockindex.o\
			$(OBJDIR)/fourier_$(FFTLIB).o\
			$(OBJDIR)/prefactors.o\
			$(OBJDIR)/gyro_average_ff.o\
			$(OBJDIR)/geometry.o\
			$(OBJDIR)/axpy.o\
			$(OBJDIR)/dchidz_term.o\
			$(OBJDIR)/all_rhs_terms.o\
			$(OBJDIR)/ff_nonlinear_term.o\
			$(OBJDIR)/df_nonlinear_term.o


$(OBJDIR)/parameters_IO.o:	$(SRCDIR)/switches.h\
			$(OBJDIR)/par_in.o\
			$(OBJDIR)/par_other.o\
			$(OBJDIR)/par_geom.o\
			$(OBJDIR)/adiabatic_response.o\
			$(OBJDIR)/file_io.o\
			$(OBJDIR)/discretization.o\
			$(OBJDIR)/discretization_adptv.o\
			$(OBJDIR)/coordinates.o\
			$(OBJDIR)/check_parameters.o\
			$(OBJDIR)/chease.o\
			$(OBJDIR)/collisions.o\
			$(OBJDIR)/geometry.o\
			$(OBJDIR)/perf_opt.o\
			$(OBJDIR)/init_cond.o\
			$(OBJDIR)/equilibrium_fields.o\
			$(OBJDIR)/external_contr.o\
			$(OBJDIR)/antenna.o\
			$(OBJDIR)/diag.o\
			$(OBJDIR)/diag_df.o\
			$(OBJDIR)/diag_vsp.o\
			$(OBJDIR)/checkpoint.o\
			$(OBJDIR)/reset_mod.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/sources_mod.o\
			$(OBJDIR)/gauss.quadrature.o\
			$(OBJDIR)/boundary_exchange_z.o\
			$(OBJDIR)/fullmatrix_aux.o\
			$(OBJDIR)/eigen_parameters.o\
			$(OBJDIR)/petsc_precond.o\
			$(OBJDIR)/gyro_average_df.o\
			$(OBJDIR)/convergence_monitoring.o\
			$(OBJDIR)/neo_equilibrium.o\
			$(OBJDIR)/time_averages.o\
			$(OBJDIR)/profiles.o\
			$(OBJDIR)/dzv_terms.o\
			$(OBJDIR)/numerical_damping.o\
			$(OBJDIR)/quasilin_model.o\
			$(OBJDIR)/miller_geometry.o\
			$(OBJDIR)/prefactors.o\
			$(OBJDIR)/profile_smoothing.o\
			$(OBJDIR)/par_poloidal_planes.o\
                        $(OBJDIR)/poloidal_planes.o

$(OBJDIR)/par_geom.o:

$(OBJDIR)/par_in.o:	$(SRCDIR)/switches.h \
			$(SRCDIR)/intrinsic_sizes.h \
			$(OBJDIR)/spectype.o

$(OBJDIR)/par_mod.o:	$(OBJDIR)/discretization.o\
			$(OBJDIR)/coordinates.o\
			$(OBJDIR)/par_in.o\
			$(OBJDIR)/par_other.o

$(OBJDIR)/par_other.o:

$(OBJDIR)/perf_opt.o:	$(OBJDIR)/par_mod.o\
			$(OBJDIR)/boundary.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/aux_fields.o\
			$(OBJDIR)/calc_rhs.o\
			$(OBJDIR)/init_physical.o\
			$(OBJDIR)/vel_space.o\
			$(OBJDIR)/time_scheme.o\
			$(OBJDIR)/compute_dt.o\
			$(OBJDIR)/check_parameters.o\
			$(OBJDIR)/file_io.o\
			$(OBJDIR)/arrays.o\
			$(OBJDIR)/field_solve_df.o\
			$(OBJDIR)/external_contr.o\
			$(OBJDIR)/diag.o\
			$(OBJDIR)/convergence_monitoring.o\
			$(OBJDIR)/petsc_precond.o\
			$(OBJDIR)/mtrandom.o\
			$(OBJDIR)/axpy.o\
			$(OBJDIR)/eigen_parameters.o\
			$(OBJDIR)/hybrid.o\
			$(OBJDIR)/geometry.o\
			$(OBJDIR)/aux_func.o\
			$(OBJDIR)/blockindex.o

$(OBJDIR)/petsc_precond.o:	$(SRCDIR)/switches.h \
			$(OBJDIR)/par_other.o\
			$(OBJDIR)/par_in.o\
			$(OBJDIR)/eigen_parameters.o\
			$(OBJDIR)/discretization.o\
			$(OBJDIR)/calc_rhs.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/geometry.o

ifeq ($(WITH_CUTILS),no)
$(OBJDIR)/gene.o:	$(OBJDIR)/dummy_c_utils.o
$(OBJDIR)/perf_opt.o:	$(OBJDIR)/dummy_c_utils.o
$(OBJDIR)/slepc_aux.o:	$(OBJDIR)/dummy_c_utils.o
else
$(OBJDIR)/gene.o:	$(OBJDIR)/c_utils.o
$(OBJDIR)/perf_opt.o:	$(OBJDIR)/c_utils.o
$(OBJDIR)/slepc_aux.o:	$(OBJDIR)/c_utils.o
endif

$(OBJDIR)/poly_roots.o:

$(OBJDIR)/prefactors.o:	$(OBJDIR)/par_mod.o\
			$(OBJDIR)/discretization.o\
			$(OBJDIR)/coordinates.o\
			$(OBJDIR)/external_contr.o\
			$(OBJDIR)/geometry.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/equilibrium_fields.o\
			$(OBJDIR)/vel_space.o

$(OBJDIR)/profile_IO.o:	$(OBJDIR)/par_mod.o\
			$(OBJDIR)/lag_interp.o\
			$(OBJDIR)/par_geom.o\
			$(OBJDIR)/file_io.o

$(OBJDIR)/profile_smoothing.o:	$(OBJDIR)/par_mod.o\
			$(OBJDIR)/lag_interp.o\
			$(OBJDIR)/geometry.o

$(OBJDIR)/profiles.o:	$(OBJDIR)/par_mod.o\
			$(OBJDIR)/geometry.o\
			$(OBJDIR)/coordinates.o\
			$(OBJDIR)/lag_interp.o\
			$(OBJDIR)/profile_IO.o\
			$(OBJDIR)/file_io.o\
			$(OBJDIR)/x_derivatives.o


$(OBJDIR)/quasilin_model.o:     $(OBJDIR)/discretization.o\
                        $(OBJDIR)/coordinates.o\
                        $(OBJDIR)/comm.o\
                        $(OBJDIR)/par_in.o\
                        $(OBJDIR)/par_other.o\
                        $(OBJDIR)/geometry.o\
			$(OBJDIR)/aux_fields.o\
			$(OBJDIR)/time_averages.o


$(OBJDIR)/reset_mod.o:	$(OBJDIR)/par_mod.o\
			$(OBJDIR)/spectype.o\
			$(OBJDIR)/vel_space.o\
			$(OBJDIR)/prefactors.o\
			$(OBJDIR)/geometry.o\
			$(OBJDIR)/field_solve_df.o\
			$(OBJDIR)/calc_rhs.o\
			$(OBJDIR)/diag.o\
			$(OBJDIR)/diag_df.o\
			$(OBJDIR)/coordinates.o\
			$(OBJDIR)/profiles.o\
			$(OBJDIR)/f0_term.o\
			$(OBJDIR)/lag_interp.o\
			$(OBJDIR)/aux_fields.o\
			$(OBJDIR)/profile_smoothing.o

$(OBJDIR)/rhs_computations.o:	$(OBJDIR)/par_mod.o\
			$(OBJDIR)/discretization.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/calc_rhs.o\
			$(OBJDIR)/full_rhs_mat.o

$(OBJDIR)/RK_coefficients.o:	$(OBJDIR)/par_mod.o\
			$(OBJDIR)/discretization.o\
			$(OBJDIR)/poly_roots.o

$(OBJDIR)/sources_mod.o:	$(OBJDIR)/comm.o\
			$(OBJDIR)/discretization.o\
			$(OBJDIR)/coordinates.o\
			$(OBJDIR)/prefactors.o\
			$(OBJDIR)/geometry.o\
			$(OBJDIR)/par_mod.o\
			$(OBJDIR)/par_in.o\
			$(OBJDIR)/profiles.o\
			$(OBJDIR)/vel_space.o\
			$(OBJDIR)/par_other.o\
			$(OBJDIR)/file_io.o\
			$(OBJDIR)/blockindex.o

$(OBJDIR)/spatial_averages.o:	$(OBJDIR)/discretization.o\
			$(OBJDIR)/geometry.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/par_other.o

$(OBJDIR)/spectype.o:

ifneq ($(COUPLE_XGC),yes)
$(OBJDIR)/stellopt_interface.o:	$(SRCDIR)/switches.h\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/parameters_IO.o
endif

$(OBJDIR)/island.o:	$(OBJDIR)/par_mod.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/coordinates.o\
			$(OBJDIR)/blockindex.o\
			$(OBJDIR)/geometry.o\
			$(OBJDIR)/file_io.o\
			$(OBJDIR)/vel_space.o

$(OBJDIR)/time_averages.o:	$(OBJDIR)/discretization.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/par_mod.o\
			$(OBJDIR)/spatial_averages.o\
			$(OBJDIR)/file_io.o

$(OBJDIR)/time_scheme.o:	$(OBJDIR)/par_mod.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/rhs_computations.o\
			$(OBJDIR)/aux_fields.o\
			$(OBJDIR)/sources_mod.o\
			$(OBJDIR)/external_contr.o\
			$(OBJDIR)/RK_coefficients.o\
			$(OBJDIR)/compute_dt.o\
			$(OBJDIR)/calc_rhs.o\
			$(OBJDIR)/mtrandom.o\
			$(OBJDIR)/antenna.o

$(OBJDIR)/tracer_aux.o:	$(OBJDIR)/discretization.o\
			$(OBJDIR)/tracer_IO.o\
			$(OBJDIR)/lag_interp.o\
			$(OBJDIR)/par_mod.o\
			$(OBJDIR)/file_io.o

$(OBJDIR)/tracer.o:	$(OBJDIR)/discretization.o\
			$(OBJDIR)/lag_interp.o\
			$(OBJDIR)/tracer_aux.o\
			$(OBJDIR)/tracer_rk5.o\
			$(OBJDIR)/tracer_IO.o\
			$(OBJDIR)/par_other.o\
			$(OBJDIR)/par_in.o\
			$(OBJDIR)/par_geom.o\
			$(OBJDIR)/comm.o


$(OBJDIR)/tracer_IO.o:	$(OBJDIR)/discretization.o\
			$(OBJDIR)/lag_interp.o\
			$(OBJDIR)/spline_interp.o\
			$(OBJDIR)/par_mod.o\
			$(OBJDIR)/mtrandom.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/par_geom.o\
			$(OBJDIR)/file_io.o

$(OBJDIR)/tracer_rk5.o:	$(OBJDIR)/tracer_rk5_util.o\
			$(OBJDIR)/tracer_IO.o

$(OBJDIR)/tracer_rk5_util.o:

$(OBJDIR)/vel_space.o:	$(OBJDIR)/par_in.o\
			$(OBJDIR)/par_other.o\
			$(OBJDIR)/discretization.o\
			$(OBJDIR)/boundary_exchange_z.o\
			$(OBJDIR)/coordinates.o\
			$(OBJDIR)/gyro_average_df.o\
			$(OBJDIR)/gyro_average_ff.o\
			$(OBJDIR)/lag_interp.o\
			$(OBJDIR)/geometry.o\
			$(OBJDIR)/axpy.o\
			$(OBJDIR)/hybrid.o\
			$(OBJDIR)/equilibrium_fields.o

$(OBJDIR)/x_derivatives.o:	$(OBJDIR)/boundary.o\
			$(OBJDIR)/discretization.o\
			$(OBJDIR)/par_other.o

$(OBJDIR)/par_poloidal_planes.o:	$(OBJDIR)/discretization.o

$(OBJDIR)/poloidal_planes.o:	$(OBJDIR)/par_poloidal_planes.o\
                                $(OBJDIR)/create_nonuni_grid.o\
				$(OBJDIR)/poloidal_planes_nonuni.o

$(OBJDIR)/poloidal_planes_aux.o:$(OBJDIR)/par_poloidal_planes.o\
                                $(OBJDIR)/discretization.o\
                                $(OBJDIR)/geometry.o\
                                $(OBJDIR)/flr_corr.o\
                                $(OBJDIR)/par_mod.o\
				$(OBJDIR)/fourier_$(FFTLIB).o

$(OBJDIR)/create_nonuni_grid.o:	$(OBJDIR)/par_other.o\
				$(OBJDIR)/par_geom.o\
				$(OBJDIR)/comm.o\
				$(OBJDIR)/coordinates.o\
				$(OBJDIR)/geometry.o\
				$(OBJDIR)/par_poloidal_planes.o

$(OBJDIR)/import_XGC_grid.o:	$(OBJDIR)/discretization.o\
				$(OBJDIR)/geometry.o\
				$(OBJDIR)/par_in.o\
				$(OBJDIR)/file_io.o\
                                $(OBJDIR)/lag_interp.o\
				$(OBJDIR)/par_poloidal_planes.o\
				$(OBJDIR)/boundary_exchange_XGC.o

$(OBJDIR)/poloidal_planes_nonuni.o:	$(OBJDIR)/par_poloidal_planes.o\
				$(OBJDIR)/create_nonuni_grid.o\
				$(OBJDIR)/poloidal_planes_aux.o\
				$(OBJDIR)/coordinates.o\
                                $(OBJDIR)/par_geom.o\
                                $(OBJDIR)/geometry.o\
                                $(OBJDIR)/lag_interp.o\
                                $(OBJDIR)/comm.o\
                                $(OBJDIR)/discretization.o\
                                $(OBJDIR)/fourier_$(FFTLIB).o\
				$(OBJDIR)/import_XGC_grid.o

$(OBJDIR)/boundary_exchange_XGC.o:	$(OBJDIR)/BoundaryDescription.o\
		                        $(OBJDIR)/boundary_exchange_general.o\
					$(OBJDIR)/discretization.o\
					$(OBJDIR)/comm.o\
					$(OBJDIR)/par_poloidal_planes.o



################################
######### Interfaces ###########
################################
$(OBJDIR)/$(LIBN).o:	$(OBJDIR)/gene_subroutine.o

$(OBJDIR)/test_$(LIBN).o:	$(SRCDIR)/test_$(LIBN).F90\
				$(OBJDIR)/$(LIBN).o

ifeq ($(LIBN),libITMgene)

$(OBJDIR)/$(LIBN)_aux.o:	$(OBJDIR)/parameters_IO.o\
				$(OBJDIR)/par_geom.o

$(OBJDIR)/$(LIBN).o:		$(OBJDIR)/libITMgene_aux.o

INCPATHS += $(ITM_INCLUDE)
LIBS += $(ITM_LIBS)

LIBELEMENTS += $(OBJDIR)/$(LIBN)_aux.o
endif

################################
######### SLEPC part ###########
################################
ifeq ($(SLEPC),yes)

$(OBJDIR)/eigen_iterative.o:	$(OBJDIR)/par_mod.o\
			$(OBJDIR)/eigen_parameters.o\
			$(OBJDIR)/par_in.o\
			$(OBJDIR)/discretization.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/rhs_computations.o\
			$(OBJDIR)/diag.o\
			$(OBJDIR)/diag_energy.o\
			$(OBJDIR)/petsc_aux.o\
			$(OBJDIR)/slepc_aux.o\
			$(OBJDIR)/file_io.o\
			$(OBJDIR)/checkpoint.o\
			$(OBJDIR)/init_cond.o\
			$(OBJDIR)/quasilin_model.o

$(OBJDIR)/eigenvalue_comp.o:	$(OBJDIR)/eigen_iterative.o

$(OBJDIR)/gene_scan.o:		$(OBJDIR)/slepc_aux.o\
				$(OBJDIR)/eigen_iterative.o

$(OBJDIR)/impl_petsc.o:	$(OBJDIR)/par_mod.o\
				$(OBJDIR)/petsc_aux.o\
				$(OBJDIR)/calc_rhs.o\
				$(OBJDIR)/convergence_monitoring.o\
				$(OBJDIR)/petsc_precond.o

$(OBJDIR)/neo_equilibrium.o:	$(OBJDIR)/petsc_aux.o\
				$(OBJDIR)/petsc_precond.o\
				$(OBJDIR)/slepc_aux.o

$(OBJDIR)/petsc_aux.o:	$(OBJDIR)/par_other.o\
			$(OBJDIR)/par_in.o\
			$(OBJDIR)/eigen_parameters.o\
			$(OBJDIR)/discretization.o\
			$(OBJDIR)/calc_rhs.o\
			$(OBJDIR)/comm.o

$(OBJDIR)/petsc_precond.o:	$(OBJDIR)/petsc_aux.o

$(OBJDIR)/rhs_computations.o:	$(OBJDIR)/impl_petsc.o

$(OBJDIR)/slepc_aux.o:		$(OBJDIR)/par_in.o\
				$(OBJDIR)/par_mod.o\
				$(OBJDIR)/par_other.o\
				$(OBJDIR)/diag.o\
				$(OBJDIR)/discretization.o\
				$(OBJDIR)/eigen_parameters.o\
				$(OBJDIR)/petsc_aux.o\
				$(OBJDIR)/petsc_precond.o\
				$(OBJDIR)/comm.o

$(OBJDIR)/compute_dt.o:	$(OBJDIR)/petsc_aux.o\
				$(OBJDIR)/slepc_aux.o

endif

################################
####### SCALAPACK part #########
################################
ifeq ($(SCALAPACK),yes)

$(OBJDIR)/adiabatic_response.o:	$(OBJDIR)/matrix.o

$(OBJDIR)/BandedMatrix.o:	$(SRCDIR)/BandedMatrix.h \
			$(OBJDIR)/processgrid.o \
			$(OBJDIR)/matrix.o \
			$(OBJDIR)/storebandedmatrix.o \
			$(OBJDIR)/Vector.o

$(OBJDIR)/comm.o:	$(OBJDIR)/matrix.o\
			$(OBJDIR)/Vector.o

$(OBJDIR)/derivative_matrix.o:	$(SRCDIR)/switches.h \
			$(OBJDIR)/matrix.o $(OBJDIR)/grid1d.o \
			$(OBJDIR)/BandedMatrix.o

$(OBJDIR)/eigen_direct.o:	$(SRCDIR)/mpl_routines.F90\
			$(OBJDIR)/eigen_parameters.o\
			$(OBJDIR)/file_io.o\
			$(OBJDIR)/fullmatrix_aux.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/diag.o\
			$(OBJDIR)/diag_energy.o\
			$(OBJDIR)/scalapack_aux.o

$(OBJDIR)/eigenvalue_comp.o: $(OBJDIR)/eigen_direct.o

$(OBJDIR)/field_solve_df.o:	$(SRCDIR)/switches.h \
			$(OBJDIR)/derivative_matrix.o \
			$(OBJDIR)/matrix.o



$(OBJDIR)/flr_corr.o:		$(SRCDIR)/switches.h\
			$(OBJDIR)/matrix.o\
			$(OBJDIR)/storefullmatrix.o\
			$(OBJDIR)/BandedMatrix.o

$(OBJDIR)/full_rhs_mat.o:	$(OBJDIR)/impl_scalapack.o

$(OBJDIR)/gyro_average_df.o:	$(SRCDIR)/switches.h\
			$(OBJDIR)/derivative_matrix.o\
			$(OBJDIR)/grid1d.o\
			$(OBJDIR)/matrix.o\
			$(OBJDIR)/Vector.o\
			$(OBJDIR)/BandedMatrix.o\
			$(OBJDIR)/boundary_exchange_x.o


$(OBJDIR)/impl_scalapack.o:	$(OBJDIR)/par_mod.o\
			$(OBJDIR)/comm.o\
			$(OBJDIR)/scalapack_aux.o\
			$(OBJDIR)/fullmatrix_aux.o\
			$(OBJDIR)/calc_rhs.o

$(OBJDIR)/ListObject.o:	$(SRCDIR)/ListObjectTemplate.f90

$(OBJDIR)/matrix.o:	$(SRCDIR)/matrix.h \
			$(OBJDIR)/processgrid.o\
			$(OBJDIR)/storefullmatrix.o

$(OBJDIR)/processgrid.o:

$(OBJDIR)/rhs_computations.o:	$(OBJDIR)/impl_scalapack.o

$(OBJDIR)/scalapack_aux.o:	$(OBJDIR)/par_mod.o

$(OBJDIR)/storebandedmatrix.o:	$(OBJDIR)/ListObject.o\
			$(OBJDIR)/storevector.o\
			$(OBJDIR)/storefullmatrix.o \
			$(OBJDIR)/processgrid.o

$(OBJDIR)/storefullmatrix.o:	$(OBJDIR)/processgrid.o\
			$(OBJDIR)/storevector.o

$(OBJDIR)/storevector.o:	$(OBJDIR)/processgrid.o


$(OBJDIR)/Vector.o:		$(OBJDIR)/processgrid.o \
				$(OBJDIR)/matrix.o \
				$(OBJDIR)/storevector.o

else
#-------- NO SCALAPACK ----------

$(OBJDIR)/adiabatic_response.o: $(OBJDIR)/dummy_matrix_module.o
$(OBJDIR)/comm.o:	$(OBJDIR)/dummy_matrix_module.o
$(OBJDIR)/derivative_matrix.o:	$(OBJDIR)/dummy_matrix_module.o
$(OBJDIR)/dummy_matrix_module.o:$(SRCDIR)/switches.h \
				$(OBJDIR)/grid1d.o
$(OBJDIR)/field_solve_df.o:	$(OBJDIR)/dummy_matrix_module.o
$(OBJDIR)/flr_corr.o:	$(OBJDIR)/dummy_matrix_module.o
$(OBJDIR)/gyro_average_df.o:	$(OBJDIR)/dummy_matrix_module.o
endif

################################
######## link MKL dfti #########
################################

ifeq ($(FFTLIB),mkl)
OBJLIST += $(OBJDIR)/mkl_dfti.o
endif

################################
##### FUTILS dep. parts ########
################################
ifeq ($(FUTILS),yes)

$(OBJDIR)/geometry.o:			$(FUTILSDIR)/libfutils.a
$(OBJDIR)/chease.o:			$(FUTILSDIR)/libfutils.a
$(OBJDIR)/convergence_monitoring.o:	$(FUTILSDIR)/libfutils.a
$(OBJDIR)/diag.o:			$(FUTILSDIR)/libfutils.a
$(OBJDIR)/parameters_IO.o:		$(FUTILSDIR)/libfutils.a
$(OBJDIR)/profiles.o:			$(FUTILSDIR)/libfutils.a
$(OBJDIR)/profile_IO.o:			$(FUTILSDIR)/libfutils.a
$(OBJDIR)/checkpoint.o:			$(FUTILSDIR)/libfutils.a
$(OBJDIR)/eigen_iterative.o:		$(FUTILSDIR)/libfutils.a
$(OBJDIR)/diag_vsp.o:			$(FUTILSDIR)/libfutils.a
$(OBJDIR)/create_nonuni_grid.o:		$(FUTILSDIR)/libfutils.a



$(OBJDIR)/poloidal_planes_nonuni.o:	$(FUTILSDIR)/libfutils.a
$(OBJDIR)/poloidal_planes.o:		$(FUTILSDIR)/libfutils.a
$(OBJDIR)/import_XGC_grid.o:            $(FUTILSDIR)/libfutils.a
ifeq ($(COUPLE_XGC),yes)
else

$(OBJDIR)/diag_energy.o:		$(FUTILSDIR)/libfutils.a
$(OBJDIR)/diag_extended.o:		$(FUTILSDIR)/libfutils.a

endif

endif



################################
### ADIOS dep. parts
################################

ifeq ($(ADIOS),yes)
    adios_deps= gene_subroutine.o \
		checkpoint.o \
		diag.o #\
		diag_coupled.o

    ADIOS_DEPS = $(addprefix $(OBJDIR)/, $(adios_deps))
    $(ADIOS_DEPS): $(OBJDIR)/adios_io.o

    $(OBJDIR)/gene_subroutine.o: $(OBJDIR)/discretization.o

    adios_uses = par_mod.o \
		 par_in.o \
		 antenna.o\
#ifneq ($(COUPLE_XGC),yes)
    adios_uses = aux_fields.o
#endif

    ADIOS_USES = $(addprefix $(OBJDIR)/, $(adios_uses))
    $(OBJDIR)/adios_io.o: $(ADIOS_USES)
endif

ifeq ($(ADIOS2),yes)
$(OBJDIR)/adios2_comm_mod.o:    $(OBJDIR)/discretization.o

$(OBJDIR)/diag.o:               $(OBJDIR)/adios2_comm_mod.o

$(OBJDIR)/gene_subroutine.o:    $(OBJDIR)/adios2_comm_mod.o
endif




ifeq ($(COUPLE_XGC),yes)
ifeq ($(ADIOS2),yes)
    $(OBJDIR)/gene_subroutine.o: $(OBJDIR)/adios2_comm_mod.o

    $(OBJDIR)/parameters_IO.o:   $(OBJDIR)/adios2_comm_mod.o
endif
ifeq ($(ADIOS),yes)
    $(OBJDIR)/gene_subroutine.o: $(OBJDIR)/adios_comm_mod.o

    $(OBJDIR)/parameters_IO.o:   $(OBJDIR)/adios_comm_mod.o
endif
$(OBJDIR)/gene_subroutine.o $(OBJDIR)/field_solver.o $(OBJDIR)/poloidal_planes.o $(OBJDIR)/poloidal_planes_nonuni.o $(OBJDIR)/import_XGC_grid.o: $(OBJDIR)/coupling_core_gene.o

    $(OBJDIR)/adios_io.o: $(OBJDIR)/coupling_core_gene.o

    $(OBJDIR)/coupling_core_gene.o:	$(OBJDIR)/coupling_base_setup.o \
					$(OBJDIR)/geometry.o \
					$(OBJDIR)/par_mod.o \
					$(OBJDIR)/par_other.o \
					$(OBJDIR)/discretization.o \
					$(OBJDIR)/par_poloidal_planes.o\
					$(OBJDIR)/geometry.o
ifeq ($(ADIOS),yes)
    $(OBJDIR)/coupling_core_gene.o:     $(OBJDIR)/adios_comm_mod.o
endif

ifeq ($(ADIOS2),yes)
    $(OBJDIR)/coupling_core_gene.o:     $(OBJDIR)/adios2_comm_mod.o
endif


    $(OBJDIR)/calc_rhs.o: $(OBJDIR)/diag.o
else

ifeq ($(INIT_XGC),yes)
    $(OBJDIR)/gene_subroutine.o: $(OBJDIR)/adios2_comm_mod.o

    $(OBJDIR)/parameters_IO.o:   $(OBJDIR)/adios2_comm_mod.o

    $(OBJDIR)/gene_subroutine.o $(OBJDIR)/field_solver.o $(OBJDIR)/poloidal_planes.o $(OBJDIR)/poloidal_planes_nonuni.o $(OBJDIR)/import_XGC_grid.o: $(OBJDIR)/coupling_core_gene.o

    $(OBJDIR)/adios_io.o: $(OBJDIR)/coupling_core_gene.o

    $(OBJDIR)/coupling_core_gene.o:	$(OBJDIR)/coupling_base_setup.o \
					$(OBJDIR)/adios2_comm_mod.o \
					$(OBJDIR)/par_mod.o \
					$(OBJDIR)/par_other.o \
					$(OBJDIR)/discretization.o \
					$(OBJDIR)/par_poloidal_planes.o\
					$(OBJDIR)/geometry.o
endif

endif

gptl_obj = gene_subroutine.o initial_value_comp.o diag.o adios_io.o checkpoint.o field_solver.o poloidal_planes_nonuni.o coupling_core_edge.o
GPTL_OBJ = $(addprefix $(OBJDIR)/, $(gptl_obj))
$(OBJDIR)/perf_monitor.o:	$(OBJDIR)/par_in.o\
				$(OBJDIR)/file_io.o
$(GPTL_OBJ): $(OBJDIR)/perf_monitor.o



####################################
##### PERFormance libraries ########
####################################
ifneq ($(strip $(USE_PERFLIB)),)
 ifeq ($(USE_PERFLIB),FR)
  $(OBJDIR)/libmyperf.a: $(MKFILE) $(OBJDIR)/FR_perf.o $(OBJDIR)/dummyperf_addon.o
	@echo "making libmyperf.a, now in the rule"
	@rm -f $(OBJDIR)/libmyperf.a
	$(ARCHIVE) $@ $(OBJDIR)/FR_perf.o $(OBJDIR)/dummyperf_addon.o
 endif

# ifeq ($(USE_PERFLIB),scalasca)
#  $(OBJDIR)/libmyperf.a: $(MKFILE) $(OBJDIR)/dummyperf.o $(OBJDIR)/dummyperf_addon.o
#	@echo "making libmyperf.a, now in the rule"
#	@rm -f $(OBJDIR)/libmyperf.a
#	$(ARCHIVE) $@ $(OBJDIR)/dummyperf_addon.o
# endif

 ifeq ($(USE_PERFLIB),scalasca)
 else ifneq ($(USE_PERFLIB),perf)
  LIBS += -L$(OBJDIR) -lmyperf
  LIBELEMENTS += $(OBJDIR)/libmyperf.a
  $(EXECDIR)/$(EXEC): $(OBJDIR)/libmyperf.a
 endif
endif

ifeq ($(COMPILE_MPI_MOD),yes)
$(MPIOBJ): $(OBJDIR)/mpimod.o
endif

####################################
OBJLIST += $(MEMCHECK_OBJ)
