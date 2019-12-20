.SUFFIXES: .o .F90

CMD=xgc2
ALL:$(CMD)

# include ${PETSC_DIR}/conf/variables
# include ${PETSC_DIR}/conf/rules

OBJ=module.o $(EXTRA_OBJ) bicub_mod.o search.o psmooth.o pol_decomp.o f0module.o charge.o \
	main.o read.o gen_perm.o sort_particles.o push.o pushe.o load.o one_d_cub_mod.o\
	setup.o efield.o interpolation.o $(MPI_OBJ) diagnosis.o limiter.o \
	bounce.o diagnosis2.o collision.o collision2.o collisionf.o collisionf2.o \
	elliptics.o diagnosis-f.o heat.o \
	turbulence.o neutral.o neutral2.o linearsolver.o fem2d.o \
	fem_ops.o poisson.o new_petsc_solver.o 
IMSL_OBJ=my_imsl.o
PORT_OBJ=bspline90_22.o taus88.o derf.o datanh.o pppack.o fmin.o
SER_OBJ=mpisingle.o
PAR_OBJ=mpi.o

