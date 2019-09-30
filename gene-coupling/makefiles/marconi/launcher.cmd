#!/bin/bash
#
## MAXWT 24
## PROCSPERNODE 48
## SUBMITCMD qsub
#PBS -N JOBNAME
#PBS -q xfuasklprod
#PBS -l walltime=WALLTIME
#PBS -l select=NODES:ncpus=PPN:mpiprocs=PPN:mem=180GB
#PBS -j oe
#PBS -A FUSIO_ALL

cd $PBS_O_WORKDIR
module load intel intelmpi mkl fftw
module load blas lapack scalapack
export OMP_NUM_THREADS=1

#determine total number of mpi procs, e.g. for scanscript applications
export NPROC=$(wc -l $PBS_NODEFILE | awk '{print $1}')

mpirun ./gene_marconi
#
