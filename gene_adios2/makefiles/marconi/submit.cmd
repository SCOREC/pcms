#!/bin/bash
### job name:
#PBS -N GENE
###see qstat -Q for available queues (xfuaprod,xfuadebug,xfuaknlprod,xfuasklprod etc):
#PBS -q xfuasklprod
## walltime depends on chosen queue:
#PBS -l walltime=01:00:00
### A3 partition: N nodes:procs per node(max 48):mpiprocs per node(max 48):
#PBS -l select=1:ncpus=48:mpiprocs=48:mem=180GB
### A2 partition: N nodes:procs per node(max 68):mpiprocs per node(max 68):
## add # to PBS -l line above and remove one # below to use A2 partition; requires xfuaknlprod queue
##PBS -l select=1:ncpus=68:mpiprocs=68:mem=86GB:mcdram=cache:numa=quadrant
### join out/err streams:
#PBS -j oe
### define account name, see "saldo -b" (A1), "saldo -b --knl" (A2), or "saldo -b --skl" (A3)
#PBS -A FUSIO_ALL

# Change to the direcotry that the job was submitted from
cd $PBS_O_WORKDIR

#set environment
module load intel intelmpi mkl fftw
module load blas lapack scalapack

export OMP_NUM_THREADS=1

#determine total number of mpi procs, e.g. for scanscript applications
export NPROC=$(wc -l $PBS_NODEFILE | awk '{print $1}')

# Launch the parallel job on the allocated compute nodes
mpirun ./gene_marconi

### to submit a parameter scan, comment the previous line and uncomment
### the following (see GENE documentation or ./scanscript --help)
#./scanscript --np $NPROC --ppn $NCPUS --mps 4

#
