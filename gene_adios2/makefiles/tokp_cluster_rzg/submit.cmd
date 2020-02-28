#!/bin/bash
### batch jobs run in /bin/bash per default
###
### check 'man qsub' for all available options
###
### job has the same environmnent variables as 
### the submission shell (currently disabled)
###$ -V
### join stdout and stderr
#$ -j y
### start the job in the current working dir
#$ -cwd
### send notification (on (a)bort, (b)egin, (e)nd, (n)ever)
#$ -m n
###
### start a parallel environment with multiples of 16 processes
#$ -pe impi_hydra.* 16
### set number of (OMP) threads per MPI task
#$ -v THREADS_PER_MPI_TASK=1
###
### wallclock limit
#$ -l h_rt=24:00:00
###
### set queue (tokp, debug)
#$ -P tokp
###
### the job's name
#$ -N GENE

#load environment (sdr-infiniband with intel compilers)
module load impi intel mkl fftw petsc-cmplx slepc-cmplx
module load hdf5-mpi

### set OMP_NUM_THREADS
export OMP_NUM_THREADS=$THREADS_PER_MPI_TASK

### set MKL to serial mode 
export MKL_SERIAL=yes

### start program
### Note: when running with OMP_NUM_THREADS>1, you may need to replace
### $NSLOTS by the actual number of mpi processes and set 
### -perhost <<MPI-PROCS-PER-HOST>> flag

mpiexec -l -n $NSLOTS ./gene_tokp_cluster_rzg

### to submit a parameter scan, comment the previous line and uncomment
### the following (see GENE documentation or ./scanscript --help)
#./scanscript --np $NSLOTS --ppn 16 --mps 4
