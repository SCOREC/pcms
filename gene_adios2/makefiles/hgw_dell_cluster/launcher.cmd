#!/bin/bash
### batch jobs run in /bin/bash per default
### MAXWT 24
### PROCSPERNODE 16
### SUBMITCMD qsub
### 
###
### check 'man qsub' for all available options
###
### start the job in the current working dir
#$ -cwd
### send notification (on (a)bort, (b)egin, (e)nd, (n)ever)
#$ -m n
###
### start a parallel environment with multiples of 8 processes
#$ -pe edge*_ompi_b NMPIPROCS
###
### wallclock limit
#$ -l h_rt=WALLTIME
###
### the job's name
#$ -N JOBNAME

#load environment
module load petsc/3.5-intel13-ompi16-cplx mkl/11.0

### set MKL to serial mode 
export MKL_SERIAL=yes

### no OpenMP
export OMP_NUM_THREADS=1

### start program
### Note: when running with OMP_NUM_THREADS>1, you may need to replace
### $NSLOTS by the actual number of mpi processes and set 
### -perhost <<MPI-PROCS-PER-HOST>> flag

mpirun -np NMPIPROCS ./gene_hgw_dell_cluster
