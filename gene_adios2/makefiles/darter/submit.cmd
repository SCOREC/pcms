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
#PBS -n GENE
#PBS -o GENE.$PBS_JOBID.oe
#PBS -j oe
#PBS -l size=16,walltime=00:10:00

#module load impi intel mkl fftw petsc slepc
module load fftw

export OMP_NUM_THREADS=1

### set MKL to serial mode
#export MKL_SERIAL=yes

### enable/disable Shared Receive Queue
### (currently we observe sporadical hang-ups if enabled)
#export MV2_USE_SRQ=0

### uncomment the following three lines if MPI memory error occur
#export MPICH_MAX_SHORT_MSG_SIZE=50000
#export MPICH_PTL_UNEX_EVENTS=80000
#export MPICH_UNEX_BUFFER_SIZE=1000M

cd $PBS_O_WORKDIR
#note: PBS_NNODES is number of cores, NOT number of nodes
aprun -n $PBS_NNODES ./gene_darter

#uncomment the block below (and comment out the above aprun) to use
#hyperthreading; CAUTION: tests saw no speedup!
#nvirtcores=`expr $PBS_NNODES \* 2`
#echo Using $nvirtcores cores in hyperthreading mode...
#aprun -j 2 -n $nvirtcores ./gene_darter

### to submit a parameter scan, comment the previous line and uncomment
### the following (see GENE documentation or ./scanscript --help)
#./scanscript --np 16 --ppn 1 --mps 4
