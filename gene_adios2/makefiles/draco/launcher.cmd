#!/bin/bash -l
#
## MAXWT 24
## PROCSPERNODE 32
## SUBMITCMD sbatch
#SBATCH --nodes=NODES
#SBATCH --ntasks-per-node=PPN
#SBATCH --cpus-per-task=1
#SBATCH --time=WALLTIME
#SBATCH -J JOBNAME
#SBATCH --partition=general
#SBATCH -o ./GENE_out.%j

#export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so
export SLURM_HINT=nomultithread
export OMP_PROC_BIND=true

##set openmp threads
export OMP_NUM_THREADS=1

##execute
srun ./gene_draco
