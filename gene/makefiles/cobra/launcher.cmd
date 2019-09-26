#!/bin/bash -l
#
## MAXWT 24
## PROCSPERNODE 40
## SUBMITCMD sbatch
#SBATCH --nodes=NODES
#SBATCH --ntasks-per-node=PPN
#SBATCH --cpus-per-task=1
#SBATCH --time=WALLTIME
#SBATCH -J JOBNAME
#SBATCH --partition=n0064
#SBATCH -o ./%x.%j.out
#SBATCH -e ./%x.%j.err

export SLURM_HINT=nomultithread
export OMP_PROC_BIND=true

##set openmp threads
export OMP_NUM_THREADS=1

##execute
srun ./gene_cobra
