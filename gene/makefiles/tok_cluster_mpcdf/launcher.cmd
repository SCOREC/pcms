#!/bin/bash -l
#
## MAXWT 48
## PROCSPERNODE 32
## SUBMITCMD sbatch
#SBATCH --nodes=NODES
#SBATCH --ntasks-per-node=PPN
#SBATCH --cpus-per-task=1
#SBATCH --time=WALLTIME
#SBATCH -J JOBNAME
#SBATCH --partition=p.tok
#SBATCH --qos=p.tok.48h
#SBATCH -o ./%x.%j.out
#SBATCH -e ./%x.%j.err

##set openmp threads
export OMP_NUM_THREADS=1

##execute
srun ./gene_tok_cluster_mpcdf
