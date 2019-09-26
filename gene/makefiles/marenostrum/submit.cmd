#!/bin/bash

#SBATCH --ntasks=768        ### Total number of MPI tasks (48 cores/node)
##SBATCH -N 16              ### Number of nodes; sometimes req'd
#SBATCH --cpus-per-task=1   ### Number of threads per task (OpenMP)
#SBATCH --time=04:00:00     ### wall clock time
#SBATCH --job-name=GENE

##list nodes for debugging
printenv SLURM_NODELIST

##set openmp threads
export OMP_NUM_THREADS=1

##execute
srun ./gene_marenostrum

### to submit a parameter scan, comment the previous line and uncomment
### the following (see GENE documentation or ./scanscript --help)
#./scanscript --np $SLURM_NTASKS --ppn 48 --mps 4
