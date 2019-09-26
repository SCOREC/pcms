#!/bin/bash -l
#SBATCH -J GENE              ### job name
#SBATCH --nodes=1            ### Total number of nodes
#SBATCH --ntasks-per-node=40 ### MPI tasks (max. 40 / hyperthreading 80)
#SBATCH --cpus-per-task=1    ### Number of threads per task (OpenMP)
#SBATCH --time=24:00:00      ### wall clock time
#SBATCH --partition=medium   ### see sinfo (express,tiny,medium,n0064,n0128,etc.)
#SBATCH -o ./%x.%j.out
#SBATCH -e ./%x.%j.err

export SLURM_HINT=nomultithread
export OMP_PROC_BIND=true

##set openmp threads
export OMP_NUM_THREADS=1

##execute
srun ./gene_cobra

### to submit a parameter scan, comment the previous line and uncomment
### the following (see GENE documentation or ./scanscript --help)
#./scanscript --np $SLURM_NTASKS --ppn 40 --mps 4 --syscall="srun ./gene_cobra"
