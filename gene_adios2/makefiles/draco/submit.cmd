#!/bin/bash -l
#SBATCH --nodes=1            ### Total number of nodes
#SBATCH --ntasks-per-node=32 ### MPI tasks (32 or 64)
#SBATCH --cpus-per-task=1    ### Number of threads per task (OpenMP)
#SBATCH --time=24:00:00      ### wall clock time
#SBATCH -J GENE              ### job name
#SBATCH --partition=general
#SBATCH -o ./GENE_out.%j

#export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so
export SLURM_HINT=nomultithread
export OMP_PROC_BIND=true

##set openmp threads
export OMP_NUM_THREADS=1

##execute
srun ./gene_draco

### to submit a parameter scan, comment the previous line and uncomment
### the following (see GENE documentation or ./scanscript --help)
#./scanscript --np $SLURM_NTASKS --ppn 32 --mps 4 --syscall="srun ./gene_draco"
