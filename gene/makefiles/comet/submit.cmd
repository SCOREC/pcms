#!/bin/sh
#SBATCH --ntasks-per-node=24  ### Number of tasks per node
#SBATCH --nodes=1           ### Number of nodes; sometimes req'd
#SBATCH --cpus-per-task=1    ### Number of threads per task (OpenMP)
#SBATCH -t 00:30:00          ### wall clock time
#SBATCH --job-name=GENE

##list nodes for debugging
printenv SLURM_NODELIST

##uncomment next line for core files
#ulimit -c unlimited

##uncomment ulimit line to increase stack size
#ulimit -s <your-desired-stack-size>

module load intel mvapich2_ib fftw mkl hdf5 scalapack

##set openmp threads
export OMP_NUM_THREADS=1
export SLURM_NODEFILE=`generate_pbs_nodefile`

##execute
mpirun_rsh -np $SLURM_NTASKS -hostfile $SLURM_NODEFILE ./gene_comet

### to submit a parameter scan, comment the previous line and uncomment
### the following (see GENE documentation or ./scanscript --help)
#./scanscript --np $SLURM_NTASKS --ppn 24 --mps 4
