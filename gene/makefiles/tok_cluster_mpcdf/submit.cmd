#!/bin/bash -l
#SBATCH --nodes=2            ### Total number of nodes
#SBATCH --ntasks-per-node=32 ### MPI tasks (32)
#SBATCH --cpus-per-task=1    ### Number of threads per task (OpenMP)
#SBATCH --time=24:00:00      ### wall clock time
#SBATCH -J GENE              ### job name
#SBATCH --partition=p.tok    ### s.tok for nodes=1, p.tok else
#SBATCH --qos=p.tok.48h      ### p.tok.2h,p.tok.48h,tok.debug
#SBATCH -o ./%x.%j.out
#SBATCH -e ./%x.%j.err

export SLURM_HINT=nomultithread

##set openmp threads
export OMP_NUM_THREADS=1

##execute
srun -l ./gene_tok_cluster_mpcdf

### to submit a parameter scan, comment the previous line and uncomment
### the following (see GENE documentation or ./scanscript --help)
#./scanscript --np $SLURM_NTASKS --ppn 32 --mps 4 --syscall="srun ./gene_tok_cluster_mpcdf"
