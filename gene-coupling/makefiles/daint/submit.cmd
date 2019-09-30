#!/bin/sh -l

#SBATCH --ntasks=576          ### Total number of MPI tasks (36 cores/node)
#SBATCH -N 16                 ### Number of nodes
#SBATCH --cpus-per-task=1     ### Number of threads per task (OpenMP)
#SBATCH --ntasks-per-node=36  ### Tasks per node, ignored if cannot match previous 3
#SBATCH --ntasks-per-core=1   ### Set to 2 for hypertheading
#SBATCH --time=24:00:00       ### wall clock time
#SBATCH --constraint=mc       ### mc=multicore, gpu=gpu
#SBATCH --job-name=GENE
##SBATCH --account=<ACCOUNT-NAME>

export NTASK=$SLURM_NTASKS
export TASKS_N=$SLURM_NTASKS_PER_NODE
export CPU_T=$SLURM_CPUS_PER_TASK
export OMP_NUM_THREADS=1

##based on 64 GB/node &36 CPU/node
export MEMORY_PER_CORE=1600
export N_PES=$SLURM_NTASKS

#set the stack size to unlimited
ulimit -s unlimited

# Launch the parallel job to the allocated compute nodes
srun -n $SLURM_NTASKS --ntasks-per-node=$SLURM_NTASKS_PER_NODE -c $SLURM_CPUS_PER_TASK  ./gene_daint

### to submit a parameter scan, comment the previous line and uncomment
### the following (see GENE documentation or ./scanscript --help)
#./scanscript --np $NTASK --ppn $SLURM_NTASKS_PER_NODE --mps 4
