#!/bin/bash -l

### Please check the queue prior to decide how many CPU to use.
### Using 24 CPU is more than OK for well resolved linear scans.
### Maximum available is 32; you can use all of them especially
### if the machine is empty (likely in weekends).
### Abuse of resource is the kind of thing that makes me upset and
### willing to do everything to prevent you from running if you do
### it systematically more than once.
### I guess everybody feels the same.
### -G. Merlo

#SBATCH --ntasks=24          ### Total number of MPI tasks (4 cores/node)
#SBATCH --cpus-per-task=1    ### Number of threads per task (OpenMP)
#SBATCH --time=72:00:00      ### wall clock time, already maximum
#SBATCH --job-name=GENE
#SBATCH --ntasks-per-node=4  ### Switch-off hyperthreading as no improvement
#SBATCH --ntasks-per-core=1  ### Switch-off hyperthreading as no improvement

##list nodes for debugging
printenv SLURM_NODELIST

##set openmp threads
export OMP_NUM_THREADS=1

##execute
srun -n $SLURM_NTASKS ./gene_ppb110

### to submit a parameter scan, comment the previous line and uncomment
### the following (see GENE documentation or ./scanscript --help)
#./scanscript --np $SLURM_NTASKS --ppn 4 --mps 1
