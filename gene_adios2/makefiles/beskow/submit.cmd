#!/bin/bash -l
# The -l above is required to get the full environment with modules

## Set the allocation to be charged for this job
## not required if you have set a default allocation
##SBATCH -A 201X-X-XX

## The name of the job
#SBATCH -J GENE

## Walltime limit, maximum is 24h
#SBATCH -t 24:00:00

## Number of MPI processes (use multiples of 32).
#SBATCH -n 32
## Number of nodes (32 cores per node)
#SBATCH --nodes=1


# Load modules
module load PrgEnv-intel fftw intel cray-petsc-complex cray-hdf5

# Run GENE:
aprun -n $SLURM_NTASKS ./gene_beskow

### to submit a parameter scan, comment the previous line and uncomment
### the following (see GENE documentation or ./scanscript --help)
#./scanscript --np $SLURM_NTASKS --ppn 32 --mps 4
