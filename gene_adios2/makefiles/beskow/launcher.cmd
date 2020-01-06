# @ shell=/bin/bash
### MAXWT 24
### PROCSPERNODE 32
### SUBMITCMD sbatch

## The name of the job
#SBATCH -J JOBNAME

## Walltime limit, maximum is 24h
#SBATCH -t WALLTIME

## Number of MPI processes (use multiples of 32).
#SBATCH -n NMPIPROCS

#SBATCH -e error_file.e
#SBATCH -o output_file.o

# Load modules
module load PrgEnv-intel fftw intel cray-petsc-complex cray-hdf5

# Run GENE:
aprun -n NMPIPROCS ./gene_beskow > gene.out 2>&1
