#!/bin/bash
#SBATCH -J GENE
#SBATCH --nodes=5
#SBATCH --ntasks=120           ### Total number of MPI tasks (24 cores/node)
#SBATCH --ntasks-per-node=24
#SBATCH --threads-per-core=1
#SBATCH --mem=64000
#SBATCH --time=24:00:00
##uncomment to use user defined stdout file name
##SBATCH --output gene.out

module purge
module load intel
#module load intelmpi
module load bullxmpi
module load hdf5 mkl fftw3

##set openmp threads
export OMP_NUM_THREADS=1

##execute
srun --mpi=pmi2 -K1 --resv-ports -n $SLURM_NTASKS ./gene_occigen

### to submit a parameter scan, comment the previous line and uncomment
### the following (see GENE documentation or ./scanscript --help)
#./scanscript --np $SLURM_NTASKS --ppn 24 --mps 4
