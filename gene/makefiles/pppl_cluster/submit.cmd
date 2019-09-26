#!/bin/tcsh
#SBATCH --partition=dawson
#SBATCH --ntasks=8
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=2GB
#SBATCH --job-name=GENE_STDTST


# --- $NPROCS is required by mpirun.
set NPROCS=`wc -l < $PBS_NODEFILE`

cd $SLURM_SUBMIT_DIR

setenv OMP_NUM_THREADS 1

module load gcc/6.1.0
module load acml/6.1.0/gfortran64
module load openmpi/3.0.0
module load blacs fftw hdf5-parallel/1.10.1 scalapack
module load petsc_complex/3.8.3
module load slepc_complex/3.8.3

mpirun -np $SLURM_NTASKS ./gene_pppl_cluster


### to submit a parameter scan, comment the previous line and uncomment
### the following (see GENE documentation or ./scanscript --help)
#./scanscript --np $SLURM_NTASKS --ppn 8 --mps 4 --syscall='srun -n $SLURM_NTASKS ./gene_pppl_cluster'
