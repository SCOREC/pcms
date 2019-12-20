#!/bin/bash -l
#SBATCH -p debug
#SBATCH -N 8
#SBATCH --ntasks-per-node 24
#SBATCH -t 30:00
#SBATCH -J GENE

## set openmp threads
export OMP_NUM_THREADS=1
#do not use file locking for hdf5
export HDF5_USE_FILE_LOCKING=FALSE

## fix formatted output in cray env
if [ "$PE_ENV" == "CRAY" ]; then
  export FILENV=my_filenenv
  assign -U on g:sf
fi

# run GENE
#srun -n $SLURM_NTASKS ./gene_edison

# run scanscript
./scanscript --np $SLURM_NTASKS --ppn 24 --mps 4 --syscall='srun -n $SLURM_NTASKS ./gene_edison' -cs
