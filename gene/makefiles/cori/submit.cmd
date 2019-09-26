#!/bin/bash -l
#SBATCH -p regular
#SBATCH -n 32
#SBATCH -t 24:00:00
#SBATCH -J GENE
#SBATCH -C haswell
#SBATCH -o ./%x.%j.out
#SBATCH -e ./%x.%j.err
##uncomment for particular account usage
##SBATCH -A <ACCOUNT>

## fix formatted output in cray env
if [ "$PE_ENV" == "CRAY" ]; then
  export FILENV=my_filenenv
  assign -U on g:sf
fi

## set openmp threads
export OMP_NUM_THREADS=1

# run GENE
srun -n $SLURM_NTASKS ./gene_cori

# run scanscript
#./scanscript --np $SLURM_NTASKS --ppn 32 --mps 4 --syscall='srun -n $SLURM_NTASKS ./gene_cori'
