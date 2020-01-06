#!/bin/bash -l
## MAXWT 24
## PROCSPERNODE 24
## SUBMITCMD sbatch
#SBATCH -p regular
#SBATCH -n NMPIPROCS
#SBATCH --ntasks-per-node PROCSPERNODE
#SBATCH -t WALLTIME
#SBATCH -J JOBNAME

## fix formatted output in cray env
if [ "$PE_ENV" == "CRAY" ]; then
  export FILENV=my_filenenv
  assign -U on g:sf
fi

## set openmp threads
export OMP_NUM_THREADS=1

# run GENE
srun -n $SLURM_NTASKS ./gene_edison
