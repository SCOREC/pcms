#!/bin/bash -l

### please make sure to run
### module swap craype-haswell craype-mic-knl
### BEFORE compiling GENE to set the correct compiler optimization and memory
#SBATCH -p debug
#SBATCH -N 2 ###number of nodes
#SBATCH -S 4 ###number of specialized cores per node (use 4 to get 64 cores per node for GENE)
#SBATCH -C knl,quad,cache  ###default setting on Cori, means 16GB fast RAM will be used as L3 cache; leaving 96GB per node
#SBATCH -t 0:30:00
#SBATCH -J GENE
#SBATCH -L SCRATCH ###need scratch access for this job
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

# $SLURM_NTASKS does not seem to work under these conditiosn, need to specify number of cores explicitly for now using -n <#cores>.
# -c 4 means that 4 logical cores (the 4 hyperthreads) will be used per MPI process, so when running 64 MPI procs per node,
# we will have one physical core per process
# --cpu_bind=cores is needed to enforce this; please read the NERSC "Example Batch Scripts for KNL" page for more info

# run GENE
srun -n 128 -c 4 --cpu_bind=cores ./gene_cori

# run scanscript
#./scanscript --np $SLURM_NTASKS --ppn 32 --mps 4 --syscall='srun -n $SLURM_NTASKS ./gene_cori'
