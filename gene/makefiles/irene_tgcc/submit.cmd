#!/bin/bash
##submit job: ccc_msub submit.cmd
##cancel job: ccc_mdel <job id>
#MSUB -r GENE              ### job name
##MSUB -N 1                ### Total number of nodes (usually not needed)
#MSUB -n 24                ### MPI tasks
#MSUB -c  1                ### Number of cores per task (OpenMP)
#MSUB -T 86000             ### wall clock time in seconds
#MSUB -q skylake
#MSUB -m work,scratch
#MSUB -o ./GENE.%j.out
#MSUB -e ./GENE.%j.err
## adjust project account in following line
#MSUB -A ra4403            ### see ccc_myproject

set -x
cd ${BRIDGE_MSUB_PWD}

### switch to project specific datadir automatically
PROJECTACC=${SLURM_JOB_ACCOUNT%"@$SLURM_JOB_PARTITION"}
module switch dfldatadir dfldatadir/$PROJECTACC

##set openmp threads
export OMP_NUM_THREADS=1

##execute
ccc_mprun ./gene_irene_tgcc

### to submit a parameter scan, comment out the previous line and uncomment
### the following (see GENE documentation or ./scanscript --help)
#./scanscript --np $SLURM_NTASKS --ppn 24 --mps 4 --syscall="ccc_mprun ./gene_irene_tgcc"
