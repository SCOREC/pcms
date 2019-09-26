#!/bin/bash
##submit with: qsub submit.cmd
##in the directory containing
##the gene executable
##
##charge to fus110 project
#PBS -A fus110

##job name:
#PBS -N GENE
##queue (regular/debug):
##PBS -q regular
##number of cores:
#PBS -l nodes=2
##walltime:
#PBS -l walltime=1:00:00
##join stdout/stderr
#PBS -j eo
##copy environment from terminal
##PBS -V
##send mail for (a)bort, (b)egin, (e)nd
#PBS -m abe

# Change to the direcotry that the job was submitted from
cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=1

# Launch the parallel job to the allocated compute nodes
aprun -n 32 -N 16 ./gene_titan

### to submit a parameter scan, comment the previous line and uncomment
### the following (see GENE documentation or ./scanscript --help)
#./scanscript --np $NPROCS --ppn $NTASKS  --mps 4
