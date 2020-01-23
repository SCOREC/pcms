#!/bin/bash
##submit with: qsub submit.cmd
##in the directory containing
##the gene executable
##
##job name:
#PBS -N GENE
##queue (regular/debug):
#PBS -q regular
##number of cores:
#PBS -l mppwidth=256
##walltime:
#PBS -l walltime=24:00:00
##join stdout/stderr
#PBS -j eo
##copy environment from terminal
#PBS -V
##send mail for (a)bort, (b)egin, (e)nd
#PBS -m abe

# Change to the direcotry that the job was submitted from
cd $PBS_O_WORKDIR

export NPROC=`qstat -f $PBS_JOBID | awk '/mppwidth/ {print $3}'`
export NTASKS=`qstat -f $PBS_JOBID | awk '/mppnppn/  {print $3}'`
export OMP_NUM_THREADS=1

# Launch the parallel job to the allocated compute nodes
aprun -n $NPROC -N $NTASKS ./gene_hopper

### to submit a parameter scan, comment the previous line and uncomment
### the following (see GENE documentation or ./scanscript --help)
#./scanscript --np $NPROCS --ppn $NTASKS  --mps 4
