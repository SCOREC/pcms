#!/bin/bash
### MAXWT 24
### SUBMITCMD llsubmit
### PROCSPERNODE 16
#@ wall_clock_limit = WALLTIME
#@ job_type = parallel
#@ job_name = JOBNAME
#@ class = general
#@ network.MPI = sn_all,not_shared,us
#@ node = NODES
#@ tasks_per_node = PPN
#@ output = $(job_name).$(jobid).out
#@ error = $(job_name).$(jobid).err
#@ queue

. /etc/profile
. /etc/profile.d/modules.sh

#setup of environment
module load intel/17.0 mkl/2017 git
module load fftw/serial petsc-complex/3.9.1c_gene
export SLEPC_DIR=/home/hpc/pr27fe/lu65nut2/tmp/slepc-3.9.1       

export MKL_SERIAL=YES
export OMP_NUM_THREADS=1

poe ./gene_supermuc
