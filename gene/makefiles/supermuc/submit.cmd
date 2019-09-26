#!/bin/bash
# DO NOT USE environment = COPY_ALL
#@ wall_clock_limit = 0:30:00
#@ job_type = parallel
###@ job_type = MPICH 
#### for intel MPI
#@ job_name = GENE
#@ class = test
#@ network.MPI = sn_all,not_shared,us
#@ energy_policy_tag = global_gene
#@ node = 1
#@ tasks_per_node = 8
#@ output = $(job_name).$(jobid).out
#@ error = $(job_name).$(jobid).err
#@ queue

. /etc/profile
. /etc/profile.d/modules.sh

export MKL_SERIAL=YES
export OMP_NUM_THREADS=1

#setup of environment
module load intel/17.0 mkl/2017 git
module load fftw/serial petsc/3.9.1c_gene
export SLEPC_DIR=/home/hpc/pr27fe/lu65nut2/tmp/slepc-3.9.1       

poe ./gene_supermuc
### to submit a parameter scan, comment the previous line and uncomment
### the following (see GENE documentation or ./scanscript --help)
#./scanscript --np 8 --ppn 8 --mps 4
