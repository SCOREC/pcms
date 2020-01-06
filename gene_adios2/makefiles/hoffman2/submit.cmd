#!/bin/csh -f
#$ -cwd
#$ -j y
#$ -pe shared 8
#$ -l h_data=4096M,h_rt=24:00:00
#
#
#  Name of application for log
#$ -v QQAPP=intelmpi
#  Email address to notify
##$ -M mark_mustermann@physics.ucla.edu
#  Notify at beginning and end of job
#$ -m bea
#  Job is not rerunable
#$ -r n
# Initialization for mpi parallel execution
#
  unalias *
  set qqversion =
  set qqapp     = "intelmpi parallel"
  set qqptasks  = 8
  set qqjob     = testsuite
  cd $SGE_O_WORKDIR
  source /u/local/bin/qq.sge/qr.runtime
  if ($status != 0) exit (1)
#
  echo "UGE job for testsuite built Fri Feb  6 14:37:23 PST 2015"
  echo ""
  echo "  testsuite directory:"
  echo "  Submitted to UGE:"
  echo "    "$qqsubmit
  echo "  'scratch' directory (on each node):"
  echo "    $qqscratch"
  echo "  8-way parallel job configuration:"
  echo "    $qqconfig" | tr "\\" "\n"
#
#
# Run the user program
#

  source /u/local/Modules/default/init/modules.csh
  module load intel/13.cs

  setenv OMP_NUM_THREADS 1

# setting I_MPI_HYDRA_BOOTSTRAP=sge ensures that mpiexec.hydra works in sge tight integration
#  not needed if using mpirun
#  setenv I_MPI_HYDRA_BOOTSTRAP sge

#  time mpirun -np 8 -env I_MPI_FABRICS ofa:ofa  \
#  time mpiexec.hydra -np 8 -env I_MPI_FABRICS ofa:ofa  \
  time mpirun -np 8 ./gene_hoffman2  >& gene.output.$JOB_ID


  echo ""
  echo "GENE finished at:  "` date `
