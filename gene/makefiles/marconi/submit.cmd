#!/bin/bash
## Submit via: sbatch submit.cmd (parameters below can be overwritten by command line options)
#SBATCH -J GENE              # job name
#SBATCH --time=24:00:00      # walltime
#SBATCH --nodes=1            # Total number of nodes
#SBATCH --ntasks-per-node=48 # MPI tasks (A3: max 48, A2: max 68)
#SBATCH --cpus-per-task=1    # Number of threads per task (OpenMP)
##SBATCH --mem=86000MB       # uncomment for recommended A2 (KNL) memory limit
##SBATCH --constraint=cache  #flat/cache
#SBATCH -e %x.%j.err         # std-error file
#SBATCH -o %x.%j.out         # std-output file
#SBATCH -p skl_fua_prod      # partition: see sinfo -d (skl_fua_pro, skl_fua_dbg, knl_fua_prod)
#SBATCH -q normal	     # QoS: see sacctmgr show -P qos (skl_qos_fuabprod, skl_qos_fualowprio, knl_qos_fuabprod etc)
#SBATCH -A FUSIO_ALL         # account: see "saldo -b" (A1), "saldo -b --knl" (A2), or "saldo -b --skl" (A3)

#set environment
module load intel intelmpi mkl fftw
module load blas lapack scalapack
export OMP_NUM_THREADS=1

# Launch the parallel job on the allocated compute nodes
srun -l -K ./gene_marconi
# If you use less than a node: srun --cpu-bind=cores -l K ./gene_marconi

### to submit a parameter scan, comment the previous line and uncomment
### the following (see GENE documentation or ./scanscript --help)
#./scanscript --np $SLURM_NTASKS --ppn $SLURM_NTASKS_PER_NODE --mps 4  --syscall="srun -l -K ./gene_marconi"

#
