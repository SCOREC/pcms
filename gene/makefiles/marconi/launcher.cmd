#!/bin/bash
#
## MAXWT 24
## PROCSPERNODE 48
## SUBMITCMD sbatch
#SBATCH -J JOBNAME
#SBATCH --time=WALLTIME
#SBATCH --nodes=NODES
#SBATCH --ntasks-per-node=PPN
#SBATCH -e %x.%j.err            # std-error file
#SBATCH -o %x.%j.out            # std-output file
#SBATCH -p skl_fua_prod # see sinfo (skl_fua_pro, skl_fua_dbg, knl_fua_prod)
#SBATCH -q normal
#SBATCH -A FUSIO_ALL      # see "saldo -b" (A1), "saldo -b --knl" (A2), or "saldo -b --skl" (A3)

module load intel intelmpi mkl fftw
module load blas lapack scalapack
export OMP_NUM_THREADS=1

srun -l -K ./gene_marconi
#
