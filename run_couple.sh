#!/bin/bash -l
#SBATCH -A m499
#SBATCH -p regular
#SBATCH -N 32
#SBATCH -C knl,quad,cache
#SBATCH --time-min=4:00:00
#SBATCH --time=4:00:00
#SBATCH --job-name=XGC_GENE_D3D
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=adesojidiekola@gmail.com
#SBATCH -o xgc_gene.%j.out
#SBATCH -e xgc_gene.%j.err

echo $(date)
## node count 
export n_nodes_XGC=24
export n_nodes_GENE=8

## XGC set up
export n_mpi_ranks_per_node_XGC=16
export n_mpi_ranks_XGC=$((${n_mpi_ranks_per_node_XGC} * ${n_nodes_XGC}))

export OMP_NUM_THREADS=8
export OMP_PLACES=threads
export OMP_PROC_BIND=spread
export OMP_MAX_ACTIVE_LEVELS=2
export OMP_NESTED=TRUE
export OMP_STACKSIZE=2G
export SWDIR=/project/projectdirs/m499/Software
export XGC_EXEC=~/wdmapp_coupling/XGC_adios2/xgc-es
export GENE_EXEC=~/wdmapp_coupling/gene_Adios2_cori/bin/gene_cori

echo ''
echo 'XGC run:'
echo 'Number of nodes: '                  $n_nodes_XGC
echo 'MPI ranks (total): '                $n_mpi_ranks_XGC
echo 'MPI ranks per node: '               $n_mpi_ranks_per_node_XGC
echo 'Number of OMP threads: '            ${OMP_NUM_THREADS}
echo 'XGC executable: '                   ${XGC_EXEC}

rm coupling/*

cd XGC/

mkdir -p restart_dir

sleep 20

OUTFILE=xgc1_${SLURM_JOB_ID}.out
ERRFILE=xgc1_${SLURM_JOB_ID}.err

srun --exclusive -N ${n_nodes_XGC} -n ${n_mpi_ranks_XGC} -c ${OMP_NUM_THREADS} --cpu_bind=cores ${XGC_EXEC} > ${OUTFILE} 2>${ERRFILE} &

## GENE set up
export HDF5_USE_FILE_LOCKING=FALSE
export n_mpi_ranks_per_node_GENE=16
export n_mpi_ranks_GENE=$((${n_mpi_ranks_per_node_GENE} * ${n_nodes_GENE}))
export OMP_NUM_THREADS=1

echo ' '
echo 'GENE run:'
echo 'Number of nodes: '                  $n_nodes_GENE
echo 'MPI ranks (total): '                $n_mpi_ranks_GENE
echo 'MPI ranks per node: '               $n_mpi_ranks_per_node_GENE
echo 'Number of OMP threads: '            ${OMP_NUM_THREADS}

sleep 15

cd ../GENE

mkdir -p out

OUTFILE=GENE_${SLURM_JOB_ID}.out
ERRFILE=GENE_${SLURM_JOB_ID}.err

srun --exclusive -N ${n_nodes_GENE} -n ${n_mpi_ranks_GENE} -l -K --cpu_bind=cores ${GENE_EXEC} > ${OUTFILE} 2>${ERRFILE} &

wait

echo $(date)
