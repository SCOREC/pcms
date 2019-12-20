#!/bin/tcsh -l
### the job's name
#$ -N GENE
### start the job in the current working dir
#$ -cwd
### start a parallel environment with multiples of 16 processes
#$ -pe edge*_ompi_b 16
#$ -R y
### wallclock limit
#$ -l h_rt=24:00:00
### mail notification (n=never)
#$ -m n
##$ -M $USER@ipp.mpg.de
#------------------------------------------

#load environment
module load petsc/3.5-intel13-ompi16-cplx mkl/11.0

### set MKL to serial mode 
setenv MKL_SERIAL yes

### no OpenMP
setenv OMP_NUM_THREADS 1

### start program
mpirun -np $NSLOTS ./gene_hgw_dell_cluster

### to submit a parameter scan, comment the previous line and uncomment
### the following (see GENE documentation or ./scanscript --help)
#./scanscript --np $NSLOTS --ppn 16 --mps 4

