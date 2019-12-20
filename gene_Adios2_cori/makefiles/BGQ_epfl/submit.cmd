#!/bin/sh

# @ job_name = GENE
# @ comment = "GENE run"
# @ job_type = bluegene
# @ error = $(pwd)$(job_name)_$(jobid).err
# @ output = $(job_name)_$(jobid).out
# @ environment = COPY_ALL;
# @ wall_clock_limit = 12:00:00
# @ notification = never
# @ class = test
# @ bg_size = 64
# @ queue


#here we state the number of processors
NP=1024

#one can run only from scratch, therefore set here diagdir
DIAGDIR=./

EXE=$(pwd)"/../bin/gene_BGQ_epfl"

mkdir -p ${DIAGDIR}
cp parameters ${DIAGDIR}

cd ${DIAGDIR}

ln -s ${EXE} .

###uncomment for single run
runjob --np ${NP} --ranks-per-node 16 : ${EXE}
