#!/bin/bash -l

export MACHINE=ubuntu_gfortran
./newprob
cd prob01
touch make_local
make -f ../makefile checkpath

sed -i "s/COMPILER *= *.*/COMPILER=gnu/" ubuntu_gfortran.mk
sed -i "s/PRECISION *= *.*/PRECISION=double/" ubuntu_gfortran.mk
sed -i 's/DEBUG *= *.*/DEBUG=yes/' ubuntu_gfortran.mk
make -f ../makefile
