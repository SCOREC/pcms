#!/bin/bash

name=$1
test=$2

cd $SCRATCH/globus-compute/$name-test
source env.sh

mkdir -p results
cd build-$name
salloc --time 00:20:00 --constraint gpu --qos=interactive --nodes=1 --ntasks-per-node=40 --cpus-per-task=1 --gpus=1 --account=m4564 ctest $test 2>&1 | tee ../results/summary.log
cat $PWD/Testing/Temporary/LastTest.log 2>&1 | tee ../results/test.log