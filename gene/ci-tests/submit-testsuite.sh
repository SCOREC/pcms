#!/bin/bash -l
compiler=$1
testset=$2
nmpi=$3
time=$4
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
if [ -x "$(command -v sbatch)" ]; then
$DIR/submit-testsuite-slurm.sh $compiler $testset $nmpi $time
elif [ -x "$(command -v qsub)" ]; then
$DIR/submit-testsuite-sge.sh $compiler $testset $nmpi $time
elif [ -x "$(command -v llsubmit)" ]; then
$DIR/submit-testsuite-ll.sh $compiler $testset $nmpi $time
else
echo "couldn't detect batch system"
exit 1
fi