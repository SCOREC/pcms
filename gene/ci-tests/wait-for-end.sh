#!/bin/bash -l

source ci-tests/function_definitions.sh

# we need a recent intel compiler
#look_for_module "intel" "18 17 16 15" "9 8 7 6 5 4 3 2 1 0"
look_for_module "impi" "2017 5.1" "4 3 2 1"
#look_for_module "mkl" "2017 12 11 10" "9 8 7 6 5 4 3 2 1 0"

module list

cd testsuite

#slurm specific part comes here
# loop until job in jobid finishes
jobid=$(cat jobid)
cat jobid
echo "jobid that is read in: $jobid"

while true; do
    status=$(sacct -j $jobid -o JobName,State |grep GENE | awk '{ print $2 }')
    echo "status=$status"
    case $status in
	COMPLETED*)
	    exitcode=$(sacct -j $jobid -o JobName,ExitCode |grep GENE | awk '{ print $2 }' | sed 's/\([0-9]*\):[0-9]*/\1/')
	    exit $exitcode
	    ;;
	RUNNING*)
	    sleep 5s
	    ;;
	PENDING*)
	    sleep 10s
	    ;;
	*)
	    echo "We have a status of $status, failing."
	    exit 1
	    ;;
    esac
done
