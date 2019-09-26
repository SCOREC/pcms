#!/bin/bash -l
#IBM LoadLeveler version for submitting testsuite jobs
#submit-testsuite-sge.sh <compiler> <testsuite set> <number of MPI processes>
compiler=$1
testset=$2
nmpi=$3
JOBNAME=GENERUNNER
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

use_all_procs=0
if [ "$nmpi" == "0" ]; then
    use_all_procs=1
fi

if [ $# -gt 3 ]; then
    runtime=$4
else
    runtime="00:30:00"
fi
if [ "$testset" == "standard" ]; then
    runtime="00:15:00"
fi
if [ "$testset" == "standard_bsg" ]; then
    runtime="00:15:00"
fi
echo "$compiler, $testset, $nmpi, $runtime"

#load $COMPILER environment from MACHINE makefile -- removing module list for security reason
MACHINE=$(make mach_wrapper 2>/dev/null)
module purge
MODULELIST=$(make COMPILER="$compiler" PRECISION=double get_required_modules | tail -n 1 | sed -r "s/module load //")
eval module load $MODULELIST git

module list

cd $DIR/../testsuite
executable=./gene_$MACHINE
case "$testset" in
    slepc|neoclassic)
	slepc_linked="$(ldd $executable 2>/dev/null | grep slepc)"
	if [ -z "$slepc_linked" ]; then
	    echo "NO EXECUTABLE WITH SLEPC LINKED IN FOUND -- DIFFERENT RUNNER THAN AT BUILD STAGE"
	    cd ..
	    echo "Sorry ... need to rebuild ... choosing $compiler-slepc-double"
	    $DIR/$compiler-slepc-double.sh
	    cd testsuite
	fi
	;;
    *)
	if [ ! -f "$executable" ]; then
	    echo "NO EXECUTABLE FOUND -- DIFFERENT RUNNER THAN AT BUILD STAGE"
	    cd ..
	    echo "Sorry ... need to rebuild ... choosing $compiler-minimal-double"
	    $DIR/$compiler-minimal-double.sh
	    cd testsuite
	fi
	;;
esac

# ------------ BATCH SYSTEM DEPENDENT PART ---------------
rm -f submit.cmd
make -f ../makefile checkpath
#echo "====================== submit.cmd ====="
#cat submit.cmd
#echo "====================== submit.cmd ====="

#LoadLeveler specific part comes here
# check if runtime is longer than 30 minutes
hour=$(sed -n "s/^\([0-9]*\):[0-9]*:[0-9]*/\1/p" <<< "$runtime")
minute=$(sed -n "s/^[0-9]*:\([0-9]*\):[0-9]*/\1/p" <<< "$runtime")
seconds=$(echo "$runtime" | awk -F ":" '{ print $3 }')
if [ -z "$seconds" ]; then
seconds="00"
fi

echo "hour=$hour, minute=$minute, seconds=$seconds"
if [ "$hour" -ne "00" ]; then
echo ""
else
echo ""
fi
sed -i 's/# @ *tasks_per_node *=.*/# @ tasks_per_node = '"$nmpi"'/' submit.cmd
sed -i 's/# @ *wall_clock_limit *=.*/# @ wall_clock_limit = '"$hour"':'"$minute"':'"$seconds"'/' submit.cmd
sed -i 's/# @ *job_name *=.*/# @ job_name = '"$JOBNAME"'/' submit.cmd
sed -i 's|module load .*|module purge \&\& module load git '"$MODULELIST"'|' submit.cmd
case "$testset" in
    slepc)
	restrict_tests="-s 1 -e 9"
	;;
    big)
	restrict_tests="-s 1 -e 6"
	;;
    *)
	restrict_tests=""
	;;
esac

if [ "${use_all_procs}" == "1" ]; then
    sed -i 's|^ *poe .*|./testsuite -t '"$testset"' '"$restrict_tests"' -syscall "poe '"${executable}"'" -o /ptmp/$USER/|' submit.cmd
else
    sed -i 's|^ *poe .*|./testsuite -t '"$testset"' '"$restrict_tests"' -syscall "poe '"${executable}"'" -o /ptmp/$USER/|' submit.cmd
fi

echo "====================== submit.cmd after modification ====="
cat submit.cmd
echo "====================== submit.cmd after modification ====="
echo ""

submitback=$(llsubmit submit.cmd)
jobid=$(echo "$submitback" | awk -F '"' '{ print $2 }')
if [ -z "$jobid" ]; then
    echo "ERROR: failed to submit job or to retrieve job id"
    exit 1
fi
jobnr=$(echo "$jobid" | awk -F "." '{ print $2 }')
echo "Submitted job ... received jobid $jobid"

# wait for some seconds to give the scheduler a chance to update the database
sleep 5s

# now wait for completion -- loop until job in jobid finishes
echo "Waiting for jobid $jobid to complete ..."
prev_status=""

while true; do
    status=$(llq -j $jobid | grep $jobid | awk '{ print $5 }')
    case $status in
	R*)
	    if [ "$prev_status" != "$status" ]; then
		prev_status=$status
		echo "$JOBNAME RUNNING"
	    fi
	    sleep 5s
	    ;;
	I*)
	    if [ "$prev_status" != "$status" ]; then
		prev_status=$status
		echo "$JOBNAME QUEUEING"
	    fi
	    sleep 10s
	    ;;
	NQ*)
	    cat $JOBNAME*$jobnr*
	    if [ "$prev_status" != "$status" ]; then
		prev_status=$status
		echo "$JOBNAME NOT QUEUED"
	    fi
	    ;;
	*)
	    cat $JOBNAME*$jobnr*
	    error_found=$(cat $JOBNAME*$jobnr* | grep 'error encountered')
	    if [ -z "$error_found" ]; then
		exit 0
	    else
		exit 1
	    fi
	    ;;
    esac
done
