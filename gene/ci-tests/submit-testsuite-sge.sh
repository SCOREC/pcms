#!/bin/bash -l
#SGE version for submitting testsuite jobs
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

if [ $nmpi -gt 32 ]; then
    nodes=$(( $nmpi/32 ))
    use_all_procs=1
else
    nodes=1
fi

if [ $# -gt 3 ]; then
    runtime=$4
else
    runtime="00:30:00"
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

#SGE specific part comes here
# check if runtime is longer than 30 minutes
hour=$(sed -n "s/^\([0-9]*\):[0-9]*:[0-9]*/\1/p" <<< "$runtime")
minute=$(sed -n "s/^[0-9]*:\([0-9]*\):[0-9]*/\1/p" <<< "$runtime")
#seconds=$(sed "s/^[0-9]*:[0-9]*:\([0-9]*\)/\1/" <<< "$runtime")

echo "hour=$hour, minute=$minute"
if [ "$hour" -ne "00" ]; then
echo ""
else
sed -i 's,^\(\#$ *-P *=\).*,\1'"debug"',' submit.cmd
sed -i 's/impi_hydra\.\*/impi_hydra_devel/g' submit.cmd
fi
sed -i 's,^\(\#$ *-pe *[^\s]* \).*,\1'"$nmpi"',' submit.cmd
sed -i 's,^\(\#$ *-l h_rt *=\).*,\1'"$runtime"',' submit.cmd
sed -i 's/#$ *-N .*/#$ -N '"$JOBNAME"'/' submit.cmd

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
    sed -i 's|^ *mpiexec .*|./testsuite -t '"$testset"' '"$restrict_tests"' -syscall "mpiexec -n '"$nmpi"' '"${executable}"'"|' submit.cmd
else
    sed -i 's|^ *mpiexec .*|./testsuite -t '"$testset"' '"$restrict_tests"' -syscall "mpiexec -n '"$nmpi"' '"${executable}"'"|' submit.cmd
fi

echo "====================== submit.cmd after modification ====="
cat submit.cmd
echo "====================== submit.cmd after modification ====="
echo ""

submitback=$(qsub submit.cmd)
jobid=$(echo "$submitback" | awk '{ print $3 }')
if [ -z "$jobid" ]; then
    echo "ERROR: failed to submit job or to retrieve job id"
    exit 1
fi
echo "Submitted job ... received jobid $jobid"

# wait for some seconds to give the scheduler a chance to update the database
sleep 5s

# now wait for completion -- loop until job in jobid finishes
echo "Waiting for jobid $jobid to complete ..."
prev_status=""

while true; do
    status=$(qstat | grep $jobid | awk '{ print $5 }')
    case $status in
	r*)
	    if [ "$prev_status" != "$status" ]; then
		prev_status=$status
		echo "$JOBNAME RUNNING"
	    fi
	    sleep 5s
	    ;;
	qw*)
	    if [ "$prev_status" != "$status" ]; then
		prev_status=$status
		echo "$JOBNAME QUEUEING"
	    fi
	    sleep 10s
	    ;;
	Eqw*)
	    cat $JOBNAME*$jobid*
	    echo "$JOBNAME IN ERROR STATE"
	    qstat -j $jobid | grep 'error reason'
	    exit 1
	    ;;
	*)
	    cat $JOBNAME*$jobid*
	    echo "waiting for batch system to collect job data ..."
	    sleep 10s
	    exitcode=$(qacct -j $jobid | grep exit_status | awk '{ print $2 }')
	    exit $exitcode
	    ;;
    esac
done
