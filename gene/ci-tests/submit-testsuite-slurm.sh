#!/bin/bash -l
#SLURM version for submitting testsuite jobs
#submit-testsuite-slurm.sh <compiler> <testsuite set> <number of MPI processes>
compiler=$1
testset=$2
nmpi=$3
JOBNAME=GENERUNNER
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

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

#slurm specific part comes here
# check if runtime is longer than 30 minutes
hour=$(sed -n "s/^\([0-9]*\):[0-9]*:[0-9]*/\1/p" <<< "$runtime")
minute=$(sed -n "s/^[0-9]*:\([0-9]*\):[0-9]*/\1/p" <<< "$runtime")
#seconds=$(sed "s/^[0-9]*:[0-9]*:\([0-9]*\)/\1/" <<< "$runtime")
echo "hour=$hour, minute=$minute"

partition=general
if [ "$hour" -ne "00" ]; then
echo ""
elif [ $minute -gt 30 ]; then
echo ""
else
    sinforeturn=$(sinfo -h -p express)
    if [ -z "$sinforeturn" ]; then
	partition=debug
    else
	partition=express
    fi
    sed -i 's,^\(\#SBATCH *--partition *=\).*,\1'"$partition"',' submit.cmd
    sed -i 's,^\(\#SBATCH *-p *\).*,\1'"$partition"',' submit.cmd
fi

#now read default tasks per node and compute required node number
ntasks_per_node=$(sed -n -e 's/^\#SBATCH *--ntasks-per-node *=* *\([0-9]*\) *.*$/\1/p' submit.cmd)
nodes=$(( ($nmpi+$ntasks_per_node-1)/$ntasks_per_node ))

#evenly distribute tasks on nodes
ntasks_per_node=$(( $nmpi/$nodes ))

sed -i 's,^\(\#SBATCH *--ntasks-per-node *=*\).*,\1'"$ntasks_per_node"',' submit.cmd
sed -i 's,^\(\#SBATCH *--nodes *=\).*,\1'"$nodes"',' submit.cmd
sed -i 's,^\(\#SBATCH *-n *\).*,\1'"$nmpi"',' submit.cmd
sed -i 's,^\(\#SBATCH *--time *=\).*,\1'"$runtime"',' submit.cmd
sed -i 's,^\(\#SBATCH *-t *\).*,\1'"$runtime"',' submit.cmd
sed -i 's,^\(\#SBATCH *-J *\).*,\1'"$JOBNAME"',' submit.cmd
sed -i 's/#SBATCH *-o .*/#SBATCH -o '"$JOBNAME"'.%j/' submit.cmd

MY_BATCH_ACCOUNT_STR=GITLAB_BATCH_ACCOUNT_${MACHINE^^}
MY_BATCH_ACCOUNT=${!MY_BATCH_ACCOUNT_STR}
if [ ! -z "$MY_BATCH_ACCOUNT" ]; then
    echo "Using $MY_BATCH_ACCOUNT_STR = $MY_BATCH_ACCOUNT"
    sed -i 's/#\+SBATCH *-A .*/#SBATCH -A '"$MY_BATCH_ACCOUNT"'/' submit.cmd
fi

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

sed -i 's|^ *srun .*|./testsuite -t '"$testset"' '"$restrict_tests"' -syscall "srun -n '"$nmpi"' -K '"${executable}"'"|' submit.cmd

echo "====================== submit.cmd after modification ====="
cat submit.cmd
echo "====================== submit.cmd after modification ====="
echo ""

submitback=$(sbatch submit.cmd)
jobid=$(echo "$submitback" |sed "s/Submitted batch job \([0-9]*\)/\1/")
if [ -z "$jobid" ]; then
    echo "ERROR: failed to submit job or to retrieve job id"
    exit 1
fi
echo "Submitted job ... received jobid $jobid"

# wait for some seconds to give the scheduler a chance to update the database
sleep 5s

# now wait for completion -- loop until job in jobid finishes

echo "Waiting for jobid $jobid to complete ..."
queue_start=$(sacct -j $jobid -o JobName,State | head -n 3)
echo "$queue_start"
prev_status=$(echo "$queue_start" | grep $JOBNAME | awk '{ print $2 }')

while true; do
    status=$(sacct -j $jobid -o JobName,State | grep $JOBNAME | awk '{ print $2 }')
    if [ -z "$status" ]; then
	sleep 5s
	continue
    fi
    case $status in
	COMPLETED*)
	    cat $JOBNAME*$jobid*
	    exitcode=$(sacct -j $jobid -o JobName,ExitCode | grep $JOBNAME | awk '{ print $2 }' | sed 's/\([0-9]*\):[0-9]*/\1/')
	    exit $exitcode
	    ;;
	RUNNING*)
	    if [ "$prev_status" != "$status" ]; then
		prev_status=$status
		sacct -j $jobid -o JobName,State | grep $JOBNAME
	    fi
	    sleep 5s
	    ;;
	PENDING*)
	    sleep 10s
	    ;;
	*)
	    echo "We have a status of $status, failing."
	    cat $JOBNAME*$jobid*
	    exit 1
	    ;;
    esac
done
