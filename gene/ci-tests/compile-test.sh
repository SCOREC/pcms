#!/bin/bash -l
#compile GENE with chosen compiler/precision for given setup (debug,miminal,slepc,all)
#compile-test.sh <compiler> <setup> <precision> <# jobs (optional)>

#sanity checks
if [ $# -lt 3 ]; then
 echo "compile-test.sh: >=3 inputs required"
 exit 1
else
 compiler=$1
 setup=$2
 precision=$3
 if [ $# -gt 3 ]; then
  jobs=$4
 else
  jobs=1
 fi
fi
echo "$compiler $setup $precision $jobs"

valid_compiler="intel gnu pgi cray"
! [[ $valid_compiler =~ (^|[[:space:]])"$compiler"($|[[:space:]]) ]] && \
    echo "compile-test.sh: invalid compiler $compiler - valid options: $valid_compiler" && exit 1
valid_setup="debug minimal slepc all"
! [[ $valid_setup =~ (^|[[:space:]])"$setup"($|[[:space:]]) ]] && \
    echo "compile-test.sh: invalid setup $setup - valid options: $valid_setup" && exit 1
valid_precision="double single"
! [[ $valid_precision =~ (^|[[:space:]])"$precision"($|[[:space:]]) ]] && \
    echo "compile-test.sh: invalid precision $precision - valid options: $valid_precision" && exit 1

#start setting up local clone for compilation
make distclean
module purge

MACHINE=$(make mach_wrapper 2>/dev/null)
if [ "$MACHINE" == "new_machine" ]; then
    export MACHINE=buildtest
fi

#get module list from MACHINE makefile -- removing module list for security reason
MODULELIST=$(make COMPILER=$compiler PRECISION=$precision get_required_modules | tail -n 1 | sed -r "s/module load //")
eval module load $MODULELIST git
module list

#make -j COMPILER=$compiler PRECISION=$precision SLEPC=no FUTILS=no SCALAPACK=yes PRODRUN=no DEBUG=yes
#use sed instead for interfacing with testsuite which may read makefile entries
sed -i "s/COMPILER *= *.*/COMPILER=$compiler/" ./bin/${MACHINE}.mk
sed -i "s/PRECISION *= *.*/PRECISION=$precision/" ./bin/${MACHINE}.mk
if [ "$setup" = "debug" ]; then
    sed -i 's/DEBUG *= *.*/DEBUG=yes/' ./bin/${MACHINE}.mk
fi
case "$setup" in
"slepc")
	sed -i "s/PRODRUN *= *.*/PRODRUN=no/" ./bin/${MACHINE}.mk
	sed -i 's/FUTILS *= *.*/FUTILS = no/' ./bin/${MACHINE}.mk
	;;
"all")
	#do nothing (keep default makefile settings)
	;;
*)
	sed -i "s/PRODRUN *= *.*/PRODRUN=no/" ./bin/${MACHINE}.mk
	sed -i 's/SLEPC *= *.*/SLEPC = no/' ./bin/${MACHINE}.mk
	sed -i 's/FUTILS *= *.*/FUTILS = no/' ./bin/${MACHINE}.mk
	;;
esac

#start compilation
make -j $jobs
