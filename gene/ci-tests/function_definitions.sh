look_for_module() {
    basename=$1
    major_range=($2)
    len_major=${#major_range[@]}
    echo "length of major_range array is ${len_major}"
    minor_range=($3)
    len_minor=${#minor_range[@]}
    echo "length of minor_range array is ${len_minor}"
    for imajor in $(seq 0 $(( len_major - 1 )) ); do
	major=${major_range[imajor]}
	# test for different minor versions
	for minor in ${minor_range[@]}; do
	    version="$basename/${major}.${minor}"
	    echo -n "Testing $version ... "
	    res=$(module av $version 2>&1)
	    if [ ! "x$res" == "x" ]; then
		echo "success"
		echo "Loading $version"
		module load $version
		break 2
	    else
		echo "failed"
	    fi
	done
	# test for module with only major version
	version="$basename/${major}"
	echo -n "Testing $version ... "
	res=$(module av $version 2>&1)
	if [ ! "x$res" == "x" ]; then
	    echo "success"
	    echo "Loading $version"
	    module load $version
	    break 1
	else
	    echo "failed"
	fi
    done
}

load_fftw_module() {
    compiler=$1
    if [ "$MACHINE" == "draco" ]; then
	case $compiler in
	    gcc)
		module load fftw/gcc/3.3.7
		;;
	    intel)
		module load fftw/3.3.7
		;;
	esac
    elif [ "$HOSTNAME" == "host-130-183-198-74" ]; then
	module load fftw-serial
    fi
}

load_petsc_slepc_modules() {
    compiler=$1
    if [ "$MACHINE" == "draco" ]; then
	case $compiler in
	    gcc)
		module load petsc-cplx slepc-cplx
		;;
	    intel)
		module load petsc-cplx slepc-cplx
		;;
	esac
    elif [ "$HOSTNAME" == "host-130-183-198-74" ]; then
	# does not exist
	echo "slepc and petsc are not yet installed"
    fi
}
