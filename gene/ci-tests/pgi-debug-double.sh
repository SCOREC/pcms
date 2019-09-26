#!/bin/bash -l
#wrapper file for easy pgi-debug-double invocation
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
$DIR/compile-test.sh pgi debug double 4
