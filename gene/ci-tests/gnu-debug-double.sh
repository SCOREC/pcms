#!/bin/bash -l
#wrapper file for easy gnu-debug-double invocation
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
$DIR/compile-test.sh gnu debug double
