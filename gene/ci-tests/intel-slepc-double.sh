#!/bin/bash -l
#wrapper file for easy intel-slepc-double invocation
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
$DIR/compile-test.sh intel slepc double 4
