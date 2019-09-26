#!/bin/bash -l
#wrapper file for easy intel-minimal-double invocation
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
$DIR/compile-test.sh intel minimal double 4
