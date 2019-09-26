#!/bin/bash -l
#wrapper file for easy pgi-minimal-double invocation
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
$DIR/compile-test.sh pgi minimal double 4
