#!/bin/bash -l
#wrapper file for easy cray-minimal-double invocation
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
$DIR/compile-test.sh cray minimal double 4
