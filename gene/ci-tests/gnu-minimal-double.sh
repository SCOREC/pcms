#!/bin/bash -l
#wrapper file for easy gnu-minimal-double invocation
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
$DIR/compile-test.sh gnu minimal double 4
