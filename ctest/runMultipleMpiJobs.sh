#!/bin/bash 
usage="Usage: <run command> <process flag> <valgrindCmd> <valgrindOptions> <name0> <procs0> <exe0> <args0> ... <name2> <procs2> <exe2> <args2>"
if [[ ($# == 0) || ("${1}" == "-h" || ${1} == "--help") ]]; then
  echo $usage
  exit 1
fi
runCmd=${1}
numProcsFlag=${2}
valgrindCmd=${3}
valgrindOptions=${4}
#clear the valgrind args if they are 'none'
if [[ ${valgrindCmd} == "none" ]]; then
  valgrindCmd=""
fi
if [[ ${valgrindOptions} == "none" ]]; then
  valgrindOptions=""
fi
shift 4

declare -a PIDS
declare -a LOGS
run() {
  local name=${1}
  local procs=${2}
  local exe=${3}
  local args=${4}
  IFS=';' read -a argsArray <<< "${args}" #split the cmake list of args
  ${runCmd} ${numProcsFlag} ${procs} ${valgrindCmd} ${valgrindOptions} ${exe} ${argsArray[@]} &> ${name}.log &
  PIDS+=($!)
  LOGS+=(${name}.log)
}

run $@
shift 4
run $@
shift 4
[[ $# == 4 ]] && run $@
for i in "${!PIDS[@]}"; do
  wait ${PIDS[$i]}
  status=$?
  if [[ $status -ne 0 ]]; then
    cat ${LOGS[$i]}
    exit $status
  fi
done
