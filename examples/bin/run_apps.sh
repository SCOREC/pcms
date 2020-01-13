export APP1_EXEC=/lore/adesoa/dev/wdmapp_coupling/examples/bin/helloSstWriter
export APP2_EXEC=/lore/adesoa/dev/wdmapp_coupling/examples/bin/helloSstReader

sleep 5

OUTFILE=app1.out
ERRFILE=app1.err


mpirun -np 4 ${APP1_EXEC} > ${OUTFILE} 2> ${ERRFILE} &


sleep 5

OUTFILE=app2.out
ERRFILE=app2.err


mpirun -np 4 ${APP2_EXEC} > ${OUTFILE} 2> ${ERRFILE} &

wait

echo $(date)
