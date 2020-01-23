#include"cpl_init.h"
#include<adios2.h>
#include <mpi.h>
// TASK_1: Execute both codes using the system commands 
// TASK_2: validate OOP implementation is necessary i.e. considering process level coordination
// TASK_3: Compare both codes


int main(int argc, char **argv){
	MPI_Init(NULL, NULL);

// This should contain the ADIOS2 calls to get data from XGC/GENE and transfer it accordingly
	return 0;

 MPI_Finalize();
}
