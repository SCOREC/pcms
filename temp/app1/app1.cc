#include<iostream>
#include<vector>

#include<adios2.h>
//#ifdef ADIOS2_HAVE_MPI
#include<mpi.h>
//#endif


int main(int argc, char *argv[])
{
   int rank = 0, nproc = 1;

//#ifdef ADIOS_HAVE_MPI
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &nproc);
//#endif

   const int NSTEP = 5;

//#ifdef ADIOS2_HAVE_MPI
   adios2::ADIOS adios(MPI_COMM_WORLD);
//#else
  // adios2::ADIOS adios;
//#endif

   std::cout << "Hello World!" << std::endl;

//#ifdef ADIOS2_HAVE_MPI
   MPI_Finalize();
//#endif
   return 0;
}
