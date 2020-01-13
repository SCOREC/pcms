#include<iostream>
#include<vector>

#include<adios2.h>
#include<mpi.h>

void send(int rank, int size);

int main(int argc, char *argv[])
{
  int rank , size;


  return 0;
}


void send(int rank, int size)
{
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  std::vector<float> myFloats = {
    (float)10.0 * rank + 0, (float)10.0 * rank + 1, (float)10.0 * rank + 2,
    (float)10.0 * rank + 3, (float)10.0 * rank + 4, (float)10.0 * rank + 5,
    (float)10.0 * rank + 6, (float)10.0 * rank + 7, (float)10.0 * rank + 8,
    (float)10.0 * rank + 9};
  const std::size_t Nx = myFloats.size();

  try
  {
    adios2::ADIOS adios(MPI_COMM_WORLD, adios2::DebugON);
    adios2::IO sstIO = adios.DeclareIO("myIO");
    sstIO.SetEngine("Sst");

    auto bpFloats = sstIO.DefineVariable<float>("bpFloats", {size * Nx}, {rank * Nx}, {Nx});

    adios2::Engine sstWriter = sstIO.Open("helloSst", adios2::Mode::Write);

    sstWriter.BeginStep();
    sstWriter.Put<float>(bpFloats, myFloats.data());
    sstWriter.EndStep();
    sstWriter.Close();
  }

  catch(std::invalid_argument &e)
  {
    std::cout << "Invalid argumant exception, STOPPING PROGRAM from rank "
      << rank << "\n";
    std::cout << e.what() << "\n";
  }
  catch(std::ios_base::failure &e)
  {
    std::cout << "IO system base failure exception, STOPPING PROGRAM from rank "
      << rank << "\n";
    std::cout << e.what() << "\n";
  }
  catch(std::exception &e)
  {
    std::cout << "Exception, STOPPING PROGRAM from rank " << rank << "\n";
    std::cout << e.what() << "\n";
  }


}
