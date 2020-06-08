#include "adios2Routines.h"

void exParFor() {
  Kokkos::parallel_for(
      4, KOKKOS_LAMBDA(const int i) {
        printf("Hello from kokkos thread i = %i\n", i);
      });
}

int main(int argc, char **argv){
  int rank, nprocs;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  Kokkos::initialize(argc, argv);
  if(!rank) {
    printf("Hello World on Kokkos execution space %s\n",
         typeid(Kokkos::DefaultExecutionSpace).name());
    exParFor();
  }

  adios2::ADIOS adios(MPI_COMM_WORLD, adios2::DebugON);
  adios2::Variable<double> send_var[2];
  const std::string dir = "../coupling";
  const int time_step = 1, RK_count = 4;

  coupler::adios2_handler gDens(adios,"gene_density");
  coupler::adios2_handler cDens(adios,"cpl_density");
  coupler::adios2_handler xFld(adios,"xgc_field");
  coupler::adios2_handler cFld(adios,"cpl_field");

  for (int i = 0; i < time_step; i++) {
    for (int j = 0; j < RK_count; j++) {
      coupler::Array2d<double>* density = coupler::receive_density(dir, gDens,MPI_COMM_WORLD,MPI_COMM_WORLD);
      coupler::printSomeDensityVals(density);
      coupler::send_density(dir, density, cDens, send_var[0]);
      coupler::destroy(density);

      coupler::Array2d<double>* field = coupler::receive_field(dir, xFld,MPI_COMM_WORLD,MPI_COMM_WORLD);
      coupler::send_field(dir, field, cFld, send_var[1]);
      coupler::destroy(field);
    }
  }

  Kokkos::finalize();
  MPI_Finalize();
  return 0;
}

