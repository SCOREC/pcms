#include "coupling.h"

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
  adios2::IO IO[6];
  adios2::Engine eng[6];
  adios2::Variable<double> send_var[2];
  const std::string dir = "../coupling";
  const int time_step = 1, RK_count = 4;

  IO[0] = adios.DeclareIO("gene_pproc");
  IO[1] = adios.DeclareIO("xgc_pproc");
  IO[2] = adios.DeclareIO("gene_density");
  IO[3] = adios.DeclareIO("cpl_density");
  IO[4] = adios.DeclareIO("xgc_field");
  IO[5] = adios.DeclareIO("cpl_field");

  //receive GENE's preproc mesh discretization values
  coupler::Array1d<int>* gene_pproc = coupler::receive_gene_pproc<int>(dir, IO[0], eng[0]);
  //coupler::Array1d*<double> xgc_pproc = coupler::receive_xgc_pproc(dir, IO[1], eng[1]);

  //Perform coupling routines
  coupler::destroy(gene_pproc);

  for (int i = 0; i < time_step; i++) {
    for (int j = 0; j < RK_count; j++) {
      coupler::Array2d<double>* density = coupler::receive_density(dir, IO[2], eng[2]);
      coupler::printSomeDensityVals(density);
      coupler::send_density(dir, density, IO[3], eng[3], send_var[0]);
      coupler::destroy(density);

      coupler::Array2d<double>* field = coupler::receive_field(dir, IO[4], eng[4]);
      coupler::send_field(dir, field, IO[5], eng[5], send_var[1]);
      coupler::destroy(field);
    }
  }

  coupler::close_engines(eng, 6);
  Kokkos::finalize();
  MPI_Finalize();
  return 0;
}

