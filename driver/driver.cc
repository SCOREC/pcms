#include "adios2Routines.h"
#include "couplingTypes.h"
#include "dataprocess.h"
#include "importpart3mesh.h"
#include "commpart1.h"
#include "BoundaryDescr3D.h"
#include "testutilities.h"
#include <string>

void exParFor() {
  Kokkos::parallel_for(
      4, KOKKOS_LAMBDA(const int i) {
        printf("Hello from kokkos thread i = %i\n", i);
      });
}

int main(int argc, char **argv){
  int rank;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  Kokkos::initialize(argc, argv);
  if(!rank) {
    printf("Hello World on Kokkos execution space %s\n",
         typeid(Kokkos::DefaultExecutionSpace).name());
    exParFor();
  }

  const int obj_count = 11; 
  adios2::ADIOS adios(MPI_COMM_WORLD, adios2::DebugON);
  adios2::IO IO[obj_count];
  adios2::Engine eng[obj_count];
  adios2::Variable<double> send_var[2];

  const std::string dir = "../coupling";

  const std::string pproc_rz = "gene_pproc_rz";
  const std::string pproc_rx = "gene_pproc_rx";
  const std::string pproc_i = "gene_pproc_i";
  const std::string pproc_c = "gene_pproc_c";

  const std::string x_pproc_v = "xgc_pproc_v";
  const std::string x_pproc_x = "xgc_pproc_x";
  const std::string x_pproc_z = "xgc_pproc_z";
  const int time_step = 1, RK_count = 4;

  IO[0] = adios.DeclareIO("gene_density");
  IO[1] = adios.DeclareIO("cpl_density");
  IO[2] = adios.DeclareIO("xgc_field");
  IO[3] = adios.DeclareIO("cpl_field");
  IO[4] = adios.DeclareIO(pproc_rz);
  IO[5] = adios.DeclareIO(pproc_rx);
  IO[6] = adios.DeclareIO(pproc_i);
  IO[7] = adios.DeclareIO(pproc_c);
  IO[8] = adios.DeclareIO(x_pproc_v);
  IO[9] = adios.DeclareIO(x_pproc_x);
  IO[10] = adios.DeclareIO(x_pproc_z);

  //receive GENE's preproc mesh discretization values
  coupler::Array1d<double>* gene_pproc_rz = coupler::receive_gene_pproc<double>(dir, IO[4], eng[4], pproc_rz);
  coupler::Array1d<double>* gene_pproc_rx = coupler::receive_gene_pproc<double>(dir, IO[5], eng[5], pproc_rx);
  coupler::Array1d<int>* gene_pproc_i = coupler::receive_gene_pproc<int>(dir, IO[6], eng[6], pproc_i);
  coupler::Array1d<std::complex<double>>* gene_pproc_c = coupler::receive_gene_pproc<std::complex<double>>(dir, IO[7], eng[7], pproc_c);

  //intialize GENE class
  const bool preproc = true;
  const bool ypar = false;
  if(!rank) std::cerr << rank << " 0.1\n"; 
  coupler::Part1ParalPar3D p1pp3d(gene_pproc_i,gene_pproc_rx,preproc);

  if(!rank) std::cerr << rank << " 0.2\n"; 
  //receive XGC's preproc mesh discretization values
  coupler::Array1d<int>* xgc_pproc_v = {0};//coupler::receive_gene_pproc<int>(dir, IO[8], eng[8], x_pproc_v);
  coupler::Array1d<double>* xgc_pproc_x = {0};//coupler::receive_gene_pproc<double>(dir, IO[9], eng[9], x_pproc_x);
  coupler::Array1d<double>* xgc_pproc_z = {0};//coupler::receive_gene_pproc<double>(dir, IO[10], eng[10], x_pproc_z);

  if(!rank) std::cerr << rank << " 0.3\n"; 
  //coupler::Part3Mesh3D p3m3d(p1pp3d, xgc_pproc_v, xgc_pproc_x, xgc_pproc_z, preproc);
  //coupler::Part3Mesh3D p3m3d(p1pp3d);
  //coupler::DatasProc3D dp3d(p1pp3d,p3m3d, preproc, ypar);
  //coupler::BoundaryDescr3D bdesc(p3m3d,p1pp3d,dp3d);

  if(!rank) std::cerr << rank << " 0.4\n"; 
  coupler::destroy(gene_pproc_rz);
  coupler::destroy(gene_pproc_rx);
  coupler::destroy(gene_pproc_i);
  coupler::destroy(gene_pproc_c);

  for (int i = 0; i < time_step; i++) {
    for (int j = 0; j < RK_count; j++) {
      coupler::Array2d<double>* density = coupler::receive_density(dir, IO[0], eng[0]);
      coupler::printSomeDensityVals(density);
      coupler::send_density(dir, density, IO[1], eng[1], send_var[0]);
      coupler::destroy(density);

      coupler::Array2d<double>* field = coupler::receive_field(dir, IO[2], eng[2]);
      coupler::send_field(dir, field, IO[3], eng[3], send_var[1]);
      coupler::destroy(field);
    }
  }

  coupler::close_engines(eng, obj_count);
  Kokkos::finalize();
  MPI_Finalize();
  return 0;
}

