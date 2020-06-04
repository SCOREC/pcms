#include "adios2Routines.h"
#include "couplingTypes.h"
#include "dataprocess.h"
#include "importpart3mesh.h"
#include "commpart1.h"
#include "BoundaryDescr3D.h"
#include "testutilities.h"
#include <string>
#include <fstream> 

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
  if(argc != 2) {
    if(!rank) printf("Usage: %s <number of timesteps>\n", argv[1]);
    exit(EXIT_FAILURE);
  }
  adios2::ADIOS adios(MPI_COMM_WORLD, adios2::DebugON);
  adios2::Variable<double> send_var[2];

  const std::string dir = "../coupling";
  const int time_step = atoi(argv[1]), RK_count = 4;

  coupler::adios2_handler gDens(adios,"gene_density");
  coupler::adios2_handler cDens(adios,"cpl_density");
  coupler::adios2_handler xFld(adios,"xgc_field");
  coupler::adios2_handler cFld(adios,"cpl_field");
  coupler::adios2_handler gQP(adios,"gene_pproc_qp");
  coupler::adios2_handler gRX(adios,"gene_pproc_rx");
  coupler::adios2_handler gInt(adios,"gene_pproc_i");
  coupler::adios2_handler gComp(adios,"gene_pproc_c");
  coupler::adios2_handler xSurf(adios,"xgc_numsurfs");
  coupler::adios2_handler xXcoord(adios,"xgc_x_coordss");
  coupler::adios2_handler xZcoord(adios,"xgc_z_coordss");
  coupler::adios2_handler xVsurf(adios,"xgc_versurfs");
  coupler::adios2_handler xCce(adios,"xgc_cce_data");

  //receive GENE's preproc mesh discretization values
  coupler::Array1d<double>* q_prof = coupler::receive_gene_pproc<double>(dir, gQP);
  coupler::Array1d<double>* gene_pproc_rx = coupler::receive_gene_pproc<double>(dir, gRX);
  coupler::Array1d<int>* gene_pproc_i = coupler::receive_gene_pproc<int>(dir, gInt);

  //intialize GENE class
  const bool preproc = true;
  const bool ypar = false;
  coupler::TestCase test_case = coupler::TestCase::off;
  coupler::Part1ParalPar3D p1pp3d(gene_pproc_i->data(),gene_pproc_rx->data(),preproc);

  //receive XGC's preproc mesh discretization values
  coupler::Array1d<int>* xgc_numsurf = coupler::receive_gene_pproc<int>(dir, xSurf);
  coupler::Array1d<double>* xgc_xcoords = coupler::receive_gene_pproc<double>(dir, xXcoord, p1pp3d.comm_x);
  coupler::Array1d<double>* xgc_zcoords = coupler::receive_gene_pproc<double>(dir, xZcoord, p1pp3d.comm_z);
  coupler::Array1d<int>* xgc_versurf = coupler::receive_gene_pproc<int>(dir, xVsurf);
  coupler::Array1d<int>* xgc_cce = coupler::receive_gene_pproc<int>(dir, xCce);

  coupler::Array1d<coupler::CV>* gene_pproc_c = coupler::receive_gene_pproc<coupler::CV>(dir, gComp);
  coupler::Part3Mesh3D p3m3d(p1pp3d, xgc_numsurf->val(0), xgc_versurf->data(), xgc_xcoords->data(), xgc_zcoords->data(), preproc);
  const int nummode = 1;
  coupler::DatasProc3D dp3d(p1pp3d, p3m3d, preproc, test_case, ypar, nummode);
  coupler::BoundaryDescr3D bdesc(p3m3d, p1pp3d, dp3d, test_case, preproc);

  coupler::destroy(q_prof);
  coupler::destroy(gene_pproc_rx);
  coupler::destroy(gene_pproc_i);
  coupler::destroy(gene_pproc_c);

  for (int i = 0; i < time_step; i++) {
    for (int j = 0; j < RK_count; j++) {
      coupler::Array2d<double>* density = coupler::receive_density(dir, gDens);
      coupler::printSomeDensityVals(density);
      /* For sends involving subcomm, the adios object must be created with the subcomm
      adios2::ADIOS adios_x(p1pp3d.comm_x, adios2::DebugON);
      coupler::adios2_handler xdens(adios_x,"cpl_dens_data");
      coupler::send_density(dir, density, xdens, send_var[0]);
      */
      coupler::send_density(dir, density, cDens, send_var[0]);
      coupler::destroy(density);

      coupler::Array2d<double>* field = coupler::receive_field(dir, xFld);
      coupler::send_field(dir, field, cFld, send_var[1]);
      coupler::destroy(field);
    }
  }

  Kokkos::finalize();
  MPI_Finalize();
  return 0;
}

