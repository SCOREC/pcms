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

  adios2::ADIOS adios(MPI_COMM_WORLD, adios2::DebugON);
  adios2::Variable<double> send_var[2];

  const std::string dir = "../coupling";
  const int time_step = 1, RK_count = 4;

  coupler::adios2_handler gDens(adios,"gene_density");
  coupler::adios2_handler cDens(adios,"cpl_density");
  coupler::adios2_handler xFld(adios,"xgc_field");
  coupler::adios2_handler cFld(adios,"cpl_field");
  coupler::adios2_handler gRZ(adios,"gene_pproc_rz");
  coupler::adios2_handler gRX(adios,"gene_pproc_rx");
  coupler::adios2_handler gInt(adios,"gene_pproc_i");
  coupler::adios2_handler gComp(adios,"gene_pproc_c");
  coupler::adios2_handler xNum(adios,"xgc_numsurf");
  coupler::adios2_handler xZNum(adios,"xgc_znum");
  coupler::adios2_handler xXcoord(adios,"xgc_xcoords");
  coupler::adios2_handler xVsurf(adios,"xgc_versurf");

  //receive GENE's preproc mesh discretization values
  coupler::Array1d<double>* gene_pproc_rz = coupler::receive_gene_pproc<double>(dir, gRZ, "gene_pproc_rz");
  coupler::Array1d<double>* gene_pproc_rx = coupler::receive_gene_pproc<double>(dir, gRX, "gene_pproc_rx");
  coupler::Array1d<int>* gene_pproc_i = coupler::receive_gene_pproc<int>(dir, gInt, "gene_pproc_i");

  //intialize GENE class
  const bool preproc = true;
  const bool ypar = false;
  coupler::Part1ParalPar3D p1pp3d(gene_pproc_i->data(),gene_pproc_rx->data(),preproc);

  //receive XGC's preproc mesh discretization values
  coupler::Array1d<int>* xgc_numsurf = coupler::receive_gene_pproc<int>(dir, xNum, "xgc_numsurfs");
  int num_surf = xgc_numsurf->val(0);
  coupler::Array1d<double>* xgc_xcoords = coupler::receive_gene_pproc<double>(dir, xXcoord, "xgc_x_coordss");
  coupler::Array1d<double>* xgc_versurf = coupler::receive_gene_pproc<double>(dir, xVsurf, "xgc_versurf");
  coupler::Array1d<int>* xgc_znum = coupler::receive_gene_pproc<int>(dir, xZNum, "xgc_pproc_zi");
  coupler::Array1d<std::complex<double>>* gene_pproc_c = coupler::receive_gene_pproc<std::complex<double>>(dir, gComp, "gene_pproc_c");

  coupler::Part3Mesh3D p3m3d(p1pp3d, 
    xgc_znum->data(), xgc_xcoords->data(),
    xgc_versurf->data(), preproc);
  //coupler::<geneclass> 
  //coupler::DatasProc3D dp3d(p1pp3d,p3m3d, preproc, ypar);
  //coupler::BoundaryDescr3D bdesc(p3m3d,p1pp3d,dp3d);

  coupler::destroy(gene_pproc_rz);
  coupler::destroy(gene_pproc_rx);
  coupler::destroy(gene_pproc_i);
  coupler::destroy(gene_pproc_c);

  for (int i = 0; i < time_step; i++) {
    for (int j = 0; j < RK_count; j++) {
      coupler::Array2d<double>* density = coupler::receive_density(dir, gDens);
      coupler::printSomeDensityVals(density);
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

