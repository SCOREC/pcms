#include "adios2Routines.h"
#include "couplingTypes.h"
#include "dataprocess.h"
#include "importpart3mesh.h"
#include "commpart1.h"
#include "BoundaryDescr3D.h"
#include "testutilities.h"
#include <string>
#include <fstream> 

struct adios2_handler{
public:
  adios2::IO IO;
  adios2::Engine eng;

  adios2_handler(adios2::ADIOS &adios, const std::string name):
	  IO(adios.DeclareIO(name))  {}
  ~adios2_handler(){
	  eng.Close();
  }
};


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

  adios2_handler h0(adios,"gene_density");
  adios2_handler h1(adios,"cpl_density");
  adios2_handler h2(adios,"xgc_field");
  adios2_handler h3(adios,"cpl_field");
  adios2_handler h4(adios,"gene_pproc_rz");
  adios2_handler h5(adios,"gene_pproc_rx");
  adios2_handler h6(adios,"gene_pproc_i");
  adios2_handler h7(adios,"gene_pproc_c");
  adios2_handler h8(adios,"xgc_numsurf");
  adios2_handler h9(adios,"xgc_znum");
  adios2_handler h10(adios,"xgc_xcoords");
  adios2_handler h11(adios,"xgc_versurf");

  //receive GENE's preproc mesh discretization values
  coupler::Array1d<double>* gene_pproc_rz = coupler::receive_gene_pproc<double>(dir, h4.IO, h4.eng, "gene_pproc_rz");
  coupler::Array1d<double>* gene_pproc_rx = coupler::receive_gene_pproc<double>(dir, h5.IO, h5.eng, "gene_pproc_rx");
  coupler::Array1d<int>* gene_pproc_i = coupler::receive_gene_pproc<int>(dir, h6.IO, h6.eng, "gene_pproc_i");

  //intialize GENE class
  const bool preproc = true;
  const bool ypar = false;
  if(!rank) std::cerr << rank << " 0.1\n"; 
  coupler::Part1ParalPar3D p1pp3d(gene_pproc_i,gene_pproc_rx,preproc);

  //receive XGC's preproc mesh discretization values
  coupler::Array1d<int>* xgc_numsurf = coupler::receive_gene_pproc<int>(dir, h8.IO, h8.eng, "xgc_numsurfs");
  int num_surf = xgc_numsurf->val(0);
  if(!rank) std::cerr << rank << " 0.2 " << num_surf <<"\n"; 
  coupler::Array1d<double>* xgc_xcoords = coupler::receive_gene_pproc<double>(dir, h10.IO, h10.eng, "xgc_x_coordss");

  //TODO: map the received xgc_xcoords correctly
  //for (int i = p1pp3d.li1; i < p1pp3d.li2; i++)
  //{
  coupler::Array1d<double>* xgc_versurf = coupler::receive_gene_pproc<double>(dir, h11.IO, h11.eng, "xgc_versurf");
  //}

  coupler::Array1d<int>* xgc_znum = coupler::receive_gene_pproc<int>(dir, h9.IO, h9.eng, "xgc_pproc_zi");
  coupler::Array1d<std::complex<double>>* gene_pproc_c = coupler::receive_gene_pproc<std::complex<double>>(dir, h7.IO, h7.eng, "gene_pproc_c");

  if(!rank) std::cerr << rank << " 0.3\n"; 
  coupler::Part3Mesh3D p3m3d(p1pp3d, 
    xgc_znum->data(), /*update*/xgc_pproc_x,
    /*update*/xgc_pproc_z, preproc);
  //coupler::<geneclass> 
  //coupler::DatasProc3D dp3d(p1pp3d,p3m3d, preproc, ypar);
  //coupler::BoundaryDescr3D bdesc(p3m3d,p1pp3d,dp3d);

  if(!rank) std::cerr << rank << " 0.4\n"; 
  coupler::destroy(gene_pproc_rz);
  coupler::destroy(gene_pproc_rx);
  coupler::destroy(gene_pproc_i);
  coupler::destroy(gene_pproc_c);

  for (int i = 0; i < time_step; i++) {
    for (int j = 0; j < RK_count; j++) {
      coupler::Array2d<double>* density = coupler::receive_density(dir, h0.IO, h0.eng);
      coupler::printSomeDensityVals(density);
      coupler::send_density(dir, density, h1.IO, h1.eng, send_var[0]);
      coupler::destroy(density);

      coupler::Array2d<double>* field = coupler::receive_field(dir, h2.IO, h2.eng);
      coupler::send_field(dir, field, h3.IO, h3.eng, send_var[1]);
      coupler::destroy(field);
    }
  }

  Kokkos::finalize();
  MPI_Finalize();
  return 0;
}

