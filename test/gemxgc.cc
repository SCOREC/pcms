#include "adios2Routines.h"
#include "couplingTypes.h"
#include "dataprocess.h"
#include "importpart3mesh.h"
#include "commpart1.h"
#include "BoundaryDescr3D.h"
#include "testutilities.h"
#include <string>
#include <fstream>

int main(int argc, char **argv){
  int rank;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  const std::string dir = "../coupling";
  const int time_step = atoi(argv[1]), RK_count = 4;

  fprintf(stderr, "%d number of time steps: %d\n", rank, time_step);

  const std::string xmlfile = "adios2cfg.xml";
  adios2::ADIOS adios(MPI_COMM_WORLD, adios2::DebugON);
  adios2::Variable<double> senddensity;
  adios2::Variable<coupler::CV> sendfield;

  coupler::adios2_handler gDens(adios,"gem_density");
  coupler::adios2_handler cDens(adios,"cpl_density");
  coupler::adios2_handler xFld(adios,"xgc_field");
  coupler::adios2_handler cFld(adios,"cpl_field");
  coupler::adios2_handler gMesh(adios,"gem_mesh");
  coupler::adios2_handler gThf(adios,"gem_thfl");
  coupler::adios2_handler gGrd(adios,"gem_grid");
  coupler::adios2_handler gQprof(adios,"gem_qprf");
  coupler::adios2_handler xCouple(adios,"xgc_couple");
  coupler::adios2_handler xRzcoords(adios,"xgc_rzcoords");  

  std::string model = "global";
  coupler::Array1d<int>* gmesh=coupler::receive_pproc<int>(dir,gMesh,model);
  if(!rank)fprintf(stderr, "ABJ done: nr %d, ntheta: %d\n",gmesh->val(5),gmesh->val(4));
  const int m =1;
  coupler::GO start[2] = {0,0};
  coupler::GO count[2] = {coupler::GO(gmesh->val(5)),coupler::GO(gmesh->val(4))};
  coupler::Array2d<double>* thflx=coupler::receive_pproc_2d<double>(dir,gThf,start,count,m);
  coupler::GO start2[2] = {0,0};
  coupler::GO count2[2] = {2,coupler::GO(gmesh->val(6))};
  coupler::Array1d<double>* qprof=coupler::receive_pproc<double>(dir,gQprof,model);
 
  //intialize GEM class
  const bool preproc = true;
  const bool ypar = false;
  coupler::TestCase test_case = coupler::TestCase::off;
  coupler::CouplingCase ccase = coupler::CouplingCase::gemxgc; 
  //coupler::Part1ParalPar3D p1pp3d(gmesh->data(),thflx->data(),qprof->data(),test_case,preproc);

  coupler::Array1d<coupler::LO>* xcouple=coupler::receive_pproc<coupler::LO>(dir,xCouple,model);
  coupler::Array1d<double>* rzcoords=coupler::receive_pproc<double>(dir,xRzcoords,model);
  coupler::Part3Mesh3D p3m3d(xcouple,rzcoords,preproc,test_case);

  return 0;
}
