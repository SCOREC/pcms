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
  adios2::Variable<double> senddensity;
  adios2::Variable<coupler::CV> sendfield;

  const std::string dir = "../coupling";
  const int time_step = atoi(argv[1]), RK_count = 4;

  coupler::adios2_handler gDens(adios,"gene_density");
  coupler::adios2_handler cDens(adios,"cpl_density");
  coupler::adios2_handler xFld(adios,"xgc_field");
  coupler::adios2_handler cFld(adios,"cpl_field");
  coupler::adios2_handler gQP(adios,"gene_pproc_qp");
  coupler::adios2_handler gRX(adios,"gene_pproc_rx");
  coupler::adios2_handler gCy(adios,"gene_cy_array");
  coupler::adios2_handler gInt(adios,"gene_pproc_i");
  coupler::adios2_handler gComp(adios,"gene_pproc_c");
  coupler::adios2_handler xXcoord(adios,"xgc_x_coordss");
  coupler::adios2_handler xSurf(adios,"xgc_numsurfs");
  coupler::adios2_handler xZcoord(adios,"xgc_z_coordss");
  coupler::adios2_handler xVsurf(adios,"xgc_versurfs");
  coupler::adios2_handler xCce(adios,"xgc_cce_data");

  //receive GENE's preproc mesh discretization values
  coupler::Array1d<double>* q_prof = coupler::receive_gene_pproc<double>(dir, gQP);
  coupler::Array1d<double>* gene_xval = coupler::receive_gene_pproc<double>(dir, gRX);//matching gene's xval arr
  coupler::Array1d<double>* gene_cy = coupler::receive_gene_pproc<double>(dir, gCy);
  double* C_y = gene_cy->data();
  double minor_r = gene_cy->val(180);
  double lx_a = gene_cy->val(181);
  double sign_phi = gene_cy->val(182);
  double dx = gene_cy->val(183);
  coupler::Array1d<int>* gene_parpar = coupler::receive_gene_pproc<int>(dir, gInt);

  //intialize GENE class
  const bool preproc = true;
  const bool ypar = false;
  coupler::TestCase test_case = coupler::TestCase::off;
  coupler::Part1ParalPar3D p1pp3d(gene_parpar->data(),gene_xval->data(),q_prof->data(),preproc);

  //receive XGC's preproc mesh discretization values
  coupler::Array1d<double>* xgc_xcoords = coupler::receive_gene_pproc<double>(dir, xXcoord);//x_XGC
  coupler::Array1d<int>* xgc_numsurf = coupler::receive_gene_pproc<int>(dir, xSurf);
  int numsurf = xgc_numsurf->val(0);
  int block_count = xgc_numsurf->val(1);

  double* xgc_zcoords = coupler::receive_gene_exact<double>(dir,xZcoord, 0, block_count);
  /*if(p1pp3d.mype_x==1) {
    std::cerr << "foo block_count " << block_count << "\n";
    for(int i=0; i < (block_count); i++)  
      std::cerr<< "foo "  << xgc_zcoords[i] << "\n";
  }*/
  int* xgc_versurf = coupler::receive_gene_exact<int>(dir,xVsurf, p1pp3d.li1, p1pp3d.nx0);
  /*if(p1pp3d.mype_x==1){
    for(int i=0;i<p1pp3d.nx0;i++){
      std::cout<<"i="<<i<<" "<<xgc_versurf[i]<<'\n';
    } 
  }*/
  coupler::Array1d<int>* xgc_cce = coupler::receive_gene_pproc<int>(dir, xCce);
  coupler::Array1d<coupler::CV>* gene_moments = coupler::receive_gene_pproc<coupler::CV>(dir, gComp);//matching gene's moments arr
  coupler::Part3Mesh3D p3m3d(p1pp3d, numsurf, block_count, xgc_versurf, xgc_cce->data(), xgc_xcoords->data(), xgc_zcoords, preproc);

  const int nummode = 1;
  coupler::DatasProc3D dp3d(p1pp3d, p3m3d, preproc, test_case, ypar, nummode);
  if(!p1pp3d.mype)std::cerr << "0.8"<< "\n";
  coupler::BoundaryDescr3D bdesc(p3m3d, p1pp3d, dp3d, test_case, preproc);
  if(!p1pp3d.mype)std::cerr << "0.9"<< "\n";

  coupler::destroy(q_prof);
  coupler::destroy(gene_xval);
  coupler::destroy(gene_parpar);
//  coupler::destroy(gene_moments);

  dp3d.InitFourierPlan3D();

  for (int i = 0; i < time_step; i++) {
    for (int j = 0; j < RK_count; j++) {
      coupler::GO start[2]={0, p1pp3d.blockstart};
      coupler::GO count[2]={coupler::GO(p1pp3d.lj0), p1pp3d.blockcount};
      coupler::Array2d<coupler::CV>* densityfromGENE = coupler::receive_density(dir, gDens,start,count,p1pp3d.comm_x);

      dp3d.DistriDensiRecvfromPart1(p3m3d,p1pp3d,densityfromGENE);

      bdesc.zDensityBoundaryBufAssign(dp3d.densin,p1pp3d);
      dp3d.InterpoDensity3D(bdesc,p3m3d,p1pp3d);
      dp3d.CmplxdataToRealdata3D();

      dp3d.AssemDensiSendtoPart3(p3m3d,p1pp3d);

      coupler::Array2d<double>* densitytoXGC = new coupler::Array2d<double>(
                                                    p3m3d.activenodes,p3m3d.lj0,p3m3d.blockcount,p3m3d.lj0, 
                                                    p3m3d.blockstart);
      double* densitytmp = densitytoXGC->data();
      for(int h=0;h<p3m3d.lj0*p3m3d.blockcount;h++){
        densitytmp[h] = dp3d.denssend[h]; 
      }
      if(p1pp3d.mype_z==0){
        coupler::send_density(dir, densitytoXGC, cDens, senddensity);
      }     
      coupler::destroy(densityfromGENE);
      coupler::destroy(densitytoXGC);

      coupler::GO start_1[2]={0,p3m3d.blockstart+p3m3d.cce_first_node-1};
      coupler::GO count_1[2]={coupler::GO(p3m3d.lj0),p3m3d.blockcount}; 
      coupler::Array2d<double>* fieldfromXGC = coupler::receive_field(dir, xFld,start_1,count_1,p1pp3d.comm_x);

      dp3d.DistriPotentRecvfromPart3(p3m3d,p1pp3d,fieldfromXGC);

      dp3d.RealdataToCmplxdata3D();
      bdesc.zPotentBoundaryBufAssign(dp3d,p3m3d,p1pp3d);
      dp3d.InterpoPotential3D(bdesc,p3m3d,p1pp3d);       

      dp3d.AssemPotentSendtoPart1(p3m3d,p1pp3d);

      coupler::Array2d<coupler::CV>* fieldtoGENE = new coupler::Array2d<coupler::CV>(
                                                   p1pp3d.totnodes,p1pp3d.lj0,p1pp3d.blockcount,
                                                   p1pp3d.lj0,p1pp3d.blockstart);
      coupler::CV* fieldtmp = fieldtoGENE->data(); 
      for(coupler::GO h=0;h<p1pp3d.lj0*p1pp3d.blockcount;h++){
        fieldtmp[h] = dp3d.potentsend[h];
      }         
      if(p1pp3d.mype_z==0){
        coupler::send_field(dir, fieldtoGENE, cFld, sendfield);
      }

      coupler::destroy(fieldtoGENE);
      coupler::destroy(fieldfromXGC);
    }
  }

  Kokkos::finalize();
  MPI_Finalize();
  return 0;
}

