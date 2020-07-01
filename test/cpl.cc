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
  const std::string dir = "../coupling";
  const int time_step = atoi(argv[1]), RK_count = 4;

  fprintf(stderr, "%d number of time steps: %d\n", rank, time_step);

  const std::string xmlfile = "adios2cfg.xml";
  adios2::ADIOS adios(MPI_COMM_WORLD, adios2::DebugON);
  adios2::Variable<double> senddensity;
  adios2::Variable<coupler::CV> sendfield;

  coupler::adios2_handler gDens(adios,"gene_density");
  coupler::adios2_handler cDens(adios,"cpl_density");
  coupler::adios2_handler xFld(adios,"xgc_field");
  coupler::adios2_handler cFld(adios,"cpl_field");
  coupler::adios2_handler gQP(adios,"gene_pproc_qp");
  coupler::adios2_handler gRX(adios,"gene_pproc_rx");
  coupler::adios2_handler gCy(adios,"gene_cy_array");
  coupler::adios2_handler gInt(adios,"gene_pproc_i");
  coupler::adios2_handler xXcoord(adios,"xgc_x_coordss");
  coupler::adios2_handler xSurf(adios,"xgc_numsurfs");
  coupler::adios2_handler xZcoord(adios,"xgc_z_coordss");
  coupler::adios2_handler xVsurf(adios,"xgc_versurfs");
  coupler::adios2_handler xCce(adios,"xgc_cce_data");

  //receive GENE's preproc mesh discretization values
  coupler::Array1d<double>* q_prof = coupler::receive_gene_pproc<double>(dir, gQP);
  coupler::Array1d<double>* gene_xval = coupler::receive_gene_pproc<double>(dir, gRX);//matching gene's xval arr
  coupler::Array1d<int>* gene_parpar = coupler::receive_gene_pproc<int>(dir, gInt);

  //intialize GENE class
  const bool preproc = true;
  const bool ypar = false;
  coupler::TestCase test_case = coupler::TestCase::off;
  coupler::Array1d<double>* gene_cy = coupler::receive_gene_pproc<double>(dir, gCy);
  coupler::Part1ParalPar3D p1pp3d(gene_parpar->data(),gene_xval->data(),q_prof->data(), gene_cy->data(), preproc);


  //receive XGC's preproc mesh discretization values
  coupler::Array1d<double>* xgc_xcoords = coupler::receive_gene_pproc<double>(dir, xXcoord);//x_XGC
  coupler::Array1d<int>* xgc_numsurf = coupler::receive_gene_pproc<int>(dir, xSurf);
  int numsurf = xgc_numsurf->val(0);
  int block_count = xgc_numsurf->val(1);

  double* xgc_zcoords = coupler::receive_gene_exact<double>(dir,xZcoord, 0, block_count);
  int* xgc_versurf = coupler::receive_gene_exact<int>(dir,xVsurf, p1pp3d.li1, p1pp3d.nx0);
  coupler::Array1d<int>* xgc_cce = coupler::receive_gene_pproc<int>(dir, xCce);
  coupler::Part3Mesh3D p3m3d(p1pp3d, numsurf, block_count, xgc_versurf, xgc_cce->data(), xgc_xcoords->data(), xgc_zcoords, preproc);
  const int nummode = 1;
  coupler::BoundaryDescr3D bdesc(p3m3d, p1pp3d, test_case, preproc);
  if(!p1pp3d.mype)std::cerr << "0.8"<< "\n"; 
  coupler::DatasProc3D dp3d(p1pp3d, p3m3d, preproc, test_case, ypar, nummode);
  if(!p1pp3d.mype)std::cerr << "0.9"<< "\n";

//  dp3d.Initmattoplane(p3m3d,p1pp3d);

  coupler::destroy(q_prof);
  coupler::destroy(gene_xval);
  coupler::destroy(gene_parpar);

  dp3d.InitFourierPlan3D();


  int m;
  for (int i = 0; i < time_step; i++) {
    for (int j = 0; j < RK_count; j++) {
      coupler::GO start[2]={0, p1pp3d.blockstart};
      coupler::GO count[2]={coupler::GO(p1pp3d.lj0), p1pp3d.blockcount};
std::cout<<"start count"<<start[0]<<" "<<start[1]<<" "<<count[0]<<" "<<count[1]<<'\n';
      m=i*RK_count+j;
      coupler::Array2d<coupler::CV>* densityfromGENE = coupler::receive_density(dir, gDens,start,count,MPI_COMM_WORLD,m);
/*
for(coupler::LO i=0;i<count[0]*count[1];i++){
std::cout<<"i="<<i<<" "<<densityfromGENE->val(i)<<'\n';
 }
*/
std::cout<<"111"<<'\n';
      dp3d.DistriDensiRecvfromPart1(p3m3d,p1pp3d,densityfromGENE);
std::cout<<"222"<<'\n';
//      coupler::destroy(densityfromGENE);
/* 
if(p1pp3d.mype==0){
for(int i=0;i<p1pp3d.li0;i++){
  for(int j=0;j<p1pp3d.lj0;j++){
    for(int k=0;k<bdesc.nzb;k++){
      std::cout<<"i,j,k="<<i<<" "<<j<<" "<<k<<" "<<dp3d.densin[i][j][p1pp3d.lk2-2+k]-bdesc.lowdenz[i][j][k]<<'\n';
    }
  }
}
}
*/

//      bdesc.zDensityBoundaryBufAssign(dp3d.densin,p1pp3d);
/*
if(p1pp3d.mype==0){
for(int i=0;i<p1pp3d.li0;i++){
  for(int j=0;j<p1pp3d.lj0;j++){
    for(int k=0;k<bdesc.nzb;k++)
     std::cout<<"i,j="<<i<<" "<<j<<" "<<bdesc.lowdenz[i][j][k]<<'\n';
  }
}
}
*/


//      dp3d.InterpoDensity3D(bdesc,p3m3d,p1pp3d);
//      dp3d.CmplxdataToRealdata3D();
//      dp3d.DensityToPart3(p3m3d);

      dp3d.AssemDensiSendtoPart3(bdesc,p3m3d,p1pp3d);
      coupler::Array2d<double>* densitytoXGC = new coupler::Array2d<double>(
                                                    p3m3d.activenodes,p3m3d.lj0,p3m3d.blockcount,p3m3d.lj0, 
                                                    p3m3d.blockstart);
      double* densitytmp = densitytoXGC->data();
      for(int h=0;h<p3m3d.lj0*p3m3d.blockcount;h++){
        densitytmp[h] = dp3d.denssend[h]; 
      }
/*
if(p1pp3d.mype==0){
for(int i=0;i<p3m3d.lj0*p3m3d.blockcount;i++){
     std::cout<<"i="<<i<<" "<<densitytmp[i]<<'\n';
}
} 
*/
 
      coupler::send_from_coupler(adios,dir,densitytoXGC,cDens.IO,cDens.eng,cDens.name,senddensity,MPI_COMM_WORLD,m);    
      coupler::destroy(densitytoXGC);
      coupler::destroy(densityfromGENE);
 
      coupler::GO start_1[2]={0,p3m3d.blockstart+p3m3d.cce_first_node-1};
      coupler::GO count_1[2]={coupler::GO(p3m3d.lj0),p3m3d.blockcount}; 

      coupler::Array2d<double>* fieldfromXGC = coupler::receive_field(dir, xFld,start_1,count_1,MPI_COMM_WORLD,m);

      dp3d.DistriPotentRecvfromPart3(p3m3d,p1pp3d,fieldfromXGC);
      coupler::destroy(fieldfromXGC);
 
      dp3d.RealdataToCmplxdata3D();
      dp3d.zPotentBoundaryBufAssign(bdesc,p3m3d,p1pp3d);
      dp3d.InterpoPotential3D(bdesc,p3m3d,p1pp3d);       

      dp3d.AssemPotentSendtoPart1(p3m3d,p1pp3d);
    
      coupler::Array2d<coupler::CV>* fieldtoGENE = new coupler::Array2d<coupler::CV>(
                                                   p1pp3d.totnodes,p1pp3d.lj0,p1pp3d.blockcount,
                                                   p1pp3d.lj0,p1pp3d.blockstart);
      coupler::CV* fieldtmp = fieldtoGENE->data(); 
      for(coupler::GO h=0;h<p1pp3d.lj0*p1pp3d.blockcount;h++){
        fieldtmp[h] = dp3d.potentsend[h];
      }         

      coupler::send_from_coupler(adios,dir,fieldtoGENE,cFld.IO,cFld.eng,cFld.name,sendfield,MPI_COMM_WORLD,m);

      MPI_Barrier(MPI_COMM_WORLD);
      std::cout<<"6666"<<'\n';
      coupler::destroy(fieldtoGENE);
      std::cerr << p1pp3d.mype << " done loop " << i << " " << j << "\n";
    }
  }

  gDens.close();
  cDens.close();
  xFld.close();
  cFld.close();
  gQP.close();
  gRX.close();
  gInt.close();
  xZcoord.close();
  xVsurf.close();
  xCce.close();
  xXcoord.close();
  xSurf.close();

  std::cerr << p1pp3d.mype << " before kokkos finalize\n";
  Kokkos::finalize();
  std::cerr << p1pp3d.mype << " done kokkos finalize\n";
  MPI_Finalize();
  return 0;
}

