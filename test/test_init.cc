#include "couplingTypes.h"
#include "dataprocess.h"
#include "importpart3mesh.h"
#include "commpart1.h"
#include "BoundaryDescr3D.h"
#include "dataprocess_impl.h"
#include "testutilities.h"
#include "inittestenv.h"
#include <string>
#include <mpi.h>

int main(int argc, char* argv[])
{
  MPI_Init(&argc,&argv);
  coupler::TestCase test_case = coupler::TestCase::t0;
  const bool preproc = true;
  const bool ypar = false;
  const int nummode=1;
//  std::string test_dir(argv[1]);
  std::string test_dir="/gpfs/u/home/MPFS/MPFSshng/barn/wdmapp_coupling_data/testdatas/";
  coupler::Part1ParalPar3D p1pp3d(preproc, test_case, test_dir);  
  coupler::Part3Mesh3D p3m3d(p1pp3d, preproc, test_case, test_dir);
  coupler::DatasProc3D dp3d(p1pp3d, p3m3d, preproc, ypar);
  coupler::BoundaryDescr3D bdesc(p3m3d,p1pp3d,dp3d,test_case,preproc);
  coupler::TestInitPotentAlongz(dp3d,p3m3d,p1pp3d.npy,nummode);
  bdesc.zPotentBoundaryBufAssign(dp3d,p3m3d,p1pp3d);
  coupler::InterpoPotential3D(bdesc,p3m3d,p1pp3d,dp3d,preproc);
  dp3d.InitFourierPlan3D(p1pp3d,p3m3d); 
  dp3d.RealdataToCmplxdata3D(p1pp3d,p3m3d);

  for(coupler::LO i=0;i<p1pp3d.li0;i++){
    for(coupler::LO j=0;j<p1pp3d.lj0;j++){
      for(coupler::LO k=0;k<p1pp3d.lk0;k++){
        dp3d.densin[i][j][k]=dp3d.potentpart1[i][j][k];
      }
    }   
  }
  dp3d.DatasProc3D::CmplxdataToRealdata3D(p1pp3d);
  coupler::zDensityBoundaryBufAssign(bdesc.nzb,p1pp3d.li0,p1pp3d.lj0*2,p1pp3d.lk0,
           bdesc.lowdenz,bdesc.updenz,dp3d.densout,p1pp3d);
  coupler::InterpoDensity3D(bdesc,p3m3d,p1pp3d,dp3d,preproc); 

if(p1pp3d.mype==2){
 for(coupler::LO i=0;i<p3m3d.li0;i++){
   for(coupler::LO k=0;k<p3m3d.mylk0[i];k++){
    for(coupler::LO j=0;j<p3m3d.lj0;j++){
       std::cout<<i<<" "<<k<<" "<<j<<" "<<dp3d.denspart3[i][j][k]-dp3d.potentin[i][j][k]<<'\n';       
      }
    }
  } 
}


/*
if(p1pp3d.mype==0){
 for(coupler::LO i=0;i<p3m3d.li0;i++){
   for(coupler::LO k=0;k<p1pp3d.lk0;k++){
    for(coupler::LO j=0;j<p3m3d.lj0;j++){
       std::cout<<i<<" "<<k<<" "<<dp3d.densout[i][j][k]-dp3d.potentinterpo[i][j][k]<<'\n';       
//std::cout<<i<<" "<<k<<" "<<dp3d.densout[i][j][k]-dp3d.potentinterpo[i][j][k]<<'\n';    
      }
    }
  } 
}
*/


/*
    for(coupler::LO j=0;j<p3m3d.lj0;j++){
       std::cout<<dp3d.potentpart3[1][j][10]<<" "<<dp3d.potentinterpo[1][j][6]<<'\n';
      }
*/
/*
if(p1pp3d.mype==0){ 
 for(coupler::LO i=0;i<p1pp3d.li0;i++){
    for(coupler::LO k=0;k<p1pp3d.lk0;k++){
      for(coupler::LO j=0;j<p1pp3d.lj0;j++){
     //   std::cout<<"check "<<dp3d.denspart3[i][j][k]-dp3d.potentin[i][j][k]<<'\n'; 
       std::cout<<i<<" "<<k<<" "<<dp3d.potentpart1[i][j][k]<<'\n';       
      }
    }
  } 
}
*/
//  coupler::FreeFourierPlan3D(dp3d); 
  MPI_Finalize(); 
  return 0;
}
