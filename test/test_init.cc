#include "couplingTypes.h"
#include "dataprocess.h"
#include "importpart3mesh.h"
#include "commpart1.h"
#include "BoundaryDescr3D.h"
#include "testutilities.h"
#include <string>
#include <mpi.h>

int main(int argc, char* argv[])
{
  MPI_Init(&argc,&argv);
  coupler::TestCase test_case = coupler::TestCase::t0;
  const bool preproc = true;
  const bool ypar = false;
  std::string test_dir(argv[1]);
  coupler::Part1ParalPar3D p1pp3d(preproc, test_case, test_dir);    
  coupler::Part3Mesh3D p3m3d(p1pp3d, preproc, test_case, test_dir);
  const int nummode = 1;
  coupler::Part1ParalPar3D* mesh1=&p1pp3d;
  coupler::Part3Mesh3D*     mesh3=&p3m3d;
  coupler::DatasProc3D dp3d(mesh1, mesh3, preproc, test_case, ypar, nummode);  

  coupler::BoundaryDescr3D bdesc(p3m3d,p1pp3d,test_case,preproc);

  dp3d.InitFourierPlan3D(); 
  dp3d.RealdataToCmplxdata3D();
  dp3d.zPotentBoundaryBufAssign(bdesc);
  dp3d.InterpoPotential3D(bdesc);

  for(coupler::LO i=0;i<p1pp3d.li0;i++){
    for(coupler::LO j=0;j<p1pp3d.lj0;j++){
      for(coupler::LO k=0;k<p1pp3d.lk0;k++){
        dp3d.densin[i][j][k]=dp3d.potentpart1[i][j][k];
      }
    }   
  }
 

  dp3d.zDensityBoundaryBufAssign(dp3d.densin,bdesc);
  dp3d.InterpoDensity3D(bdesc); 
  dp3d.CmplxdataToRealdata3D();


  
  if(p1pp3d.mype==2){
   for(coupler::LO i=0;i<p3m3d.li0;i++){
     for(coupler::LO k=0;k<p3m3d.mylk0[i];k++){
      for(coupler::LO j=0;j<p3m3d.lj0;j++){
	 std::cout<<i<<" "<<k<<" "<<j<<" "<<dp3d.denspart3[i][j][k]-dp3d.potentin[i][j][k]<<'\n';       
	}
      }
    } 
  }


if(p1pp3d.mype==2){
  for(coupler::LO i=0;i<p1pp3d.li0;i++){
    for(coupler::LO j=0;j<p1pp3d.lj0;j++){
      for(coupler::LO k=0;k<p1pp3d.lk0;k++){
         std::cout<<i<<" "<<k<<" "<<j<<" "<<dp3d.potentpart1[i][j][k]<<'\n';
      }
    }
  }
}
  MPI_Finalize(); 
  return 0;
}
