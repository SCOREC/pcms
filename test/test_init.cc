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
  const int nummode=3;
//  std::string test_dir(argv[1]);
  std::string test_dir="/gpfs/u/home/MPFS/MPFSshng/barn/wdmapp_coupling_data/testdatas/";
  coupler::Part1ParalPar3D p1pp3d(preproc, test_case, test_dir);  
  coupler::Part3Mesh3D p3m3d(p1pp3d, preproc, test_case, test_dir);
  coupler::DatasProc3D dp3d(p1pp3d, p3m3d, preproc, ypar);
  coupler::BoundaryDescr3D bdesc(p3m3d,p1pp3d,dp3d,test_case);
  coupler::TestInitPotentAlongz(dp3d,p3m3d,p1pp3d.npy,nummode);
  bdesc.zPotentBoundaryBufAssign(dp3d,p3m3d,p1pp3d);
  coupler::InterpoPotential3D(bdesc,p3m3d,p1pp3d,dp3d,preproc);
  dp3d.InitFourierPlan3D(p1pp3d,p3m3d); 
  dp3d.RealdataToCmplxdata3D(p1pp3d,p3m3d);
//std::cout<<"3"<<'\n'; 
//MPI_Barrier(MPI_COMM_WORLD);
//  MpiFreeComm(p1pp3d);

//  coupler::FreeFourierPlan3D(dp3d); 
  MPI_Finalize(); 
  return 0;
}
