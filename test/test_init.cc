#include "classes.h"
#include "commpart1.h"
#include "importpart3mesh.h"
#include "dataprocess.h"
#include "testutilities.h"

int main(int argc, char* argv[])
{
  MPI_Init(&argc,&argv);
  coupler::Part1ParalPar3D p1pp3d;
  coupler::Part3Mesh3D p3m3d;
  coupler::DatasProc3D dp3d;
  coupler::BoundaryDescr3D bdesc;  
  coupler::InitPart1ParalPar3D(p1pp3d);
  coupler::ImportPart3Mesh3D(p3m3d,p1pp3d);
  coupler::InitDatasProc3Dparameters(dp3d,p1pp3d,p3m3d);
  coupler::AllocDatasProc3dDensityArrays(dp3d,p1pp3d,p3m3d);
  coupler::AllocDatasProc3dPotentArrays(dp3d,p1pp3d,p3m3d);
  coupler::InitBoundaryDescr3D(bdesc,p3m3d,p1pp3d,dp3d);
  coupler::TestInitPotentAlongz(dp3d,p3m3d,p1pp3d.npy,3);
  coupler::zPotentBoundaryBufAssign(MPI_DOUBLE,bdesc,dp3d,p3m3d,p1pp3d);
  coupler::InterpoPotential3D(bdesc,p3m3d,p1pp3d,dp3d);

//  coupler::InitFourierPlan3D(dp3d,p1pp3d,p3m3d);  
//  coupler::PotentRealdataToCmplxdata3D(dp3d,p1pp3d,p3m3d); 

//std::cout<<"3"<<'\n'; 
//MPI_Barrier(MPI_COMM_WORLD);
//  MpiFreeComm(p1pp3d);

//  coupler::FreeFourierPlan3D(dp3d); 
  MPI_Finalize(); 
  return 0;
}
