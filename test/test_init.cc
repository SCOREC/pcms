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
//std::cout<<"3"<<'\n'; 
//MPI_Barrier(MPI_COMM_WORLD);
//  MpiFreeComm(p1pp3d);
 
  MPI_Finalize(); 
  return 0;
}
