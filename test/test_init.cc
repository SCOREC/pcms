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

  coupler::DatasProc3D dp3d(p1pp3d,p3m3d, preproc, ypar);
  coupler::BoundaryDescr3D bdesc(p3m3d,p1pp3d,dp3d);
//std::cout<<"3"<<'\n'; 
//MPI_Barrier(MPI_COMM_WORLD);
//  MpiFreeComm(p1pp3d);
 
  MPI_Finalize(); 
  return 0;
}
