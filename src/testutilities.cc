//#include "importpart3mesh.h"
#include "testutilities.h"

namespace coupler {

void InitPart1paral3DInCoupler(Part1ParalPar3D  &p1pp3d)
{
/*
  LO size;
  LO rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
*/
  LO* data=new LO[12];
  std::string fname=test_dir+"parpart1.nml";
  InputfromFile(data,12,fname); 
 
  p1pp3d.npx=data[0];
  p1pp3d.nx0=data[1];
  p1pp3d.nxb=data[2];
  p1pp3d.li0=data[3];

  p1pp3d.npy=data[4];
  p1pp3d.ny0=data[5];
  p1pp3d.nyb=data[6];
  p1pp3d.lj0=data[7];

  p1pp3d.npz=data[8];
  p1pp3d.nz0=data[9];
  p1pp3d.nzb=data[10];
  p1pp3d.lk0=data[11];

  p1pp3d.NP=p1pp3d.npx*p1pp3d.npy*p1pp3d.npz;
  CreateSubCommunicators(p1pp3d);
  p1pp3d.li1=p1pp3d.mype_x*p1pp3d.li0;
  p1pp3d.li2=p1pp3d.li1+p1pp3d.li0-1;
  p1pp3d.lj1=p1pp3d.mype_y*p1pp3d.lj0;
  p1pp3d.lj2=p1pp3d.lj1+p1pp3d.lj0-1;
  p1pp3d.lk1=p1pp3d.mype_z*p1pp3d.lk0;
  p1pp3d.lk2=p1pp3d.lk1+p1pp3d.lk0-1;
  
  delete[] data;
}






}
