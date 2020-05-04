//#include "importpart3mesh.h"
#include "testutilities.h"
#include "commpart1.h"

namespace coupler {

void InitPart1paral3DInCoupler(Part1ParalPar3D  &p1pp3d)
{
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

void TestInitPotentAlongz(DatasProc3D& dp3d,Part3Mesh3D& p3m3d,LO npy,LO n)
{
  if(npy==1){
    LO li0,lj0,lk0;
    li0=p3m3d.li0;
    lj0=p3m3d.lj0;
    double ylen;
    double sum;
    double dy=2.0*cplPI/double(lj0);
    for(LO i=0;i<li0;i++){
      lk0=p3m3d.mylk0[i];
      for(LO k=0;k<lk0;k++){
	ylen=0.0;
	for(LO j=0;j<lj0;j++){
          ylen=double(j)*dy;
          sum=0.0;
          for(LO h=0;h<n;h++){
            sum+=cos(double(h+1)*ylen);
          }
          dp3d.potentin[i][j][k]=sum;
	}
      }
    }  
  } 
}


}
