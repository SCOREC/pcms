#include "inittestenv.h"

namespace coupler{

void TestInitPotentAlongz(DatasProc3D& dp3d, const Part3Mesh3D& p3m3d, const LO npy, const LO n)
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
            sum+=cos(double(h+1)*ylen-cplPI);
          }
          dp3d.potentin[i][j][k]=sum;
        }
      }
    }  
  } 
}

//void testInitFourierPlan()
   
}
