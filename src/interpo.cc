#include "dataprocess.h"

namespace coupler {

void InterpoDensity3D(BoundaryDescr3D &bdesc,Part3Mesh3D& p3m3d, Part1ParalPar3D &p1pp3d,DatasProc3D& dp3d)
{
  double* yin;
  double* xin;
  double* yout;
  double* xout;
  LO nzb=bdesc.nzb;
  xin=new double[p1pp3d.lk0+2*nzb];
  yin=new double[p1pp3d.lk0+2*nzb];  
  for(LO l=0;l<nzb;l++){
    xin[l]=p1pp3d.pzp[0]-double(nzb-l)*p1pp3d.dz;
    xin[p1pp3d.lk0+nzb+l]=p1pp3d.pzp[p1pp3d.lk0-1]+double(l+1)*p1pp3d.dz;
  }
  for(LO k=0;k<p1pp3d.lk0-1;k++){  
    xin[nzb+k]=p1pp3d.pzp[k];
  }   
  if(p1pp3d.preproc==true){
    for(LO i=0;i<p3m3d.li0;i++){
      for(LO j=0;j<p3m3d.lj0;j++){
        for(LO l=0;l<nzb-1;l++){
          yin[l]=bdesc.lowdenz[i][j][l];
          yin[p1pp3d.lk0+nzb+l]=bdesc.updenz[i][j][l];
        }
        for(LO k=0;k<p1pp3d.lk0-1;k++){  
          yin[nzb+k]=dp3d.densout[i][j][k];
          
        }
        xout=p3m3d.pzcoords[i];
        yout=dp3d.denspart3[i][j]; 
        Lag3dArray(yin,xin,p1pp3d.lk0+2*nzb,yout,xout,p3m3d.mylk0[i]);
      }
    }   
  }
  delete[] xin,yin,yout,xout;
}


void InterpoPotential3D(BoundaryDescr3D &bdesc,Part3Mesh3D& p3m3d, Part1ParalPar3D &p1pp3d,DatasProc3D& dp3d)
{
  std::complex<double>* yin;
  std::complex<double>* yout;
  double* xin;
  double* xout;
  yout=new std::complex<double>[p1pp3d.lk0];
  LO nzb=bdesc.nzb;
  if(p1pp3d.preproc==true){
    for(LO i=0;i<p3m3d.li0;i++){
      yin=new std::complex<double>[p3m3d.mylk0[i]+2*nzb];
      xin=new double[p3m3d.mylk0[i]+2*nzb];
      for(LO l=0;l<nzb-1;l++){  
        xin[l]=bdesc.lowzpart3[i][l];
        xin[p3m3d.mylk0[i]+nzb+l]=bdesc.upzpart3[i][l];      
      }
      for(LO j=0;j<p3m3d.lj0;j++){
        for(LO l=0;l<nzb-1;l++){
          yin[l]=bdesc.lowpotentz[i][j][l];
          yin[p3m3d.mylk0[i]+nzb+l]=bdesc.uppotentz[i][j][l];
        }
        for(LO k=0;k<p3m3d.mylk0[i]-1;k++){
          xin[k]=p3m3d.pzcoords[i][k];
          yin[k]=dp3d.potentin[i][j][k];
        }
        xout=p1pp3d.pzp;
        yout=dp3d.potentout[i][j];
        Lag3dArray(yin,xin,p3m3d.mylk0[i]+2*nzb,yout,xout,p1pp3d.lk0);
     }
     delete[] xin,yin; 
    }
  }
 delete[] yout;
}


}
