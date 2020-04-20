#include <coupling1.h>
#include <dataprocess.h>

namespace coupler {

void InterpoDensity3D(BoundaryDescr3D &bdesc,Part3Mesh3D& p3m3d, Part1ParalPar3D &p1pp3d,DatasProc3D& dp3d)
{
  double* yin,xin,yout,xout;
  xin=new double[p1pp3d.lz0+2*nzb];
  yin=new double[p1pp2s.lz0+2*nzb];
  GO nzb=bdesc.nzb;  
  for(int l=0;l<nzb-1;l++){
    xin[l]=p1pp3d.pzp[0]-double(nzb-l)*p1pp3d.dz;
    xin[p1pp3d.lz0+nzb+l]=p1pp3d.pzp[p1pp3d.lz0-1]+double(l+1)*p1pp3d.dz;
  }
  for(int k=0;k<p1pp3d.lz0-1;k++){  
    xin[nzb+k]=p1pp3d.pzp[k];
  }   
  if(prepro==true){
    for(int i=0;i<p3m3d.li0;i++){
      for(int j=0;j<p3m3d.lj0;j++){
        for(int l=0;l<nzb-1;l++){
          yin[l]=bdesc.lowdenz[i][j][l];
          yin[p1pp3d.lz0+nzb+l]=bdesc.updenz[i][j][l];
        }
        for(int k=0;k<p1pp3d.lz0-1;k++){  
          yin[nzb+k]=dp3d.densout[i][j][k];
          
        }
        xout=p3m3d.pzcoords[i];
        yout=dp3ddp3d.denspart3[i][j]; 
        Lag3dArray(yin,xin,p1pp2s.lz0+2*nzb,yout,xout,p3m3d.mylk0[i]);
      }
    }   
  }
  delete[] xin,yin,yout,xout;
}


void InterpoPotential3D(BoundaryDescr3D &bdesc,Part3Mesh3D& p3m3d, Part1ParalPar3D &p1pp3d,DatasProc3D& dp3d)
{
  std::complex<double>* yin,yout;
  double* xin,xout;
  yout=std::complex<double>[p1pp3d.lz0];
  GO nzb=bdesc.nzb;
  if(prepro==true){
    for(int i=0;i<p3m3d.li0;i++){
      yin=new std::complex<double>[p3m3d.mylk0[i]+2*nzb];
      xin=new double[p3mp3d.mylk0[i]+2*nzb];
      for(int l=0;l<nzb-1;l++){  
        xin[l]=bdesc.lowzpart3[l];
        xin[p3m3d.mylk0[i]+nzb+l]=bdesc.upzpart3[l];      
      }
      for(int j=0;j<p3m3d.lj0;j++){
        for(int l=0;l<nzb-1;l++){
          yin[l]=bdesc.lowpotentz[i][j][l];
          yin[p3m3d.mylk0[i]+nzb+l]=bdesc.uppotentz[i][j][l];
        }
        for(int k=0;k<p3m3d.mylk0[i]-1;k++){
          xin[k]=p3m3d.pzcoords[i][k];
          yin[k]=dp3d.potentin[i][j][k];
        }
        xout=p1pp3d.pzp;
        yout=dp3d.potentout[i][j];
        Lag3dArray(yin,xin,p3m3d.mylk0[i]+2*nzb,yout,xout,pipp3d.lk0);
     }
  
    }
  }
 delete[] xin,yin,xout,yout;
}


}
