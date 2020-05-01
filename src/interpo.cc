#include "couplingTypes.h"
#include "dataprocess.h"
#include "importpart3mesh.h"
#include "commpart1.h"
#include "BoundaryDescr3D.h"

namespace coupler {

// routines for interpolation
template<class T>
T Lag3dInterpo1D(const T yin[4],const double xin[4],const double x)
{
  double l0,l1,l2,l3;
  T yout;
  l0=(x-xin[1])*(x-xin[2])*(x-xin[3])/(xin[0]-xin[1])/(xin[0]-xin[2])/(xin[0]-xin[3]);
  l1=(x-xin[0])*(x-xin[2])*(x-xin[3])/(xin[1]-xin[0])/(xin[1]-xin[2])/(xin[1]-xin[3]);
  l2=(x-xin[0])*(x-xin[1])*(x-xin[3])/(xin[2]-xin[0])/(xin[2]-xin[1])/(xin[2]-xin[3]);
  l3=(x-xin[0])*(x-xin[1])*(x-xin[2])/(xin[3]-xin[0])/(xin[3]-xin[1])/(xin[3]-xin[2]);
  yout=yin[0]*l0+yin[1]*l1+yin[2]*l2+yin[3]*l3;
  return yout;
}

//central 3rd order Lagrange polynormal interpolation
template<class T>
void Lag3dArray(T* yin,double* xin,LO nin,T* yout,double* xout,LO nout){
  LO jstart=2;
  LO j1=jstart;
  LO j2,j0,jm;
  double x;
  T func[4];
  double coords[4];
  for(LO j=0;j<nout;j++){
    x=xout[j];
    while(x>=xin[j1] && j1<nin-2 && j1>1){
      j1=+1;
    }
    j2=j1+1;
    j0=j1-1;
    jm=j1-2;
    coords[0]=xin[jm];
    coords[1]=xin[j0];
    coords[2]=xin[j1];
    coords[3]=xin[j2];
    func[0]=yin[jm];
    func[1]=yin[j0];
    func[2]=yin[j1];
    func[3]=yin[j2];
    yout[j]=Lag3dInterpo1D(func,coords,x);
  }
}

//FIXME What does this function output - none of the input
//FIXME class member variables are modified.
//FIXME Added optional arg for preproc - once this is in
//FIXME a class that can be removed.
void InterpoDensity3D(const BoundaryDescr3D &bdesc,
    const Part3Mesh3D& p3m3d, 
    const Part1ParalPar3D &p1pp3d,
    const DatasProc3D& dp3d,
    bool preproc = true)
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
  if(preproc==true){
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


//FIXME What does this function output - none of the input
//FIXME class member variables are modified.
//FIXME Added optional arg for preproc - once this is in
//FIXME a class that can be removed.
void InterpoPotential3D(const BoundaryDescr3D &bdesc,
    const Part3Mesh3D& p3m3d,
    const Part1ParalPar3D &p1pp3d,
    const DatasProc3D& dp3d,
    bool preproc = true)
{
  std::complex<double>* yin;
  std::complex<double>* yout;
  double* xin;
  double* xout;
  yout=new std::complex<double>[p1pp3d.lk0];
  LO nzb=bdesc.nzb;
  if(preproc==true){
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
