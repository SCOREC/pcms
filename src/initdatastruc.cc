#include <coupling1.h>
#include <complex>
#include <dataprocessing.h>

namespace coupler {

void InitDatasProc3Dparameters(DatasProc3D& dp3d,Part1ParalPar3D& p1pp3d,Part3Mesh3D &p3m3d )
{
  if(para.prepro==true){
   if(dp3d.yparal==true){
     if(p1pp3d.li0%p1pp3d.npy==0){
       dp3d.part1li0=p1pp3d.li0/p1pp3d.npy;
       dp3d.part3li0=dp3d.part1li0;
     } else{
        if(p1pp3d.my_pey==p1pp3d.npy-1){
          dp3d.part1li0=p1pp3d.li0%p1pp3d.npy;
          dp3d.part3li0=dp3d.part1li0;
        }  else{
            dp3d.part1li0=p1pp3d.li0%p1pp3d.npy;
            dp3d.part3li0=dp3d.part1li0;
        }
     }
    dp3d.part1lj0=2*p1pp3d.nj0;
    dp3d.part3lj0=p1pp3d.nj0;   // here may need rethinking.
   } else{
     dp3d.part1li0=p1pp3d.li0;
     dp3d.part3li0=p1pp3d.li0;
     dp3d.part1lj0=2*p1pp3d.nj0;
     dp3d.part3lj0=p1pp3d.nj0;   // here may need rethinking.
   }
 }
  for(int i=0;i<p3m3d.li0;i++)  dp3d.sum=+p3m3d.myli0(i);
}

void AllocDatasProc3dDensityArraies(DatasProc3D& dp3d,Part1ParalPar3D& p1pp3d,Part3Mesh3D& p3m3d)
{
  dp3d.densin=new std::complex<double>**[p1pp3d.li0];
  for(int i=0;i<p1pp3d.li0;i++){
    dp3d.densin[i]=new std::complex<double>*[p1pp3d.lj0];
      for(int j=0;j<p1pp3d.lj0;j++)
        dp3d.densin[i][j]=new std::complex<double>[p1pp3d.lk0];
  }
  dp3d.densintmp=new std::complex<double>*[dp3d.part1li0];
  for(int k=0;k<dp3d.part1li0;k++){
      dp3d.densintmp[k]=new std::complex<double>[p1pp2d.nj0];
    }
  dp3d.densouttmp=new double*[dp3d.part1li0];
  for(int k=0;k<dp3d.part1li0;k++)
      dp3d.densouttmp[k]=new double[dp3d.part1lj0];

  dp3d.densout=new double**[p3m3d.li0];
  for(int i=0;i<p3m3d.li0;i++){
    dp3d.densout[i]=new double*[p3m3d.lj0]
      for(GO j=0;k<p3m3d.lj0;k++){
        dp3d.densout[i][j]=new double*[p1pp3d.lk0];
      }
  }

  dp3d.denspart3=new double**[p3m3d.li0];
  for(int i=0;i<p3m3d.li0;i++){
    dp3d.denspart3[i]=new double*[p3m3d.lj0]
    for(int j=0; j<p3m3d.lj0; j++)
      dp3d.denspart3[i][j]=new double*[p3m3d.mylk0[i]];
  }
} 

void AllocDatasProc3dPotentArraies(DatasProc3D& dp3d,Part1ParalPar3D& p1pp3d,Part3Mesh3D &p3m3d)
{
  dp3d.potentin=new double**[p3m3d.li0];
  for(int i=0;i<p3m3d.li0;i++){
    dp3d.potentin[i]=new double*[p3m3d.lj0];
      for(int j=0;j<p3m3d.lj0;j++)
        dp3d.densin[i][j]=new double[p3m3d.mylk0[i]];
  }
  dp3d.potenttmp=new double*[dp3d.sum];
  for(int k=0;k<dp3d.sum;k++){
      dp3d.potenttmp[k]=new double[p3m3d.lj0];
    }
  dp3d.potentout=new std:complex<double>**[p3d3d.li0];
  for(int i=0;i<p3m3d.li0;i++){
    dp3d.potentout[i]=new std:complex<double>*[dp3d.part3lj0];
      for(int j=0;j<p3d3d.lj0;j++)
         dp3d.potentout[i][j]=new std:complex<double>[p3m3d.mylk0[i]];
  }

  dp3d.potentpart1=new std:complex<double>**[p1pp3d.li0];
  for(int i=0;i<p1pp3d.li0;i++){
    dp3d.potentpart1[i]=new std::complex<double>*[p1pp3d.lj0];
    for(int j=0;j<p1pp3d.lj0;j++){
      dp3d.potentpart1[i][j]=new std::complex<double>[p1pp3d.lk0];
    }
  }

}

}


