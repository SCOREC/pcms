#include "dataprocess.h"

namespace coupler {

void InitDatasProc3Dparameters(DatasProc3D& dp3d,Part1ParalPar3D& p1pp3d,Part3Mesh3D &p3m3d )
{
  if(preproc==true){
    if(dp3d.yparal==true){
      if(p1pp3d.li0%p1pp3d.npy==0){
        dp3d.part1li0=p1pp3d.li0/p1pp3d.npy;
        dp3d.part3li0=dp3d.part1li0;
      } else{
        if(p1pp3d.mype_y==p1pp3d.npy-1){
          dp3d.part1li0=p1pp3d.li0%p1pp3d.npy;
          dp3d.part3li0=dp3d.part1li0;
        }  else{
            dp3d.part1li0=p1pp3d.li0%p1pp3d.npy;
            dp3d.part3li0=dp3d.part1li0;
        }
     }
    dp3d.part1lj0=2*p1pp3d.ny0;
    dp3d.part3lj0=p1pp3d.ny0;   // here may need rethinking.
   } else{
     dp3d.part1li0=p1pp3d.li0;
     dp3d.part3li0=p1pp3d.li0;
     dp3d.part1lj0=2*p1pp3d.ny0;
     dp3d.part3lj0=p1pp3d.ny0;   // here may need rethinking.
   }
 }
  for(LO i=0;i<p3m3d.li0;i++)  dp3d.sum+=p3m3d.mylk0[i];
}

void AllocDatasProc3dDensityArrays(DatasProc3D& dp3d,Part1ParalPar3D& p1pp3d,Part3Mesh3D& p3m3d)
{
 if(dp3d.yparal==false){
  dp3d.densin=new std::complex<double>**[p1pp3d.li0];
  for(LO i=0;i<p1pp3d.li0;i++){
    dp3d.densin[i]=new std::complex<double>*[p1pp3d.lj0];
    for(LO j=0;j<p1pp3d.lj0;j++)
      dp3d.densin[i][j]=new std::complex<double>[p1pp3d.lk0];
  }
  GO num=p1pp3d.li0*p1pp3d.lj0*p1pp3d.lk0;
  dp3d.densintmp=new std::complex<double>[num];

  num=p3m3d.li0*p3m3d.lj0*p1pp3d.lk0;
  dp3d.densouttmp=new double[num];

  dp3d.densout=new double**[p1pp3d.li0];
  for(LO i=0;i<p3m3d.li0;i++){
    dp3d.densout[i]=new double*[p1pp3d.lj0*2];
    for(GO j=0;j<p3m3d.lj0;j++){
      dp3d.densout[i][j]=new double[p1pp3d.lk0];
      }
  }

  dp3d.denspart3=new double**[p3m3d.li0];
  for(LO i=0;i<p3m3d.li0;i++){
    dp3d.denspart3[i]=new double*[p3m3d.lj0];
    for(LO j=0; j<p3m3d.lj0; j++)
      dp3d.denspart3[i][j]=new double[p3m3d.mylk0[i]];
  }
} 
}

void AllocDatasProc3dPotentArrays(DatasProc3D& dp3d,Part1ParalPar3D& p1pp3d,Part3Mesh3D& p3m3d)
{
 if(dp3d.yparal==false){
  dp3d.potentin=new double**[p3m3d.li0];
  for(LO i=0;i<p3m3d.li0;i++){
    dp3d.potentin[i]=new double*[p3m3d.lj0];
    for(LO j=0;j<p3m3d.lj0;j++)
      dp3d.potentin[i][j]=new double[p3m3d.mylk0[i]];
  }
  dp3d.potentintmp=new double[dp3d.sum*p3m3d.lj0];
/*  for(LO k=0;k<dp3d.sum;k++){
      dp3d.potenttmp[k]=new double[p3m3d.lj0];
    }
*/
  dp3d.potentouttmp=new std::complex<double>[dp3d.sum*p3m3d.lj0/2];
  dp3d.potentout=new std::complex<double>**[p3m3d.li0];
  for(LO i=0;i<p3m3d.li0;i++){
    dp3d.potentout[i]=new std::complex<double>*[dp3d.part3lj0];
    for(LO j=0;j<p3m3d.lj0;j++)
      dp3d.potentout[i][j]=new std::complex<double>[p3m3d.mylk0[i]];
  }

  dp3d.potentpart1=new std::complex<double>**[p1pp3d.li0];
  for(LO i=0;i<p1pp3d.li0;i++){
    dp3d.potentpart1[i]=new std::complex<double>*[p1pp3d.lj0];
    for(LO j=0;j<p1pp3d.lj0;j++){
      dp3d.potentpart1[i][j]=new std::complex<double>[p1pp3d.lk0];
    }
  }
 }
}

DatasProc3D::~DatasProc3D()
{
  if(densin!=NULL) delete[] densin;
  if(densintmp!=NULL) delete[] densintmp;
  if(densouttmp!=NULL) delete[] densouttmp;
  if(densout!=NULL) delete[] densout;
  if(denspart3!=NULL) delete[] denspart3;
  if(potentin!=NULL) delete[] potentin;
  if(potentintmp!=NULL) delete[] potentintmp;
  if(potentouttmp!=NULL) delete[] potentouttmp;
  if(potentout!=NULL) delete[] potentout;
  if(potentpart1!=NULL) delete[] potentpart1;       
}


}


