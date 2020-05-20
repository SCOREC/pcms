#include "dataprocess.h"
#include "importpart3mesh.h"
#include "commpart1.h"

namespace coupler {

DatasProc3D::DatasProc3D(const Part1ParalPar3D& p1pp3d,
    const Part3Mesh3D &p3m3d,
    bool pproc,
    bool ypar)
  : preproc(pproc),
    yparal(ypar),
    p1(p1pp3d.li0,p1pp3d.lj0,p1pp3d.lk0,
	 p1pp3d.ny0, p1pp3d.npy, p1pp3d.mype_y, 
         p1pp3d.res_fact),
    p3(p3m3d.li0,p3m3d.lj0,p3m3d.mylk0)
  {
    init();
    AllocDensityArrays();
    AllocPotentArrays();
}


void DatasProc3D::init()
{
  if(preproc==true){
    if(yparal==true){
      if(p1.li0%p1.npy==0){
        part1li0=p1.li0/p1.npy;
        part3li0=part1li0;
      } else{
        if(p1.mype_y==p1.npy-1){
          part1li0=p1.li0%p1.npy;
          part3li0=part1li0;
        }  else{
          part1li0=p1.li0%p1.npy;
          part3li0=part1li0;
        }
      }
      part1lj0=2*p1.ny0;
      part3lj0=p1.ny0;   // here may need rethinking.
    } else{
      part1li0=p1.li0;
      part3li0=p1.li0;
      part1lj0=2*p1.ny0;
      part3lj0=p1.ny0;   // here may need rethinking.
    }
  }
  sum=0;
  for(LO i=0;i<p3.li0;i++)  sum+=p3.mylk0[i];
}

void DatasProc3D::AllocDensityArrays()
{
  if(yparal==false){
    densin=new CV**[p1.li0];
    for(LO i=0;i<p1.li0;i++){
      densin[i]=new CV*[p1.lj0];
      for(LO j=0;j<p1.lj0;j++)
        densin[i][j]=new CV[p1.lk0];
    }
    densintmp=new CV[p1.lj0];
    densouttmp=new double[p1.lj0*2];

    densout=new double**[p1.li0];
    for(LO i=0;i<p3.li0;i++){
      densout[i]=new double*[p1.lj0*2];
      for(GO j=0;j<p3.lj0;j++){
        densout[i][j]=new double[p1.lk0];
      }
    }

    denspart3=new double**[p3.li0];
    for(LO i=0;i<p3.li0;i++){
      denspart3[i]=new double*[p3.lj0];
      for(LO j=0; j<p3.lj0; j++)
        denspart3[i][j]=new double[p3.mylk0[i]];
    }
  } 
}

void DatasProc3D::AllocPotentArrays()
{ 
  if(yparal==false){
    potentin=new double**[p3.li0];
    for(LO i=0;i<p3.li0;i++){
      potentin[i]=new double*[p3.lj0];
      for(LO j=0;j<p3.lj0;j++)
        potentin[i][j]=new double[p3.mylk0[i]];
    }

    potentinterpo=new double**[p3.lj0];
    for(LO i=0;i<p3.li0;i++){
      potentinterpo[i]=new double*[p3.lj0];
      for(LO j=0;j<p3.lj0;j++)
        potentinterpo[i][j]=new double[p1.lk0];
    }

    potentintmp=new double[p3.lj0];
    potentouttmp=new CV[p3.lj0/2+1];
   
    potentpart1=new CV**[p1.li0];
    for(LO i=0;i<p1.li0;i++){
      potentpart1[i]=new CV*[p1.lj0];
      for(LO j=0;j<p1.lj0;j++){
        potentpart1[i][j]=new CV[p1.lk0];
      }
    }
  }
}


DatasProc3D::~DatasProc3D()
{
  FreeFourierPlan3D();
  if(densin!=NULL) delete[] densin;
  if(densintmp!=NULL) delete[] densintmp;
  if(densouttmp!=NULL) delete[] densouttmp;
  if(densout!=NULL) delete[] densout;
  if(denspart3!=NULL) delete[] denspart3;
  if(potentin!=NULL) delete[] potentin;
  if(potentouttmp!=NULL) delete[] potentouttmp;
  if(potentinterpo!=NULL) delete[] potentinterpo;
  if(potentpart1!=NULL) delete[] potentpart1;       
}

}

