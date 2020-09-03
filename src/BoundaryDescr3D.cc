#include "BoundaryDescr3D.h"
#include "couplingTypes.h"
#include "importpart3mesh.h"
#include "dataprocess.h"
#include "commpart1.h"
#include "testutilities.h"
#include "sendrecv_impl.h"
#include <mpi.h>
#include <math.h>

namespace coupler {

BoundaryDescr3D::BoundaryDescr3D(
    const Part3Mesh3D& p3m3d,
    const Part1ParalPar3D &p1pp3d,
    const TestCase tcase,
    bool pproc):test_case(tcase), preproc(pproc)
{
  if(preproc==true){
    nzb=p1pp3d.nzb;
    updenz=new CV**[p1pp3d.li0[p1pp3d.mype_x]];
    lowdenz=new CV**[p1pp3d.li0[p1pp3d.mype_x]];
    for(LO i=0;i<p1pp3d.li0[p1pp3d.mype_x];i++){
      updenz[i]=new CV*[p1pp3d.lj0];
      lowdenz[i]=new CV*[p1pp3d.lj0];
      for(LO j=0;j<p1pp3d.lj0;j++){
	updenz[i][j]=new CV[nzb];
	lowdenz[i][j]=new CV[nzb];
      }
    }

    if(p1pp3d.mype_z==0){
      if(lowpbmat==NULL){
        lowpbmat=new CV*[p1pp3d.li0[p1pp3d.mype_x]];
        for(LO i=0;i<p1pp3d.li0[p1pp3d.mype_x];i++){
          lowpbmat[i]=new CV[p1pp3d.lj0];      
        } 
      }
     }else if(p1pp3d.mype_z==p1pp3d.npz-1){
      if(uppbmat==NULL){
        uppbmat=new CV*[p1pp3d.li0[p1pp3d.mype_x]];
        for(LO i=0;i<p1pp3d.li0[p1pp3d.mype_x];i++){
          uppbmat[i]=new CV[p1pp3d.lj0]; 
        }
       }
     }

    uppotentz=new CV**[p3m3d.xboxinds[p1pp3d.mype_x][0]];
    lowpotentz=new CV**[p3m3d.xboxinds[p1pp3d.mype_x][0]];
    upzpart3=new double*[p3m3d.xboxinds[p1pp3d.mype_x][0]];
    lowzpart3=new double*[p3m3d.xboxinds[p1pp3d.mype_x][0]];
    for(LO i=0;i<p3m3d.xboxinds[p1pp3d.mype_x][0];i++){
      uppotentz[i]=new CV*[p3m3d.lj0/2];
      lowpotentz[i]=new CV*[p3m3d.lj0/2]; 
      upzpart3[i]=new double[nzb];
      lowzpart3[i]=new double[nzb];
      for(LO j=0;j<p3m3d.lj0/2;j++){
	uppotentz[i][j]=new CV[nzb];
	lowpotentz[i][j]=new CV[nzb];
      }
    }    

   initpbmat(p1pp3d);
  } 
}

void BoundaryDescr3D::initpbmat(const Part1ParalPar3D &p1pp3d)
{
   LO num;
   if(p1pp3d.mype_z==0){
     for(LO i=0;i<p1pp3d.li0[p1pp3d.mype_x];i++){
       num=p1pp3d.li1+i;
       for(LO j=0;j<p1pp3d.lj0;j++){  
	 lowpbmat[i][j]=exp(CV(0.0,1.0)*2.0*cplPI*double(p1pp3d.n0_global)
                        *double(j+p1pp3d.ky0_ind)*p1pp3d.q_prof[num]);
       }
     } 
   } else if(p1pp3d.mype_z==p1pp3d.npz-1){
     for(LO i=0;i<p1pp3d.li0[p1pp3d.mype_x];i++){
       num=p1pp3d.li1+i;
       for(LO j=0;j<p1pp3d.lj0;j++){
         uppbmat[i][j]=exp(-CV(0.0,1.0)*2.0*cplPI*double(p1pp3d.n0_global)
                        *double(j+p1pp3d.ky0_ind)*p1pp3d.q_prof[num]);
       }
     }
   }
}


BoundaryDescr3D::~BoundaryDescr3D()
{
  if(upzpart3!=NULL) delete[] upzpart3;
  if(lowzpart3!=NULL) delete[] lowzpart3;
  if(updenz!=NULL) delete[] updenz;
  if(lowdenz!=NULL) delete[] lowdenz;
  if(uppotentz!=NULL) delete[] uppotentz;
  if(lowpotentz!=NULL) delete[] lowpotentz;
}



}
 
