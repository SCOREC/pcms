#include "BoundaryDescr3D.h"
#include "couplingTypes.h"
#include "importpart3mesh.h"
#include "dataprocess.h"
#include "commpart1.h"
#include "sendrecv_impl.h"
#include <mpi.h>

namespace coupler {

BoundaryDescr3D::BoundaryDescr3D(
    const Part3Mesh3D& p3m3d,
    const Part1ParalPar3D &p1pp3d,
    const DatasProc3D& dp3d)
{
  nzb=p1pp3d.nzb;
  updenz=new double**[p1pp3d.li0];
  lowdenz=new double**[p1pp3d.li0];
  for(LO i=0;i<p1pp3d.li0;i++){
    updenz[i]=new double*[dp3d.part1lj0];
    lowdenz[i]=new double*[dp3d.part1lj0];
    for(LO j=0;j<dp3d.part1lj0;j++){
      updenz[i][j]=new double[nzb];
      lowdenz[i][j]=new double[nzb];
    }
  } 
  uppotentz=new CV**[p3m3d.xboxinds[0][p1pp3d.mype_x]];
  lowpotentz=new CV**[p3m3d.xboxinds[0][p1pp3d.mype_x]];
  upzpart3=new double*[p3m3d.xboxinds[0][p1pp3d.mype_x]];
  lowzpart3=new double*[p3m3d.xboxinds[0][p1pp3d.mype_x]];
  for(LO i=0;i<p3m3d.xboxinds[0][p1pp3d.mype_x];i++){
    uppotentz[i]=new CV*[p3m3d.lj0];
    lowpotentz[i]=new CV*[p3m3d.lj0]; 
    upzpart3[i]=new double[nzb];
    lowzpart3[i]=new double[nzb];
    for(LO j=0;j<dp3d.part3lj0;j++){
      uppotentz[i][j]=new CV[nzb];
      lowpotentz[i][j]=new CV[nzb];
    }
  }
}

void BoundaryDescr3D::zPotentBoundaryBufAssign(
    const DatasProc3D& dp3d, 
    const Part3Mesh3D& p3m3d,
    const Part1ParalPar3D &p1pp3d)
{
  if(lowpotentz==NULL||uppotentz==NULL){
    std::cout<<"ERROR:the boundary buffer of the potential must be allocated beforing invoking this routine.";
    std::exit(EXIT_FAILURE);
  }
  LO li0,lj0,lk0;
  li0=p3m3d.xboxinds[0][p1pp3d.mype_x];
  lj0=p3m3d.lj0;
  if(p1pp3d.npz>1){
    for(LO i=0;i<li0;i++){
      lk0=p3m3d.mylk0[i];
      if(lk0<nzb){
         std::cout<<"ERROR: the interpolation order is larger than the box count along z dimension.";
         std::exit(EXIT_FAILURE);
       } 
      if(p1pp3d.periods[2]==1){ 
        mpisendrecv_aux1D(p1pp3d.comm_z,nzb,li0,lj0,lk0,lowzpart3[i],upzpart3[i],
          p3m3d.pzcoords[i]); //double
        if(p1pp3d.comm_z==0) 
        for(LO j=0;j<lj0;j++){
          mpisendrecv_aux1D(p1pp3d.comm_z,nzb,li0,lj0,lk0,lowpotentz[i][j],uppotentz[i][j],
              dp3d.potentout[i][j]); //complex
        }
      } else {
         std::cout<<"The topology is not right for the parallel domain."<<'\n';
         std::exit(EXIT_FAILURE);
      }
    }  
  } else{
      if(p1pp3d.periods[2]==1){
        for(LO i=0;i<li0;i++){
          lk0=p3m3d.mylk0[i];
          if(lk0<nzb){
            std::cout<<"ERROR: the LOerpolation order is larger than the box count along z dimension.";
            std::exit(EXIT_FAILURE);
          }  
          for(LO k=0;k<nzb-1;k++){
            lowzpart3[i][k]=p3m3d.pzcoords[i][lk0-nzb+k];
            upzpart3[i][k]=p3m3d.pzcoords[i][k];
          }
          for(LO j=0;j<lj0;j++){
            for(LO k=0;k<nzb;k++){
              lowpotentz[i][j][k]=dp3d.potentout[i][j][k];
              lowpotentz[i][j][k]=dp3d.potentout[i][j][lk0-nzb+k];
            }  
         }     
       }
     } else {
         std::cout<<"The topology is not right for serial y domain."<<'\n';
         std::exit(EXIT_FAILURE);
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
 
