#include "dataprocess.h"

namespace coupler {

void InitBoundaryDescr3D(BoundaryDescr3D &bdesc,Part3Mesh3D& p3m3d, Part1ParalPar3D &p1pp3d,DatasProc3D& dp3d)  
{
  bdesc.nzb=p1pp3d.nzb;
  bdesc.updenz=new double**[p1pp3d.li0];
  bdesc.lowdenz=new double**[p1pp3d.li0];
  for(LO i=0;i<p1pp3d.li0;i++){
    bdesc.updenz[i]=new double*[dp3d.part1lj0];
    bdesc.lowdenz[i]=new double*[dp3d.part1lj0];
    for(LO j=0;j<dp3d.part1lj0;j++){
      bdesc.updenz[i][j]=new double[bdesc.nzb];
      bdesc.lowdenz[i][j]=new double[bdesc.nzb];
    }
  } 
  bdesc.uppotentz=new std::complex<double>**[p3m3d.xboxinds[0][p1pp3d.mype_x]];
  bdesc.lowpotentz=new std::complex<double>**[p3m3d.xboxinds[0][p1pp3d.mype_x]];
  bdesc.upzpart3=new double*[p3m3d.xboxinds[0][p1pp3d.mype_x]];
  bdesc.lowzpart3=new double*[p3m3d.xboxinds[0][p1pp3d.mype_x]];
  for(LO i=0;i<p3m3d.xboxinds[0][p1pp3d.mype_x];i++){
    bdesc.uppotentz[i]=new std::complex<double>*[p3m3d.lj0];
    bdesc.lowpotentz[i]=new std::complex<double>*[p3m3d.lj0]; 
    bdesc.upzpart3[i]=new double[bdesc.nzb];
    bdesc.lowzpart3[i]=new double[bdesc.nzb];
    for(LO j=0;j<dp3d.part3lj0;j++){
      bdesc.uppotentz[i][j]=new std::complex<double>[bdesc.nzb];
      bdesc.lowpotentz[i][j]=new std::complex<double>[bdesc.nzb];
    }
  }
}

void zPotentBoundaryBufAssign(MPI_Datatype mpitype,BoundaryDescr3D &bdesc,DatasProc3D& dp3d, Part3Mesh3D& p3m3d,\
     Part1ParalPar3D &p1pp3d)
{
  if(bdesc.lowpotentz==NULL||bdesc.uppotentz==NULL){
    std::cout<<"ERROR:the boundary buffer of the potential must be allocated beforing invoking this routine.";
    std::exit(EXIT_FAILURE);
  }
  LO li0,lj0,lk0,nzb;
  li0=p3m3d.xboxinds[0][p1pp3d.mype_x];
  lj0=p3m3d.lj0;
  nzb=bdesc.nzb;
  if(p1pp3d.npz>1){
    for(LO i=0;i<li0;i++){
      lk0=p3m3d.mylk0[i];
      if(lk0<nzb){
         std::cout<<"ERROR: the interpolation order is larger than the box count along z dimension.";
         std::exit(EXIT_FAILURE);
       } 
      if(p1pp3d.periods[2]==1){ 
        mpisendrecv_aux1D(mpitype,p1pp3d.comm_z,nzb,li0,lj0,lk0,bdesc.lowzpart3[i],bdesc.upzpart3[i], \
          p3m3d.pzcoords[i]); 
        if(p1pp3d.comm_z==0) 
        for(LO j=0;j<lj0;j++){
          mpisendrecv_aux1D(mpitype,p1pp3d.comm_z,nzb,li0,lj0,lk0,bdesc.lowpotentz[i][j],bdesc.uppotentz[i][j],\
              dp3d.potentout[i][j]); 
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
            bdesc.lowzpart3[i][k]=p3m3d.pzcoords[i][lk0-nzb+k];
            bdesc.upzpart3[i][k]=p3m3d.pzcoords[i][k];
          }
          for(LO j=0;j<lj0;j++){
            for(LO k=0;k<nzb;k++){
              bdesc.lowpotentz[i][j][k]=dp3d.potentout[i][j][k];
              bdesc.lowpotentz[i][j][k]=dp3d.potentout[i][j][lk0-nzb+k];
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
 
