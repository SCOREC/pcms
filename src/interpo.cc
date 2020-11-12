#include "couplingTypes.h"
#include "dataprocess.h"
#include "interpoutil.h"
#include "importpart3mesh.h"
#include "commpart1.h"
#include "BoundaryDescr3D.h"

namespace coupler {

void DatasProc3D::mesh1dforDensityInterpo()
{
   LO nzb=bdesc->nzb;
   if(p1->periods[2]==1){
     if(p1->mype_z==0){
       for(LO l=0;l<nzb;l++){
	 mesh1ddens[l]=p1->pzcoords[p1->nz0-nzb+l]-2.0*cplPI;
	 mesh1ddens[p1->lk0+nzb+l]=p1->pzcoords[p1->lk2+l+1];          
       }
     } else if(p1->mype_z==p1->npz-1){
	 for(LO l=0;l<nzb;l++){
	   mesh1ddens[l]=p1->pzcoords[p1->lk1-nzb+l];
	   mesh1ddens[p1->lk0+nzb+l]=p1->pzcoords[l]+2.0*cplPI;
	 }        
     }else{
	 for(LO l=0;l<nzb;l++){
	   mesh1ddens[l]=p1->pzcoords[p1->lk1-nzb+l];
	   mesh1ddens[p1->lk0+nzb+l]=p1->pzcoords[p1->lk2+l+1];
	 }
     }
     for(LO k=0;k<p1->lk0;k++){  
       mesh1ddens[nzb+k]=p1->pzcoords[p1->lk1+k];
     }      
   }

}

void DatasProc3D::InterpoDensity3D()
{
  CV* yin;
  CV* yout;
  double* xout;
  LO nzb=bdesc->nzb;
  yin=new CV[p1->lk0+2*nzb];  
  if(preproc==true){
   for(LO i=0;i<p1->li0;i++){
      for(LO j=0;j<p1->lj0;j++){
        for(LO l=0;l<nzb;l++){
          yin[l]=bdesc->lowdenz[i][j][l];
          yin[p1->lk0+nzb+l]=bdesc->updenz[i][j][l];
        }
        for(LO k=0;k<p1->lk0;k++){  
          yin[nzb+k]=densin[i][j][k];
        }
        xout=p3->pzcoords[i];
        yout=densinterpo[i][j];
        Lag3dArray(yin,mesh1ddens,p1->lk0+2*nzb,yout,xout,p3->mylk0[i]);
      }
    }   
  }
  delete[] yin;
  yin=NULL;
}

void DatasProc3D::mesh1dforPotentialInterpo()
{
  LO nzb=bdesc->nzb;
  for(LO i=0;i<p3->li0;i++){
    for(LO l=0;l<nzb;l++){
      mesh1dpotent[i][l]=bdesc->lowzpart3[i][l];
      mesh1dpotent[i][p3->mylk0[i]+nzb+l]=bdesc->upzpart3[i][l];
    }

    for(LO j=0;j<p3->lj0/2;j++){
      for(LO k=0;k<p3->mylk0[i];k++){
	mesh1dpotent[i][k+nzb]=p3->pzcoords[i][k];
      }
    }
  }
}


void DatasProc3D::InterpoPotential3D()
{
  CV* yin;
  CV* yout;
  LO nzb=bdesc->nzb;
  if(preproc==true){
    for(LO i=0;i<p3->li0;i++){
      yin=new CV[p3->mylk0[i]+2*nzb];
      for(LO j=0;j<p3->lj0/2;j++){
        for(LO l=0;l<nzb;l++){
          yin[l]=bdesc->lowpotentz[i][j][l];
          yin[p3->mylk0[i]+nzb+l]=bdesc->uppotentz[i][j][l];
        }
        for(LO k=0;k<p3->mylk0[i];k++){
          yin[k+nzb]=potentinterpo[i][j][k];
        }
        yout=potentpart1[i][j];
        Lag3dArray(yin,mesh1dpotent[i],p3->mylk0[i]+2*nzb,yout,p1->pzp,p1->lk0);
     }

     delete[] yin; 
    }
  }
}


void gemXgcDatasProc3D::interpoDensityAlongZ(double*** box)
{
  double* yin;
  double* yout;
  double* xout;
  LO nzb=bdesc->nzb;
  yin=new double[p1->lk0+2*nzb];
  if(preproc==true){
   printf("reach interpo\n");
   MPI_Barrier(MPI_COMM_WORLD);  
   for(LO i=0;i<p1->li0;i++){
      for(LO j=0;j<p1->lj0;j++){
//        if(p1->mype==0) printf("i1: %d, j1: %d \n", i, j);
//        MPI_Barrier(MPI_COMM_WORLD);
        for(LO l=0;l<nzb;l++){
          yin[l]=bdesc->lowdenzgemxgc[i][j][l];
          yin[p1->lk0+nzb+l]=bdesc->updenzgemxgc[i][j][l];
        }
	
//        if(p1->mype==0) printf("i2: %d, j2: %d \n", i, j);
//        MPI_Barrier(MPI_COMM_WORLD);
        for(LO k=0;k<p1->lk0;k++){
          yin[nzb+k]=densin[i][j][k];
        }
//        if(p1->mype==0) printf("i4: %d, j4: %d \n", i, j);
//        MPI_Barrier(MPI_COMM_WORLD);
        xout=p3->theta_geo[i];
        yout=box[i][j];
        Lag3dArray(yin,bdesc->thetameshgem,p1->lk0+2*nzb,yout,xout,p3->mylk0[i]);
//        if(p1->mype==0) printf("i5: %d, j5: %d \n", i, j);
//        MPI_Barrier(MPI_COMM_WORLD);
  	}
    }
  }
  delete[] yin;
  yin=NULL;
}

void gemXgcDatasProc3D::interpoDensityAlongY()
{
  double* yin;
  double* xin;
  double* yout = new double[p1->lj0];
  double* xout;
  LO nzb=bdesc->nzb;
  yin=new double[p1->lj0+2*nzb];
  xin=new double[p1->lj0+2*nzb];
  for(LO i=0;i<p1->li0;i++){
    for(LO k=0;k<p1->lk0;k++){
      for(LO j=0;j<p1->lj0;j++){
        yin[nzb+j]=densinterone[i][j][k]; 
      }
      for(LO b=0;b<nzb;b++){
        yin[b]=densinterone[i][p1->lj0+b-nzb][k];
        yin[p1->lj0+nzb-1+b]=densinterone[i][b][k];
      }
      xout=bdesc->ymeshxgc[i][k];
      Lag3dArray(yin,bdesc->ymeshgem,p1->lj0+2*nzb,yout,xout,p1->lj0); 
      for(LO j=0;j<p1->lj0;j++){
        densintertwo[i][j][k]=yout[j];
      }
    }
  }
}

void gemXgcDatasProc3D::InterpoPotential3DAlongZ(double*** boxin, double*** boxout)
{
  double* yin;
  double* yout;
  LO nzb=bdesc->nzb;
  if(preproc==true){
    for(LO i=0;i<p1->li0;i++){
      yin=new double[p3->mylk0[i]+2*nzb];
      for(LO j=0;j<p1->lj0;j++){
        for(LO l=0;l<nzb;l++){
          yin[l]=bdesc->lowpotentzgemxgc[i][j][l];
          yin[p3->mylk0[i]+nzb+l]=bdesc->uppotentzgemxgc[i][j][l];
        }
        for(LO k=0;k<p3->mylk0[i];k++){
          yin[k+nzb]=boxin[i][j][k];
        }
        yout=boxout[i][j];
        Lag3dArray(yin,bdesc->thflxmeshxgc[i],p3->mylk0[i]+2*nzb,yout,p1->thflx[i],p1->lk0);
     }
    }
  }

}


}
