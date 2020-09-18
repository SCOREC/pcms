#include "dataprocess.h"
#include "importpart3mesh.h"
#include "commpart1.h"
#include "BoundaryDescr3D.h"
#include "sendrecv_impl.h"
#include <cassert>
#include <cmath>

namespace coupler {
 void  gemXgcDatasProc3D::gemXgcDatasProc3D(bool pproc,TestCase test_case, bool ypar)
 : preproc(pproc),
   testcase(test_case),
   yparal(ypar)  
 {    
   AllocDensityArrays(); 

   AllocPotentArrays();
 }


 void gemXgcDatasProc3D::AllocDensityArrays()
 {
   densin=new double**[p1->imx+1];
   for(LO i=0;i<p1->imx+;i++){
     densin[i]=new double*[p1->jmx+1];
     for(LO j=0;j<p1->jmx+1;j++)
       densin[i][j]=new double[2];       
   }

   densgem=new double**[p1->li0];
   for(LO i=0;i<p1->li0;i++){ 
     densgem[i]=new double*[p1->lj0]{
       for(LO j=0;j<p1->lj0;j++)
         densgem[i][j]=new double[p1->lk0];
     }
   }

   densxgc=new double**[p3->li0];
   for(LO i=0;i<p3->li0;i++){
     densxgc=new double*[p3->lj0];
     for(LO j=0;j<p3->lj0;j++)
       densxgc=new double[p3->mylk0[i]];
   }    
 }

 void gemXgcDatasProc3D::AllocPotentArrays()
 {
   double**** pot_gem_fl=new double***[p1->li0];
   double*** pot_ygem=new double**[p1->li0];
   for(LO i=0;i<p1->li0;i++){
     pot_gem_fl[i]=new double**[p1->lj0];
     pot_ygem[i]=new double*[p1->lj0];
     for(hj=0;j<p1->lj0;j++){
       pot_gem_fl[i][j]=new double*[>mylk0[i]];
       pot_ygem[i][j]=new double[mylk0[i]];
	 for(LO k=0;k<mylk0[i];k++)
	   pot_gem_fl[i][j][k]=new double[4];
     }
   }

 } 

 void allocSendRecvbuff()
 {
   numsend=new LO[p1->NP];
   numrecv=new LO[p1->NP];
   sdispls=new GO[p1->NP];
   rdispls=new GO[p1->NP];

   // sending side is in tube_comm X grid_comm collective
   for(LO i=0;i<NP;i++){
     numsend[i]=0;
     numrecv[i]=0;
     sdispls[i]=0;
     rdispls[i]=0;
   }
   LO rank;
   for(LO i=0;i<sendOverlap_x.size();i++){
     for(LO j=0;j<sendOverlap_y.size();j++){
        rank=sendOverlap_x[i][0]
 
     }
   }
    
  

 }

 void gemXgcDatasProc3D::DistriDensiRecvfromGem(const Array3d<double>* densityfromGEM)
 {    
   double* array=densityfromGEM->data();
   for(LO i=0;i<p1->imx+1;i++){ 
     for(LO j=0;j<p1->jmx+1;j++){
       for(LO k=0;k<2;k++){
         densin[i][j][k]=array[(jmx+1)*i+2*j+k];
       }
     }
   }

   densityfromGemToCoupler();   
   zDensityBoundaryBufAssign(densgem,bdesc);
   interpoDensityAlongZ(densinterone);
   interpoDensityAlongY();
    
//fixme: distribute and assemble the datas to get the matrix to be sent to xgc by the adios2 routine     

 }


void gemXgcDatasProc3D::DistriPotentRecvfromXGC(const Array3d<double>* potentfromXGC)
{
  double*** tmp=new double**[p1->li0];
  for(LO i=0;i<p1->li0;i++){
    tmp[i]=new double*[nphi];
    for(LO j=0;j<nphi;j++){
      tmp[i][j]=new double[p3->versurf[i+p1->li1]]; //The boundary buffer is assigned. 
    }
  }    

  //FIXME: distribute potentfromXGC->datas() to tmp;  

  // First interpolation along the field line; 
  double* tmppotent;
  double* tmpflx=new double[4];
  double* tmplength;
  for(LO i=0;i<p1->li0;i++){
    for(LO j=0;j<p1->lj0;j++){ 
/*
      double* tmpthflx=new double[2.0*bdesc->nzb+p3->mylk0[i]];
      for(LO h=0;h<bdesc->nzb;h++) tmpthflx[h]=bdesc->lowpotentgemxgc[h];
      for(LO l=0;l<bdesc->) 
*/ 
      for(LO k=0;k<p3->mylk0[i];k++){ 
        for(LO h=0;h<4;h++){
          tmppotent=thetaflx_pot[i][j][k][h];
          for(LO l=0;l<4;l++) tmpflx[l]=tmp[i][j][thetaflx_ind_pot[i][j][k][h][l]];
          pot_gem_fl[i][j][k][h]=Lag3dInterpo1D(tmppotent,tmpflx,thetaflx_pot[i][j][k][h][4]); 
      }
      tmppotent=pot_gem_fl[i][j][k];
      tmplength=nodesdist_fl[i][j][k];
      pot_ygem[i][j][k]=Lag3dInterpo1D(tmppotent,tmplength,nodesdist_fl[i][j][k][4]);     
     }
   } 
 }  

 // The 2nd interpolation along theta
 InterpoPotential3DAlongZ(pot_ygem,pot_ythgem);  
 // FIXME: distribute pot_ythgem
}




void gemXgcDatasProc3D::densityfromGemToCoupler(const double*** tmpdensity)
{
   numsend=new LO[p1->NP];
   numrecv=new LO[p1->NP];
   sdispls=new GO[p1->NP];
   rdispls=new GO[p1->NP];
 
// sending side is in tube_comm X grid_comm collective    
   

}

void prepareSendRecv()
 {
   


 }

void gemXgcDatasProc3D::densityfromGemToCoupler(const double*** tmpdensity)
 {
   LO* tmplk1=new double[p1->npz];
   LO* tmplk2=new double[p1->npz];
   MPI_Allgather(&p1->lk1,1,MPI_INT,tmplk1,1,MPI_INT,p1->comm_z);
   MPI_Allgather(&p1->lk2,1,MPI_INT,tmplk1,1,MPI_INT,p1->comm_z);
   LO i=0;
   while(tmplk1[i]-mype_grid!=0) i+=1;
   LO minrank=i;
   i=0
   while(tmplk2[i]-mype_grid!=0) i+=1;
   LO maxrank=i;
 
   LO* numsend=new LO[p1->npz];
   LO* numrecv=new LO[p1->npz];
   LO* sdispls=new GO[p1->npz];
   LO* rdispls=new GO[p1->npz];
   GO  scounts=0;
   GO  scounts=0;

   for(LO i=0;i<p1->npz;i++){
     if(i<minrank){
       numsend[i]=0;
     }else if(i>maxrank){
       numsend[i]=0;
     }else {
       numsend[i]=li0*lj0;
     }
     
     if(i<p1->lk1){
       numrecv[i]=0;
     }else if(i>p1->lk2){
       numrecv[i]=0;
     }else{
       numrecv[i]=li0*lj0;
     }
   } 

   sdispls[0]=0;
   rdispls[0]=0;
   for(LO i=1;i<p1->npz;i++){
     sdispls[i]=sdispls[i-1]+numsend[i-1];
     rdispls[i]=rdispls[i-1]+numrecv[i-1];
   }

   for(LO i=0;i<p1->npz;i++){
     scounts+=numsend[i];
     rcounts+=numrecv[i];
   }
   
   double* sendbuf=new double[scounts];
   double* recvbuf=new double[rcounts];

   for(LO l=0;l<p1->npz;l++){
     if(numsend[l]!=0){
       for(LO i=0;i<p1->li0;i++){
         for(LO j=0;j<p1->lj0;j++){
           sendbuf[sdispls[l]+i*lj0+j]=den_tmp[p1->li1+i][j];
         } 
       }
     }
   }

   MPI_Alltoallv(sendbuf,scounts,sdispls,MPI_DOUBLE,recvbuf,rcounts,rdispls,MPI_DOUBLE,p1->comm_z);

   for(LO l=0;l<p1->npz;l++){
     if(numrecv[l]!=0){
       for(LO i=0;i<p1->li0;i++){
         for(LO j=0;j<p1->lj0;j++){
           densgem[i][j][l-p1->lk0]=recvbuf[rdispls[l]+i*p1->lj0+j];
         }
       }
     }
   }      
 }


}
