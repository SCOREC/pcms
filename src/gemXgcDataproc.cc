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
   allocDensityArrays(); 

   allocPotentArrays();

   allocSendRecvbuff();
 }


 void gemXgcDatasProc3D::allocDensityArrays()
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

 void gemXgcDatasProc3D::allocPotentArrays()
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
/*
   numsend=new LO[p1->NP];
   numrecv=new LO[p1->NP];
   sdispls=new GO[p1->NP];
   rdispls=new GO[p1->NP];
*/
   // sending side is in tube_comm X grid_comm collective
   sendnum=0;
   recvnum=0;
   for(LO i=0;i<NP;i++){
     numsend[i]=0;
     numrecv[i]=0;
     sdispls[i]=0;
     rdispls[i]=0;
   }
   LO rank;
   for(LO i=0;i<p1->sendOverlap_x.size();i++){
     for(LO j=0;j<p1->sendOverlap_y.size();j++){
        rank=p1->sendOverlap_x[i][0]*p1->npz+p1->sendOverlap_th[j][0];
        numsend[rank]=(p1->sendOverlap_x[i][2]-p1->sendOverlap_x[i][1]+1)
         *(p1->sendOverlap_th[j][2]-p1->sendOverlap_th[j][1]+1)*p1->lj0;
     }
   } 
   for(LO i=1;i<p1->NP;i++){ 
     sendnum=numsend[i];
     sdispls[i]=sdispls[i-1]+numsend[i-1];
   } 
   for(LO i=0;i<p1->recvOverlap_x.size();i++){
     for(LO j=0;j<p1->recvOverlap_y.size();j++){
        rank=p1->recvOverlap_x[i][0]*p1->kmx+p1->recvOverlap_th[j][0];
        numrecv[rank]=(p1->recvOverlap_x[i][2]-p1->recvOverlap_x[i][1]+1)
         *(p1->recvOverlap_th[j][2]-p1->recvOverlap_th[j][1]+1)*p1->lj0;
     }
   }  
   for(LO i=1;i<p1->NP;i++){
     rdispls[i]=rdispls[i-1]+numrecv[i-1];
     recvnum=numrecv[i];
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

   distriDensityFromGemToCoupler();  
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




void gemXgcDatasProc3D::densityFromGemToCoupler()
{
// sending happens in MPI_COMM_WORLD, but the process in grid_comm and tube_comm is mappped to comm_x and comm_y    
   double* sendbuf=new double[sendnum];
   double* recvbuf=new double[recvnum];
   LO xnum=0;
   LO thnum=0;
   LO nrank=0;

   //for initialize the sendbuf  
   for(LO i=0;i<p1->npx;i++){
     for(LO k=0;k<p1->npz;k++)
     nrank=i*p1->npz+k;
     if(numsend[nrank]!=0){
       xnum=0;
       while(assert(i!=sendOverlap_x[xnum][0])){
         xnum+=1;
         assert(xnum<=sendOverlap_x.size());
       }
       thnum=0;
       while(assert(k!=sendOverlap_th[thnum][0])){
         thnum+=1;
         assert(thnum<=sendOverlap_th.size());
       }              
       for(LO i1=p1->sendOverlap_x[xnum][1];i1<p1->sendOverlap_x[xnum][2];i1++){
         for(LO j1=0;j1<p1->lj0;j1++){
           for(LO k1=p1->sendOverlap_th[thnum][1];k1++<p1->sendOverlap_th[thnum][2];k1++){
             sendbuf[sdispls[n]+(i1-p1->sendOverlap_x[xnum][1])*p1->lj0
               *p1->sendOverlap_th[thnum][3]+j1*p1->sendOverlap_th[thnum][3]
               +k1-p1->sendOverlap_th[thnum][1]]
               =den_tmp[i1-p1->tli0][j1][k1-p1->glk0];
           }
         }
       }     
     }
   }    
   
   MPI_Alltoallv(sendbuf,sendnum,sdispls,MPI_DOUBLE,recvbuf,recvnum,rdispls,MPI_DOUBLE,MPI_COMM_WORLD);

   LO tubenum=0;
   LO gridnum=0;
   LO nrank=0;
   for(LO i=0;i<p1->ntude;i++){
     for(LO k=0;k<p1->gnpz;k++)
     nrank=i*p1->gnpz+k;
     if(numrecv[nrank]!=0){
       tubenum=0;
       while(assert(i!=recvOverlap_x[tubenum][0])){
         tubenum+=1;
         assert(tubenum<=recvOverlap_x.size());
       }
       gridnum=0;
       while(assert(k!=recvOverlap_th[gridnum][0])){
         gridnum+=1;
         assert(gridnum<=recvOverlap_th.size());
       }              

       for(LO i1=p1->recvOverlap_x[tubenum][1];i1<p1->recvOverlap_x[tubenum][2];i1++){
         for(LO j1=0;j1<p1->lj0;j1++){
           for(LO k1=p1->recvOverlap_th[gridnum][1];k1++<p1->recvOverlap_th[gridnum][2];k1++){
             densgem[i1-p1->li0][j1][k1-p1->lk0]=recvbuf[rdispls[nrank]+(i1-p1->recvOverlap_x[tubenum][1])
               *p1->lj0*recvOverlap_th[gridnum][3]+j1*recvOverlap_th[gridnum][3]+k1-p1->recvOverlap_th[gridnum][1]];
           }
         }
       }          
     }
   } 
}




}
