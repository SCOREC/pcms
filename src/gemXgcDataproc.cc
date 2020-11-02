#include "dataprocess.h"
#include "importpart3mesh.h"
#include "commpart1.h"
#include "BoundaryDescr3D.h"
#include "sendrecv_impl.h"
#include "interpoutil.h"
#include <cassert>
#include <cmath>

namespace coupler {

  gemXgcDatasProc3D::gemXgcDatasProc3D(const Part1ParalPar3D* p1pp3d,
    const Part3Mesh3D* p3m3d,
    const BoundaryDescr3D* bdesc_,
    bool pproc,
    TestCase test_case, 
    bool ypar)
 :  preproc(pproc),
    testcase(test_case),
    yparal(ypar)  
 {  
   p1 = p1pp3d;
   p3 = p3m3d;
   bdesc = bdesc_;  

   allocDensityArrays(); 

   allocPotentArrays();

   allocSendRecvbuff();

 }


 void gemXgcDatasProc3D::allocDensityArrays()
 {
   densin=new double**[p1->imx+1];
   for(LO i=0;i<p1->imx+1;i++){
     densin[i]=new double*[p1->jmx+1];
     for(LO j=0;j<p1->jmx+1;j++)
       densin[i][j]=new double[2];       
   }

   densCpl=new double**[p1->li0];
   for(LO i=0;i<p1->li0;i++){ 
     densCpl[i]=new double*[p1->lj0];
       for(LO j=0;j<p1->lj0;j++)
         densCpl[i][j]=new double[p1->lk0];     
   }

   densXgc = new double**[p3->li0];
   for(LO i=0;i<p3->li0;i++){
     densXgc[i]=new double*[p3->lj0];
     for(LO j=0;j<p3->lj0;j++)
       densXgc[i][j]=new double[p3->mylk0[i]];
   }    
 }

 void gemXgcDatasProc3D::allocPotentArrays()
 {
   pot_gem_fl=new double***[p1->li0];
   potyCpl=new double**[p1->li0];
   for(LO i=0;i<p1->li0;i++){
     pot_gem_fl[i]=new double**[p1->lj0];
     potyCpl[i]=new double*[p1->lj0];
     for(LO j=0;j<p1->lj0;j++){
       pot_gem_fl[i][j]=new double*[p3->mylk0[i]];
       potyCpl[i][j]=new double[p3->mylk0[i]];
	 for(LO k=0;k<p3->mylk0[i];k++)
	   pot_gem_fl[i][j][k]=new double[4];
     }
   }
  
   potythCpl=new double**[p1->li0];
   for(LO i=0;i<p1->li0;i++){
     potythCpl[i]=new double*[p1->lj0];
     for(LO j=0;j<p1->lj0;j++) potythCpl[i][j]=new double[p1->lk0];
   }

   potGem=new double[p1->tli0*p1->lj0*p1->glk0];
 } 

 void gemXgcDatasProc3D::allocSendRecvbuff()
 {
   numsend=new LO[p1->NP];
   numrecv=new LO[p1->NP];
   sdispls=new LO[p1->NP];
   rdispls=new LO[p1->NP];

   // sending side is in tube_comm X grid_comm collective
   sendnum=0;
   recvnum=0;
   for(LO i=0;i<p1->NP;i++){
     numsend[i]=0;
     numrecv[i]=0;
     sdispls[i]=0;
     rdispls[i]=0;
   }
   LO rank;
   for(LO i=0;i<p1->sendOverlap_x.size();i++){
     for(LO j=0;j<p1->sendOverlap_th.size();j++){
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
     for(LO j=0;j<p1->recvOverlap_th.size();j++){
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
   /*
   double* array=densityfromGEM->data();
   LO n=0;
   for(LO i=0;i<p1->tli0;i++){ 
     for(LO j=0;j<p1->lj0;j++){
       for(LO k=0;k<glk0;k++){
         n++;
         densin[i][j][k]=array[n];
       }
     }
   }
   */
   densityFromGemToCoupler(densityfromGEM);  
   zDensityBoundaryBufAssign(densCpl);
   interpoDensityAlongZ(densinterone);
   interpoDensityAlongY();
   
   //fixme: distribute and assemble the datas to get the matrix to be sent to xgc by the adios2 routine     
   
 }


void gemXgcDatasProc3D::DistriPotentRecvfromXGC(const Array3d<double>* potentfromXGC)
{
  double*** tmp=new double**[p1->li0];
  for(LO i=0;i<p1->li0;i++){
    tmp[i]=new double*[p3->nphi];
    for(LO j=0;j<p3->nphi;j++){
      tmp[i][j]=new double[p3->versurf[i+p1->li1]]; //The boundary buffer is assigned. 
    }
  }    

  //FIXME: distribute potentfromXGC->datas() to tmp;
  double* potent=potentfromXGC->data();  
  LO n=0;
  for(LO i=0;i<p1->li0;i++){
    for(LO j=0;j<p3->nphi;j++){
      for(LO k=0;k<p3->versurf[i+p1->li1];j++){
        n++;
        tmp[i][j][k]=potent[n];
      }
    }
  }  


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
          tmppotent=p3->thetaflx_pot[i][j][k][h];
          for(LO l=0;l<4;l++) tmpflx[l]=tmp[i][j][p3->thetaflx_ind_pot[i][j][k][h][l]];
          pot_gem_fl[i][j][k][h]=Lag3dInterpo1D(tmppotent,tmpflx,p3->thetaflx_pot[i][j][k][h][4]); 
        }
        tmppotent=pot_gem_fl[i][j][k];
        tmplength=p3->nodesdist_fl[i][j][k];
        potyCpl[i][j][k]=Lag3dInterpo1D(tmppotent,tmplength,p3->nodesdist_fl[i][j][k][4]);     
      }
    } 
  }  

  // The 2nd interpolation along theta
  InterpoPotential3DAlongZ(potyCpl,potythCpl);  
  potentFromCouplerToGem();
}


void gemXgcDatasProc3D::densityFromGemToCoupler(const Array3d<double>* densityfromGEM)
{
   double* array=densityfromGEM->data();
// sending happens in MPI_COMM_WORLD, but the process in grid_comm and tube_comm is mappped to comm_x and comm_y    
   double* sendbuf=new double[sendnum];
   double* recvbuf=new double[recvnum];
   LO xnum=0;
   LO thnum=0;
   LO nrank=0;
   LO n;

   //for initialize the sendbuf  
   for(LO i=0;i<p1->npx;i++){
     for(LO k=0;k<p1->npz;k++){
       nrank=i*p1->npz+k;
       if(numsend[nrank]!=0){
	 xnum=0;
	 while(i!=p1->sendOverlap_x[xnum][0]){
	   xnum+=1;
	   assert(xnum<=p1->sendOverlap_x.size());
	 }
	 thnum=0;
	 while(k!=p1->sendOverlap_th[thnum][0]){
	   thnum+=1;
	   assert(thnum<=p1->sendOverlap_th.size());
	 }      
		
	 for(LO i1=p1->sendOverlap_x[xnum][1];i1<p1->sendOverlap_x[xnum][2];i1++){
	   for(LO j1=0;j1<p1->lj0;j1++){
	     for(LO k1=p1->sendOverlap_th[thnum][1];k1++<p1->sendOverlap_th[thnum][2];k1++){
	       n=(i1-p1->tli0)*p1->lj0*p1->glk0+j1*p1->glk0+k1-p1->glk1;
	       sendbuf[sdispls[n]+(i1-p1->sendOverlap_x[xnum][1])*p1->lj0
		 *p1->sendOverlap_th[thnum][3]+j1*p1->sendOverlap_th[thnum][3]
		 +k1-p1->sendOverlap_th[thnum][1]]=array[n];  //=densin[i1-p1->tli0][j1][k1-p1->glk0];
	     }
	   }
	 }     
       }
     }    
   }
   MPI_Alltoallv(sendbuf,numsend,sdispls,MPI_DOUBLE,recvbuf,numrecv,rdispls,MPI_DOUBLE,MPI_COMM_WORLD);

   LO tubenum=0;
   LO gridnum=0;
   nrank=0;
   for(LO i=0;i<p1->ntube;i++){
     for(LO k=0;k<p1->gnpz;k++){
       nrank=i*p1->gnpz+k;
       if(numrecv[nrank]!=0){
	 tubenum=0;
	 while(i!=p1->recvOverlap_x[tubenum][0]){
	   tubenum+=1;
	   assert(tubenum<=p1->recvOverlap_x.size());
	 }
	 gridnum=0;
	 while(k!=p1->recvOverlap_th[gridnum][0]){
	   gridnum+=1;
	   assert(gridnum<=p1->recvOverlap_th.size());
	 }              

	 for(LO i1=p1->recvOverlap_x[tubenum][1];i1<p1->recvOverlap_x[tubenum][2];i1++){
	   for(LO j1=0;j1<p1->lj0;j1++){
	     for(LO k1=p1->recvOverlap_th[gridnum][1];k1++<p1->recvOverlap_th[gridnum][2];k1++){
	       densCpl[i1-p1->li0][j1][k1-p1->lk0]=recvbuf[rdispls[nrank]+(i1-p1->recvOverlap_x[tubenum][1])
		 *p1->lj0*p1->recvOverlap_th[gridnum][3]+j1*p1->recvOverlap_th[gridnum][3]+k1-p1->recvOverlap_th[gridnum][1]];
	     }
	   }
	 }          
       }
     } 
   }  
}

// It's a reverse procedure of densityFromGemToCoupler
void gemXgcDatasProc3D::potentFromCouplerToGem()
{
   GO scounts=recvnum;
   GO rcounts=sendnum;
   double* sendbuf=new double[scounts];
   double* recvbuf=new double[rcounts];

   //for initialize the sendbuf  
   LO tubenum=0;
   LO gridnum=0;
   LO nrank=0;
   for(LO i=0;i<p1->ntube;i++){
     for(LO k=0;k<p1->gnpz;k++){
       nrank=i*p1->gnpz+k;
       if(numrecv[nrank]!=0){
	 tubenum=0;
	 while(i!=p1->recvOverlap_x[tubenum][0]){
	   tubenum+=1;
	   assert(tubenum<=p1->recvOverlap_x.size());
	 }
	 gridnum=0;
	 while(k!=p1->recvOverlap_th[gridnum][0]){
	   gridnum+=1;
	   assert(gridnum<=p1->recvOverlap_th.size());
	 }              

	 for(LO i1 = p1->recvOverlap_x[tubenum][1]; i1 < p1->recvOverlap_x[tubenum][2]; i1++){
	   for(LO j1 = 0; j1<p1->lj0; j1++){
	     for(LO k1 = p1->recvOverlap_th[gridnum][1]; k1 < p1->recvOverlap_th[gridnum][2]; k1++){
	       sendbuf[rdispls[nrank]+(i1-p1->recvOverlap_x[tubenum][1])
		 *p1->lj0*p1->recvOverlap_th[gridnum][3]+j1*p1->recvOverlap_th[gridnum][3]+k1-p1->recvOverlap_th[gridnum][1]]
	       =potythCpl[i1-p1->li0][j1][k1-p1->lk0];
	     }
	   }
	 }          
       }
     }
   } 

   MPI_Alltoallv(sendbuf,numrecv,rdispls,MPI_DOUBLE,recvbuf,numsend,sdispls,MPI_DOUBLE,MPI_COMM_WORLD);

   LO xnum=0;
   LO thnum=0;
   nrank=0;
   LO n;

   //for initialize the sendbuf  
   for(LO i=0;i<p1->npx;i++){
     for(LO k=0;k<p1->npz;k++){
       nrank=i*p1->npz+k;
       if(numsend[nrank]!=0){
	 xnum=0;
	 while(i!=p1->sendOverlap_x[xnum][0]){
	   xnum+=1;
	   assert(xnum<=p1->sendOverlap_x.size());
	 }
	 thnum=0;
	 while(k!=p1->sendOverlap_th[thnum][0]){
	   thnum+=1;
	   assert(thnum<=p1->sendOverlap_th.size());
	 }      
		
	 for(LO i1=p1->sendOverlap_x[xnum][1];i1<p1->sendOverlap_x[xnum][2];i1++){
	   for(LO j1=0;j1<p1->lj0;j1++){
	     for(LO k1=p1->sendOverlap_th[thnum][1];k1++<p1->sendOverlap_th[thnum][2];k1++){
		potGem[(k1-p1->sendOverlap_th[thnum][1])*p1->lj0*(p1->imx+1)+j1*(p1->imx+1)+i1]
		=recvbuf[sdispls[n]+(i1-p1->sendOverlap_x[xnum][1])*p1->lj0
		 *p1->sendOverlap_th[thnum][3]+j1*p1->sendOverlap_th[thnum][3]
		 +k1-p1->sendOverlap_th[thnum][1]];
	     }
	   }
	 }     
       }
     }
   } 
   MPI_Reduce(MPI_IN_PLACE,potGem,p1->tli0*p1->lj0*p1->glk0,MPI_DOUBLE,MPI_SUM,0,p1->tube_comm);

}




}
