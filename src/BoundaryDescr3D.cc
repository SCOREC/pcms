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
    const Part1ParalPar3D& p1pp3d,
    const DatasProc3D& dp3d,
    const TestCase tcase,
    bool pproc):test_case(tcase), preproc(pproc)
{
  if(preproc==true){
    nzb=p1pp3d.nzb;
    updenz=new CV**[p1pp3d.li0];
    lowdenz=new CV**[p1pp3d.li0];

    for(LO i=0;i<p1pp3d.li0;i++){
      updenz[i]=new CV*[p1pp3d.lj0];
      lowdenz[i]=new CV*[p1pp3d.lj0];
      for(LO j=0;j<p1pp3d.lj0;j++){
	updenz[i][j]=new CV[nzb];
	lowdenz[i][j]=new CV[nzb];
      }
    }

    if(p1pp3d.mype_z==0){
      if(lowpbmat==NULL){
        lowpbmat=new CV*[p1pp3d.li0];
        for(LO i=0;i<p1pp3d.li0;i++){
          lowpbmat[i]=new CV[p1pp3d.lj0];      
        } 
      }
     }else if(p1pp3d.mype_z==p1pp3d.npz-1){
      if(uppbmat==NULL){
        uppbmat=new CV*[p1pp3d.li0];
        for(LO i=0;i<p1pp3d.li0;i++){
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
      for(LO j=0;j<p3m3d.lj0;j++){
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
     for(LO i=0;i<p1pp3d.li0;i++){
       num=p1pp3d.li1+i;
       for(LO j=0;j<p1pp3d.lj0;j++){  
	 lowpbmat[i][j]=exp(CV(0.0,1.0)*2.0*cplPI*double(p1pp3d.n0_global)
                        *double(j+p1pp3d.ky0_ind)*p1pp3d.q_prof[num]);
       }
     } 
   } else if(p1pp3d.mype_z==p1pp3d.npz-1){
     for(LO i=0;i<p1pp3d.li0;i++){
       num=p1pp3d.li1+i;
       for(LO j=0;j<p1pp3d.lj0;j++){
         uppbmat[i][j]=exp(-CV(0.0,1.0)*2.0*cplPI*double(p1pp3d.n0_global)
                        *double(j+p1pp3d.ky0_ind)*p1pp3d.q_prof[num]);
       }
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
  li0=p3m3d.xboxinds[p1pp3d.mype_x][0];
  lj0=p3m3d.lj0/2;
  if(p1pp3d.npz>1){
    if(p1pp3d.periods[2]==1){  
      for(LO i=0;i<li0;i++){
        lk0=p3m3d.mylk0[i];
        if(lk0<nzb){
          std::cout<<"ERROR: the interpolation order is larger than the box count along z dimension.";
          std::exit(EXIT_FAILURE);
        } 
        mpisendrecv_aux1D(p1pp3d.comm_z,nzb,li0,lj0,lk0,lowzpart3[i],upzpart3[i],
          p3m3d.pzcoords[i]);   
     
       //for test debugging
       if(test_case==TestCase::t0){
         if(p1pp3d.mype_z==0){
             std::cout<<"lowzpart3="<<lowzpart3[i][0]<<" "<<lowzpart3[i][1]<<'\n';
             std::cout<<"lowpzcoords="<<p3m3d.pzcoords[i][0]<<" "<<p3m3d.pzcoords[i][1]<<'\n';
           
         } else if(p1pp3d.mype_z==p1pp3d.npz-1){
            std::cout<<"upzpart3="<<upzpart3[i][0]<<" "<<upzpart3[i][1]<<'\n';
            std::cout<<"upzcoords="<<p3m3d.pzcoords[i][lk0-2]<<" "<<p3m3d.pzcoords[i][lk0-1]<<'\n'; 
         }
       }

       for(LO j=0;j<lj0;j++){
          mpisendrecv_aux1D(p1pp3d.comm_z,nzb,li0,lj0,lk0,lowpotentz[i][j],uppotentz[i][j],
              dp3d.potentinterpo[i][j]); 
        //enforce the parallel boundary condition
          if(p1pp3d.mype_z==0){
            for(LO k=0;k<nzb;k++)
              lowpotentz[i][j][k]=lowpotentz[i][j][k]*lowpbmat[i][j];
          } else if(p1pp3d.mype_z==p1pp3d.npz-1){
             for(LO k=0;k<nzb;k++)
               uppotentz[i][j][k]=uppotentz[i][j][k]*uppbmat[i][j];
          }
       }
     }
  
     if(p1pp3d.mype_z==0){
         for(LO h=0;h<li0;h++){
           for(LO k=0;k<nzb;k++){  
             lowzpart3[h][k]=lowzpart3[h][k]-2.0*cplPI;
           }
         }
     }else if(p1pp3d.mype_z==p1pp3d.npz-1){
         for(LO h=0;h<li0;h++){
           for(LO k=0;k<nzb;k++){
              upzpart3[h][k]=upzpart3[h][k]+2.0*cplPI;
           }
         }          
     }
     //for test debugging
     if(test_case==TestCase::t0){
         if(p1pp3d.mype_z==0){
           for(LO k=0;k<li0;k++){
             std::cout<<"lowzpart3["<<k<<"][1]="<<lowzpart3[k][1]<<'\n';
           }
         }else if(p1pp3d.mype_z==p1pp3d.npz-1){
            for(LO k=0;k<li0;k++){ 
              std::cout<<"upzpart3["<<k<<"][1]="<<upzpart3[k][1]<<'\n'; 
            }  
         }
      } 
    } else if(p1pp3d.periods[2]==0){
         std::cout<<"The parallelization of 'z' domain is not down with unperiodic boundary condiiton"
                  <<" and npz>1"<<'\n';
         std::exit(EXIT_FAILURE);
    }
  } else {
      if(p1pp3d.periods[2]==1){ 
	for(LO i=0;i<li0;i++){
	   lk0=p3m3d.mylk0[i];
	   if(lk0<nzb){
	     std::cout<<"ERROR: the interpolation order is larger than the box count along z dimension.";
	     std::exit(EXIT_FAILURE);
	   }  
	   for(LO k=0;k<nzb;k++){
	     lowzpart3[i][k]=p3m3d.pzcoords[i][lk0-nzb+k];
	     upzpart3[i][k]=p3m3d.pzcoords[i][k];
	   }
	   for(LO j=0;j<lj0;j++){
	     for(LO k=0;k<nzb;k++){
	       lowpotentz[i][j][k]=dp3d.potentin[i][j][k];
	       lowpotentz[i][j][k]=dp3d.potentin[i][j][lk0-nzb+k];
	     }  
	   }     
	}
     } else if(p1pp3d.periods[2]==0) {
         std::cout<<"The parallelization of 'z' domain is not down with unperiodic boundary condiiton"
                  <<" and npz=1"<<'\n';
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

void BoundaryDescr3D::zDensityBoundaryBufAssign(CV*** box,
    const Part1ParalPar3D& p1pp3d) {
  if (lowdenz == NULL || updenz == NULL) {
    std::cout << "ERROR:the boundary buffer must be alloctted before "
                 "calling this routine.";
    std::exit(EXIT_FAILURE);
  }
  const LO lx = p1pp3d.li0;
  const LO ly = p1pp3d.lj0;
  const LO lz = p1pp3d.lk0;
  if (p1pp3d.npz > 1) {
    if (lz >= nzb) {
      mpisendrecv_aux2D(p1pp3d.comm_z, nzb, lx, ly, lz, lowdenz, updenz, box);
    } else {
      std::cout << "ERROR: nzb is larger than lz. A larger lz is required.";
      std::exit(EXIT_FAILURE);
    }
  } else {
    if (p1pp3d.periods[2] == 1) {
      for (LO i = 0; i < lx ; i++) {
        for (LO j = 0; j < ly; j++) {
          for (LO k = 0; k < nzb; k++) {
            lowdenz[i][j][k] = box[i][j][lz - nzb + k];
            updenz[i][j][k] = box[i][j][k];
            if(p1pp3d.mype_z==0){
              lowdenz[i][j][k]=lowdenz[i][j][k]*lowpbmat[i][j];
            } else if(p1pp3d.mype_z==p1pp3d.npz-1){
              updenz[i][j][k]=updenz[i][j][k]*uppbmat[i][j];
            }
          }
        }
      }
    } else {
      std::cout << "The topology is not right." << '\n';
      std::exit(EXIT_FAILURE);
    }
  }
}

}
 
