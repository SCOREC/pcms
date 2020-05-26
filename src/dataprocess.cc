#include "dataprocess.h"
#include "importpart3mesh.h"
#include "commpart1.h"

namespace coupler {

DatasProc3D::DatasProc3D(const Part1ParalPar3D& p1pp3d,
    const Part3Mesh3D &p3m3d,
    bool pproc,
    TestCase test_case,
    bool ypar,
    int nummode)
  : preproc(pproc),
    testcase(test_case),
    yparal(ypar),
    p1(p1pp3d.li0,p1pp3d.lj0,p1pp3d.lk0,
	 p1pp3d.ny0, p1pp3d.npy, p1pp3d.mype_y, 
         p1pp3d.res_fact),
    p3(p3m3d.li0,p3m3d.lj0,p3m3d.mylk0)
  {
    init();
    AllocDensityArrays();
    AllocPotentArrays();
    if(testcase==TestCase::t0) {
      TestInitPotentAlongz(p3m3d, p1pp3d.npy, nummode);
    }
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
    densinterpo=new CV**[p1.li0];
    for(LO i=0;i<p1.li0;i++){
      densinterpo[i]=new CV*[p1.lj0];
      for(GO j=0;j<p1.lj0;j++){
        densinterpo[i][j]=new CV[p3.mylk0[i]];
      }
    }

   densintmp=new CV[p1.lj0];
   densouttmp=new double[p1.lj0*2];

   denspart3=new double**[p3.li0];
   for(LO i=0;i<p3.li0;i++){
     denspart3[i]=new double*[p3.lj0];
     for(LO j=0; j<p3.lj0; j++)
        denspart3[i][j]=new double[p3.mylk0[i]];
   }

//notice: the dimension may be gave more detaild describtion based on 
//the first coupling surface and last coupling surface
   densrecv = new CV*[p1.ny0];
   for(LO i=0;i<p1.ny0;i++){
     densrecv[i] = new CV*[p3.sum];
   }   

   denssend = new double*[p1.ny0*2];
   for(LO i=0;i<p1.ny0; i++){
     denssend = new double[p3.sum];
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

    potentintmp=new double[p3.lj0];
    potentouttmp=new CV[p3.lj0/2+1];
 
    potentinterpo=new CV**[p3.li0];
    for(LO i=0;i<p3.li0;i++){
      potentinterpo[i]=new CV*[p3.lj0/2];
      for(LO j=0;j<p3.lj0/2;j++)
        potentinterpo[i][j]=new CV[p3.mylk0[i]];
    }
  
    potentpart1=new CV**[p1.li0];
    for(LO i=0;i<p1.li0;i++){
      potentpart1[i]=new CV*[p1.lj0];
      for(LO j=0;j<p1.lj0;j++){
        potentpart1[i][j]=new CV[p1.lk0];
      }
    }
 
//notice: the dimension may be gave more detaild describtion based on
////the first coupling surface and last coupling surface 
    potentrecv = new double*[p1.ny0*2];
    for(LO i=0;i<p1.ny0;i++){
      potentrecv[i] = new double[sum];
    }

    potentsend = new double*[p1.ny0];
    for(LO i=0;i<p1.ny0;i++){
      potentsend[i] = new double[sum];
    }
  }
}


void DatasProc3D:DistriPotentRecvfromPart3(const LO* nstart,const LO* versurf, double** potentrecv)
{ 
  double** tmp
  for(LO j=0;j<p3.lj0;j++){
    tmp = new double*[p3.li0];
    LO xl=0
    for(LO i=0;i<p3.li0;i++){
      xl=p3.li1+i;
      tmp[i]=new double[versurf[xl]];
      GO sumbegin=0;
      for(LO h=0;h<xl;h++){
        sumbegin+=(GO)versurf[h];
      }
      for(LO m=0;m<versurf[xl];m++){
        tmp[i][m]=potentrecv[i][sumbegin+m];
      }
      reshuffleforward(tmp[xl],nstart[xl],vesurf[xl]);
      for(LO k=0;k<versurf[xl];k++){
        potentin[i][j][k]=tmp[i][lk1+k];
      }
      delete[] tmp[i];
    }
  }

 } 


void DatasProc3D:AssemPotentSendtoPart1(const LO* nstart,const LO* versurf, double** potentsend)
{
  LO* recvcount = new LO[p1.npz];
  LO* rdispls = new LO[p1.npz];
  MPI_Datatype mpitype = getDtype(LO)      
  MPI_Allgather(&p1.lk0,1,mpitype,recvcount,1,mpitype,p1.comm_z); 
  rdispls[0]=0;
  for(LO i=1;i<p1.npz;i++){
    rdispls[i]=ridspls[0]+recvcout[i];
  }
  for(LO j=0;j<p1.lj0;j++){
    LO xl=0;
    for(LO i=0;i<p1.li0;i++){
      xl=p1.li1+i;
      GO sumbegin=0;
      for(LO h=0;h<xl;h++){
        sumbegin+=(GO)p1.nz0
      }      
      CV* tmp = new CV[versurf[xl]];
      MPI_Allgatherv(potentpart1[i][j],p1.lk0,MPI_CXX_DOUBLE_COMPLEX,tmp,recvcount,rdispls,
                    MPI_CXX_DOUBLE_COMPLEX,p1.comm_z);    
      for(LO m=0;m<p1.nz0;m++){
        potentsend[j][sumbegin+m]=tmp[m];
      }      
      delete[] tmp; 
    }
  } 
  delete[] recvcount,rdispls;
}


void DatasProc3D:DistriDensiRecvfromPart1(const LO* versurf, CV** densrecv)
{
  CV** tmp
  tmp = new CV*[p1.li0];
  for(LO i=0;i<p1.li0;i++)
    tmp[i]=new CV[p1.nz0];
 
  for(LO j=0;j<p1.lj0;j++){
    LO xl=0
    for(LO i=0;i<p1.li0;i++){
      xl=p1.li1+i;
      GO sumbegin=0;
      for(LO h=0;h<xl;h++){
        sumbegin+=(GO)p1.nz0;
      }
      for(LO m=0;m<p1.nz0;m++){
        tmp[i][m]=densrecv[i][sumbegin+m];
      }
      for(LO k=0;k<lk0;k++){
        potentin[i][j][k]=tmp[i][lk1+k];
      }
     }
   }    
   for(LO i=0;i<p1.li0;i++){
     delete[] tmp[i];
   } 
}


void DatasProc3D:AssemDensSendtoPart3(const LO* nstart,const LO* versurf, CV** denssend)
{
  for(LO j=0;j<p1.lj0;j++){
    LO xl=0;
    for(LO i=0;i<p1.li0;i++){
      LO* recvcount = new LO[p1.npz];
      LO* rdispls = new LO[p1.npz];
      MPI_Datatype mpitype = getDtype(LO)      
      MPI_Allgather(&p3.mylk0[i],1,mpitype,recvcount,1,mpitype,p1.comm_z); 
      rdispls[0]=0;
      for(LO i=1;i<p1.npz;i++){
	rdispls[i]=ridspls[0]+recvcout[i];
      }

      xl=p3.li1+i;     
      CV* tmp = new CV[versurf[xl]];
      MPI_Allgatherv(denspart3[i][j],p3.mylk0[i],MPI_CXX_DOUBLE_COMPLEX,tmp,recvcount,rdispls,
                    MPI_CXX_DOUBLE_COMPLEX,p1.comm_z);    
      reshufflebackward(tmp,p3.nstart,versurf[xl]);
      GO sumbegin=0;
      for(LO h=0;h<xl;h++){
        sumbegin+=(GO)versurf[h];
      } 
      for(LO m=0;m<versurf[xl];m++){
        potentsend[j][sumbegin+m]=tmp[m];
      }      
      delete[] tmp; 
    }
    delete[] recvcount,rdispls;
  }  
 


}


void DatasProc3D::TestInitPotentAlongz(const Part3Mesh3D& p3m3d,
    const LO npy, const LO n) {
  if(npy==1){
    LO li0,lj0,lk0;
    li0=p3m3d.li0;
    lj0=p3m3d.lj0;
    double ylen;
    double sum;
    double dy=2.0*cplPI/double(lj0);
    for(LO i=0;i<li0;i++){
      lk0=p3m3d.mylk0[i];
      for(LO k=0;k<lk0;k++){
        ylen=0.0;
        for(LO j=0;j<lj0;j++){
          ylen=double(j)*dy;
          sum=0.0;
          for(LO h=0;h<n;h++){
            sum+=cos(double(h+1)*ylen-cplPI);
          }
          potentin[i][j][k]=sum;
        }
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
  if(densinterpo!=NULL) delete[] densinterpo;
  if(denspart3!=NULL) delete[] denspart3;
  if(potentin!=NULL) delete[] potentin;
  if(potentouttmp!=NULL) delete[] potentouttmp;
  if(potentinterpo!=NULL) delete[] potentinterpo;
  if(potentpart1!=NULL) delete[] potentpart1;       
}

}

