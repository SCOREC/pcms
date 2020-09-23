#include "commpart1.h"
#include "interpoutil.h"
#include <cassert>
#include <cmath>
#include <algorithm>

namespace coupler {

// Input GENE's mesh
void Part1ParalPar3D::initTest0(std::string test_dir)
{
  assert(!test_dir.empty());
  LO* data=new LO[12];
  std::string fname=test_dir+"parpart1.nml";
  InputfromFile(data,12,fname); 
 
  npx=data[0];
  nx0=data[1];
  nxb=data[2];
  li0=data[3];

  npy=data[4];
  ny0=data[5];
  nyb=data[6];
  lj0=data[7];

  npz=data[8];
  nz0=data[9];
  nzb=data[10];
  lk0=data[11];

  NP=npx*npy*npz;
  CreateSubCommunicators();
  li1=mype_x*li0;
  li2=li1+li0-1;
  lj1=mype_y*lj0;
  lj2=lj1+lj0-1;
  lk1=mype_z*lk0;
  lk2=lk1+lk0-1;
  
  delete[] data;

  n0_global=8.0;
  ky0_ind=0;
 
  q_prof=new double[nx0]; 
  fname=test_dir+"q_prof.nml";
  InputfromFile(q_prof,nx0,fname);
  
}

//read the paralllization parameters
void Part1ParalPar3D::init(LO* parpar, double* xzcoords, double* q_prof, double* gene_cy, std::string test_dir)
{
 if(preproc==true){ 
   if(test_case==TestCase::t0){
     assert(!parpar && !xzcoords && !q_prof);//when not testing, arrays come from ADIOS2
     initTest0(test_dir);
//  The following two lines looks weird here. 
     xzcoords = new double[nx0];
     parpar = new LO[1];// This is unused but it is deleted
   }else{
     assert(parpar);
     assert(q_prof);
     assert(gene_cy);
     npx=parpar[0];
     nx0=parpar[1];
     nxb=parpar[2];
     li0=parpar[3];
     li1=parpar[4];
     li2=parpar[5];
     lg0=parpar[6];
     lg1=parpar[7];
     lg2=parpar[8];

     npy=parpar[9];
     ny0=parpar[10];
     nyb=parpar[11];
     llj0=parpar[12];
     llj1=parpar[13];
     llj2=parpar[14];
     lm0=parpar[15];
     lm1=parpar[16];
     lm2=parpar[17];

     npz=parpar[18];
     nz0=parpar[19];
     nzb=parpar[20];
     lk0=parpar[21];
     lk1=parpar[22];
     lk2=parpar[23];
     ln0=parpar[24];
     ln1=parpar[25];
     ln2=parpar[26];

     n0_global=parpar[27];
     ky0_ind=parpar[28];    
     NP=npx*npy*npz; 
     lj0=llj0;
     lj1=llj1;
     lj2=llj2;    

     CreateSubCommunicators(); 
   }

   std::cout<<"mype,mype_x,mype_z,li1,li2,lj1,lj2,lk1,lk2="<<mype<<" "<<mype_x<<" "<<mype_z<<" "<<li1<<" "
   <<li2<<" "<<lj1<<" "<<lj2<<" "<<lk1<<" "<<lk2<<'\n';

   // initialize the radial locations of the flux surface and poloidal angles
   pzcoords=new double[nz0];
   xcoords=new double[nx0];

   n_cuts = 2*lj0;

   C_y = gene_cy;
   minor_r = gene_cy[nx0];
   lx_a = gene_cy[nx0+1];
   sign_phi = gene_cy[nx0+2];
   dx = gene_cy[nx0+3];
   rhostar=gene_cy[nx0+4];
   Tref=gene_cy[nx0+5];
   nref=gene_cy[nx0+6];
   mref=gene_cy[nx0+7];
   Bref=gene_cy[nx0+8];
   Lref=gene_cy[nx0+9];

   L_tor=sign_phi*2.0*cplPI/double(n0_global*lj0*2);
   if(mype==0){
     std::cout<<"n_cuts="<<n_cuts<<'\n';
     std::cout<<"minor_r="<<minor_r<<'\n';
     std::cout<<"n0_global="<<n0_global<<'\n';
     std::cout<<"sign_phi="<<sign_phi<<'\n';
     std::cout<<"lj0="<<lj0<<'\n';
     std::cout<<"L_tor="<<L_tor<<'\n'; 
     std::cout<<"rhostar="<<rhostar<<'\n';
     std::cout<<"Tref="<<Tref<<'\n';
     std::cout<<"nref="<<nref<<'\n';
     std::cout<<"mref="<<mref<<'\n';
     std::cout<<"Bref="<<Bref<<'\n';
     std::cout<<"Lref="<<Lref<<'\n';
     std::cout<<"ky0_ind="<<ky0_ind<<'\n';
   }   
   phi_cut = new double[lj0*2];
   for(LO i=0;i<lj0*2;i++){
     phi_cut[i] = L_tor*double(i+1); // need more check
   } 
   res_fact=n0_global;
   dy=dx*rhostar*minor_r/res_fact;
   y_res=2.0*lj0;
   y_res_back=y_res*res_fact;
   norm_fact_dens=nref * 1.0e+19 * rhostar* minor_r;
   norm_fact_field=Tref * 1.0e+3* rhostar * minor_r;   
   if(mype==0){
     std::cout<<"norm_fact_dens="<<norm_fact_dens<<'\n';
     std::cout<<"norm_fact_field="<<norm_fact_field<<'\n';
     std::cout<<"y_red_back="<<y_res_back<<'\n';
   }

   totnodes=nx0*nz0;
   if(test_case==TestCase::t0){
      assert(!test_dir.empty());
      std::string fname=test_dir+"xcoords.nml";
      InputfromFile(xzcoords,nx0,fname);
   }else{
      assert(xzcoords);
   }

   for(LO i=0;i<nx0;i++){
     xcoords[i]=xzcoords[i];
   }
   if(test_case==TestCase::t0){
     dz=2.0*cplPI/nz0;
   }else{
     dz=xzcoords[nx0];
   }

   for(LO i=0;i<nz0;i++){
     pzcoords[i]=-1.0*cplPI+(double)i*dz;
     if(mype==0) printf("i=: %3d,pz=:%19.17f \n",i,pzcoords[i]);
//std::cout<<"i="<<i<<" "<<"pz="<<pzcoords[i]<<'\n';
   }
   pzp=new double[lk0];
   for(LO i=0;i<lk0;i++){
     pzp[i]=double(lk1+i)*dz-1.0*cplPI;
   }
   blockindice();  

   if(test_case==TestCase::t0){ 
     delete[] parpar;
     delete[] xzcoords;
   }
 }
}

void Part1ParalPar3D::CreateSubCommunicators()
{
   // create 3D parallel cart with z being periodic
   int rorder = 0;
   int dim[3]={(int)npx,(int)npy,(int)npz};
   MPI_Cart_create(MPI_COMM_WORLD,3,dim,periods,rorder,&comm_cart);

   MPI_Comm subcomuni[3];
   for(int i=0;i<3;i++){
     int remain[3]={0,0,0};
     remain[i]=1;
     MPI_Cart_sub(comm_cart,remain,&subcomuni[i]);
   }

   comm_x=subcomuni[0];
   comm_y=subcomuni[1];
   comm_z=subcomuni[2];

   MPI_Comm_rank(MPI_COMM_WORLD,&mype);
   MPI_Comm_rank(comm_x,&mype_x);
   MPI_Comm_rank(comm_y,&mype_y);
   MPI_Comm_rank(comm_z,&mype_z);
}

void Part1ParalPar3D::MpiFreeComm()
{
  MPI_Comm_free(&comm_x);
  MPI_Comm_free(&comm_y);
  MPI_Comm_free(&comm_z);
  MPI_Comm_free(&comm_cart);
}

void Part1ParalPar3D::blockindice()
{
   blockcount = GO(nz0*li0);
   blockstart = GO(nz0*li1);
   blockend = GO(nz0*(li2+1)-1);

}

//Input GEM's mesh

void Part1ParalPar3D::initGem(const Array1d<int>* gemmesh, const Array2d<double>* thfnz)
{  
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  LO* tmp=gemmesh->data();
  ntube=tmp[0];
  imx=tmp[1];
  jmx=tmp[2];
  kmx=tmp[3]; 
  ntheta=tmp[4];

  thetagrideq=new double[ntheta];
  thetagrideq[ntheta/2]=0.0;
  for(LO i=ntheta/2+1;i<ntheta+1;i++){
    thetagrideq[i]=double(i-ntheta/2)*dth;
    thetagrideq[ntheta-i]=-thetagrideq[i];
  }
  
  CreateGemsubcommunicators();
  if(npx>imx) std::cout<<"Error: npx>imx; radial mesh is not dense enough"<<'\n';
  std::exit(1); 
  decomposeGemMeshforCoupling();
  CreateSubCommunicators();
  double* tmpreal;
  tmpreal=thfnz->data();
  lz=tmpreal[0];
  ly=tmpreal[1];
  dz=lz/double(kmx);
  delz=lz/double(ntheta);
  dy=ly/double(jmx);
  dth=2.0*cplPI/double(ntheta);

  thflxeq=new double*[li0];
  for(LO i=0;i<li0;i++) thflxeq[i]=new double[ntheta+1];  
//  LO surfx=0;
//  for(i=0;i<mype_x;i++) surfx+=li0[i];
  for(LO i=0;i<li0;i++){    
    for(LO k=0;k<ntheta+1;k++)
      thflxeq[i][k]=tmpreal[(li1+i)*(ntheta+1)+k];
  }   
 
  thflx=new double*[li0];
  for(LO i=0;i<li0;i++) thflx[i]=new double[lk0];
 
  //interpolation for obtaining the flux theta of mesh for the perturbation 
  double* tmpth=new double[kmx+1];
  double tmpdth=2.0*cplPI/double(kmx); 
  for(LO i=kmx/2+1;i<kmx+1;i++){      // Here, another way is to minus cplPI
    theta[i]=double(i-kmx/2)*tmpdth;
    theta[kmx-i]=-theta[i];
  }

  double* tmpthetaeq=new double[ntheta+5];
  tmpthetaeq[0]=thetagrideq[0]-2.0*dth;
  tmpthetaeq[1]=thetagrideq[1]-dth;
  tmpthetaeq[ntheta+2]=thetagrideq[ntheta]+dth;
  tmpthetaeq[ntheta+3]=thetagrideq[ntheta]+2.0*dth;
  for(LO k=2;k<ntheta+2;k++) tmpthetaeq[k]=thetagrideq[k-2];
  
  //Here, the continusous boundary condition is used for the 3rd-order Lagrangain interpolaiton; It's better to replace it with the cubic spline interpolation
  double* tmpflxeq=new double[ntheta+5];  
  double* tmpflx=new double[kmx+1]; 
  for(LO i=0;i<li0;i++){
    tmpflxeq[0]=thflxeq[li1+i][0];
    tmpflxeq[1]=thflxeq[li1+i][0];
    tmpflxeq[ntheta+2]=thflxeq[li1+i][ntheta-1];   
    tmpflxeq[ntheta+3]=thflxeq[li1+i][ntheta-1];
    for(LO k=2;k<ntheta+2;i++) tmpflxeq[k]=thflxeq[li1+i][k-2];
    Lag3dArray(tmpflxeq,tmpthetaeq,ntheta+5,tmpflx,theta,kmx+1); 

    //Then, the initialization of theflx
    for(LO k=0;k<lk0;k++) thflx[i][k]=tmpflx[lk1+k];
  }  

  double* y_gem = new double[jmx+1];
  for(LO j=0;j<jmx+1;j++){
    y_gem[j]=double(j)*ly/dy;
  }
  delete[] tmpthetaeq;
  delete[] tmpflxeq;
  delete[] tmpflx;

}


void Part1ParalPar3D::CreateGemsubcommunicators()
{
//split the communicator for GEM's domain decomposition
  LO gclr,tclr;
  MPI_Comm_size(MPI_COMM_WORLD, &NP);
  MPI_Comm_rank(MPI_COMM_WORLD, &mype);   
  
  //Here: NP=ntube*(NP/ntube)
  gclr=int(mype/ntube);
  tclr=mype%ntube;
  MPI_Comm_split(MPI_COMM_WORLD,gclr,tclr,&grid_comm);
  MPI_Comm_split(MPI_COMM_WORLD,tclr,gclr,&tube_comm);
  MPI_Comm_rank(grid_comm,&mype_g);
  MPI_Comm_rank(tube_comm,&mype_t); 
 
  if(mype_g<kmx-1){
    glk0=mype_g;
    glk1=mype_g;
    glk2=mype_g;
  }else{
    glk0=mype_g;
    glk1=mype_g;
    glk2=mype_g+1;
  }
  //x domain decomposition of GEM 
  if(mype_t!=ntube-1){
    tli0=(imx+1)/ntube;
    tli1=tli0*mype_t;
    tli2=tli1+tli0-1;
  }else{
    tli0=(imx+1)%ntube;
    tli2=imx;
    tli1=imx-tli0+1;
  }  

  //split MPI_COMM_WORLD for GEM-XGC mapping
  npz=int(sqrt(kmx+1));
  while(floor((kmx+1)/npz)<(kmx+1)/npz){
    if((kmx+1)/npz>2) npz++;     
  }
  if((kmx+1)/npz<2) std::cout<<"Error: the number of processes is not chosen right."<<'\n';
  std::exit(1);
  npx=numprocs/npz; 
  npy=1;
 } 
     
void Part1ParalPar3D::decomposeGemMeshforCoupling()
{
  LO n=int((imx+1)/npx);
  if(mype_x<npx-1){
    lk0=n;
    li1=mype_x*n;
    li2=li1+n-1;
  }else{
    li1=mype_x*n;
    li0=imx-li1+2;
    li2=imx;
  }
  n=int((kmx+1)/npz);
  if(mype_z<npz-1){
    lk0=n;
    lk1=mype_z*n;
    lk2=lk1+n-1;
  }else{
    lk1=mype_z*n;
    lk0=kmx-lk1+1;
    lk2=kmx;
  }
  lj0=jmx+1;
} 

//Mapping the rank (mype_x,mype_z) and (mype_g,mype_t)
void Part1ParalPar3D::rankMapping()
{
  mype_xztg=new LO*[NP];
  for(LO i=0;i<NP;i++){ 
    mype_xztg[i]=new LO[4];  
    mype_xztg[i][0]=mype%npx;
    mype_xztg[i][1]=LO(mype/npx);
    mype_xztg[i][2]=mype%ntube;
    mype_xztg[i][3]=LO(mype/ntube);
  }
}

void Part1ParalPar3D::overlapBox()
{
  LO xsendlow[npx],xsendup[npx]; 
  MPI_Allgather(&li1,1,MPI_INT,&xsendlow,1,MPI_INT,comm_x);
  MPI_Allgather(&li2,1,MPI_INT,&xsendup,1,MPI_INT,comm_x);
  getOverlapBox(sendOverlap_x,xsendlow,xsendup,npx,tli1,tli2);
  
  LO thsendlow[npz],thsendup[npz];
  MPI_Allgather(&lk1,1,MPI_INT,&thsendlow,1,MPI_INT,comm_z);
  MPI_Allgather(&lk2,1,MPI_INT,&thsendup,1,MPI_INT,comm_z);  
  getOverlapBox(sendOverlap_th,thsendlow,thsendup,npz,glk1,glk2); 

  LO xrecvlow[ntube],xrecvup[ntube];
  MPI_Allgather(&tli1,1,MPI_INT,&xrecvlow,1,MPI_INT,tube_comm);
  MPI_Allgather(&tli2,1,MPI_INT,&xrecvup,1,MPI_INT,tube_comm);
  getOverlapBox(recvOverlap_x,xrecvlow,xrecvup,ntube,li1,li2);  
   
  LO threcvlow[kmx+1],threcvup[kmx+1];
  MPI_Allgather(&glk1,1,MPI_INT,&threcvlow,1,MPI_INT,grid_comm);
  MPI_Allgather(&glk2,1,MPI_INT,&threcvup,1,MPI_INT,grid_comm);
  getOverlapBox(recvOverlap_th,threcvlow,threcvup,kmx+1,lk1,lk2);  
}

void Part1ParalPar3D::getOverlapBox(vecint2d vec2d,LO* lowind,LO* upind,LO numproc2,LO low,LO up)
{
  LO min,max;
  bool overlap;
//   vecint2d tmp2d;  
  for(LO j=0;j<numproc2;j++){
    vecint1d tmp1d={0,0,0,0};
    overlap=false;
    if(low>upind[j]){
      break;
    }else{
      if(low>lowind[j]){
	overlap=true;
	min=lowind[j];
	if(up>upind[j]){
	   max=upind[j];
	}else{
	   max=low;
	}
      }else{
	if(lowind[j]>up){
	  break;
	}else{
	  overlap=true;
	  min=lowind[j];
	  if(up>upind[j]){
	     max=upind[j];
	  }else{
	     max=up;
	  }
	}          
      }
    }     
    if(overlap==true){ 
      tmp1d[0]=j;
      tmp1d[1]=min;
      tmp1d[2]=max;
      tmp1d[3]=max-min+1;
    }
    vec2d.push_back(tmp1d); 
  } 
}






 /*
  void Part1ParalPar3D::CreateGroupComm()
  {   
    MPI_Group z_group;
    MPI_Comm_group(p1->comm_z,&z_group,);
 
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
    LO* ranksend=new LO[maxrank-minrank+1];
    for(i=0;i<maxrank-minrank+1;i++) ranksend[i]=i+minrank; 
    
    MPI_Group zsend_group;
    MPI_Comm  zsend_comm;
    MPI_Group_incl(z_group,maxrank-minrank+1,ranksend,&zsend_group);
    MPI_Comm_create_group(p1->comm_z,zsend_group,0,&zsend_comm);

    LO* rankrecv=new LO[lk0];
    for(i=0;i<lk0;i++) rankrecv[i]=lk1+i;
    MPI_Group zrecv_group;
    MPI_Comm  zrecv_comm;
    MPI_Group_incl(z_group,,rankredv,&zrecv_group);
    MPI_Comm_create_group(p1->comm_z,zrecv_group,1,&zrecv_comm); 
  } 
*/
}
