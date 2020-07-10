#include "commpart1.h"
#include <cassert>

namespace coupler {

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
     if(mype==0) std::cout<<"i="<<i<<" "<<"pz="<<pzcoords[i]<<'\n';
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


}
