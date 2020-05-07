#include "commpart1.h"
#include <cassert>

namespace coupler {

void Part1ParalPar3D::initTest0(std::string test_dir)
{
  int pe;
  MPI_Comm_rank(MPI_COMM_WORLD,&pe);
  assert(!test_dir.empty());
/*
  LO size;
  LO rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
*/
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
  if(pe==0) fprintf(stderr, " 0.15, npx: %d, npy: %d, npz: %d \n", npx,npy,npz);
  CreateSubCommunicators();
  if(pe==0) fprintf(stderr, " 0.16 \n");
  li1=mype_x*li0;
  if(pe==0) fprintf(stderr, " 0.17 \n");
  li2=li1+li0-1;
  if(pe==0) fprintf(stderr, " 0.18 \n");
  lj1=mype_y*lj0;
  if(pe==0) fprintf(stderr, " 0.19 \n");
  lj2=lj1+lj0-1;
  if(pe==0) fprintf(stderr, " 0.20 \n");
  lk1=mype_z*lk0;
  if(pe==0) fprintf(stderr, " 0.21 \n");
  lk2=lk1+lk0-1;
  if(pe==0) fprintf(stderr, " 0.22 \n");
  
  delete[] data;
}

//read the paralllization parameters
void Part1ParalPar3D::init(std::string test_dir)
{
  int pe;
  MPI_Comm_rank(MPI_COMM_WORLD,&pe);
  if(pe==0) fprintf(stderr, " 0.0 ");

 if(preproc==true){ 
//   LO* parpar=new LO[26];   
   if(test_case==TestCase::t0){
  if(pe==0) fprintf(stderr, " 0.1 ");
     initTest0(test_dir);
  if(pe==0) fprintf(stderr, " 0.3 ");
   }else{
  if(pe==0) fprintf(stderr, " 0.2 ");
//     receive_field1D(GO &parpar, "../coupling","para_parameters",9,MPI_COMM_WORLD);
     npx=parpar->val(0);
     nx0=parpar->val(1);
     nxb=parpar->val(2);
     li0=parpar->val(3);
     li1=parpar->val(4);
     li2=parpar->val(5);
     lg0=parpar->val(6);
     lg1=parpar->val(7);
     lg2=parpar->val(8);
   
     npy=parpar->val(9);
     ny0=parpar->val(10);
     nyb=parpar->val(11);
     lj0=parpar->val(12);
     lj1=parpar->val(13);
     lj2=parpar->val(14);
     lm0=parpar->val(15);
     lm1=parpar->val(16);
     lm2=parpar->val(17);
   
     npz=parpar->val(18);
     nz0=parpar->val(19);
     nzb=parpar->val(20);
     lk0=parpar->val(21);
     lk1=parpar->val(22);
     lk2=parpar->val(23);
     ln0=parpar->val(24);
     ln1=parpar->val(25);
     ln2=parpar->val(26);
  
     NP=npx*npy*npz;  
     CreateSubCommunicators();
   }

   // initialize the radial locations of the flux surface and poloidal angles
   pzcoords=new double[nz0];
   xcoords=new double[nx0];
   if(test_case==TestCase::t0){
      assert(!test_dir.empty());
      std::string fname=test_dir+"xcoords.nml";
      InputfromFile_(xzcoords,nx0,fname);
   }else{
   }

   for(LO i=0;i<nx0;i++){
     xcoords[i]=xzcoords->val(i);
   }
   if(test_case==TestCase::t0){
     dz=2.0*cplPI/nz0;
   }else{
     dz=xzcoords->val(nx0-1); 
   }
 
   for(LO i=0;i<nz0;i++){
     pzcoords[i]=-1.0*cplPI+(double)i*dz;
   }
  delete[] parpar;
  delete[] xzcoords;
 }
}

void Part1ParalPar3D::CreateSubCommunicators()
{
int nprocs;
  MPI_Comm_rank(MPI_COMM_WORLD,&mype);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
MPI_Barrier(MPI_COMM_WORLD);
  npw = 4;//TODO: this should be sent in via ADIOS
   // create 3D parallel cart with z being periodic
   //int rorder = 0;//TODO: Does 0 mean true or false
   int dim[4]={(int)npx,(int)npy,(int)npz, (int)npw};
   if(mype==0) fprintf(stderr," npx: %d, npy: %d, npz: %d, npw: %d \n", npx, npy, npz, npw); 
   MPI_Cart_create(MPI_COMM_WORLD,4,dim,periods,false,&comm_cart);
  if(mype==0) fprintf(stderr, " 0.151 ,nprocs:%d \n", nprocs);
  if(mype==0) fprintf(stderr, " 0.1511 comm_cart:%p \n", comm_cart);

   MPI_Comm subcomuni[4];
  if(mype==0) fprintf(stderr, " 0.152 ");
   for(int i=0;i<4;i++){
     int remain[4]={0,0,0,0};
     remain[i]=1;
     MPI_Cart_sub(comm_cart,remain,&subcomuni[i]);
     if(mype==0) fprintf(stderr, " 0.153 subcomm[i]: %p \n",&subcomuni[i]);
   }
   if(mype==0) fprintf(stderr," subcomun[1]: %p, subcomun[2]: %p, subcomun[3]: %p, subcomun[4]: %p \n",subcomuni[0],subcomuni[1],subcomuni[2],subcomuni[3]); 
  if(mype==0) fprintf(stderr, " 0.154 ");

   comm_x=subcomuni[0];
   comm_y=subcomuni[1];
   comm_z=subcomuni[2];
   comm_w=subcomuni[3];

   if(mype==0) fprintf(stderr," comm_x: %p, comm_y: %p, comm_z: %p, comm_w: %p \n ", comm_x ,comm_y , comm_z,comm_w);
   MPI_Comm_rank(MPI_COMM_WORLD,&mype);
   MPI_Comm_rank(comm_x,&mype_x);
   MPI_Comm_rank(comm_y,&mype_y);
   MPI_Comm_rank(comm_z,&mype_z);
   MPI_Comm_rank(comm_w,&mype_w);
   if(mype==0) fprintf(stderr," mype_x : %d, mype_y : %d, mype_z: %d, mype_w: %d \n ", mype_x, mype_y , mype_z,mype_w);

}

void Part1ParalPar3D::MpiFreeComm()
{
  MPI_Comm_free(&comm_x);
  MPI_Comm_free(&comm_y);
  MPI_Comm_free(&comm_z);
  MPI_Comm_free(&comm_w);
  MPI_Comm_free(&comm_cart);
}


}
