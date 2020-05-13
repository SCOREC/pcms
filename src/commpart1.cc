#include "commpart1.h"
#include <cassert>

namespace coupler {

void Part1ParalPar3D::initTest0(std::string test_dir)
{
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
  CreateSubCommunicators();
  li1=mype_x*li0;
  li2=li1+li0-1;
  lj1=mype_y*lj0;
  lj2=lj1+lj0-1;
  lk1=mype_z*lk0;
  lk2=lk1+lk0-1;
  
  delete[] data;
}

//read the paralllization parameters
void Part1ParalPar3D::init(std::string test_dir)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

 if(preproc==true){ 
   if(test_case==TestCase::t0){
     initTest0(test_dir);
   }else{
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
if(!rank) fprintf(stderr," npx: %d, npy: %d, npz: %d \n", npx, npy, npz);  
  
     NP=npx*npy*npz;  
     CreateSubCommunicators();
   }
   if(!rank) fprintf(stderr,"0.11 \n");

   // initialize the radial locations of the flux surface and poloidal angles
   pzcoords=new double[nz0];
   xcoords=new double[nx0];

   if(!rank) fprintf(stderr,"0.12 \n");
   // If running the test case, take this route
   if(test_case==TestCase::t0){
      double* xzcoord; // this lie may be deleted
      xzcoord=new double[nx0];
      assert(!test_dir.empty());
      std::string fname=test_dir+"xcoords.nml";
      InputfromFile(xzcoord,nx0,fname);
      for(LO i=0;i<nx0;i++){
        xcoords[i]=xzcoord[i];
      }
      if(test_case==TestCase::t0){
        dz=2.0*cplPI/nz0;
      }else{
        dz=xzcoord[nx0-1]; 
      }
                                             
      for(LO i=0;i<nz0;i++){
        pzcoords[i]=-1.0*cplPI+(double)i*dz;
      }
   // if running the GENE_cuth test, take this route
   }else{
   if(!rank) fprintf(stderr,"0.13\n");
      for(LO i=0;i<nx0;i++){
        xcoords[i]=xzcoords->val(i);
      }
      if(test_case==TestCase::t0){
        dz=2.0*cplPI/nz0;
      }else{
        dz=xzcoords->val(nx0-1);
      }
                                                
   if(!rank) fprintf(stderr,"0.14\n");
      for(LO i=0;i<nz0;i++){
        pzcoords[i]=-1.0*cplPI+(double)i*dz;
      }
   }
   
   if(!rank) fprintf(stderr,"0.15\n");
  destroy(parpar);
   if(!rank) fprintf(stderr,"0.16\n");
  destroy(xzcoords);
   if(!rank) fprintf(stderr,"0.17\n");
 }
}

void Part1ParalPar3D::CreateSubCommunicators()
{
   // create 3D parallel cart with z being periodic
   int rorder = 0;
   int dim[3]={(int)npx,(int)npy,(int)npz};
   MPI_Cart_create(MPI_COMM_WORLD,3,dim,periods,rorder,&comm_cart);

   MPI_Comm subcomuni[3];
   for (int i=0;i<3;i++){
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


}
