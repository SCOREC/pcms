#include "commpart1.h"
#include "testutilities.h"

namespace coupler {

//read the paralllization parameters
void InitPart1ParalPar3D (Part1ParalPar3D  &p1pp3d)
{
 if(preproc==true){ 
   LO* parpar=new LO[26];   
   if(test_case==0){
     InitPart1paral3DInCoupler(p1pp3d);
   }else{
//     receive_field1D(GO &parpar, "../coupling","para_parameters",9,MPI_COMM_WORLD);
     p1pp3d.npx=parpar[0];
     p1pp3d.nx0=parpar[1];
     p1pp3d.nxb=parpar[2];
     p1pp3d.li0=parpar[3];
     p1pp3d.li1=parpar[4];
     p1pp3d.li2=parpar[5];
     p1pp3d.lg0=parpar[6];
     p1pp3d.lg1=parpar[7];
     p1pp3d.lg2=parpar[8];

     p1pp3d.npy=parpar[9];
     p1pp3d.ny0=parpar[10];
     p1pp3d.nyb=parpar[11];
     p1pp3d.lj0=parpar[12];
     p1pp3d.lj1=parpar[13];
     p1pp3d.lj2=parpar[14];
     p1pp3d.lm0=parpar[15];
     p1pp3d.lm1=parpar[16];
     p1pp3d.lm2=parpar[17];

     p1pp3d.npz=parpar[18];
     p1pp3d.nz0=parpar[19];
     p1pp3d.nzb=parpar[20];
     p1pp3d.lk0=parpar[21];
     p1pp3d.lk1=parpar[22];
     p1pp3d.lk2=parpar[23];
     p1pp3d.ln0=parpar[24];
     p1pp3d.ln1=parpar[25];
     p1pp3d.ln2=parpar[26];
    
     p1pp3d.NP=p1pp3d.npx*p1pp3d.npy*p1pp3d.npz;  
     CreateSubCommunicators(p1pp3d);
   }

   // initialize the radial locations of the flux surface and poloidal angles
   p1pp3d.pzcoords=new double[p1pp3d.nz0];
   p1pp3d.xcoords=new double[p1pp3d.nx0];
   double* xzcoords; // this lie may be deleted
   xzcoords=new double[p1pp3d.nx0];
   if(test_case==0){
      std::string fname=test_dir+"xcoords.nml";
      InputfromFile(xzcoords,p1pp3d.nx0,fname);
   }else{
   //receive_field1D(double& xzcoord, "../coupling","xcoords_dz",p1pp3d.nx0+1,MPI_COMM_WORLD);
   }

   for(LO i=0;i<p1pp3d.nx0;i++){
     p1pp3d.xcoords[i]=xzcoords[i];
   }
   if(test_case==0){
     p1pp3d.dz=2.0*cplPI/p1pp3d.nz0;
   }else{
     p1pp3d.dz=xzcoords[p1pp3d.nx0]; 
   }
 
   for(LO i=0;i<p1pp3d.nz0;i++){
     p1pp3d.pzcoords[i]=-1.0*cplPI+(double)i*p1pp3d.dz;
   }
  delete[] parpar;
  delete[] xzcoords;
 }
}

void CreateSubCommunicators(Part1ParalPar3D  &p1pp3d)
{
   // create 3D parallel cart with z being periodic
   int rorder = 0;
   int dim[3]={(int)p1pp3d.npx,(int)p1pp3d.npy,(int)p1pp3d.npz};
//   MPI_Comm comm_cart;
   MPI_Cart_create(MPI_COMM_WORLD,3,dim,p1pp3d.periods,rorder,&p1pp3d.comm_cart);

   MPI_Comm subcomuni[3];
   for(int i=0;i<3;i++){
     int remain[3]={0,0,0};
     remain[i]=1;
     MPI_Cart_sub(p1pp3d.comm_cart,remain,&subcomuni[i]);
   }

   p1pp3d.comm_x=subcomuni[0];
   p1pp3d.comm_y=subcomuni[1];
   p1pp3d.comm_z=subcomuni[2];

   MPI_Comm_rank(MPI_COMM_WORLD,&p1pp3d.mype);
   MPI_Comm_rank(p1pp3d.comm_x,&p1pp3d.mype_x);
   MPI_Comm_rank(p1pp3d.comm_y,&p1pp3d.mype_y);
   MPI_Comm_rank(p1pp3d.comm_z,&p1pp3d.mype_z);

}

void MpiFreeComm(Part1ParalPar3D  &p1pp3d)
{
  MPI_Comm_free(&p1pp3d.comm_x);
  MPI_Comm_free(&p1pp3d.comm_y);
  MPI_Comm_free(&p1pp3d.comm_z);
  MPI_Comm_free(&p1pp3d.comm_cart);
}


}
