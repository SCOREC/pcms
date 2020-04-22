#include "commpart1.h"
namespace coupler {

//read the paralllization parameters
void InitPart1ParalPar3D (Part1ParalPar3D  &p1pp3d)
{
   if(p1pp3d.preproc==true){ 
   GO parpar[26];
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
   // create 3D parallel cart with z being periodic
//   int period[3]={1,1,0};
   int rorder = 1;
   int dim[3]={(int)p1pp3d.npx,(int)p1pp3d.npy,(int)p1pp3d.npz};
   MPI_Comm comm_cart;
   MPI_Cart_create(MPI_COMM_WORLD,3,dim,p1pp3d.periods,rorder,&comm_cart);

   int remain[3]={0,0,0};
   MPI_Comm subcomuni[3];
   for(int i=0;i<2;i++){
     remain[i]=1;
     MPI_Cart_sub(comm_cart,remain,&subcomuni[i]);
   }

   p1pp3d.comm_x=subcomuni[0];
   p1pp3d.comm_y=subcomuni[1];
   p1pp3d.comm_z=subcomuni[2];

   MPI_Comm_rank(MPI_COMM_WORLD,&p1pp3d.mype);
   MPI_Comm_rank(p1pp3d.comm_x,&p1pp3d.mype_x);
   MPI_Comm_rank(p1pp3d.comm_y,&p1pp3d.mype_y);
   MPI_Comm_rank(p1pp3d.comm_z,&p1pp3d.mype_z);

   // initialize the radial locations of the flux surface and poloidal angles
   double* xzcoord; // this lie may be deleted
   //receive_field1D(double& xzcoord, "../coupling","xcoords_dz",p1pp3d.nx0+1,MPI_COMM_WORLD);
   for(GO i=0;i<p1pp3d.nx0;i++){
     p1pp3d.xcoords[i]=xzcoord[i];
   }
   p1pp3d.dz=xzcoord[p1pp3d.nx0]; 
     
   for(int i=0;i<p1pp3d.nz0-1;i++){
     p1pp3d.pzcoords[i]=cplPI+(double)i*p1pp3d.dz;
   }
 } 
}
}
