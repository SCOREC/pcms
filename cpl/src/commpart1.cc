#include<coupling1.h>
#include<algorithm>
namespace coupler {

class Part1ParalPar3D {
  public:
    GO mype ;
    GO mype_x;
    GO mype_y;
    GO mype_z;
    GO comm_x;
    GO comm_y;
    GO comm_z;
    GO NP;
    GO npx,nx0,nxb,li0,li1,li2,lg0,lg1,lg2;
    GO npy,ny0,nyb,lh0,lh1,lh2,lm0,lm1,lm2;
    GO npz,nz0 nzb,lk0,lk1,lk2,ln0,ln1,ln2;
    double* xcoords;
    double* pzcoords;
    double dz;
       
}


//read the paralllization parameters
void InitPart1ParalPar3D (Part1ParalPar3D  &p1pp3d){
  if(prepro==true){ 
   receive_field1D(Array1D<GO> &parpar, "../coupling","para_parameters",9,MPI_COMM_WORLD);
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
   p1pp3d.lh0=parpar[12];
   p1pp3d.lh1=parpar[13];
   p1pp3d.lh2=parpar[14];
   p1pp3d.lm0=parpar[15];
   p1pp3d.lm1=parpar[16];
   p1pp3d.lm2=parpar[17];

   p1pp3d.npz=parpar[13];
   p1pp3d.nz0=parpar[14];
   p1pp3d.nzb=parpar[15];
   p1pp3d.lk0=parpar[12];
   p1pp3d.lk1=parpar[13];
   p1pp3d.lk2=parpar[14];
   p1pp3d.ln0=parpar[15];
   p1pp3d.ln1=parpar[16];
   p1pp3d.ln2=parpar[17];
  
   p1pp3d.NP=p1pp3d.npx*p1pp3d.npy*p1pp3d.npz  
   // create 3D parallel cart with z being periodic
   int period[3]={1,1,0};
   int rorder = 1;
   int dim[3]={npx,npy,npz};
   MPI_Comm comm_cart;
   MPI_Cart_create(MPI_COMM_WORLD,3,dim,period,order,&comm_cart);
   
   int remain[3]=0;
   int subcomuni[3];
   for(int i=0;i<2;i++){
     remain[i]=1;
     MPI_Cart_sub(comm_cart,remain,&subcomuni[i]);
   }

   p1pp3d.comm_x=subcomuni[0];
   p1pp3d.comm_y=subcomuni[1];
   p1pp3d.comm_z=subcomuni[2]; 
    
   MPI_Comm_rank(MPI_COMM_WORLD,p1pp3d.mype);
   MPI_Comm_rank(p1pp3d.comm_x,&p1pp3d.mype_x);
   MPI_Comm_rank(p1pp3d.comm_y,&p1pp3d.mype_y);
   MPI_Comm_rank(p1pp3d.comm_z,&p1pp3d.mype_z);
    
   // initialize the radial locations of the flux surface and poloidal angles
   receive_field1D(Array1D<double> &xzcoord, "../coupling","xcoords_dz",p1pp3d.nx0+1,MPI_COMM_WORLD);  
   for(GO i==0;i<p1pp3d.nx0;i++)
     p1pp3d.xcoords[i]=xzcoord[i];
   p1pp3d.dz=xzcoord[nx0]; 
     
   for(int i=0;i<nz0-1;i++){
     p1pp3d.pzcoords[i]=pi_+(double)i*p1pp3d.dz;
   }

}

 
}
