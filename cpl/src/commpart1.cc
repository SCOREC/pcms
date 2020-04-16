#include"coupling.h"
#include<algorithm>
namespace coupler {

  class Part1ParalPar3D {
    public:
      Part1ParalPar3D(GO num_){
      }
      int mype ;
      int mype_x;
      int mype_y;
      int mype_z;
      MPI_Comm comm_x;
      MPI_Comm comm_y;
      MPI_Comm comm_z;
      GO NP;
      GO npx,nx0,nxb,li0,li1,li2,lg0,lg1,lg2;
      GO npy,ny0,nyb,lh0,lh1,lh2,lm0,lm1,lm2;
      GO npz,nz0,nzb,lk0,lk1,lk2,ln0,ln1,ln2;
      double* xcoords;
      double* pzcoords;
      double dz;
  };
  
  
  //read the paralllization parameters
  void InitPart1ParalPar3D (Part1ParalPar3D  &p1pp3d, Array1d<int>* parpar, Array1d<double>* xzcoord){
  	//where is it defined prepro??/
   // if(true){
//     Array1d<int>* parpar = receive_gene_pproc("../coupling","para_parameters",9,MPI_COMM_WORLD);
     p1pp3d.npx=parpar->val(0);
     p1pp3d.nx0=parpar->val(1);
     p1pp3d.nxb=parpar->val(2);
     p1pp3d.li0=parpar->val(3);
     p1pp3d.li1=parpar->val(4);
     p1pp3d.li2=parpar->val(5);
     p1pp3d.lg0=parpar->val(6);
     p1pp3d.lg1=parpar->val(7);
     p1pp3d.lg2=parpar->val(8);
  
     p1pp3d.npy=parpar->val(9);
     p1pp3d.ny0=parpar->val(10);
     p1pp3d.nyb=parpar->val(11);
     p1pp3d.lh0=parpar->val(12);
     p1pp3d.lh1=parpar->val(13);
     p1pp3d.lh2=parpar->val(14);
     p1pp3d.lm0=parpar->val(15);
     p1pp3d.lm1=parpar->val(16);
     p1pp3d.lm2=parpar->val(17);
  
     p1pp3d.npz=parpar->val(18);
     p1pp3d.nz0=parpar->val(19);
     p1pp3d.nzb=parpar->val(20);
     p1pp3d.lk0=parpar->val(21);
     p1pp3d.lk1=parpar->val(22);
     p1pp3d.lk2=parpar->val(23);
     p1pp3d.ln0=parpar->val(24);
     p1pp3d.ln1=parpar->val(25);
     p1pp3d.ln2=parpar->val(26);
  
     p1pp3d.NP=p1pp3d.npx*p1pp3d.npy*p1pp3d.npz;
     // create 3D parallel cart with z being periodic
     int period[3]={1,1,0};
     int order = 1;
     int dim[3]={(int)p1pp3d.npx,(int)p1pp3d.npy,(int)p1pp3d.npz};
     MPI_Comm comm_cart;
     MPI_Cart_create(MPI_COMM_WORLD,3,dim,period,order,&comm_cart);
  
     int remain[3]={0};
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
     // std::vector<double> xzcoord = {0.0};
     // receive_field1D(&xzcoord, "../coupling","xcoords_dz",p1pp3d.nx0+1,MPI_COMM_WORLD);
     for(GO i=0;i<p1pp3d.nx0;i++)
       p1pp3d.xcoords[i]=xzcoord->val(i);
     p1pp3d.dz=xzcoord->val(p1pp3d.nx0);
  
     for(int i=0;i<p1pp3d.nz0-1;i++){
       p1pp3d.pzcoords[i]=cplPI+(double)i*p1pp3d.dz;
     }
  
   // }
  }
}
