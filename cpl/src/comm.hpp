#include<adios2.h>
#include<iostream>
#include<mpi.h>
#include<cassert>
#include<typeinfo>

class part1par3D {
  public:
    int mype ;
    int mype_x;
    int mype_y;
    int mype_z;
    int Np;
    int npx,nx0,nxb,li0,li1,li2,lg0,lg1,lg2;
    int npy,ny0,nyb,lh0,lh1,lh2,lm0,lm1,lm2;
    int npz,nz0 nzb,lk0,lk1,lk2,ln0,ln1,ln2;
    int* x_part1 
    int* z_part1
    
    void InitParalBox();
    


}

   void InitParaBox(){
       

}
