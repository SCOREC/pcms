#ifndef CLASSES_H
#define CLASSES_H

#include "coupling.h"

namespace coupler {

class Part1ParalPar3D {
  public:
    LO mype ; // my process rank in mpi_comm_world
    LO mype_x; // my process rank in comm_x
    LO mype_y; // my process rank in comm_y
    LO mype_z; // my process rank in comm_z
    MPI_Comm comm_x;
    MPI_Comm comm_y;
    MPI_Comm comm_z;
    MPI_Comm comm_cart;
    LO  NP; // The total number of processes
    LO npx,npy,npz;
    LO nx0,nxb,li0,li1,li2,lg0,lg1,lg2;
    LO ny0,nyb,lj0,lj1,lj2,lm0,lm1,lm2;
    LO nz0,nzb,lk0,lk1,lk2,ln0,ln1,ln2;
//    LO myli0,mylj1,myl12;  // The indexes of box y after Fourier transform
    int periods[3]={0,1,1};
    double* xcoords=NULL; // The 1d array storing the radial position of all flux surfaces
    double* pzcoords=NULL; // The 1d array storing the poloidal angle of all vertices along the poloidal surface curve.
    double* pzp=NULL; // The 1d array storing the poloial on each process.
    double dz;  // The equal step length along the poloidal flux curve.
 // parameters for creating the magnetic field, density and temperature ground.
    LO res_fact;
    ~Part1ParalPar3D()
    {
      if(xcoords!=NULL)  delete[] xcoords;
      if(pzcoords!=NULL) delete[] pzcoords;
      if(pzp!=NULL)      delete[] pzp;
    }
};

class Part3Mesh3D{
  public:
    LO  nsurf;    // number of flux surfaces
    LO* versurf=NULL; // numbers of vertice on the flux surfaces
    double* xcoords=NULL;
    LO  li0,li1,li2;
    LO** xboxinds=NULL;  //The indexes of all boxes on the radial dimension
    LO lj0;
    LO* mylk0=NULL;
    LO* mylk1=NULL;
    LO* mylk2=NULL; // The indexes of box along z dimension
    GO boxstar,boxend,boxcount; // The  indexes of the 2d box

    double** Rcoords=NULL;  // The R coordinate of all vertices within the 2d box
    double** Zcoords=NULL;  // The Z coordinate of all vertices within the 2d box
    double** pzcoords=NULL;  // The z coordinates of all points with the 2d box.
    ~Part3Mesh3D()
   {
     if(versurf!=NULL) delete[] versurf;
     if(xboxinds!=NULL) delete[] xboxinds;
     if(xcoords!=NULL) delete[] xcoords;
     if(mylk0!=NULL) delete[] mylk0;
     if(mylk1!=NULL) delete[] mylk1;
     if(mylk2!=NULL) delete[] mylk2;
     if(Rcoords!=NULL) delete[] Rcoords;
     if(Zcoords!=NULL) delete[] Zcoords;
     if(pzcoords!=NULL) delete[] pzcoords;
   }
};

 class DatasProc3D{
   public:
     LO part1li0; // part1li0 is the element number of subdomain of x on each process belonging to comm_y. This x subdomain
                // belongs to 2d box on the x-y plane storing the input charged ion density for fourier transform.
     LO part3li0; // part3li0 is the element number of subdomain of x on each process belonging to comm_y. This x subdomain
                // belongs to 2d box on the x-y plane storing the input electrostatic potential.
     LO part1lj0; //the count of elements on y domain on each process after backward Fourier transform
     LO part3lj0; ////the count of elements on y domain on each process after forward Fourier transform


     LO sum;
// here, pointers must be assigned a NULL;
     std::complex<double>*** densin=NULL;  // input 3d density in complex number
     std::complex<double>*  densintmp=NULL;  // temporary 2d density array prepared for backward fourier transform
     double* densouttmp=NULL; // store the x-y 2d real density after backward fourier transform
     double*** densout=NULL;   // store xyz 3d real density
     double*** denspart3=NULL; // storing the density being sent to the part3
     double*** potentin=NULL;   // the input real electrostatic potential in 3d xyz
     double*** potentinterpo=NULL; // 3d temporary array stroring real electrostatic potential after interpolation
     double* potentintmp=NULL;  // temporary xy 2d potential array for forward fourier transform
     std::complex<double>* potentouttmp=NULL; //
     std::complex<double>*** potentpart1=NULL; // storing the electrostatic potential being sent to the part1.
     bool yparal=false;
     fftw_plan plan_forward, plan_backward;
//The following parameters for yparal=true;
     LO myli0;
     ~DatasProc3D();
};

class BoundaryDescr3D{
  public:
    LO nzb;
    double** upzpart3=NULL;
    double** lowzpart3=NULL;
    double*** updenz=NULL; // The upper  boundary buffer on z domain for interpolation and storing the real quantiies resulted from the backward Fourier transform of complex charged density.
    double*** lowdenz=NULL;
    double*** uppotentz=NULL; //The upper  boundary buffer on z domain for interpolation and storing the complex  quantiies resulted from the forward Fourier transform of electrosttic potential.
    double*** lowpotentz=NULL;
    ~BoundaryDescr3D();
};



}

#endif
