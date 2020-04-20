#include<commpart1.h>
#include<importpart3mesh.h>
#include<BoundExchange.h>

namespace coupler {

 class DatasProc3D{
   public:
     GO part1li0; // part1li0 is the element number of subdomain of x on each process belonging to comm_y. This x subdomain
                // belongs to 2d box on the x-y plane storing the input charged ion density for fourier transform.
     GO part3li0; // part3li0 is the element number of subdomain of x on each process belonging to comm_y. This x subdomain
                // belongs to 2d box on the x-y plane storing the input electrostatic potential.
     GO part1lj0; //the count of elements on y domain on each process after backward Fourier transform
     GO part3lj0; ////the count of elements on y domain on each process after forward Fourier transform
     GO sum;
// here, pointers must be assigned a NULL;
     std::complex<double>*** densin=NULL;  // input 3d density in complex number
     std::complex<double>**  densintmp=NULL;  // temporary 2d density array prepared for backward fourier transform
     double** densouttmp=NULL; // store the x-y 2d real density after backward fourier transform
     double*** densout=NULL;   // store xyz 3d real density
     double*** potentin=NULL;   // the input real electrostatic potential in 3d xyz
     double** potenttmp=NULL;  // temporary xy 2d potential array for forward fourier transform
     std::complex<double>*** potentout=NULL; // 3d temporary array stroring complex electrostatic potential
     bool yparal;
     fftw_plan plan_forward, plan_backward;
// The following is for the interpolation
     double*** denspart3=NULL; // storing the density being sent to the part3
     std::complex<double>*** potentpart1=NULL; // storing the electrostatic potential being sent to the part1.
}

void InitDatasProc3Dparameters(DatasProc3D& dp3d,Part1ParalPar3D& p1pp3d,Part3Mesh3D &p3m3d );

void AllocDatasProc3dDensityArraies(DatasProc3D& dp3d,Part1ParalPar3D& p1pp3d,Part3Mesh3D& p3m3d);

void AllocDatasProc3dPotentArraies(DatasProc3D& dp3d,Part1ParalPar3D& p1pp3d,Part3Mesh3D &p3m3d);


//routines for Fourier transform

void CmplxdataToRealdata3D(DatasProc3D& dp3d, Part1ParalPar3D& p1pp3d);

void RealdataToCmplxdata3D(DatasProc3D& dp3d, Part1ParalPar3D& p1pp3d,Part3Mesh3D &p3m3d);

void TransposeComplex(InMatrix,OutMatrix, DatasProc3D& dp3d, Part1ParalPar3D& p1pp3d);

void ExecuteCmplToReal(DatasProc3D& dp3d, Part1ParalPar3D& p1pp3d);

void ExecuteRealToCmplx(DatasProc3D& dp3d, Part1ParalPar3D& p1pp3d);

void InitFourierPlan3D( DatasProc3D& dp3d);


// routines for interpolation
template<class T>
T Lag3dInterpo1D(const T yin[4],const double xin[4],const double x){
  double l0,l1,l2,l3;
  T yout;
  l0=(x-xin[1])*(x-xin[2])*(x-xin[3])/(xin[0]-xin[1])/(xin[0]-xin[2])/(xin[0]-xin[3]);
  l1=(x-xin[0])*(x-xin[2])*(x-xin[3])/(xin[1]-xin[0])/(xin[1]-xin[2])/(xin[1]-xin[3]);
  l2=(x-xin[0])*(x-xin[1])*(x-xin[3])/(xin[2]-xin[0])/(xin[2]-xin[1])/(xin[2]-xin[3]);
  l3=(x-xin[0])*(x-xin[1])*(x-xin[2])/(xin[3]-xin[0])/(xin[3]-xin[1])/(xin[3]-xin[2]);
  yout=yin[0]*l0+yin[1]*l1+yin[2]*l2+yin[3]*l3;
  return yout;
}

//central 3rd order Lagrange polynormal interpolation
template<class T>
void Lag3dAarray(T* yin,T* xin,GO nin,T* yout,T* xout,GO nout){
       GO jstart=2;
       GO j1=jstart;
       GO j2,j0,jm;
       T x;
       T func[4];
       double coords[4];
       for(GO j=0;j<nout;j++){
         x=xout[j];
         while(x>=xin[j1] && j1<nin-2 && j1>1){
           j1=+1;
         }
         j2=j1+1;
         j0=j1-1;
         jm=j1-2;
         coords[0]=xin[jm];
         coords[1]=xin[j0];
         coords[2]=xin[j1];
         coords[3]=xin[j2];
         func[0]=yin[jm];
         func[1]=yin[j0];
         func[2]=yin[j1];
         func[3]=yin[j3];
         yout[j]=Lag3dInterpo1D(func[4],coords[4],x);
       }

     }

void InterpoDensity3D(BoundaryDescr3D &bdesc,Part3Mesh3D& p3m3d, \  
     Part1ParalPar3D &p1pp3d,DatasProc3D& dp3d);

void InterpoPotential3D(BoundaryDescr3D &bdesc,Part3Mesh3D& p3m3d, \
     Part1ParalPar3D &p1pp3d,DatasProc3D& dp3d);

}
