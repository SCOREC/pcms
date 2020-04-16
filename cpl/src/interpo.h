#include<coupling1.h>

namespace coupler {
  
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










}
