#include <coupling1.h>
#include <complex>

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

 void InitDatasProc3Dparameters(DatasProc3D& dp3d,Part1ParalPar3D& p1pp3d,Part3Mesh3D &p3m3d ){
  if(para.prepro==true){
   if(dp3d.yparal==true){
     if(p1pp3d.li0%p1pp3d.npy==0){
       dp3d.part1li0=p1pp3d.li0/p1pp3d.npy;
       dp3d.part3li0=dp3d.part1li0;
     } else{
        if(p1pp3d.my_pey==p1pp3d.npy-1){
          dp3d.part1li0=p1pp3d.li0%p1pp3d.npy;
          dp3d.part3li0=dp3d.part1li0;
        }  else{
	    dp3d.part1li0=p1pp3d.li0%p1pp3d.npy;
            dp3d.part3li0=dp3d.part1li0;
        }
     }
    dp3d.part1lj0=2*p1pp3d.nj0;
    dp3d.part3lj0=p1pp3d.nj0;   // here may need rethinking. 
   } else{
     dp3d.part1li0=p1pp3d.li0;
     dp3d.part3li0=p1pp3d.li0;
     dp3d.part1lj0=2*p1pp3d.nj0;
     dp3d.part3lj0=p1pp3d.nj0;   // here may need rethinking.
   }
 }
  for(int i=0;i<p3m3d.li0;i++)  dp3d.sum=+p3m3d.myli0(i); 
}

void DensityAllocDatasProc3dArraies(DatasProc3D& dp3d,Part1ParalPar3D& p1pp3d,Part3Mesh3D& p3m3d){
  dp3d.densin=new std::complex<double>**[p1pp3d.li0];
  for(int i=0;i<p1pp3d.li0;i++){
    dp3d.densin[i]=new std::complex<double>*[p1pp3d.lj0];
      for(int j=0;j<p1pp3d.lj0;j++)
        dp3d.densin[i][j]=new std::complex<double>[p1pp3d.lk0];
  }
  dp3d.densintmp=new std::complex<double>*[dp3d.part1li0];
  for(int k=0;k<dp3d.part1li0;k++){
      dp3d.densintmp[k]=new std::complex<double>[p1pp2d.nj0]; 
    }
  dp3d.densouttmp=new double*[dp3d.part1li0];
  for(int k=0;k<dp3d.part1li0;k++)
      dp3d.densouttmp[k]=new double[dp3d.part1lj0];
 
  dp3d.denspart3=new double**[p3m3d.li0];
  for(int i=0;i<p3m3d.li0;i++){
    dp3d.denspart3[i]=new double*[p3m3d.lj0]
    for(int j=0; j<p3m3d.lj0; j++)
      dp3d.denspart3[i][j]=new double*[p3m3d.mylk0[i]];
  }  

}

void AllocDatasProc3dPotentArraies(DatasProc3D& dp3d,Part1ParalPar3D& p1pp3d,Part3Mesh3D &p3m3d){
  dp3d.potentin=new double**[p3m3d.li0];
  for(int i=0;i<p3m3d.li0;i++){
    dp3d.potentin[i]=new double*[p3m3d.lj0];
      for(int j=0;j<p3m3d.lj0;j++)
        dp3d.densin[i][j]=new double[p3m3d.mylk0[i]];
  }
  dp3d.potenttmp=new double*[dp3d.sum];
  for(int k=0;k<dp3d.sum;k++){
      dp3d.potenttmp[k]=new double[p3m3d.lj0]; 
    }
  dp3d.potentout=new std:complex<double>*[p3d3d.li0];
  for(int i=0;i<p3m3d.li0;i++){
    dp3d.potentout[i]=new std:complex<double>[dp3d.part3lj0]; 
      for(int j=0;j<p3d3d.lj0;j++)
         dp3d.potentout[i][j]=new std:complex<double>[p3m3d.mylk0[i]];
  }   

  dp3d.potentpart1=new std:complex<double>**[p1pp3d.li0];
  for(int i=0;i<p1pp3d.li0;i++){
    dp3d.potentpart1[i]=new std::complex<double>*[p1pp3d.lj0];
    for(int j=0;j<p1pp3d.lj0;j++){
      dp3d.potentpart1[i][j]=new std::complex<double>[p1pp3d.lk0];
    }
  }

} 
// read densin and potentin routines will be written here. 

// Beforing invoking this routine, the FFtw plan routine must be invoked. 
void CmplxdataToRealdata3D(DatasProc3D& dp3d, Part1ParalPar3D& p1pp3d){
  if(dp3d.densin==NULL){
    std::cout<<"ERROR: before invoking CmplxdataToRealdata3D, DatasProc3.densin must be allocated"<<'\n';
    std::exit();
  }
  if(dp3d.densout==NULL){
    std::cout<<"ERROR: before invoking CmplxdataToRealdata3D, DatasProc3.densout must be allocated"<<'\n';
    std::exit();
  }
  for(int k=0; k<p1pp3d.lk0;k++){
    for(int i=0;i<p1pp3d.li0;i++){
      for(int j=0;j<p1pp3d.lj0;j++)
        dp3d.densintmp[i][j]=dp3d.densin[i][j][k];      
    }
    ExecuteCmplToReal(dp3d,p1pp3d)
    for(int i=0;i<dp3d.part1li0;i++){
      for(int j=0;j<dp3d.part1lj0;j++)   
        dp3d.densout[i][j][k]=dp3d.densouttmp[i][j];
    }
  }
//    delete[] dp3d.denstin;
//    dp3d.denstin=NULL; 
}

void CmplxdataToRealdata3D(DatasProc3D& dp3d, Part1ParalPar3D& p1pp3d,Part3Mesh3D &p3m3d){
  if(dp3d.potentin==NULL){
    std::cout<<"ERROR: before invoking CmplxdataToRealdata3, DatasProc3D.potentin must be allocated"<<'\n';
    std::exit();
  }
  if(dp3d.potentout==NULL){
    std::cout<<"ERROR: before invoking CmplxdataToRealdata3,DatasProc3D.potentout must be allocated"<<'\n';
    std::exit();
  }
  GO sum=0;
  for(GO i=0;i<p3m3d.li0){
    for(GO k=0;k<p3m3d.mylk0(i);k++){
      sum=+1;
      for(GO j=0;j<p3m3d.lj0-1;j++)
        dp3d.potenttmp[sum,j]=dp3d.potentin[i][j][k];  
    }
  }
  ExecuteRealToCmplx(dp3d,p1pp3d);
  sum=0;
  for(GO i=0;i<p3m3d.li0){
    for(GO k=0;k<p3m3d.mylk0(i);k++){
      sum=+1;
      for(GO j=0;j<p3m3d.lj0-1;j++)
        dp3d.potentout[i][j][k]=dp3d.potenttmp[sum,j];      
    }
  }
}

// This routine is not required in first verion of coupler, but would be modifed 
// for in the 2nd version. Here, the indexes may need exchange. 
 void TransposeComplex(InMatrix,OutMatrix, DatasProc3D& dp3d, Part1ParalPar3D& p1pp3d){
  std::complex<double>*** sbuf,rbufc;
  sbuf=new std::complex<double>**[p1pp3d.nj0];
  rbufc=new std::complex<double>**[p1pp3d.nj0];
  for(int i=0;i<p1pp3d.nj0;i++){
    sbuf[i]=new std::complex<double>*[dp3d.part1li0];
    rbufc[i]=new std::complex<double>*[dp3d.part1li0];
    for(int j=0;j<dp3d.part1li0;j++){
      sbuf[i][j]=new std::complex<double>*[p1pp3d.npy];
      rbufc[i][j]=new std::complex<double>*[p1pp3d.npy];
    }
  }

// It's not finished for dp3d.yparal==true here.
 void ExecuteCmplToReal(DatasProc3D& dp3d, Part1ParalPar3D& p1pp3d){
   if(dp3d.yparal==true){
     std::complex<double>** tmp_cmplx;
     tmp_cmpl=new std::complex<double>[p1pp3d.nj0];
     double** tmp_re 
     tmp_re=new double*[2*p1pp3d.nj0*para.res_fact]; 
     for(GO i=0;i<dp3d.nj0-1;i++)
       tmp_cmpl[i]=new std::complex<double>[dp3d.myli0]; 
     for(GO j=0;j<2*p1pp3d.nj0*para.res_fact-1;j++)
       tmp_re[j]=new double[dp3d.myli0];
     // Here is not finished for yparal=true

     } else{
     fftw_execute(dp3d.plan_backward);
    } 
 }

 void ExecuteRealToCmplx(DatasProc3D& dp3d, Part1ParalPar3D& p1pp3d){
   if(dp3d.yparal==true){
     // Here is not finished for yparal=true
   } else{
     fftw_execute(dp3d.plan_forward);
   }

 } 

 void InitFourierForwardPlan( DatasProc3D& dp3d){

    fftpa.backward=fft_plan_dft_c2r_2d(fftpa.par1li0, part1nj0, fftpa.densintmp, fftpa.densouttmp,FFTW_BACKWARD);

    fftpa.forward=fft_plan_dft_r2c_2d(fftpa.sum, part3nj0, fftpa.potentmp, fftpa.potentout,FFTW_ESTIMATE);
}

}  
