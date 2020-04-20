#include <coupling1.h>
#include <dataprocess.h>
#include <complex>

namespace coupler {

// read densin and potentin routines will be written here. 

// Beforing invoking this routine, the FFtw plan routine must be invoked. 
void CmplxdataToRealdata3D(DatasProc3D& dp3d, Part1ParalPar3D& p1pp3d)
{
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
    ExecuteCmplToReal(&dp3d,&p1pp3d)
    for(int i=0;i<dp3d.part1li0;i++){
      for(int j=0;j<dp3d.part1lj0;j++)   
        dp3d.densout[i][j][k]=dp3d.densouttmp[i][j];
    }
  }
//    delete[] dp3d.denstin;
//    dp3d.denstin=NULL; 
}

void RealdataToCmplxdata3D(DatasProc3D& dp3d, Part1ParalPar3D& p1pp3d,Part3Mesh3D &p3m3d)
{
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
  ExecuteRealToCmplx(&dp3d,&p1pp3d);
  sum=0;
  for(GO i=0;i<p3m3d.li0,i++){
    for(GO k=0;k<p3m3d.mylk0(i);k++){
      sum=+1;
      for(GO j=0;j<p3m3d.lj0-1;j++)
        dp3d.potentout[i][j][k]=dp3d.potenttmp[sum,j];      
    }
  }
}

// This routine is not required in first verion of coupler, but would be modifed 
// for in the 2nd version. Here, the indexes may need exchange. 
void TransposeComplex(InMatrix,OutMatrix, DatasProc3D& dp3d, Part1ParalPar3D& p1pp3d)
{
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
void ExecuteCmplToReal(DatasProc3D& dp3d, Part1ParalPar3D& p1pp3d)
{
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

void ExecuteRealToCmplx(DatasProc3D& dp3d, Part1ParalPar3D& p1pp3d)
{
   if(dp3d.yparal==true){
     // Here is not finished for yparal=true
   } else{
     fftw_execute(dp3d.plan_forward);
   }

 } 

void InitFourierPlan3D( DatasProc3D& dp3d)
{

    fftpa.backward=fft_plan_dft_c2r_2d(dp3d.par1li0, dp3d.part1nj0, reinterpret_cast<fftw_complex*>(dp3d.densintmp),\
                    reinterpret_cast<fftw_complex*>(dp3d.densouttmp),FFTW_BACKWARD);

    fftpa.forward=fft_plan_dft_r2c_2d(dp3d.sum, dp3d.part3nj0, reinterpret_cast<fftw_complex*>(dp3d.potenttmp), \
                   reinterpret_cast<fftw_complex*>(dp3d.potentout),FFTW_ESTIMATE);
}

void FreeFourierPlan3D(DatasProc3D& dp3d)
{
  fft_destroy_plan(dp3d.plan_forward);
  fft_destroy_plan(dp3d.plan_backward);  
  delete[] dp3d.densintmp;
  delete[] dp3d.densouttmp;
  delete[] dp3d.potenttmp;
  delete[] dp3d.potentout;

}  
