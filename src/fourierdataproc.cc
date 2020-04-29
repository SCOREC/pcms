#include "dataprocess.h"

namespace coupler {

// read densin and potentin routines will be written here. 

// Beforing invoking this routine, the FFtw plan routine must be invoked. 
void CmplxdataToRealdata3D(DatasProc3D& dp3d, Part1ParalPar3D& p1pp3d)
{
  if(dp3d.densin==NULL){
    std::cout<<"ERROR: before invoking CmplxdataToRealdata3D, DatasProc3.densin must be allocated"<<'\n';
    std::exit(EXIT_FAILURE);
  }
  if(dp3d.densout==NULL){
    std::cout<<"ERROR: before invoking CmplxdataToRealdata3D, DatasProc3.densout must be allocated"<<'\n';
    std::exit(EXIT_FAILURE);
  }
  GO num;
  if(dp3d.yparal==false){
    for(LO i=0; i<p1pp3d.li0;i++){
      for(LO j=0;j<p1pp3d.lj0;j++){
	for(LO k=0;k<p1pp3d.lk0;k++){
	  num=(GO)(i*p1pp3d.lj0*p1pp3d.lk0+j*p1pp3d.lk0+k+1);
	  dp3d.densintmp[num]=dp3d.densin[i][j][k];
	}      
      }
      ExecuteCmplxToReal(dp3d,p1pp3d);
      for(LO i=0;i<p1pp3d.li0;i++){
	for(LO j=0;j<dp3d.part1lj0;j++){ 
	   for(LO k=0;k<p1pp3d.lk0;k++){ 
	     num=(GO)(i*dp3d.part1lj0*p1pp3d.lk0+j*p1pp3d.lk0+k+1);
	     dp3d.densout[i][j][k]=dp3d.densouttmp[num];          
	   }
	}
  //    delete[] dp3d.denstin;
  //    dp3d.denstin=NULL; 
     }
    }
   }
 }

void RealdataToCmplxdata3D(DatasProc3D& dp3d, Part1ParalPar3D& p1pp3d,Part3Mesh3D& p3m3d)
{
  if(dp3d.potentin==NULL){
    std::cout<<"ERROR: before invoking CmplxdataToRealdata3, DatasProc3D.potentin must be allocated"<<'\n';
    std::exit(EXIT_FAILURE);
  }
  if(dp3d.potentout==NULL){
    std::cout<<"ERROR: before invoking CmplxdataToRealdata3,DatasProc3D.potentout must be allocated"<<'\n';
    std::exit(EXIT_FAILURE);
  }
  GO sum=0;
  if(dp3d.yparal==false){
    for(LO i=0;i<p3m3d.li0;i++){
      for(LO k=0;k<p3m3d.mylk0[i];k++){
	sum+=1;
	for(LO j=0;j<p3m3d.lj0;j++)
	  dp3d.potentintmp[sum*p3m3d.lj0+j+1]=dp3d.potentin[i][j][k];  
      }
    }
    ExecuteRealToCmplx(dp3d,p1pp3d);
    sum=0;
    for(LO i=0;i<p3m3d.li0;i++){
      for(LO k=0;k<p3m3d.mylk0[i];k++){
	sum+=1;
	for(LO j=0;j<dp3d.part3lj0;j++)
	  dp3d.potentout[i][j][k]=dp3d.potentouttmp[sum*dp3d.part1lj0+j+1];      
      }
    }
   }
}

// This routine is not required in first verion of coupler, but would be modifed 
// for in the 2nd version. Here, the indexes may need exchange. 
void TransposeComplex(std::complex<double>** InMatrix,std::complex<double>** OutMatrix, DatasProc3D& dp3d, \
     Part1ParalPar3D& p1pp3d)
{
  std::complex<double>*** sbuf;
  std::complex<double>*** rbuf;
  sbuf=new std::complex<double>**[p1pp3d.ny0];
  rbuf=new std::complex<double>**[p1pp3d.ny0];
  for(int i=0;i<p1pp3d.ny0;i++){
    sbuf[i]=new std::complex<double>*[dp3d.part1li0];
    rbuf[i]=new std::complex<double>*[dp3d.part1li0];
    for(int j=0;j<dp3d.part1li0;j++){
      sbuf[i][j]=new std::complex<double>[p1pp3d.npy];
      rbuf[i][j]=new std::complex<double>[p1pp3d.npy];
    }
  }
}
// It's not finished for dp3d.yparal==true here.
void ExecuteCmplxToReal(DatasProc3D& dp3d, Part1ParalPar3D& p1pp3d)
{
   if(dp3d.yparal==true){
     std::complex<double>** tmp_cmplx;
     tmp_cmplx=new std::complex<double>*[p1pp3d.ny0];
     double** tmp_re; 
     tmp_re=new double*[2*p1pp3d.ny0*p1pp3d.res_fact]; 
     for(LO i=0;i<p1pp3d.ny0;i++)
       tmp_cmplx[i]=new std::complex<double>[dp3d.myli0]; 
     for(LO j=0;j<2*p1pp3d.ny0*p1pp3d.res_fact;j++)
       tmp_re[j]=new double[dp3d.myli0];
     // Here is not finished for yparal=true

    }else{
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

void InitFourierPlan3D( DatasProc3D& dp3d,Part1ParalPar3D& p1pp3d, Part3Mesh3D &p3m3d)
{ 
  if(dp3d.yparal==true){
    dp3d.plan_backward=fftw_plan_dft_c2r_2d(p1pp3d.li0*p1pp3d.lk0, p1pp3d.lj0,\
                       reinterpret_cast<fftw_complex*>(dp3d.densintmp),dp3d.densouttmp,FFTW_BACKWARD); 
    dp3d.plan_forward=fftw_plan_dft_r2c_2d(dp3d.sum, p3m3d.lj0,dp3d.potentintmp, \
                      reinterpret_cast<fftw_complex*>(dp3d.potentouttmp),FFTW_ESTIMATE);
  }
}

void FreeFourierPlan3D(DatasProc3D& dp3d)
{
  fftw_destroy_plan(dp3d.plan_forward);
  fftw_destroy_plan(dp3d.plan_backward);  
  delete[] dp3d.densintmp;
  dp3d.densintmp=NULL;
  delete[] dp3d.densouttmp;
  dp3d.densintmp=NULL;
  delete[] dp3d.potentintmp;
  dp3d.potentintmp=NULL;
  delete[] dp3d.potentouttmp;
  dp3d.potentouttmp=NULL;
}

}  
