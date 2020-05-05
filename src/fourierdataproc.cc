#include "dataprocess.h"
#include "commpart1.h"
#include "importpart3mesh.h"
#include "couplingTypes.h"
#include <iostream>

namespace coupler {

// read densin and potentin routines will be written here. 

// Beforing invoking this routine, the FFtw plan routine must be invoked. 
void DatasProc3D::CmplxdataToRealdata3D(const Part1ParalPar3D& p1pp3d)
{
  if(densin==NULL){
    std::cout<<"ERROR: before invoking CmplxdataToRealdata3D, DatasProc3.densin must be allocated"<<'\n';
    std::exit(EXIT_FAILURE);
  }
  if(densout==NULL){
    std::cout<<"ERROR: before invoking CmplxdataToRealdata3D, DatasProc3.densout must be allocated"<<'\n';
    std::exit(EXIT_FAILURE);
  }
  GO num;
  if(yparal==false){
    for(LO i=0; i<p1pp3d.li0;i++){
      for(LO j=0;j<p1pp3d.lj0;j++){
	for(LO k=0;k<p1pp3d.lk0;k++){
	  num=(GO)(i*p1pp3d.lj0*p1pp3d.lk0+j*p1pp3d.lk0+k+1);
	  densintmp[num]=densin[i][j][k];
	}      
      }
      ExecuteCmplxToReal(p1pp3d);
      for(LO i=0;i<p1pp3d.li0;i++){
	for(LO j=0;j<part1lj0;j++){ 
	   for(LO k=0;k<p1pp3d.lk0;k++){ 
	     num=(GO)(i*part1lj0*p1pp3d.lk0+j*p1pp3d.lk0+k+1);
	     densout[i][j][k]=densouttmp[num];          
	   }
	}
  //    delete[] denstin;
  //    denstin=NULL; 
     }
    }
   }
 }

void DatasProc3D::RealdataToCmplxdata3D(const Part1ParalPar3D& p1pp3d, const Part3Mesh3D& p3m3d)
{
  if(potentin==NULL){
    std::cout<<"ERROR: before invoking CmplxdataToRealdata3, DatasProc3D.potentin must be allocated"<<'\n';
    std::exit(EXIT_FAILURE);
  }
  if(potentout==NULL){
    std::cout<<"ERROR: before invoking CmplxdataToRealdata3,DatasProc3D.potentout must be allocated"<<'\n';
    std::exit(EXIT_FAILURE);
  }
  GO sum=0;
  if(yparal==false){
    for(LO i=0;i<p3m3d.li0;i++){
      for(LO k=0;k<p3m3d.mylk0[i];k++){
	sum+=1;
	for(LO j=0;j<p3m3d.lj0;j++)
	  potentintmp[sum*p3m3d.lj0+j+1]=potentin[i][j][k];  
      }
    }
    ExecuteRealToCmplx(p1pp3d);
    sum=0;
    for(LO i=0;i<p3m3d.li0;i++){
      for(LO k=0;k<p3m3d.mylk0[i];k++){
	sum+=1;
	for(LO j=0;j<part3lj0;j++)
	  potentout[i][j][k]=potentouttmp[sum*part1lj0+j+1];      
      }
    }
   }
}

//TODO Incomplete, Not used, put into the class it supports, I assume DatasProc3D
// This routine is not required in first verion of coupler, but would be modifed 
// for in the 2nd version. Here, the indexes may need exchange. 
void TransposeComplex(CV** InMatrix,CV** OutMatrix, DatasProc3D& dp3d,
     Part1ParalPar3D& p1pp3d)
{
  CV*** sbuf;
  CV*** rbuf;
  sbuf=new CV**[p1pp3d.ny0];
  rbuf=new CV**[p1pp3d.ny0];
  for(int i=0;i<p1pp3d.ny0;i++){
    sbuf[i]=new CV*[dp3d.part1li0];
    rbuf[i]=new CV*[dp3d.part1li0];
    for(int j=0;j<dp3d.part1li0;j++){
      sbuf[i][j]=new CV[p1pp3d.npy];
      rbuf[i][j]=new CV[p1pp3d.npy];
    }
  }
}
// It's not finished for yparal==true here.
void DatasProc3D::ExecuteCmplxToReal(const Part1ParalPar3D& p1pp3d)
{
   if(yparal==true){
     CV** tmp_cmplx;
     tmp_cmplx=new CV*[p1pp3d.ny0];
     double** tmp_re; 
     tmp_re=new double*[2*p1pp3d.ny0*p1pp3d.res_fact]; 
     for(LO i=0;i<p1pp3d.ny0;i++)
       tmp_cmplx[i]=new CV[myli0]; 
     for(LO j=0;j<2*p1pp3d.ny0*p1pp3d.res_fact;j++)
       tmp_re[j]=new double[myli0];
     // Here is not finished for yparal=true

    }else{
     fftw_execute(plan_backward);
    } 
 }

void DatasProc3D::ExecuteRealToCmplx(const Part1ParalPar3D&)
{
   if(yparal==true){
     // Here is not finished for yparal=true
   } else{
     fftw_execute(plan_forward);
   }

 } 

void DatasProc3D::InitFourierPlan3D(const Part1ParalPar3D& p1pp3d, const Part3Mesh3D &p3m3d)
{ 
  if(yparal==true){
    plan_backward=fftw_plan_dft_c2r_2d(p1pp3d.li0*p1pp3d.lk0, p1pp3d.lj0,
                       reinterpret_cast<fftw_complex*>(densintmp),densouttmp,FFTW_BACKWARD); 
    plan_forward=fftw_plan_dft_r2c_2d(sum, p3m3d.lj0,potentintmp,
                      reinterpret_cast<fftw_complex*>(potentouttmp),FFTW_ESTIMATE);
  }
}

void DatasProc3D::FreeFourierPlan3D()
{
  if(plan_forward)
    fftw_destroy_plan(plan_forward);
  if(plan_backward)
    fftw_destroy_plan(plan_backward);  
}

}  
