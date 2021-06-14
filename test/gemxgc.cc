#include "adios2Routines.h"
#include "couplingTypes.h"
#include "dataprocess.h"
#include "importpart3mesh.h"
#include "commpart1.h"
#include "BoundaryDescr3D.h"
#include "testutilities.h"
#include <string>
#include <fstream>

int main(int argc, char **argv){
  int rank;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  const std::string dir = "../coupling";
  const int time_step = atoi(argv[1]), RK_count = 4;

  fprintf(stderr, "%d number of time steps: %d\n", rank, time_step);

  const std::string xmlfile = "adios2cfg.xml";
  adios2::ADIOS adios(MPI_COMM_WORLD, adios2::DebugON);
  adios2::Variable<double> senddensity;
  adios2::Variable<double> sendfield;

  coupler::adios2_handler gDens(adios,"gem_density");
  coupler::adios2_handler cDens(adios,"cpl_density");
  coupler::adios2_handler xFld(adios,"xgc_field");
  coupler::adios2_handler cFld(adios,"cpl_field");
  coupler::adios2_handler gMesh(adios,"gem_mesh");
  coupler::adios2_handler gThf(adios,"gem_thfl");
  coupler::adios2_handler gGrd(adios,"gem_grid");
  coupler::adios2_handler gQprof(adios,"gem_qprf");
  coupler::adios2_handler xCouple(adios,"xgc_mesh");
  coupler::adios2_handler xRzcoords(adios,"xgc_rzcoords");  

  std::string model="global";
  int m =0;
  coupler::GO start[2] = {0,0};
  const coupler::Array1d<int>* gmesh=coupler::receive_pproc<int>(dir,gMesh,model);
  coupler::GO ntheta = (coupler::GO)gmesh->val(4);
  coupler::GO nr = (coupler::GO)gmesh->val(1);
  coupler::GO nnode = (coupler::GO)gmesh->val(5);
  printf("%d,%d,%d,%d,%d,%d \n", gmesh->val(0),gmesh->val(1),gmesh->val(2),
            gmesh->val(3), gmesh->val(4), gmesh->val(5));
  printf("mype: %d \n", rank);
  MPI_Barrier(MPI_COMM_WORLD);
/*
  coupler::GO count[2] = {2,(coupler::GO)gmesh->val(5)};
  coupler::GO count2[2] = {0,10};
*/
  const coupler::Array1d<double>* thfl_qprof=coupler::receive_pproc<double>(dir,gQprof,model);
  for(coupler::LO i=0; i<2; i++) if(!rank) fprintf(stderr,"array[%d]: %f\n",i,thfl_qprof->val(i));

  coupler::Array1d<double>* rzcoords = coupler::receive_pproc<double>(dir,gGrd,model);
  if(!rank) fprintf(stderr, "rzcoords[2*nnode-2]: %f, rzcoords[2*nnode-1]: %f \n", 
          rzcoords->val(2*nnode-2),rzcoords->val(2*nnode-1));

  coupler::Array1d<coupler::LO>* xcouple = coupler::receive_pproc<coupler::LO>(dir, xCouple, model);
  if(!rank) fprintf(stderr, "xcouple[5]: %d, xcouple[6]: %d \n", xcouple->val(5), xcouple->val(6));

  //intialize GEM class
  const bool preproc = true;
  const bool ypar = false;
  coupler::TestCase test_case = coupler::TestCase::off;
  coupler::CouplingCase ccase = coupler::CouplingCase::gemxgc; 
  
  coupler::Part1ParalPar3D p1pp3d(gmesh, thfl_qprof, test_case, preproc);
  MPI_Barrier(MPI_COMM_WORLD);
  coupler::Part1ParalPar3D* p1 = &p1pp3d;  

  coupler::Part3Mesh3D p3m3d(p1, xcouple, rzcoords, preproc, test_case);
  MPI_Barrier(MPI_COMM_WORLD);
//  fprintf(stderr, "sxz 00 \n");
  
  coupler::BoundaryDescr3D bdesc(p3m3d, p1pp3d, ccase, test_case, preproc);  
 
  coupler::Part3Mesh3D*     p3 = &p3m3d;
  coupler::BoundaryDescr3D* bound = &bdesc;

  coupler::gemXgcDatasProc3D  gxdp3d(p1,p3,bound,preproc,test_case,ypar);
  coupler::gemXgcDatasProc3D  gx3d_ = gxdp3d; 

  coupler::GO* start_adios = new coupler::GO[3];
  coupler::GO* count = new coupler::GO[3];
  coupler::GO dim0, dim1, dim2;
      
  coupler::GO INDS1d[2]={0,p3m3d.lj0*p3m3d.blockcount};

  coupler::GO globcount[3] = {p1->imx+1, p1->jmx+1, p1->kmx+1};
  coupler::GO loccount[3] = {p1->tli0, p1->lj0, p1->glk0};
  coupler::GO start3d[3] = {p1->tli1, 0, p1->glk1};

  coupler::Array3d<double>* densityfromGEM;
  coupler::Array2d<double>* densitytoXGC; 
  coupler::Array2d<double>* fieldfromXGC; 
  coupler::Array3d<double>* fieldtoGEM;
 
  double* densitytmp;
  double* fieldgem;  
  bool debug = false;
  start_adios[0] = p1->glk1;
  start_adios[1] = 0;
  start_adios[2] = p1->tli1;
  count[0] = p1->glk0;
  count[1] = p1->lj0;
  count[2] = p1->tli0;      

  coupler::GO start_1[2]={0, p3m3d.blockstart+p3m3d.cce_first_node-1};
  coupler::GO count_1[2]={coupler::GO(p3m3d.nphi), p3m3d.blockcount};
 

  for (int i = 0; i < time_step; i++) {
    for (int j = 0; j < RK_count; j++) {
      m = i*RK_count+j;
      if(p1->mype == 0) printf("m equals %d \n", m);
     // receive density from GEM to coupler
      densityfromGEM = coupler::receive_pproc_3d<double>(dir, gDens, start_adios, count, m, MPI_COMM_WORLD); 
      
      debug = true;
      if( debug ) {
         coupler::printminmax1d(densityfromGEM->data(), count[0]*count[1]*count[2], p1->mype, 
        "densityfromGEM", m, false);
/*
          double* tmp = densityfromGEM->data();
          if (p1->mype == 4) {
            for (coupler::LO k=0; k< count[0]*count[1]*count[2]; k++) {
              printf("k: %d, receve density: %f \n", k, tmp[k]);
            }
          }
*/
     }

      gxdp3d.DistriDensiRecvfromGem(densityfromGEM);
      MPI_Barrier(MPI_COMM_WORLD);
      densitytoXGC = new coupler::Array2d<double>(
                                  p3m3d.activenodes,p3m3d.lj0,p3m3d.blockcount,p3m3d.lj0,
                                  p3m3d.blockstart);
      densitytmp = densitytoXGC->data();
      for(int h=0; h<p3m3d.nphi*p3m3d.blockcount; h++){
        densitytmp[h] = gxdp3d.denssend[h];
       if(isnan(densitytmp[h])) {
         printf("h: %d, densitytmp[h] is nan. \n", h);
         exit(1);
       }
      }
      debug = true;
      if(debug) {
        coupler::printminmax1d(densitytmp, p3m3d.nphi*p3m3d.blockcount, p1->mype, "densitytoXGC", m, false);
      }
//      realsum=0.0;
      // send density from coupler to xgc
      MPI_Barrier(MPI_COMM_WORLD);
      coupler::send_from_coupler(adios,dir,densitytoXGC,cDens.IO,cDens.eng,cDens.name,
      senddensity,MPI_COMM_WORLD,m);

      // receive field from xgc to coupler
      fieldfromXGC = coupler::receive_field(dir, xFld,start_1, count_1, MPI_COMM_WORLD,m);
      debug = true;
      if (debug) {
        coupler::printminmax1d(fieldfromXGC->data(), count_1[0]*count_1[1], p1->mype,
	        "fieldfromXGC", m, false);
      }
      gxdp3d.DistriPotentRecvfromXGC(fieldfromXGC);

      MPI_Barrier(MPI_COMM_WORLD);
//    fprintf(stderr, "sxz 666 \n");
      fieldtoGEM = new coupler::Array3d<double>(p1->nx0, p1->lj0, p1->nz0, p1->tli0, p1->lj0,
                   p1->glk0, start3d);
      fieldgem = fieldtoGEM->data();
      MPI_Barrier(MPI_COMM_WORLD);
//    fprintf(stderr, "sxz 777 \n");
      for (int h=0; h<p1->tli0*p1->lj0*p1->glk0; h++) {      
        fieldgem[h] = gxdp3d.potGem[h];
      }      
      debug = false;
      if( debug && p1->mype == 1) {
	for(int h=0; h<p1->tli0*p1->lj0*p1->glk0; h++){
          fprintf(stderr, "h: %d, fieldgem: %19.18f \n", h, fieldgem[h]);
	  if(isnan(fieldgem[h])) {
	    printf("h: %d, fieldgem[h] is nan. \n", h);
	    exit(1);
	  }
	}
      } 
     
      debug = true;
      if (debug) {
        coupler::printminmax1d(fieldtoGEM->data(), p1->tli0*p1->lj0*p1->glk0, p1->mype,
                   "fieldtoGEM", m, false);
      }

      // send field from coupler to gem
  
      MPI_Barrier(MPI_COMM_WORLD);
//      fprintf(stderr, "sxz 444 \n");
      coupler::send3D_from_coupler(adios, dir, fieldtoGEM, cFld.IO, cFld.eng, cFld.name,
      sendfield, MPI_COMM_WORLD, m);
    }
  }

  coupler::destroy(densitytoXGC); 
  coupler::destroy(densityfromGEM); 
  coupler::destroy(fieldfromXGC);
  coupler::destroy(fieldtoGEM);         


  gDens.close();
  cDens.close();

  xFld.close();

//  cFld.close();

  gMesh.close();
//  gThf.close();
  gGrd.close();
  gQprof.close();
  xCouple.close();
//  xRzcoords.close();

  MPI_Finalize();
  printf("MPI is finalized.\n");
  return 0;
}
/*gemxgc.cc*/
