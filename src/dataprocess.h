#ifndef DATAPROCESS_H
#define DATAPROCESS_H

#include "couplingTypes.h"
#include "testutilities.h"
#include <fftw3.h>
#include "importpart3mesh.h"
#include "commpart1.h"
#include "BoundaryDescr3D.h"

namespace coupler {

// forward declare
template<class T>
class Array2d;

class Part3Mesh3D;
class Part1ParalPar3D;
class BoundaryDescr3D;

class DatasProc3D {
public:
  LO part1li0; // part1li0 is the element number of subdomain of x on each
               // process belonging to comm_y. This x subdomain belongs to 2d
               // box on the x-y plane storing the input charged ion density for
               // fourier transform.
  LO part3li0; // part3li0 is the element number of subdomain of x on each
               // process belonging to comm_y. This x subdomain belongs to 2d 
               // box on the x-y plane storing the input electrostatic 
               // potential.
  LO part1lj0; // the count of elements on y domain on each process after
               // backward Fourier transform
  LO part3lj0; ////the count of elements on y domain on each process after
               ///forward Fourier transform  
  LO sum;
  // here, pointers must be assigned a NULL;
  CV** densrecv = NULL; // the 2d array density recived by the coupler.
                           // and sent by part1. This is received in comm_x and comm_y communicators.
  CV*** densin = NULL; // input 3d density in complex number
  CV*** densinterpo = NULL; 
  CV* densintmp = NULL; // temporary 2d density array prepared for backward
                        // fourier transform
  double* densouttmp = NULL; // store the x-y 2d real density after backward 
                             // fourier transform
  double*** denspart3 = NULL; // storing the density being sent to the part3
  double*** densTOpart3 = NULL;

  double*  denssend = NULL; // the 1d array density  sent to part3
                             //  would be sent in comm_x and comm_y communicators


  double** potentrecv = NULL; // the 2d array potential received by the coupler 
                                 // and sent by part3. This is received in comm_x and comm_y communicators.
  double*** potentin = NULL; // the input real electrostatic potential in 3d xyz
  double* potentintmp = NULL;
  CV* potentouttmp = NULL; // 
  CV*** potentinterpo = NULL; // temporary xy 2d potential array for forward 
                              // fourier transform
  CV*** potentpart1 = NULL; // storing the electrostatic potential being sent
                            // to the part1.
  CV*  potentsend = NULL; // the 1d array complex potential sent  to part1, and would be sent 
                           // in comm_x and comm_y communicators. 

// matrix for the transformation between planes and xyz
   double**** mattoplane=NULL;
 
   CV****      mat_to_plane=NULL; 
   double**** mat_from_weight=NULL;
   int****    mat_from_ind_plane=NULL;
   int****    mat_from_ind_n=NULL; 


  fftw_plan plan_forward = NULL, plan_backward = NULL;
  // The following parameters for yparal=true;
  LO myli0;


  /* constructor
   * optional argument supports setting
   * the prepoc and yparal modes
   */
  DatasProc3D(const Part1ParalPar3D& p1pp3d,
      const Part3Mesh3D &p3m3d,
      bool pproc = true,
      TestCase test_case = TestCase::off,
      bool ypar = false,
      int nummode = 1);
  ~DatasProc3D();
  void InterpoDensity3D(const BoundaryDescr3D& bdesc, const Part3Mesh3D& p3m3d,
                        const Part1ParalPar3D& p1pp3d);
  void InterpoPotential3D(const BoundaryDescr3D& bdesc, const Part3Mesh3D& p3m3d,
                        const Part1ParalPar3D& p1pp3d);
  //routines for Fourier transform
  void CmplxdataToRealdata3D();
  void RealdataToCmplxdata3D();
  void InitFourierPlan3D();
  void AllocMatXYZtoPlane();
  void Prepare_mats_from_planes();
  void Initmattoplane();

  void oldAssemDensiSendtoPart3(BoundaryDescr3D& bdesc,const Part3Mesh3D &p3m3d, const Part1ParalPar3D& p1pp3d);
  void AssemDensiSendtoPart3(BoundaryDescr3D& bdesc,const Part3Mesh3D &p3m3d, const Part1ParalPar3D& p1pp3d); 
  void DistriDensiRecvfromPart1(const Part3Mesh3D &p3m3d, const Part1ParalPar3D& p1pp3d,
       const Array2d<CV>* densityfromGENE);
  void AssemPotentSendtoPart1(const Part3Mesh3D &p3m3d, const Part1ParalPar3D& p1pp3d);
  void DensityToPart3(const Part3Mesh3D& p3m3d);
  void oldDistriPotentRecvfromPart3(const Part3Mesh3D& p3m3d, const Part1ParalPar3D& p1pp3d,
     const Array2d<double>* fieldfromXGC);
  void DistriPotentRecvfromPart3(const Part3Mesh3D &p3m3d, const Part1ParalPar3D& p1pp3d,
       const Array2d<double>* fieldfromXGC); 
  void oldInitmattoplane(const Part3Mesh3D& p3m3d,const Part1ParalPar3D& p1pp3d);  

  LO getP1li0() { return p1.li0; };
  LO getP1ny0() { return p1.ny0; };
  LO getP1npy() { return p1.npy; };


private:
  const bool preproc;
  const TestCase testcase;
  const bool yparal;
  Part1ParalPar3D p1;
  Part3Mesh3D p3;
/*
  // this struct contains the read-only values from Part1ParalPar3D class
  const struct P1Data {
    P1Data(LO li, LO lj, LO lk, LO ny, LO np, LO pe_y,GO blockcount_, LO res_, LO n0_global_,const double* q_prof_,
           double L_tor_) : li0(li), lj0(lj), lk0(lk), ny0(ny), npy(np), mype_y(pe_y), blockcount(blockcount_),
           res_fact(res_),n0_global(n0_global_),q_prof(q_prof_),L_tor(L_tor_){};
    const LO li0;
    const LO lj0;
    const LO lk0;
    const LO ny0;
    const LO npy;
    const LO mype_y;
    const LO res_fact;
    const GO blockcount;
    const LO n0_global;
    const double* q_prof;
    const double L_tor; 
  } p1;

  // this struct contains the read-only values from Part3Mesh3D class
  const struct P3Data {
    P3Data(LO li, LO lj,GO blockcount_,LO* mylk,LO y_res_b,const double** pzcoords_) : li0(li), lj0(lj), 
          blockcount(blockcount_), mylk0(mylk),y_res_back(y_res_b),pzcoord(pzcoords_) {};
      const LO li0;
      const LO lj0;
      const LO y_res_back;
      const GO blockcount;
      LO const* const mylk0;
      const double** pzcoords;
  } p3;
*/
  /* helper function for destructor */
  void FreeFourierPlan3D(); // called from the destructor - does that make sense?
  /* helper functions for constructor */
  void init();
  void AllocDensityArrays();
  void AllocPotentArrays();
  void TestInitPotentAlongz(const Part3Mesh3D& p3m3d, const LO npy, const LO n);
  /* helper functions for CmplxdataToRealdata3D and RealdataToCmplxdata3D */
  void ExecuteRealToCmplx();
  void ExecuteCmplxToReal();

  };
 
void TransposeComplex(CV** InMatrix,CV** OutMatrix, DatasProc3D& dp3d,
     Part1ParalPar3D& p1pp3d);


} // namespace coupler

#endif
