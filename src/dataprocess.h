#ifndef DATAPROCESS_H
#define DATAPROCESS_H

#include "couplingTypes.h"
#include <fftw3.h>

namespace coupler {

// forward declare
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
  CV*** densin = NULL; // input 3d density in complex number
  CV* densintmp = NULL; // temporary 2d density array prepared for backward
                        // fourier transform
  double* densouttmp = NULL; // store the x-y 2d real density after backward 
                             // fourier transform
  double*** densout = NULL; // store xyz 3d real density
  double*** denspart3 = NULL; // storing the density being sent to the part3
  double*** potentin = NULL; // the input real electrostatic potential in 3d xyz
  double* potentintmp = NULL; // temporary xy 2d potential array for forward 
                              // fourier transform
  CV* potentouttmp = NULL;
  CV*** potentout = NULL; // 3d temporary array stroring complex
                          // electrostatic potential
  CV*** potentpart1 = NULL; // storing the electrostatic potential being sent
                            // to the part1.
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
      bool ypar = false);
  ~DatasProc3D();
  //routines for Fourier transform
  void CmplxdataToRealdata3D();
  LO getP1li0() { return p1.li0; };
  LO getP1ny0() { return p1.ny0; };
  LO getP1npy() { return p1.npy; };

private:
  const bool preproc;
  const bool yparal;

  const struct P1Data {
    P1Data(LO li, LO lj, LO lk, LO ny, LO np, LO pe_y, LO res) : 
	    li0(li), lj0(lj), lk0(lk), 
	  ny0(ny), npy(np), mype_y(pe_y), res_fact(res)
	  {};
    const LO li0;
    const LO lj0;
    const LO lk0;
    const LO ny0;
    const LO npy;
    const LO mype_y;
    const LO res_fact;
  } p1;

  const struct P3Data {
    P3Data(LO li, LO lj, LO* mylk) : li0(li), lj0(lj), mylk0(mylk) {};
    const LO li0;
    const LO lj0;
    const LO* mylk0;
  } p3;

  /* helper function for destructor */
  void FreeFourierPlan3D(); // called from the destructor - does that make sense?
  /* helper functions for constructor */
  void init();
  void AllocDensityArrays();
  void AllocPotentArrays();
  /* helper functions for CmplxdataToRealdata3D and RealdataToCmplxdata3D */
  void ExecuteRealToCmplx(P1Data &p1);
  void ExecuteCmplxToReal(P1Data &p1);

public:
  
  void RealdataToCmplxdata3D(P1Data &p1, P3Data &p3);
  void InitFourierPlan3D(P1Data &p1,P3Data &p3);
};

void TransposeComplex(CV** InMatrix,CV** OutMatrix, DatasProc3D& dp3d,
     LO* p1);

void InterpoDensity3D(const BoundaryDescr3D& bdesc, const DatasProc3D& dp3d);

void InterpoPotential3D(const BoundaryDescr3D& bdesc, const DatasProc3D& dp3d);

} // namespace coupler

#endif
