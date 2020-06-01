#ifndef BOUNDARY_DESCR_3D_H
#define BOUNDARY_DESCR_3D_H

#include "couplingTypes.h"
#include "testutilities.h"

namespace coupler {

//forward declare
class DatasProc3D;
class Part3Mesh3D;
class Part1ParalPar3D;

//boundary treatment
class BoundaryDescr3D{
  public:
    LO nzb;
    double** upzpart3=NULL;
    double** lowzpart3=NULL;
    CV*** updenz=NULL; // The upper  boundary buffer on z domain for interpolation and storing the real quantiies resulted from the backward Fourier transform of complex charged density.
    CV*** lowdenz=NULL;
    CV*** uppotentz=NULL; //The upper  boundary buffer on z domain for interpolation and storing the complex  quantiies resulted from the forward Fourier transform of electrosttic potential.
    CV*** lowpotentz=NULL;
    CV** uppbmat=NULL; //the upper matrix for the parallel boundary condition
    CV** lowpbmat=NULL; // the lower matrix for the parallel boundary condition
   /* constructor */
    BoundaryDescr3D(const Part3Mesh3D& p3m3d,
        const Part1ParalPar3D &p1pp3d,
        const DatasProc3D& dp3d,
        TestCase tcase = TestCase::off,
        bool pproc = true);
    /* destructor */
    ~BoundaryDescr3D();
    void zPotentBoundaryBufAssign(const DatasProc3D& dp3d, 
        const Part3Mesh3D& p3m3d,
        const Part1ParalPar3D &p1pp3d);
    void zDensityBoundaryBufAssign(CV*** box, const Part1ParalPar3D& p1pp3d);

  private:
    /* prevent users from calling this */  
    BoundaryDescr3D() : test_case(TestCase::off), preproc(false) {};
    void initpbmat(const Part1ParalPar3D &p1pp3d);
};
    const TestCase test_case;
    const bool preproc;

}

#endif
