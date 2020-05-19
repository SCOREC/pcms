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
    double*** updenz=NULL; // The upper  boundary buffer on z domain for interpolation and storing the real quantiies resulted from the backward Fourier transform of complex charged density.
    double*** lowdenz=NULL;
    double*** uppotentz=NULL; //The upper  boundary buffer on z domain for interpolation and storing the complex  quantiies resulted from the forward Fourier transform of electrosttic potential.
    double*** lowpotentz=NULL;
   /* constructor */
    BoundaryDescr3D(const Part3Mesh3D& p3m3d,
        const Part1ParalPar3D &p1pp3d,
        const DatasProc3D& dp3d,
        const TestCase tcase = TestCase::off,
        const bool pproc = true);
    /* destructor */
    ~BoundaryDescr3D();
    //Not used
    void zPotentBoundaryBufAssign(const DatasProc3D& dp3d, 
        const Part3Mesh3D& p3m3d,
        const Part1ParalPar3D &p1pp3d);
  private:
    /* prevent users from calling this */  
    BoundaryDescr3D() : test_case(TestCase::off), preproc(false) {};
    const TestCase test_case;
    const bool preproc;
};

}

#endif
