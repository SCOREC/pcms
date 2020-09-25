#ifndef BOUNDARY_DESCR_3D_H
#define BOUNDARY_DESCR_3D_H

#include "couplingTypes.h"
#include "testutilities.h"

namespace coupler {

//forward declare
class Part3Mesh3D;
class Part1ParalPar3D;

//boundary treatment
class BoundaryDescr3D{
  public:
    LO nzb;
    double** upzpart3=NULL;
    double** lowzpart3=NULL;
    CV*** updenz=NULL; // The upper  boundary buffer on z domain for interpolation and storing the 
                       //real quantiies resulted from the backward Fourier transform of complex charged density.
    CV*** lowdenz=NULL;
    CV*** uppotentz=NULL; //The upper  boundary buffer on z domain for interpolation and storing the complex  
                          //quantiies resulted from the forward Fourier transform of electrosttic potential.
    CV*** lowpotentz=NULL;
    CV** uppbmat=NULL; //the upper matrix for the parallel boundary condition
    CV** lowpbmat=NULL; // the lower matrix for the parallel boundary condition

    // specially for GEM-XGC coupling
    double*** lowdenzgemxgc=NULL;
    double*** updenzgemxgc=NULL;
    double*** uppotentzgemxgc=NULL;
    double*** lowpotentzgemxgc=NULL;

    double*** ymeshxgc=NULL;
    double* ymeshgem=NULL;
    double** thflxmeshxgc=NULL;
    double* thetameshgem=NULL;   


   /* constructor */
    BoundaryDescr3D(const Part3Mesh3D& p3m3d,
        const Part1ParalPar3D &p1pp3d,
        const CouplingCase c_case,
        const TestCase tcase = TestCase::off,
        const bool pproc = true)
      : ccase(c_case),test_case(tcase),preproc(pproc){
        if(c_case==CouplingCase::genexgc){
          initGeneXgc(p3m3d,p1pp3d);
        }else if(c_case==CouplingCase::gemxgc){
          initGemXgc(p3m3d,p1pp3d);
        }       
      }
    /* destructor */
    ~BoundaryDescr3D();
  private:
    /* prevent users from calling this */  
//    BoundaryDescr3D() : test_case(TestCase::off), preproc(false) {};
    const TestCase test_case;
    const CouplingCase ccase;
    const bool preproc;
    const Part1ParalPar3D* p1;
    const Part3Mesh3D* p3;
    void initpbmat(const Part1ParalPar3D &p1pp3d);
    void initGeneXgc(const Part3Mesh3D& p3m3d,const Part1ParalPar3D &p1pp3d);
    void initGemXgc(const Part3Mesh3D& p3m3d,const Part1ParalPar3D &p1pp3d);
};

}

#endif
