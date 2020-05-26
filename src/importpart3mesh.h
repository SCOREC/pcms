#ifndef IMPORTPART3MESH_H
#define IMPORTPART3MESH_H

#include "couplingTypes.h"
#include "testutilities.h"

namespace coupler {

//forward declare the Part1ParalPar3D class; we don't care how 
// the class is defined in this header since we don't call
// any of its functions or access its member variables
class Part1ParalPar3D;

class Part3Mesh3D{
  public:
    LO nsurf;    // number of flux surfaces
    LO* versurf = NULL; // numbers of vertice on the flux surfaces //
    double* xcoords = NULL;
    double* zcoordall = NULL;
    LO  li0,li1,li2;
    LO** xboxinds=NULL;  //The indexes of all boxes on the radial dimension
    LO lj0;
    LO* mylk0=NULL;
    LO* mylk1=NULL;
    LO* mylk2=NULL; // The indexes of box along z dimension
    GO boxstar,boxend,boxcount; // The  indexes of the 2d box

    double** Rcoords=NULL;  // The R coordinate of all vertices within the 2d box
    double** Zcoords=NULL;  // The Z coordinate of all vertices within the 2d box
    double** pzcoords=NULL;  // The z coordinates of all points within the 2d box.
    /* constructor - versurf has length = numsurf & versurf[i] = the number of nodes surface[i]
     * xcoords saves the radial coordinate of each surface.
     * zcoords saves the poloidal angle of each node on each surface.
     */
    Part3Mesh3D(Part1ParalPar3D &p1pp3d,
        LO nsurf_,
        LO* versurf_,
        double* xcoords_,
        double* zcoord_,
        bool pproc = true)
      : nsurf(nsurf_),
        versurf(versurf_),
        xcoords(xcoords_),
        zcoordall(zcoord_),
	preproc(pproc), test_case(TestCase::off) {
      init(p1pp3d);
    }
    /* testing constructor
     * arguments support reading
     * the test case number and directory
     * from the user
     */
    Part3Mesh3D(Part1ParalPar3D &p1pp3d,
        bool pproc,
        TestCase tcase,
        std::string tdir)
      : preproc(pproc), test_case(tcase) {
      assert(tcase < TestCase::invalid);
      init(p1pp3d,tdir);
    }
    ~Part3Mesh3D()
    {
     if(versurf!=NULL) delete[] versurf;
     if(xboxinds!=NULL) delete[] xboxinds;
     if(xcoords!=NULL) delete[] xcoords;
     if(mylk0!=NULL) delete[] mylk0;
     if(mylk1!=NULL) delete[] mylk1;
     if(mylk2!=NULL) delete[] mylk2;
     if(Rcoords!=NULL) delete[] Rcoords;
     if(Zcoords!=NULL) delete[] Zcoords;
     if(pzcoords!=NULL) delete[] pzcoords;
   }
  private:
    const bool preproc;
    const TestCase test_case;
    void init(const Part1ParalPar3D &p1pp3d, const std::string test_dir="");
    /* helper function called by init */
    void DistriPart3zcoords(const Part1ParalPar3D  &p1pp3d,
        const std::string test_dir="");
    /* helper function called by DistriPart3zcoords */
    void DistributePoints(double* exterarr, LO gstart, LO li,
        double* interarr, const Part1ParalPar3D  &p1pp3d);
    /* default constructor 
     * put this in private to prevent users from calling it
     */
    Part3Mesh3D() : preproc(true), test_case(TestCase::off) {};
};


// The following utility functions should go into a header for utilities if they
// need to be called by other other classes or code.
// If they are not needed elsewhere, the prototypes/definitions below should
// be removed.
LO  minloc(const double* zcoords, const LO n); 

void reshuffle_nodes(double* zcoords,const LO nstart,const LO vertnum);

double minimalvalue(const double* array, const LO n);

}

#endif
