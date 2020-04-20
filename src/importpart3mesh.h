#include<commpart1.h>
#include<importpart3mesh.cc>

namespace coupler {
class Part3Mesh3D{
  public:
    GO  nsurf;    // number of flux surfaces
    GO* versurf; // numbers of vertice on the flux surfaces
    double* xcoords;
    GO  li0,li1,li2;
    GO** xboxinds;  //The indexes of all boxes on the radial dimension
    GO lj0;
    GO* mylk0,mylk1,mylk2; // The indexes of box along z dimension
    GO boxstar,boxend,boxcount; // The  indexes of the 2d box

    double** Rcoords;  // The R coordinate of all vertices within the 2d box
    double** Zcoords;  // The Z coordinate of all vertices within the 2d box
    double** pzcoords;  // The z coordinates of all points with the 2d box.

}

void ImportPart3Mesh3D(Part3Mesh3D &p3m3d, Part1ParalPar3D  &p1pp3d);

void DistriPart3zcoords(Part3Mesh3D &p3m3d, Part1ParalPar3D  &p1pp3d);


GO  minloc(const double* zcoords, const GO n); 

void reshuffle_nodes(double* zcoords,const GO nstart,const GO vertnum);

void DistributePoints(double* exterarr,GO gstart,GO li, double* interarr,Part3Mesh3D &p3m3d,  \
     Part1ParalPar3D  &p1pp3d);


}
