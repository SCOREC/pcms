#ifndef IMPORTPART3MESH_H
#define IMPORTPART3MESH_H

#include "classes.h"
#include "coupling.h"

namespace coupler {

void ImportPart3Mesh3D(Part3Mesh3D &p3m3d, Part1ParalPar3D  &p1pp3d);

void DistriPart3zcoords(Part3Mesh3D &p3m3d, Part1ParalPar3D  &p1pp3d);


LO  minloc(const double* zcoords, const LO n); 

void reshuffle_nodes(double* zcoords,const LO nstart,const LO vertnum);

void DistributePoints(double* exterarr,LO gstart,LO li, double* interarr,Part3Mesh3D &p3m3d,
     Part1ParalPar3D  &p1pp3d);

double minimalvalue(const double* array, const LO n);

void InitzcoordsInCoupler(double* zcoords,LO* versurf,LO nsurf);

}

#endif
