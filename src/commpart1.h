#ifndef COMMPART1_H
#define COMMPART1_H	

#include "classes.h"
#include "coupling.h"
namespace coupler {

void InitPart1ParalPar3D(Part1ParalPar3D& p1pp3d);

void InitPart1paral3DInCoupler(Part1ParalPar3D  &p1pp3d);

void CreateSubCommunicators(Part1ParalPar3D  &p1pp3d);

void MpiFreeComm(Part1ParalPar3D  &p1pp3d);

}



#endif







