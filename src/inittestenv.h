#ifndef INITTESTENV_H
#define INITTESTENV_H

#include "commpart1.h"
#include "importpart3mesh.h"
#include "dataprocess.h"
#include "BoundaryDescr3D.h"

namespace coupler { 

void TestInitPotentAlongz(DatasProc3D& dp3d, const Part3Mesh3D& p3m3d, const LO npy, const LO n);


}

#endif
