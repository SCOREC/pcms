#ifndef DATAPROCESS_IMPL_H
#define DATAPROCESS_IMPL_H

#include "sendrecv_impl.h"
#include "commpart1.h"

namespace coupler {

template <class T>
void zDensityBoundaryBufAssign(LO nzb, LO lx, LO ly, LO lz, T*** lowbuf,
                               T*** upbuf, T*** box, Part1ParalPar3D& p1pp3d) {
  if (lowbuf == NULL || upbuf == NULL) {
    std::cout << "ERROR:the boundary buffer must be alloctted beforing "
                 "involing this routine.";
    std::exit(EXIT_FAILURE);
  }
  if (p1pp3d.npz > 1) {
    if (lz >= nzb) {
      mpisendrecv_aux2D(p1pp3d.comm_z, nzb, lx, ly, lz, lowbuf, upbuf, box);
    } else {
      std::cout << "ERROR: nzb is larger than lz. A larger lz is required.";
      std::exit(EXIT_FAILURE);
    }
  } else {
    if (p1pp3d.periods[2] == 1) {
      for (LO i = 0; i < lx ; i++) {
        for (LO j = 0; j < ly; j++) {
          for (LO k = 0; k < nzb; k++) {
            lowbuf[i][j][k] = box[i][j][lz - nzb + k];
            upbuf[i][j][k] = box[i][j][k];
          }
        }
      }
    } else {
      std::cout << "The topology is not right." << '\n';
      std::exit(EXIT_FAILURE);
    }
  }
}

} // namespace coupler

#endif
