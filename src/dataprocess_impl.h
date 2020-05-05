#ifndef DATAPROCESS_IMPL_H
#define DATAPROCESS_IMPL_H

#include "sendrecv_impl.h"
#include "commpart1.h"

namespace coupler {

// FIXME Not used; was not compiled; based on the
// FIXME arguments this should go into DatasProc3D.
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
      // SZ The call to mpisendrecv_aux2D was missing the
      // SZ mpi communicator argument, nzb, and was passing
      // SZ p1pp3d which did not match the function definition.
      MPI_Comm FIXME = MPI_COMM_WORLD; // HACK to get this to compile
      mpisendrecv_aux2D(FIXME, nzb, lx, ly, lz, lowbuf, upbuf, box);
    } else {
      std::cout << "ERROR: nzb is larger than lzb. Large lzb is required.";
      std::exit(EXIT_FAILURE);
    }
  } else {
    if (p1pp3d.periods[2] == 1) {
      for (LO i = 0; i < lx - 1; i++) {
        for (LO j = 0; j < ly - 1; j++) {
          for (LO k = 0; k < nzb - 1; k++) {
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
