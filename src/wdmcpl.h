#ifndef WDMCPL_H_
#define WDMCPL_H_
#include <mpi.h>
#include <redev.h>

namespace wdmcpl {
  void init(int argc, char** argv, MPI_Comm comm);
}

#endif
