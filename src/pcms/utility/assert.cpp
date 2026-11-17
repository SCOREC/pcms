#include "pcms/utility/assert.hpp"
#include "pcms/utility/print.hpp"

#include <mpi.h>

#include <cstdlib>

namespace pcms
{
namespace detail
{
void Pcms_Assert_Fail(const char* msg)
{
  printError("%s", msg);
  abort();
}

void AssertFail(const char* cond, const char* file, int line, const char* msg)
{
  int rank = -1;
  int initialized = 0;
  int finalized = 0;
  MPI_Initialized(&initialized);
  MPI_Finalized(&finalized);
  if (initialized != 0 && finalized == 0) {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  }
  AssertFailOnRank(cond, file, line, rank, msg);
}

void AssertFailOnRank(const char* cond, const char* file, int line, int rank,
                      const char* msg)
{
  if (msg != nullptr) {
    printError("%s failed at %s + %d on rank %d: %s\n", cond, file, line, rank,
               msg);
  } else {
    printError("%s failed at %s + %d on rank %d\n", cond, file, line, rank);
  }
  abort();
}
} // namespace detail
} // namespace pcms
