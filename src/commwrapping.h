#include <coupling.h>

namespace coupler {

template<class T1, class T2>
void mpiallgatherwrapper(T1* in,T1* out,T2 n,MPI_Datatype mpitype,MPI_Comm comm)
{
  assert(n>0);
  MPI_Allgather(in,n,mpitype,out,n,mpitype,comm);  
}






}
