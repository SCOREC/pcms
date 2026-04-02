// petsc requires these headers be included before petscvec
#if defined(PETSC_HAVE_KOKKOS)
#include <petscmat_kokkos.hpp>
#include <petscvec_kokkos.hpp>
#endif

#include "pcms/transfer/petsc_utils.hpp"

namespace pcms
{
PetscErrorCode createSeqAIJMat(MPI_Comm comm, PetscInt m, PetscInt n,
                               PetscInt nz, const PetscInt nnz[], Mat* mat)
{
  PetscFunctionBeginUser;
#if defined(PETSC_HAVE_KOKKOS)
  PetscCall(MatCreateSeqAIJKokkos(comm, m, n, nz, nnz, mat));
#else
  PetscCall(MatCreateSeqAIJ(comm, m, n, nz, nnz, mat));
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode createSeqVec(MPI_Comm comm, PetscInt n, Vec* vec)
{
  PetscFunctionBeginUser;
#if defined(PETSC_HAVE_KOKKOS)
  PetscCall(VecCreateSeqKokkos(comm, n, vec));
#else
  PetscCall(VecCreateSeq(comm, n, vec));
#endif
  PetscFunctionReturn(PETSC_SUCCESS);
}
} // namespace pcms
