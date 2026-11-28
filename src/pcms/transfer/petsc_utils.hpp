#ifndef PCMS_TRANSFER_PETSC_UTILS_HPP
#define PCMS_TRANSFER_PETSC_UTILS_HPP

#include <petscconf.h>
#include <petscmat.h>
#include <petscvec.h>
#include <Kokkos_Macros.hpp>

// GPU builds require PETSc with Kokkos support for device-resident
// matrix/vector types
#if (defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)) &&             \
  !defined(PETSC_HAVE_KOKKOS)
#error "GPU build (CUDA/HIP) requires PETSc compiled with Kokkos support \
(configure with --with-kokkos)"
#endif

namespace pcms::transfer
{
// Creates a sequential AIJ matrix: Kokkos AIJ on GPU builds, AIJ otherwise.
// Pass nz=0 and nnz=nullptr when using COO preallocation afterwards.
PetscErrorCode createSeqAIJMat(MPI_Comm comm, PetscInt m, PetscInt n,
                               PetscInt nz, const PetscInt nnz[], Mat* mat);

// Creates a sequential vector: a Kokkos vector on GPU builds, Seq otherwise.
PetscErrorCode createSeqVec(MPI_Comm comm, PetscInt n, Vec* vec);
} // namespace pcms::transfer

#endif // PCMS_TRANSFER_PETSC_UTILS_HPP
