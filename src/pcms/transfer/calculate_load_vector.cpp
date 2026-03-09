#include "pcms/transfer/calculate_load_vector.hpp"
#include "pcms/transfer/coo_assembly_utils.hpp"

namespace pcms {
PetscErrorCode calculateLoadVector(Omega_h::Mesh &target_mesh,
  Omega_h::Mesh &source_mesh,
  const IntersectionResults &intersection,
  const Omega_h::Reals &source_values,
  Vec*loadVec_out) {

  PetscFunctionBeginUser;
  PetscInt nnz = 0;
  PetscInt* coo_i = nullptr;
  PetscCall(build_linear_triangle_coo_rows(target_mesh, &coo_i, &nnz));
  PetscScalar* coo_vals = nullptr;
  PetscCall(PetscMalloc1(nnz, &coo_vals));

  // Fill COO values
  auto elmLoadVector =
      buildLoadVector(target_mesh, source_mesh, intersection, source_values);

  auto hostElmLoadVector = Kokkos::create_mirror_view(elmLoadVector);
  Kokkos::deep_copy(hostElmLoadVector, elmLoadVector);
  PetscCheck(static_cast<PetscInt>(hostElmLoadVector.extent(0)) == nnz,
             PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ,
             "Element load vector size (%d) does not match COO nnz (%d)",
             static_cast<PetscInt>(hostElmLoadVector.extent(0)), nnz);

  for (PetscInt e = 0; e < nnz; ++e) {
    coo_vals[e] = hostElmLoadVector(e);
  }

  // create vector with preallocated COO structure
  Vec vec;
  PetscCall(VecCreate(PETSC_COMM_WORLD, &vec));
  PetscCall(VecSetSizes(vec, target_mesh.nverts(), PETSC_DECIDE));
  PetscCall(VecSetFromOptions(vec));
  PetscCall(VecSetPreallocationCOO(vec, nnz, coo_i));
  PetscCall(VecSetValuesCOO(vec, coo_vals, ADD_VALUES));
  PetscCall(PetscFree(coo_i));
  PetscCall(PetscFree(coo_vals));

  if (target_mesh.nelems() < 10) {
    PetscCall(VecView(vec, PETSC_VIEWER_STDOUT_WORLD));
  }

  *loadVec_out = vec;
  PetscFunctionReturn(PETSC_SUCCESS);
}
}
