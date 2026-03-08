#include "pcms/transfer/calculate_load_vector.hpp"

namespace pcms {
PetscErrorCode calculateLoadVector(Omega_h::Mesh &target_mesh,
  Omega_h::Mesh &source_mesh,
  const IntersectionResults &intersection,
  const Omega_h::Reals &source_values,
  Vec*loadVec_out) {

  PetscFunctionBeginUser;
  const int numNodesPerTri = 3;

  const int nnz = target_mesh.nelems() * numNodesPerTri;

  // Allocate COO indices and values
  PetscInt* coo_i;
  PetscScalar* coo_vals;
  PetscCall(PetscMalloc2(nnz, &coo_i, nnz, &coo_vals));

  // Fill COO global indices and values
  auto elmVerts = Omega_h::HostRead(target_mesh.ask_elem_verts());
  auto elmLoadVector =
      buildLoadVector(target_mesh, source_mesh, intersection, source_values);

  auto hostElmLoadVector = Kokkos::create_mirror_view(elmLoadVector);
  Kokkos::deep_copy(hostElmLoadVector, elmLoadVector);

  PetscInt idx = 0;
  for (PetscInt e = 0; e < target_mesh.nelems(); ++e) {
    for (PetscInt vi = 0; vi < numNodesPerTri; ++vi) {
      coo_i[idx] = elmVerts[numNodesPerTri * e + vi];
      coo_vals[idx] = hostElmLoadVector(numNodesPerTri * e + vi);
      ++idx;
    }
  }

  // create vector with preallocated COO structure
  Vec vec;
  PetscCall(VecCreate(PETSC_COMM_WORLD, &vec));
  PetscCall(VecSetSizes(vec, target_mesh.nverts(), PETSC_DECIDE));
  PetscCall(VecSetFromOptions(vec));
  PetscCall(VecSetPreallocationCOO(vec, nnz, coo_i));
  PetscCall(VecSetValuesCOO(vec, coo_vals, ADD_VALUES));
  PetscCall(PetscFree2(coo_i, coo_vals));

  if (target_mesh.nelems() < 10) {
    PetscCall(VecView(vec, PETSC_VIEWER_STDOUT_WORLD));
  }

  *loadVec_out = vec;
  PetscFunctionReturn(PETSC_SUCCESS);
}
}
