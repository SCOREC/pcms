/**
 * @file calculate_load_vector.hpp
 * @brief Routines for assembling global load vector in conservative field
 * projection.
 *
 * Provides functionality to compute and assemble the global load vector
 * used in Galerkin-based conservative field transfer between non-conforming
 * meshes.
 *
 * @created by Abhiyan Paudel
 * @date August 2025
 */

#ifndef PCMS_INTERPOLATOR_CALCULATE_LOAD_VECTOR_HPP
#define PCMS_INTERPOLATOR_CALCULATE_LOAD_VECTOR_HPP
#include <Omega_h_adapt.hpp>
#include <Omega_h_array_ops.hpp>
#include <Omega_h_atomics.hpp> //Omega_h::atomic_fetch_add
#include <Omega_h_build.hpp>
#include <Omega_h_class.hpp>
#include <Omega_h_compare.hpp>
#include <Omega_h_dbg.hpp>
#include <Omega_h_file.hpp> //Omega_h::binary
#include <Omega_h_for.hpp>
#include <Omega_h_recover.hpp> //project_by_fit
#include <Omega_h_shape.hpp>
#include <Omega_h_timer.hpp>
#include <iomanip> //precision
#include <iostream>
#include <petscvec_kokkos.hpp>
#include <sstream> //ostringstream

#include <pcms/interpolator/mesh_interscetion/load_vector_integrator.hpp>
#include <MeshField.hpp>

#include <petscmat.h>

// detect floating point exceptions
#include <fenv.h>

using ExecutionSpace = Kokkos::DefaultExecutionSpace;
using MemorySpace = Kokkos::DefaultExecutionSpace::memory_space;

/**
 * @brief Assembles the global load vector.
 *
 * This function computes the unassembled local load vector contributions for
 * each triangular element in the target mesh using `buildLoadVector()` and then
 * assembles them into a global PETSc vector in COO format.
 *
 *
 * @param target_mesh The target Omega_h mesh to which the scalar field is being
 * projected.
 * @param source_mesh The source Omega_h mesh containing the original scalar
 * field values.
 * @param intersection Precomputed intersection data for each target element.
 *                     Includes the number and indices of intersecting source
 * elements.
 * @param source_values Nodal scalar field values defined on the source mesh.
 * @param[out] loadVec_out Pointer to a PETSc Vec where the assembled load
 * vector will be stored.
 *
 * @return PetscErrorCode Returns PETSC_SUCCESS if successful, or an appropriate
 * PETSc error code otherwise.
 *
 * @note
 * - Works for 2D linear triangular elements.
 * - Uses COO-style preallocation and insertion into the PETSc vector.
 * - Internally calls `buildLoadVector()` to compute per-element contributions.
 * - The resulting vector is used as the right-hand side (RHS) in a projection
 * solve.
 *
 * @see buildLoadVector,IntersectionResults
 */

namespace pcms
{

inline PetscErrorCode calculateLoadVector(
  Omega_h::Mesh& target_mesh, Omega_h::Mesh& source_mesh,
  const IntersectionResults& intersection, const Omega_h::Reals& source_values,
  Vec* loadVec_out)
{

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
} // namespace pcms
#endif
