/**
 * @file calculate_load_vector.hpp
 * @brief Routines for assembling global load vector in conservative field
 * projection.
 *
 * Provides functionality to compute and assemble the global load vector
 * used in Galerkin-based conservative field transfer between non-conforming
 * meshes.
 *
 */

#ifndef PCMS_INTERPOLATOR_CALCULATE_LOAD_VECTOR_HPP
#define PCMS_INTERPOLATOR_CALCULATE_LOAD_VECTOR_HPP
#include <Omega_h_mesh.hpp>
#include <pcms/transfer/load_vector_integrator.hpp>
#include <petscmat.h>

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
// FIXME use PCMS error handling rather than returning a PETSC error code
PetscErrorCode calculateLoadVector(
  Omega_h::Mesh& target_mesh, Omega_h::Mesh& source_mesh,
  const IntersectionResults& intersection, const Omega_h::Reals& source_values,
  Vec* loadVec_out);
} // namespace pcms
#endif
