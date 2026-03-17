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

#ifndef PCMS_TRANSFER_CALCULATE_LOAD_VECTOR_HPP
#define PCMS_TRANSFER_CALCULATE_LOAD_VECTOR_HPP
#include <Omega_h_mesh.hpp>
#include <pcms/transfer/mesh_intersection.hpp>
#include <petscvec.h>

namespace pcms
{

/**
 * @brief Assembles the global load vector using mesh intersection method.
 *
 * This function computes the unassembled local load vector contributions for
 * each triangular element in the target mesh using `buildLoadVectorMI()` and
 * then assembles them into a global PETSc vector in COO format.
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
 * - Internally calls `buildLoadVectorMI()` to compute per-element
 * contributions.
 * - The resulting vector is used as the right-hand side (RHS) in a projection
 * solve.
 *
 * @see buildLoadVectorMI,IntersectionResults
 */
// FIXME use PCMS error handling rather than returning a PETSC error code
PetscErrorCode calculateLoadVectorMI(Omega_h::Mesh& target_mesh,
                                     Omega_h::Mesh& source_mesh,
                                     const IntersectionResults& intersection,
                                     const Omega_h::Reals& source_values,
                                     Vec* loadVec_out);

/**
 * @brief Assembles the global load vector using Monte Carlo integration.
 *
 * This function computes the unassembled local load vector contributions for
 * each triangular element in the target mesh using Monte Carlo integration and
 * then assembles them into a global PETSc vector in COO format.
 *
 * The element-local contributions are computed from source-field values already
 * evaluated at the sampled physical points in each target element. The sampling
 * pattern is defined on the reference triangle and may be generated either by
 * uniform random sampling or by precomputed Sobol barycentric samples.
 *
 * @param target_mesh The target Omega_h mesh to which the scalar field is being
 * projected.
 * @param field_values_at_points Source-field values evaluated at the sampled
 * physical points in the target mesh. The data is expected to be stored
 * element-by-element, with total size
 * `target_mesh.nelems() * npoints_each_tri`.
 * @param npoints_each_tri Number of sample points used in each target triangle
 * for the Monte Carlo integration.
 * @param method Sampling method used to define the reference barycentric sample
 * coordinates. Supported options include uniform random sampling and Sobol
 * sampling.
 * @param sobol_filename Path to the file containing precomputed Sobol
 * barycentric samples when `method` is `SamplingMethod::SOBOL`.
 * @param[out] loadVec_out Pointer to a PETSc Vec where the assembled load
 * vector will be stored.
 *
 * @return PetscErrorCode Returns PETSC_SUCCESS if successful, or an appropriate
 * PETSc error code otherwise.
 *
 * @note
 * - Works for 2D linear triangular elements.
 * - Uses COO-style preallocation and insertion into the PETSc vector.
 * - Internally computes per-element Monte Carlo load vector contributions
 *   before assembling the global vector.
 * - For linear triangular elements, the barycentric coordinates are equal to
 *   the local shape-function values used in the Monte Carlo estimator.
 * - The resulting vector is used as the right-hand side (RHS) in a projection
 *   solve.
 *
 * @see buildLoadVectorMC, SamplingMethod
 */
PetscErrorCode calculateLoadVectorMC(
  Omega_h::Mesh& target_mesh, const Omega_h::Reals& field_values_at_points,
  const int npoints_each_tri, SamplingMethod method,
  const std::string sobol_filename = "", Vec* loadVec_out);
} // namespace pcms
#endif // PCMS_TRANSFER_CALCULATE_LOAD_VECTOR_HPP
