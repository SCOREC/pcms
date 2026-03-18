/**
 * @file conservative_projection_solver.hpp
 * @brief Solves the conservative projection of scalar fields between
 * non-matching meshes.
 *
 * Provides the main interface to perform Galerkin projection of scalar fields
 * from a source mesh to a target mesh using conservative transfer using a
 * supermesh generated from mesh intersections.
 *
 * The solver computes the right-hand side (load vector), assembles the mass
 * matrix, and solves the resulting linear system to obtain projected nodal
 * values.
 *
 */

#ifndef PCMS_TRANSFER_CONSERVATIVE_PROJECTION_SOLVER_HPP
#define PCMS_TRANSFER_CONSERVATIVE_PROJECTION_SOLVER_HPP

#include <Omega_h_array.hpp>
#include <Omega_h_mesh.hpp>

#include <pcms/transfer/mesh_intersection.hpp>

namespace pcms
{

/**
 * @brief Solves a conservative Galerkin projection problem to transfer a scalar
 * field onto a target mesh.
 *
 * This routine assembles and solves a linear system of the form
 * \f[
 *   M \cdot x = b
 * \f]
 * where
 * - \f$M\f$ is the mass matrix on the target mesh based on linear (\f$P_1\f$)
 *   finite elements,
 * - \f$b\f$ is the load vector associated with the projection,
 * - \f$x\f$ is the unknown nodal field on the target mesh.
 *
 * The load vector is computed by integrating the source field against the
 * target basis functions. Depending on the solver variant, these integrals are
 * evaluated either exactly using mesh-intersection method or approximately
 * using Monte Carlo sampling.
 *
 * ### Algorithm steps
 * 1. Assemble the target mass matrix.
 * 2. Assemble the load vector corresponding to the chosen projection method.
 * 3. Solve the resulting linear system using PETSc.
 * 4. Return the projected nodal field on the target mesh.
 *
 * @return A vector of nodal values on the target mesh after projection
 *         (`Omega_h::Reals`).
 *
 * @note The specific inputs required to assemble the load vector depend on the
 *       projection method used by the solver variant.
 */

/**
 * @brief Solves the conservative Galerkin projection using exact integration
 * over mesh-intersection regions.
 *
 * @param target_mesh The target Omega_h mesh where the field is projected.
 * @param source_mesh The source Omega_h mesh containing the original field.
 * @param intersection Precomputed intersection information between source and
 *                     target meshes.
 * @param source_values Nodal scalar field values on the source mesh.
 * @return Projected nodal values on the target mesh.
 */
Omega_h::Reals solveGalerkinProjectionMI(
  Omega_h::Mesh& target_mesh, Omega_h::Mesh& source_mesh,
  const IntersectionResults& intersection, const Omega_h::Reals& source_values);

/**
 * @brief Solves the conservative Galerkin projection using Monte Carlo
 * integration of the load vector.
 *
 * @param target_mesh The target Omega_h mesh where the field is projected.
 * @param field_values_at_points Source-field values evaluated at the sampled
 *                               physical points in the target mesh.
 * @param npoints_each_tri Number of sample points used in each target triangle.
 * @param method Sampling method used to define the reference sample pattern.
 * @param loadVec_out PETSc load vector associated with the sampled projection
 *                    data, if required by the workflow.
 * @param sobol_filename Path to the file containing precomputed Sobol samples
 *                       when Sobol sampling is selected.
 * @return Projected nodal values on the target mesh.
 */
Omega_h::Reals solveGalerkinProjectionMC(
  Omega_h::Mesh& target_mesh, const Omega_h::Reals& field_values_at_points,
  const int npoints_each_tri, SamplingMethod method, Vec* loadVec_out,
  const std::string sobol_filename = "");

/**
 * @brief Computes the Galerkin right-hand-side vector using  mesh-intersection
 * method.
 *
 * This routine computes the global load vector associated with the conservative
 * Galerkin projection of a scalar source field onto the target mesh. The load
 * vector entries are evaluated by exact integration over the intersection
 * regions between the source and target meshes.
 *
 * The returned vector corresponds to the right-hand side
 * \f[
 *   b_i = \int_{\Omega} f^s(\mathbf{x}) \, \psi_i(\mathbf{x}) \,
 * d\mathbf{\Omega},
 * \f]
 * where \f$f^s\f$ is the source field and \f$\psi_i\f$ are the target basis
 * functions.
 *
 * This function computes only the RHS/load vector. It does not assemble or
 * solve the full projection system.
 *
 * @param[in] target_mesh Target Omega_h mesh on which the projected field is
 *                        represented.
 * @param[in] source_mesh Source Omega_h mesh containing the original scalar
 *                        field.
 * @param[in] intersection Precomputed intersection data between source and
 *                         target meshes.
 * @param[in] source_values Nodal scalar field values defined on the source
 *                          mesh.
 *
 * @return Global Galerkin RHS/load vector as `Omega_h::Reals`.
 *
 * @note This routine is intended for exact conservative transfer based on mesh
 *       intersections.
 * @note The returned vector may be used independently for diagnostics,
 *       verification, comparison against approximate methods, or as input to a
 *       subsequent projection solve.
 *
 * @see calculateLoadVectorMI
 */
Omega_h::Reals computeRhsVectorMI(Omega_h::Mesh& target_mesh,
                                  Omega_h::Mesh& source_mesh,
                                  const IntersectionResults& intersection,
                                  const Omega_h::Reals& source_values);

/**
 * @brief Computes the Galerkin right-hand-side vector using Monte Carlo
 * integration.
 *
 * This routine computes the global load vector associated with the conservative
 * Galerkin projection of a scalar field onto the target mesh using Monte Carlo
 * integration on each target element.
 *
 * The returned vector corresponds to the right-hand side
 * \f[
 *   b_i = \int_{\Omega} f^s(\mathbf{x}) \, \psi_i(\mathbf{x}) \, d\mathbf{x},
 * \f]
 * approximated from sampled source-field values at physical points in the
 * target mesh.
 *
 * The sampling pattern is defined on the reference triangle and may be
 * generated either by uniform random sampling or by Sobol-based sampling.
 *
 * This function computes only the RHS/load vector. It does not assemble or
 * solve the full projection system.
 *
 * @param[in] target_mesh Target Omega_h mesh on which the projected field is
 *                        represented.
 * @param[in] field_values_at_points Source-field values evaluated at the sample
 *                                   points in the target elements. The expected
 *                                   layout is element-wise flattened, with
 *                                   total size
 *                                   `target_mesh.nelems() * npoints_each_tri`.
 * @param[in] npoints_each_tri Number of sample points used in each target
 *                             triangle.
 * @param[in] method Sampling method used to define the reference barycentric
 *                   sampling method.
 * @param[in] sobol_filename Path to the file containing precomputed Sobol
 *                           barycentric samples when Sobol sampling is used.
 *
 * @return Global Galerkin RHS/load vector as `Omega_h::Reals`.
 *
 * @note This routine is intended for approximate conservative transfer using
 *       Monte Carlo integration.
 * @note The returned vector may be used independently for diagnostics,
 *       verification, comparison against exact mesh-intersection integration,
 *       or as input to a subsequent projection solve.
 *
 * @see calculateLoadVectorMC, SamplingMethod
 */
Omega_h::Reals computeRhsVectorMC(Omega_h::Mesh& target_mesh,
                                  const Omega_h::Reals& field_values_at_points,
                                  const int npoints_each_tri,
                                  SamplingMethod method,
                                  const std::string& sobol_filename = "");
} // namespace pcms

#endif // PCMS_TRANSFER_CONSERVATIVE_PROJECTION_SOLVER_HPP
