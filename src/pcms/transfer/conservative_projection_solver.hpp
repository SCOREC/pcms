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

#ifndef PCMS_TRANSFER_GALERKIN_PROJECTION_SOLVER_HPP
#define PCMS_TRANSFER_GALERKIN_PROJECTION_SOLVER_HPP

#include <Omega_h_array.hpp>
#include <Omega_h_mesh.hpp>

#include <pcms/transfer/calculate_load_vector.hpp>

namespace pcms
{

/**
 * @brief Solves a conservative galerkin projection problem to transfer scalar
 * field values onto a target mesh.
 *
 * This function assembles and solves a linear system of the form:
 * \f[
 *   M \cdot x = f
 * \f]
 * where:
 * - \f$M\f$ is the mass matrix on the target mesh (based on P1 finite
 * elements),
 * - \f$f\f$ is the load vector computed on the supermesh,
 * - \f$x\f$ is the unknown nodal field on the target mesh (solution).
 *
 * The method computes the conservative field transfer between two non-matching
 * meshes using mesh  intersections (supermesh).
 *
 * ### Algorithm Steps:
 * 1. Compute and assemble mass matrix and load vector
 * 2. Solve the linear system using PETSc.
 * 3. Return the solution as a nodal field on the target mesh.
 *
 * @param target_mesh The Omega_h mesh where the field is projected.
 * @param source_mesh The Omega_h mesh containing the original field data.
 * @param intersection Precomputed intersection information between source and
 * target meshes.
 * @param source_values Nodal scalar field values on the source mesh.
 *
 * @return A vector of nodal values on the target mesh after projection
 * (Omega_h::Reals).
 *
 *
 */

Omega_h::Reals solveGalerkinProjection(Omega_h::Mesh& target_mesh,
                                       Omega_h::Mesh& source_mesh,
                                       const IntersectionResults& intersection,
                                       const Omega_h::Reals& source_values);

Omega_h::Reals rhsVectorMI(Omega_h::Mesh& target_mesh,
                           Omega_h::Mesh& source_mesh,
                           const IntersectionResults& intersection,
                           const Omega_h::Reals& source_values);
} // namespace pcms

#endif // PCMS_TRANSFER_GALERKIN_PROJECTION_SOLVER_HPP
