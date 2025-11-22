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
 * @created by Abhiyan Paudel
 * @date August 2025
 */

#ifndef PCMS_INTERPOLATOR_GALERKIN_PROJECTION_SOLVER_HPP
#ifndef PCMS_INTERPOLATOR_GALERKIN_PROJECTION_SOLVER_HPP

#include <Kokkos_Core.hpp>
#include <MeshField_Shape.hpp>
#include <Omega_h_array.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_shape.hpp>
#include <pcms/point_search.h>
#include <petscvec_kokkos.hpp> // Must be included before other PETSc headers
#include <petscksp.h>
#include <petscmat_kokkos.hpp>
#include <petscsys.h>
#include <petscvec.h>

#include <pcms/interpolator/mesh_intersection/calculate_load_vector.hpp>
#include <pcms/interpolator/mesh_intersection/calculate_mass_matrix.hpp>

/**
 * @brief Solves a linear system Ax = b using PETSc's KSP solvers
 *
 * Uses PETSc's Krylov Subspace solvers to find x in Ax = b.
 * The solver can be configured through PETSc runtime options.
 *
 * @param A The system matrix
 * @param b The right-hand side vector
 * @return Vec Solution vector x
 */
static Vec solveLinearSystem(Mat A, Vec b)
{
  PetscInt m, n;
  PetscErrorCode ierr;

  ierr = MatGetSize(A, &m, &n);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  Vec x;
  ierr = VecCreateSeqKokkos(PETSC_COMM_SELF, n, &x);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  KSP ksp;
  ierr = KSPCreate(PETSC_COMM_SELF, &ksp);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  ierr = KSPSetOperators(ksp, A, A);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  ierr = KSPSetComputeSingularValues(ksp, PETSC_TRUE);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  ierr = KSPSetFromOptions(ksp);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  ierr = KSPSetUp(ksp);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  ierr = KSPSolve(ksp, b, x);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  /// compute and print condition number estimate
  PetscReal smax = 0.0, smin = 0.0;
  ierr = KSPComputeExtremeSingularValues(ksp, &smax, &smin);
  if (!ierr && smin > 0.0) {
    PetscPrintf(PETSC_COMM_WORLD,
                "Estimated condition number of matrix A: %.6e\n", smax / smin);
  } else {
    PetscPrintf(PETSC_COMM_WORLD,
                "Condition number estimate unavailable (smin <= 0 or error)\n");
  }

  ierr = KSPDestroy(&ksp);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  return x;
}

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

namespace pcms
{
Omega_h::Reals solveGalerkinProjection(Omega_h::Mesh& target_mesh,
                                       Omega_h::Mesh& source_mesh,
                                       const IntersectionResults& intersection,
                                       const Omega_h::Reals& source_values)
{

  if ((PetscInt)source_values.size() !=
      source_mesh.coords().size() / source_mesh.dim()) {
    std::cerr << "ERROR: source_values size (" << source_values.size()
              << ") doesn't match expected size ("
              << source_mesh.coords().size() / source_mesh.dim() << ")"
              << std::endl;
    throw std::runtime_error("source_values length mismatch");
  }

  Mat mass;
  PetscErrorCode ierr = calculateMassMatrix(target_mesh, &mass);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  Vec vec;
  ierr = calculateLoadVector(target_mesh, source_mesh, intersection,
                             source_values, &vec);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  Vec x = solveLinearSystem(mass, vec);

  Omega_h::Write<Omega_h::Real> solution_vector(
    target_mesh.nverts(), 0, "stores the solution coefficients");

  PetscScalar* array;

  ierr = VecGetArray(x, &array);

  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  auto solution_host = Omega_h::HostWrite<Omega_h::Real>(target_mesh.nverts());

  for (PetscInt i = 0; i < target_mesh.nverts(); ++i) {
    solution_host[i] = array[i];
  }

  solution_vector = Omega_h::Write<Omega_h::Real>(solution_host);
  ierr = VecRestoreArray(x, &array);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  ierr = VecDestroy(&x);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  ierr = MatDestroy(&mass);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  ierr = VecDestroy(&vec);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  return Omega_h::read(solution_vector);
}

Omega_h::Reals rhsVectorMI(Omega_h::Mesh& target_mesh,
                           Omega_h::Mesh& source_mesh,
                           const IntersectionResults& intersection,
                           const Omega_h::Reals& source_values)
{
  Vec vec;
  PetscErrorCode ierr;
  ierr = calculateLoadVector(target_mesh, source_mesh, intersection,
                             source_values, &vec);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  Omega_h::Write<Omega_h::Real> rhsvector(target_mesh.nverts(), 0,
                                          "stores the rhs vector");
  {
    PetscErrorCode ierr;
    PetscScalar* array;
    ierr = VecGetArray(vec, &array);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);

    auto rhsvec_host = Omega_h::HostWrite<Omega_h::Real>(target_mesh.nverts());

    for (PetscInt i = 0; i < target_mesh.nverts(); ++i) {
      rhsvec_host[i] = array[i];
    }

    rhsvector = Omega_h::Write<Omega_h::Real>(rhsvec_host);
    ierr = VecRestoreArray(vec, &array);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
  }

  ierr = VecDestroy(&vec);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  return Omega_h::read(rhsvector);
}
} // namespace pcms

#endif // PCMS_INTERPOLATOR_GALERKIN_PROJECTION_SOLVER_HPP
