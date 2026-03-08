#include <petscvec_kokkos.hpp>
#include <petscmat_kokkos.hpp>
#include <petscksp.h>

#include "pcms/transfer/conservative_projection_solver.hpp"
#include "pcms/transfer/calculate_mass_matrix.hpp"

namespace pcms {
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
static Vec solveLinearSystem(Mat A, Vec b) {
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
                "Estimated condition number of matrix A: %.6e\n",
                smax / smin);
  } else {
    PetscPrintf(PETSC_COMM_WORLD,
                "Condition number estimate unavailable (smin <= 0 or error)\n");
  }

  ierr = KSPDestroy(&ksp);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  return x;
}
Omega_h::Reals solveGalerkinProjection(Omega_h::Mesh &target_mesh,
                                       Omega_h::Mesh &source_mesh,
                                       const IntersectionResults &intersection,
                                       const Omega_h::Reals &source_values) {
  if ((PetscInt) source_values.size() !=
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
  ierr = calculateLoadVector(target_mesh,
                             source_mesh,
                             intersection,
                             source_values,
                             &vec);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  Vec x = solveLinearSystem(mass, vec);

  Omega_h::Write<Omega_h::Real> solution_vector(
    target_mesh.nverts(),
    0,
    "stores the solution coefficients");

  PetscScalar *array;

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
Omega_h::Reals rhsVectorMI(Omega_h::Mesh &target_mesh,
                           Omega_h::Mesh &source_mesh,
                           const IntersectionResults &intersection,
                           const Omega_h::Reals &source_values) {
  Vec vec;
  PetscErrorCode ierr;
  ierr = calculateLoadVector(target_mesh,
                             source_mesh,
                             intersection,
                             source_values,
                             &vec);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  Omega_h::Write<Omega_h::Real> rhsvector(target_mesh.nverts(),
                                          0,
                                          "stores the rhs vector"); {
    PetscErrorCode ierr;
    PetscScalar *array;
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
}
