#include <petscksp.h>

#include "pcms/transfer/conservative_projection_solver.hpp"
#include "pcms/utility/arrays.h"
#include "pcms/transfer/petsc_utils.hpp"
#include <Kokkos_Core.hpp>

namespace pcms
{

static Omega_h::Reals vecToOmegaHReals(Vec vec)
{
  PetscInt n = 0;
  PetscErrorCode ierr = VecGetSize(vec, &n);
  CHKERRABORT(PETSC_COMM_SELF, ierr);

  const PetscScalar* array = nullptr;
  ierr = VecGetArrayRead(vec, &array);
  CHKERRABORT(PETSC_COMM_SELF, ierr);

  auto values_host = Omega_h::HostWrite<Omega_h::Real>(n);
  for (PetscInt i = 0; i < n; ++i) {
    values_host[i] = array[i];
  }

  ierr = VecRestoreArrayRead(vec, &array);
  CHKERRABORT(PETSC_COMM_SELF, ierr);

  return Omega_h::Reals(values_host);
}

// ---------------------------------------------------------------------------
// GalerkinProjectionSolver
// ---------------------------------------------------------------------------

GalerkinProjectionSolver::GalerkinProjectionSolver(
  const BilinearFormIntegrator& mass_integrator,
  LinearFormIntegrator& rhs_integrator)
  : rhs_integrator_(&rhs_integrator)
{
  Mat A = mass_integrator.GetMatrix();

  PetscInt m = 0, n = 0;
  PetscErrorCode ierr = MatGetSize(A, &m, &n);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  nverts_ = m;

  const std::size_t num_pts =
    rhs_integrator_->GetIntegrationPoints().GetValues().extent(0);
  sampled_values_ =
    Kokkos::View<Real**, DeviceMemorySpace>("rhs_sampled", num_pts, 1);

  ierr = KSPCreate(PETSC_COMM_SELF, &ksp_);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  ierr = KSPSetOperators(ksp_, A, A);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  ierr = KSPSetFromOptions(ksp_);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  ierr = KSPSetUp(ksp_);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
}

GalerkinProjectionSolver::~GalerkinProjectionSolver()
{
  if (ksp_) {
    KSPDestroy(&ksp_);
  }
  // mat_ is owned by the BilinearFormIntegrator; KSP holds its own reference.
}

Omega_h::Reals GalerkinProjectionSolver::Solve(
  const PointEvaluator<Real>& evaluator, const Field<Real>& source_field) const
{
  evaluator.Evaluate(source_field, MakeRank2View(sampled_values_));
  return Solve(MakeConstRank2View(sampled_values_));
}

Omega_h::Reals GalerkinProjectionSolver::Solve(
  Rank2View<const Real, DeviceMemorySpace> sampled_values) const
{
  rhs_integrator_->Assemble(sampled_values);
  Vec rhs_vector = rhs_integrator_->GetVector();

  Vec solution = nullptr;
  PetscErrorCode ierr = createSeqVec(PETSC_COMM_SELF, nverts_, &solution);
  CHKERRABORT(PETSC_COMM_SELF, ierr);

  ierr = KSPSolve(ksp_, rhs_vector, solution);
  CHKERRABORT(PETSC_COMM_SELF, ierr);

  auto result = vecToOmegaHReals(solution);

  ierr = VecDestroy(&solution);
  CHKERRABORT(PETSC_COMM_SELF, ierr);

  return result;
}

} // namespace pcms
