#include "pcms/transfer/mass_smoother.hpp"

#include <petscmat.h>

#include "pcms/field/field.h"
#include "pcms/field/layout/omega_h_lagrange.h"
#include "pcms/transfer/omega_h_mass_integrator.hpp"
#include "pcms/transfer/petsc_utils.hpp"
#include "pcms/utility/arrays.h"
#include "pcms/utility/assert.h"
#include "pcms/utility/profile.h"

namespace pcms
{

namespace
{

void Check(PetscErrorCode ierr)
{
  CHKERRABORT(PETSC_COMM_SELF, ierr);
}

void CheckSmootherMatrix(const BilinearFormIntegrator& mass, PetscInt num_rows,
                         Vec lumped_mass)
{
  if (mass.IsDiagonal()) {
    throw pcms_error("MassSmoother: the mass matrix is diagonal, so the "
                     "smoother would be the identity");
  }
  Mat M = mass.GetMatrix();
  PetscInt m = 0;
  PetscInt n = 0;
  Check(MatGetSize(M, &m, &n));
  if (m != num_rows || n != num_rows) {
    throw pcms_error("MassSmoother: mass matrix size does not match the "
                     "layout's number of DOF holders");
  }

  Check(MatGetRowSum(M, lumped_mass));
  PetscReal min_row_sum = 0.0;
  Check(VecMin(lumped_mass, nullptr, &min_row_sum));
  if (!(min_row_sum > 0.0)) {
    throw pcms_error("MassSmoother: the mass matrix has a non-positive row "
                     "sum; inv(M_L) is undefined");
  }

  for (PetscInt row = 0; row < m; ++row) {
    PetscInt ncols = 0;
    const PetscScalar* vals = nullptr;
    Check(MatGetRow(M, row, &ncols, nullptr, &vals));
    bool negative = false;
    for (PetscInt k = 0; k < ncols; ++k) {
      negative = negative || (PetscRealPart(vals[k]) < 0.0);
    }
    Check(MatRestoreRow(M, row, &ncols, nullptr, &vals));
    if (negative) {
      throw pcms_error("MassSmoother: the mass matrix has a negative entry, "
                       "so a sweep would not be a convex combination");
    }
  }
}

} // namespace

MassSmoother::MassSmoother(std::shared_ptr<const FieldLayout> layout,
                           std::unique_ptr<BilinearFormIntegrator> mass)
  : layout_(std::move(layout)), mass_(std::move(mass))
{
  if (!layout_ || !mass_) {
    throw pcms_error("MassSmoother: layout and mass integrator must be set");
  }
  const auto num_rows = static_cast<PetscInt>(layout_->GetNumOwnedDofHolder());
  Check(createSeqVec(PETSC_COMM_SELF, num_rows, &lumped_mass_));
  CheckSmootherMatrix(*mass_, num_rows, lumped_mass_);
  Check(createSeqVec(PETSC_COMM_SELF, num_rows, &q_));
  Check(createSeqVec(PETSC_COMM_SELF, num_rows, &tmp_));
  host_values_.resize(static_cast<std::size_t>(layout_->OwnedSize()));
}

MassSmoother::MassSmoother(const FunctionSpace& space,
                           std::unique_ptr<BilinearFormIntegrator> mass)
  : MassSmoother(space.GetLayout(), std::move(mass))
{
}

MassSmoother::~MassSmoother()
{
  for (Vec* v : {&lumped_mass_, &q_, &tmp_}) {
    if (*v != nullptr) {
      VecDestroy(v);
    }
  }
}

void MassSmoother::Apply(Field<Real>& field) const
{
  PCMS_FUNCTION_TIMER;
  if (&field.GetLayout() != layout_.get()) {
    throw pcms_error("MassSmoother::Apply: field layout mismatch");
  }
  const LO num_holders = layout_->GetNumOwnedDofHolder();
  const int num_components = layout_->GetNumComponents();
  const auto values = field.GetDOFHolderDataHost();
  const auto row_of_holder = layout_->GetGlobalToLocalPermutationHost();
  Mat M = mass_->GetMatrix();

  for (int c = 0; c < num_components; ++c) {
    PetscScalar* q = nullptr;
    Check(VecGetArray(q_, &q));
    for (LO i = 0; i < num_holders; ++i) {
      q[row_of_holder(i)] = values(i, c);
    }
    Check(VecRestoreArray(q_, &q));

    Check(MatMult(M, q_, tmp_));
    Check(VecPointwiseDivide(q_, tmp_, lumped_mass_));

    const PetscScalar* smoothed = nullptr;
    Check(VecGetArrayRead(q_, &smoothed));
    for (LO i = 0; i < num_holders; ++i) {
      host_values_[static_cast<std::size_t>(i) * num_components + c] =
        PetscRealPart(smoothed[row_of_holder(i)]);
    }
    Check(VecRestoreArrayRead(q_, &smoothed));
  }

  field.SetDOFHolderDataHost(Rank2View<const Real, HostMemorySpace>(
    host_values_.data(), num_holders, num_components));
}

std::unique_ptr<BilinearFormIntegrator> BuildMassIntegrator(
  const FunctionSpace& space, MassMatrixType mass_type)
{
  if (std::dynamic_pointer_cast<const OmegaHLagrangeLayout>(
        space.GetLayout())) {
    return BuildOmegaHMassIntegrator(space, mass_type);
  }
  throw pcms_error(
    "BuildMassIntegrator: no mass integrator for this layout type");
}

std::unique_ptr<MassSmoother> BuildMassSmoother(const FunctionSpace& space)
{
  return std::make_unique<MassSmoother>(space.GetLayout(),
                                        BuildMassIntegrator(space));
}

} // namespace pcms
