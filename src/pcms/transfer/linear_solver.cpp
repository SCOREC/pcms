#include "pcms/transfer/linear_solver.hpp"
#include "pcms/configuration.hpp"
#include "pcms/field/field_ops.hpp"
#include "pcms/utility/assert.hpp"
#include <Kokkos_Core.hpp>
#ifdef PCMS_ENABLE_PETSC
#include "pcms/transfer/petsc_utils.hpp"
#include <petscksp.h>
#endif
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

namespace pcms::transfer::detail
{

namespace
{

#ifdef PCMS_ENABLE_PETSC
// Owns a PETSc object and destroys it.
template <typename Handle, PetscErrorCode (*Destroy)(Handle*)>
class PetscOwner
{
public:
  PetscOwner() = default;
  PetscOwner(const PetscOwner&) = delete;
  PetscOwner& operator=(const PetscOwner&) = delete;
  ~PetscOwner()
  {
    if (handle_ != nullptr) {
      Destroy(&handle_);
    }
  }
  Handle* operator&() noexcept { return &handle_; }
  operator Handle() const noexcept { return handle_; }

private:
  Handle handle_ = nullptr;
};

using OwnedMat = PetscOwner<Mat, MatDestroy>;
using OwnedVec = PetscOwner<Vec, VecDestroy>;
using OwnedKsp = PetscOwner<KSP, KSPDestroy>;
#endif

} // namespace

void Divide(const Field<Real>& numerator, const Field<Real>& denominator,
            Field<Real>& out)
{
  out.SetDOFHolderDataUnchecked(MakeConstRank2View(pcms::detail::MapOwned<Real>(
    numerator.GetDOFHolderData().GetValues(),
    denominator.GetDOFHolderData().GetValues(),
    KOKKOS_LAMBDA(const Real n, const Real d) { return n / d; })));
  out.Synchronize();
}

struct LinearSolver::Petsc
{
#ifdef PCMS_ENABLE_PETSC
  OwnedMat matrix;
  OwnedKsp ksp;
  OwnedVec rhs;
  OwnedVec solution;
  PetscInt size = 0;
  // the row of the matrix of each owned DOF holder
  std::vector<PetscInt> row_of_owned;
#endif
};

LinearSolver::LinearSolver(const BilinearForm& form, LinearSolveOptions options)
  : form_(form), options_(std::move(options))
{
  const auto space = form_.Form().Space();
  const auto& layout = *space->GetLayout();
  if (options_.method == LinearSolveOptions::Method::MatrixFree) {
    diagonal_.emplace(ScratchFunction(*space));
    cg_ = std::make_unique<CgScratch>(
      CgScratch{ScratchFunction(*space), ScratchFunction(*space),
                ScratchFunction(*space), ScratchFunction(*space)});
    if (form_.SupportsMatrix() && !layout.IsDistributed()) {
      form_.AssembleDiagonal(*diagonal_);
    } else {
      diagonal_.reset();
    }
    return;
  }

#ifndef PCMS_ENABLE_PETSC
  throw pcms_error("LinearSolver: the assembled solve needs PETSc; configure "
                   "PCMS with -DPCMS_ENABLE_PETSC=ON or use the MatrixFree "
                   "method");
#else
  if (layout.IsDistributed()) {
    throw pcms_error("LinearSolver: the assembled solve is serial; the "
                     "layout is distributed");
  }
  petsc_ = std::make_unique<Petsc>();
  const auto coo = form_.AssembleCoo();
  const auto gids = layout.GetGidsHost();
  petsc_->size = static_cast<PetscInt>(gids.extent(0));
  // the rows of the matrix are the DOF holders in the order of their gids
  std::vector<GO> sorted(gids.extent(0));
  for (std::size_t i = 0; i < sorted.size(); ++i) {
    sorted[i] = gids(i);
  }
  std::sort(sorted.begin(), sorted.end());
  const auto row_of = [&](GO gid) {
    return static_cast<PetscInt>(
      std::lower_bound(sorted.begin(), sorted.end(), gid) - sorted.begin());
  };
  const auto rows_host =
    Kokkos::create_mirror_view_and_copy(HostMemorySpace{}, coo.rows);
  const auto columns_host =
    Kokkos::create_mirror_view_and_copy(HostMemorySpace{}, coo.columns);
  const auto values_host =
    Kokkos::create_mirror_view_and_copy(HostMemorySpace{}, coo.values);
  const auto nnz = static_cast<PetscCount>(rows_host.extent(0));
  std::vector<PetscInt> rows(nnz), columns(nnz);
  std::vector<PetscScalar> values(nnz);
  for (PetscCount k = 0; k < nnz; ++k) {
    rows[k] = row_of(rows_host(k));
    columns[k] = row_of(columns_host(k));
    values[k] = values_host(k);
  }
  petsc_->row_of_owned.resize(gids.extent(0));
  for (std::size_t o = 0; o < sorted.size(); ++o) {
    petsc_->row_of_owned[o] = row_of(gids(o));
  }

  const auto n = petsc_->size;
  PetscErrorCode ierr =
    createSeqAIJMat(PETSC_COMM_SELF, n, n, 0, nullptr, &petsc_->matrix);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  ierr =
    MatSetPreallocationCOO(petsc_->matrix, nnz, rows.data(), columns.data());
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  ierr = MatSetValuesCOO(petsc_->matrix, values.data(), INSERT_VALUES);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  ierr = createSeqVec(PETSC_COMM_SELF, n, &petsc_->rhs);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  ierr = createSeqVec(PETSC_COMM_SELF, n, &petsc_->solution);
  CHKERRABORT(PETSC_COMM_SELF, ierr);

  ierr = KSPCreate(PETSC_COMM_SELF, &petsc_->ksp);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  ierr = KSPSetOperators(petsc_->ksp, petsc_->matrix, petsc_->matrix);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  // a direct solve unless the options say otherwise
  ierr = KSPSetType(petsc_->ksp, KSPPREONLY);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  PC pc = nullptr;
  ierr = KSPGetPC(petsc_->ksp, &pc);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  ierr = PCSetType(pc, PCLU);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  if (!options_.petsc_options_prefix.empty()) {
    ierr =
      KSPSetOptionsPrefix(petsc_->ksp, options_.petsc_options_prefix.c_str());
    CHKERRABORT(PETSC_COMM_SELF, ierr);
  }
  ierr = KSPSetFromOptions(petsc_->ksp);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  ierr = KSPSetUp(petsc_->ksp);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
#endif
}

LinearSolver::~LinearSolver() = default;

int LinearSolver::Solve(const Field<Real>& functional,
                        Field<Real>& solution) const
{
  if (&functional.GetLayout() != &solution.GetLayout()) {
    throw pcms_error("LinearSolver::Solve: the functional and the solution "
                     "must have the same layout");
  }
#ifdef PCMS_ENABLE_PETSC
  if (petsc_ != nullptr) {
    // one KSP solve per component
    const auto b = functional.GetDOFHolderDataHost().GetValues();
    Kokkos::View<Real**, HostMemorySpace> x("linear_solver_solution",
                                            b.extent(0), b.extent(1));
    for (std::size_t c = 0; c < b.extent(1); ++c) {
      PetscScalar* array = nullptr;
      PetscErrorCode ierr = VecGetArray(petsc_->rhs, &array);
      CHKERRABORT(PETSC_COMM_SELF, ierr);
      for (std::size_t o = 0; o < b.extent(0); ++o) {
        array[petsc_->row_of_owned[o]] = b(o, c);
      }
      ierr = VecRestoreArray(petsc_->rhs, &array);
      CHKERRABORT(PETSC_COMM_SELF, ierr);
      ierr = KSPSolve(petsc_->ksp, petsc_->rhs, petsc_->solution);
      CHKERRABORT(PETSC_COMM_SELF, ierr);
      const PetscScalar* result = nullptr;
      ierr = VecGetArrayRead(petsc_->solution, &result);
      CHKERRABORT(PETSC_COMM_SELF, ierr);
      for (std::size_t o = 0; o < b.extent(0); ++o) {
        x(o, c) = result[petsc_->row_of_owned[o]];
      }
      ierr = VecRestoreArrayRead(petsc_->solution, &result);
      CHKERRABORT(PETSC_COMM_SELF, ierr);
    }
    solution.SetDOFHolderDataUncheckedHost(MakeConstRank2View(x));
    solution.Synchronize();
    return 0;
  }
#endif

  // Preconditioned conjugate gradients, all components jointly.
  auto& r = cg_->residual;
  auto& p = cg_->direction;
  auto& ap = cg_->image;
  auto& z = cg_->preconditioned;
  const auto precondition = [&]() {
    if (diagonal_) {
      Divide(r, *diagonal_, z);
    } else {
      pcms::detail::CopyValues(r, z);
    }
  };
  pcms::detail::SetZero(solution);
  pcms::detail::CopyValues(functional, r);
  const Real target = options_.relative_tolerance * Norm2(functional);
  if (Norm2(r) <= target) {
    return 0;
  }
  precondition();
  pcms::detail::CopyValues(z, p);
  Real rz = Apply(r, z);
  const int max_iterations =
    options_.max_iterations > 0
      ? options_.max_iterations
      : static_cast<int>(form_.Form().Space()->GetLayout()->GlobalSize());
  for (int iteration = 1; iteration <= max_iterations; ++iteration) {
    form_.Apply(p, ap);
    const Real curvature = Apply(ap, p);
    if (!(curvature > 0)) {
      break;
    }
    const Real alpha = rz / curvature;
    Axpy(alpha, p, solution);
    Axpy(-alpha, ap, r);
    if (Norm2(r) <= target) {
      return iteration;
    }
    precondition();
    const Real rz_next = Apply(r, z);
    Scale(rz_next / rz, p);
    Axpy(Real(1), z, p);
    rz = rz_next;
  }
  throw pcms_error("LinearSolver::Solve: the conjugate gradients did not "
                   "converge");
}

} // namespace pcms::transfer::detail
