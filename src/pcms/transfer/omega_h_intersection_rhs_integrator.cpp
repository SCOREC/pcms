#include "pcms/transfer/omega_h_intersection_rhs_integrator.hpp"
#include "pcms/transfer/petsc_utils.hpp"
#include "pcms/utility/arrays.h"
#include "pcms/utility/assert.h"
#include <Kokkos_Core.hpp>
#include <petscksp.h>

namespace pcms
{

OmegaHIntersectionRHSIntegrator::OmegaHIntersectionRHSIntegrator(
  const FunctionSpace& source_space, const FunctionSpace& target_space)
  : OmegaHIntersectionRHSIntegrator(
      std::make_shared<OmegaHIntersectionQuadrature>(source_space,
                                                     target_space))
{
}

OmegaHIntersectionRHSIntegrator::OmegaHIntersectionRHSIntegrator(
  const FunctionSpace& source_space, const FunctionSpace& target_space,
  std::shared_ptr<const MeshIntersection> intersection)
  : OmegaHIntersectionRHSIntegrator(
      std::make_shared<OmegaHIntersectionQuadrature>(source_space, target_space,
                                                     std::move(intersection)))
{
}

OmegaHIntersectionRHSIntegrator::OmegaHIntersectionRHSIntegrator(
  std::shared_ptr<const OmegaHLagrangeLayout> source_layout,
  CoordinateSystem source_coordinate_system,
  std::shared_ptr<const OmegaHLagrangeLayout> target_layout,
  CoordinateSystem target_coordinate_system)
  : OmegaHIntersectionRHSIntegrator(
      std::make_shared<OmegaHIntersectionQuadrature>(
        std::move(source_layout), source_coordinate_system,
        std::move(target_layout), target_coordinate_system))
{
}

OmegaHIntersectionRHSIntegrator::OmegaHIntersectionRHSIntegrator(
  std::shared_ptr<const OmegaHIntersectionQuadrature> quadrature)
  : quadrature_(std::move(quadrature))
{
  if (!quadrature_) {
    throw pcms_error("OmegaHIntersectionRHSIntegrator: quadrature must be set");
  }
  const auto target_dofs = quadrature_->GetTargetDofs();
  const PetscInt nnz = static_cast<PetscInt>(target_dofs.extent(0));
  PetscErrorCode ierr =
    createSeqVec(PETSC_COMM_SELF, quadrature_->GetNumTargetDofs(), &vec_);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  // VecSetPreallocationCOO takes the COO indices on the host
  // TODO: ask Todd/PETSc folks if there is a better way to do this for GPU
  // support
  auto target_dofs_host =
    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, target_dofs);
  ierr = VecSetPreallocationCOO(vec_, nnz, target_dofs_host.data());
  CHKERRABORT(PETSC_COMM_SELF, ierr);
}

OmegaHIntersectionRHSIntegrator::~OmegaHIntersectionRHSIntegrator()
{
  if (vec_) {
    VecDestroy(&vec_);
  }
}

const OmegaHIntersectionQuadrature&
OmegaHIntersectionRHSIntegrator::GetQuadrature() const noexcept
{
  return *quadrature_;
}

CoordinateView<DeviceMemorySpace>
OmegaHIntersectionRHSIntegrator::GetIntegrationPoints() const noexcept
{
  return quadrature_->GetIntegrationPoints();
}

Vec OmegaHIntersectionRHSIntegrator::GetVector() const noexcept
{
  return vec_;
}

void OmegaHIntersectionRHSIntegrator::Assemble(
  Rank2View<const Real, DeviceMemorySpace> sampled_values)
{
  const int ndof = quadrature_->GetDofsPerElement();
  const std::size_t num_pts =
    static_cast<std::size_t>(quadrature_->GetNumPoints());
  PCMS_ALWAYS_ASSERT(static_cast<std::size_t>(sampled_values.extent(0)) ==
                     num_pts);
  PCMS_ALWAYS_ASSERT(sampled_values.extent(1) >= 1);

  PetscErrorCode ierr = VecZeroEntries(vec_);
  CHKERRABORT(PETSC_COMM_SELF, ierr);

  auto sv = Kokkos::View<const Real**, Kokkos::LayoutRight, DeviceMemorySpace,
                         Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
    sampled_values.data_handle(), sampled_values.extent(0),
    sampled_values.extent(1));
  Kokkos::View<PetscScalar*, DeviceMemorySpace> coo_vals("rhs_coo_vals",
                                                         num_pts * ndof);
  auto coeffs = quadrature_->GetWeights();
  Kokkos::parallel_for(
    "rhs_coo_vals", static_cast<int>(num_pts), KOKKOS_LAMBDA(int i) {
      const PetscScalar f = static_cast<PetscScalar>(sv(i, 0));
      for (int j = 0; j < ndof; ++j) {
        coo_vals(i * ndof + j) =
          static_cast<PetscScalar>(coeffs(i * ndof + j)) * f;
      }
    });

  ierr = VecSetValuesCOO(vec_, coo_vals.data(), ADD_VALUES);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
}

// ---------------------------------------------------------------------------
// Builder
// ---------------------------------------------------------------------------

std::unique_ptr<LinearFormIntegrator> BuildOmegaHConservativeRHSIntegrator(
  const FunctionSpace& source_space, const FunctionSpace& target_space)
{
  return std::make_unique<OmegaHIntersectionRHSIntegrator>(source_space,
                                                           target_space);
}

} // namespace pcms
