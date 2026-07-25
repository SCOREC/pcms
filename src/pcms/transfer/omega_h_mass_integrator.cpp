#include "pcms/transfer/omega_h_mass_integrator.hpp"
#include "pcms/transfer/mass_matrix_integrator.hpp"
#include "pcms/transfer/omega_h_form_integrator_utils.hpp"
#include "pcms/transfer/petsc_utils.hpp"
#include "pcms/utility/assert.h"
#include "pcms/utility/memory_spaces.h"
#include <KokkosController.hpp>
#include <MeshField.hpp>
#include <MeshField_Element.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_shape.hpp>
#include <petscksp.h>

namespace pcms
{

namespace
{

// Fills the diagonal COO entries of a P0 (piecewise-constant) mass matrix: one
// entry per element on its own DOF, valued at the element area.
void FillP0MassCoo(
  int nelems, const Omega_h::Reals& coords, const Omega_h::LOs& faces2nodes,
  const Kokkos::View<const LO*, DeviceMemorySpace>& global_to_local,
  const Kokkos::View<PetscInt*, DeviceMemorySpace>& coo_rows,
  const Kokkos::View<PetscInt*, DeviceMemorySpace>& coo_cols,
  const Kokkos::View<PetscScalar*, DeviceMemorySpace>& vals)
{
  Kokkos::parallel_for(
    "mass_p0_diag", nelems, KOKKOS_LAMBDA(int e) {
      const auto verts = Omega_h::gather_verts<3>(faces2nodes, e);
      const Omega_h::Matrix<2, 3> vm =
        Omega_h::gather_vectors<3, 2>(coords, verts);
      Omega_h::Few<Omega_h::Vector<2>, 2> basis;
      basis[0] = vm[1] - vm[0];
      basis[1] = vm[2] - vm[0];
      const Omega_h::Real area =
        Kokkos::fabs(Omega_h::triangle_area_from_basis(basis));
      const PetscInt g = static_cast<PetscInt>(global_to_local(e));
      coo_rows(e) = g;
      coo_cols(e) = g;
      vals(e) = static_cast<PetscScalar>(area);
    });
}

// Fills the 3x3-block COO sparsity pattern of a P1 (linear) mass matrix: each
// element contributes a dense block coupling its three vertex DOFs.
void FillP1MassCooPattern(
  int nelems, const Omega_h::LOs& faces2nodes,
  const Kokkos::View<const LO*, DeviceMemorySpace>& global_to_local,
  const Kokkos::View<PetscInt*, DeviceMemorySpace>& coo_rows,
  const Kokkos::View<PetscInt*, DeviceMemorySpace>& coo_cols)
{
  Kokkos::parallel_for(
    "mass_coo_pattern", nelems, KOKKOS_LAMBDA(int e) {
      const auto verts = Omega_h::gather_verts<3>(faces2nodes, e);
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          const int idx = e * 9 + i * 3 + j;
          coo_rows(idx) = static_cast<PetscInt>(global_to_local(verts[i]));
          coo_cols(idx) = static_cast<PetscInt>(global_to_local(verts[j]));
        }
      }
    });
}

} // namespace

OmegaHMassIntegrator::OmegaHMassIntegrator(const FunctionSpace& target_space)
  : OmegaHMassIntegrator(std::dynamic_pointer_cast<const OmegaHLagrangeLayout>(
                           target_space.GetLayout()),
                         target_space.GetCoordinateSystem())
{
}

OmegaHMassIntegrator::OmegaHMassIntegrator(
  std::shared_ptr<const OmegaHLagrangeLayout> target_layout,
  CoordinateSystem coordinate_system)
{
  detail::CheckOmegaHScalarLagrangeLayout(coordinate_system, target_layout,
                                          "OmegaHMassIntegrator", "target");

  Omega_h::Mesh& mesh = target_layout->GetMesh();
  const auto global_to_local = target_layout->GetGlobalToLocalPermutation();
  const PetscInt num_dofs =
    static_cast<PetscInt>(target_layout->GetNumOwnedDofHolder());
  const int nelems = mesh.nelems();

  if (target_layout->GetOrder() == 0) {
    // P0 target: piecewise-constant basis functions have disjoint support, so
    // the mass matrix is diagonal with M_ee = area(e). One COO entry per
    // element on its own (element-id) diagonal.
    const auto& coords = mesh.coords();
    const auto& faces2nodes = mesh.ask_down(Omega_h::FACE, Omega_h::VERT).ab2b;

    const PetscInt nnz = static_cast<PetscInt>(nelems);
    Kokkos::View<PetscInt*, DeviceMemorySpace> coo_rows("mass_coo_rows", nnz);
    Kokkos::View<PetscInt*, DeviceMemorySpace> coo_cols("mass_coo_cols", nnz);
    Kokkos::View<PetscScalar*, DeviceMemorySpace> vals("mass_vals", nnz);
    FillP0MassCoo(nelems, coords, faces2nodes, global_to_local, coo_rows,
                  coo_cols, vals);

    PetscErrorCode ierr =
      createSeqAIJMat(PETSC_COMM_SELF, num_dofs, num_dofs, 0, nullptr, &mat_);
    CHKERRABORT(PETSC_COMM_SELF, ierr);
    auto coo_rows_host =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, coo_rows);
    auto coo_cols_host =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, coo_cols);
    ierr = MatSetPreallocationCOO(mat_, nnz, coo_rows_host.data(),
                                  coo_cols_host.data());
    CHKERRABORT(PETSC_COMM_SELF, ierr);
    ierr = MatSetValuesCOO(mat_, vals.data(), INSERT_VALUES);
    CHKERRABORT(PETSC_COMM_SELF, ierr);
    return;
  }

  // P1 target: consistent mass matrix assembled from MeshField per-element 3x3
  // blocks. (Higher MeshField orders extend this branch via
  // getTriangleElement.)
  MeshField::OmegahMeshField<DefaultExecutionSpace, 2,
                             MeshField::KokkosController>
    omf(mesh);
  auto coordField = omf.getCoordField();
  const auto [shp, map] = MeshField::Omegah::getTriangleElement<1>(mesh);
  MeshField::FieldElement coordFe(mesh.nelems(), coordField, shp, map);
  auto elm_mass_dev = buildElementMassMatrix(mesh, coordFe);

  // Build COO sparsity pattern on device: each element contributes a 3x3 block.
  const auto& faces2nodes = mesh.ask_down(Omega_h::FACE, Omega_h::VERT).ab2b;

  const PetscInt nnz = static_cast<PetscInt>(nelems) * 9;
  Kokkos::View<PetscInt*, DeviceMemorySpace> coo_rows("mass_coo_rows", nnz);
  Kokkos::View<PetscInt*, DeviceMemorySpace> coo_cols("mass_coo_cols", nnz);
  FillP1MassCooPattern(nelems, faces2nodes, global_to_local, coo_rows,
                       coo_cols);

  // Create sparse matrix, preallocate with COO pattern, then bulk-set values.
  // elm_mass_dev is in the same element-major order as coo_rows/coo_cols, so
  // it can be passed directly to MatSetValuesCOO — no host copy needed.
  PetscErrorCode ierr =
    createSeqAIJMat(PETSC_COMM_SELF, num_dofs, num_dofs, 0, nullptr, &mat_);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  auto coo_rows_host =
    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, coo_rows);
  auto coo_cols_host =
    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, coo_cols);
  ierr = MatSetPreallocationCOO(mat_, nnz, coo_rows_host.data(),
                                coo_cols_host.data());
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  ierr = MatSetValuesCOO(mat_, elm_mass_dev.data(), INSERT_VALUES);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
}

OmegaHMassIntegrator::~OmegaHMassIntegrator()
{
  if (mat_) {
    MatDestroy(&mat_);
  }
}

Mat OmegaHMassIntegrator::GetMatrix() const noexcept
{
  return mat_;
}

std::unique_ptr<BilinearFormIntegrator> BuildOmegaHMassIntegrator(
  const FunctionSpace& target_space)
{
  return std::make_unique<OmegaHMassIntegrator>(target_space);
}

} // namespace pcms
