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

// Measure (area in 2D, volume in 3D) of a simplex from its vertex-difference
// basis.
template <int Dim>
KOKKOS_INLINE_FUNCTION Omega_h::Real SimplexMeasure(
  const Omega_h::Few<Omega_h::Vector<Dim>, Dim>& basis)
{
  if constexpr (Dim == 3) {
    return Kokkos::fabs(Omega_h::tet_volume_from_basis(basis));
  } else {
    return Kokkos::fabs(Omega_h::triangle_area_from_basis(basis));
  }
}

// Fills the diagonal COO entries of a P0 (piecewise-constant) mass matrix: one
// entry per element on its own DOF, valued at the element measure.
template <int Dim>
void FillP0MassCoo(
  int nelems, const Omega_h::Reals& coords, const Omega_h::LOs& elems2nodes,
  const Kokkos::View<const LO*, DeviceMemorySpace>& global_to_local,
  const Kokkos::View<PetscInt*, DeviceMemorySpace>& coo_rows,
  const Kokkos::View<PetscInt*, DeviceMemorySpace>& coo_cols,
  const Kokkos::View<PetscScalar*, DeviceMemorySpace>& vals)
{
  Kokkos::parallel_for(
    "mass_p0_diag", nelems, KOKKOS_LAMBDA(int e) {
      const auto verts = Omega_h::gather_verts<Dim + 1>(elems2nodes, e);
      const Omega_h::Matrix<Dim, Dim + 1> vm =
        Omega_h::gather_vectors<Dim + 1, Dim>(coords, verts);
      Omega_h::Few<Omega_h::Vector<Dim>, Dim> basis;
      for (int d = 0; d < Dim; ++d) {
        basis[d] = vm[d + 1] - vm[0];
      }
      const Omega_h::Real measure = SimplexMeasure<Dim>(basis);
      const PetscInt g = static_cast<PetscInt>(global_to_local(e));
      coo_rows(e) = g;
      coo_cols(e) = g;
      vals(e) = static_cast<PetscScalar>(measure);
    });
}

// Fills the (Dim+1)x(Dim+1)-block COO sparsity pattern of a P1 (linear) mass
// matrix: each element contributes a dense block coupling its vertex DOFs. When
// lumped, every block entry is mapped onto the row's diagonal instead; PETSc's
// COO assembly sums repeated indices, so the unmodified element values
// accumulate into the row-sum lumped diagonal.
template <int Dim>
void FillP1MassCooPattern(
  int nelems, const Omega_h::LOs& elems2nodes,
  const Kokkos::View<const LO*, DeviceMemorySpace>& global_to_local,
  bool lumped, const Kokkos::View<PetscInt*, DeviceMemorySpace>& coo_rows,
  const Kokkos::View<PetscInt*, DeviceMemorySpace>& coo_cols)
{
  constexpr int nv = Dim + 1;
  Kokkos::parallel_for(
    "mass_coo_pattern", nelems, KOKKOS_LAMBDA(int e) {
      const auto verts = Omega_h::gather_verts<nv>(elems2nodes, e);
      for (int i = 0; i < nv; ++i) {
        for (int j = 0; j < nv; ++j) {
          const int idx = e * (nv * nv) + i * nv + j;
          const auto row = static_cast<PetscInt>(global_to_local(verts[i]));
          coo_rows(idx) = row;
          coo_cols(idx) =
            lumped ? row : static_cast<PetscInt>(global_to_local(verts[j]));
        }
      }
    });
}

// Assembles the target-space mass matrix for spatial dimension Dim (triangles
// for Dim==2, tetrahedra for Dim==3) and returns the owned PETSc matrix.
template <int Dim>
Mat BuildOmegaHMassMatrixImpl(Omega_h::Mesh& mesh,
                              const OmegaHLagrangeLayout& target_layout,
                              MassMatrixType mass_type)
{
  const auto global_to_local = target_layout.GetGlobalToLocalPermutation();
  const PetscInt num_dofs =
    static_cast<PetscInt>(target_layout.GetNumLocalDofHolder());
  const int nelems = mesh.nelems();
  Mat mat = nullptr;

  if (target_layout.GetOrder() == 0) {
    // P0 target: piecewise-constant basis functions have disjoint support, so
    // the mass matrix is diagonal with M_ee = measure(e). One COO entry per
    // element on its own (element-id) diagonal.
    const auto& coords = mesh.coords();
    const auto& elems2nodes = mesh.ask_down(Dim, Omega_h::VERT).ab2b;

    const PetscInt nnz = static_cast<PetscInt>(nelems);
    Kokkos::View<PetscInt*, DeviceMemorySpace> coo_rows("mass_coo_rows", nnz);
    Kokkos::View<PetscInt*, DeviceMemorySpace> coo_cols("mass_coo_cols", nnz);
    Kokkos::View<PetscScalar*, DeviceMemorySpace> vals("mass_vals", nnz);
    FillP0MassCoo<Dim>(nelems, coords, elems2nodes, global_to_local, coo_rows,
                       coo_cols, vals);

    PetscErrorCode ierr =
      createSeqAIJMat(PETSC_COMM_SELF, num_dofs, num_dofs, 0, nullptr, &mat);
    CHKERRABORT(PETSC_COMM_SELF, ierr);
    auto coo_rows_host =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, coo_rows);
    auto coo_cols_host =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, coo_cols);
    ierr = MatSetPreallocationCOO(mat, nnz, coo_rows_host.data(),
                                  coo_cols_host.data());
    CHKERRABORT(PETSC_COMM_SELF, ierr);
    ierr = MatSetValuesCOO(mat, vals.data(), INSERT_VALUES);
    CHKERRABORT(PETSC_COMM_SELF, ierr);
    return mat;
  }

  // P1 target: consistent mass matrix assembled from MeshField per-element
  // (Dim+1)x(Dim+1) blocks.
  MeshField::OmegahMeshField<DefaultExecutionSpace, Dim,
                             MeshField::KokkosController>
    omf(mesh);
  auto coordField = omf.getCoordField();
  Kokkos::View<MeshField::Real*> elm_mass_dev;
  if constexpr (Dim == 3) {
    const auto [shp, map] = MeshField::Omegah::getTetrahedronElement<1>(mesh);
    MeshField::FieldElement coordFe(mesh.nelems(), coordField.field, shp, map);
    elm_mass_dev = buildElementMassMatrix(mesh, coordFe);
  } else {
    const auto [shp, map] = MeshField::Omegah::getTriangleElement<1>(mesh);
    MeshField::FieldElement coordFe(mesh.nelems(), coordField.field, shp, map);
    elm_mass_dev = buildElementMassMatrix(mesh, coordFe);
  }

  // Build COO sparsity pattern on device: each element contributes a
  // (Dim+1)x(Dim+1) block.
  const auto& elems2nodes = mesh.ask_down(Dim, Omega_h::VERT).ab2b;

  constexpr int nv = Dim + 1;
  const PetscInt nnz = static_cast<PetscInt>(nelems) * (nv * nv);
  Kokkos::View<PetscInt*, DeviceMemorySpace> coo_rows("mass_coo_rows", nnz);
  Kokkos::View<PetscInt*, DeviceMemorySpace> coo_cols("mass_coo_cols", nnz);
  FillP1MassCooPattern<Dim>(nelems, elems2nodes, global_to_local,
                            mass_type == MassMatrixType::Lumped, coo_rows,
                            coo_cols);

  // Create sparse matrix, preallocate with COO pattern, then bulk-set values.
  // elm_mass_dev is in the same element-major order as coo_rows/coo_cols, so
  // it can be passed directly to MatSetValuesCOO — no host copy needed.
  PetscErrorCode ierr =
    createSeqAIJMat(PETSC_COMM_SELF, num_dofs, num_dofs, 0, nullptr, &mat);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  auto coo_rows_host =
    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, coo_rows);
  auto coo_cols_host =
    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, coo_cols);
  ierr = MatSetPreallocationCOO(mat, nnz, coo_rows_host.data(),
                                coo_cols_host.data());
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  ierr = MatSetValuesCOO(mat, elm_mass_dev.data(), INSERT_VALUES);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  return mat;
}

} // namespace

OmegaHMassIntegrator::OmegaHMassIntegrator(const FunctionSpace& target_space,
                                           MassMatrixType mass_type)
  : OmegaHMassIntegrator(std::dynamic_pointer_cast<const OmegaHLagrangeLayout>(
                           target_space.GetLayout()),
                         target_space.GetCoordinateSystem(), mass_type)
{
}

OmegaHMassIntegrator::OmegaHMassIntegrator(
  std::shared_ptr<const OmegaHLagrangeLayout> target_layout,
  CoordinateSystem coordinate_system, MassMatrixType mass_type)
{
  detail::CheckOmegaHScalarLagrangeLayout(coordinate_system, target_layout,
                                          "OmegaHMassIntegrator", "target");

  Omega_h::Mesh& mesh = target_layout->GetMesh();
  diagonal_ =
    target_layout->GetOrder() == 0 || mass_type == MassMatrixType::Lumped;
  if (mesh.dim() == 3) {
    mat_ = BuildOmegaHMassMatrixImpl<3>(mesh, *target_layout, mass_type);
  } else {
    mat_ = BuildOmegaHMassMatrixImpl<2>(mesh, *target_layout, mass_type);
  }
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

bool OmegaHMassIntegrator::IsDiagonal() const noexcept
{
  return diagonal_;
}

std::unique_ptr<BilinearFormIntegrator> BuildOmegaHMassIntegrator(
  const FunctionSpace& target_space, MassMatrixType mass_type)
{
  return std::make_unique<OmegaHMassIntegrator>(target_space, mass_type);
}

} // namespace pcms
