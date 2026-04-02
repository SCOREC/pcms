#include "pcms/transfer/calculate_mass_matrix.hpp"
#include "pcms/transfer/coo_assembly_utils.hpp"
#include "pcms/transfer/mass_matrix_integrator.hpp"
#include "pcms/transfer/petsc_utils.hpp"
#include "pcms/utility/memory_spaces.h"

namespace pcms
{
/**
 * @brief Creates a PETSc matrix based on mesh connectivity
 *
 * This function creates a sparse matrix with the proper sparsity pattern
 * according to the mesh connectivity. The matrix size corresponds to the
 * number of vertices in the mesh.
 *
 * @param mesh The Omega_h mesh to create the matrix from
 * @param[out] A Pointer to the PETSc matrix to be created
 * @return PetscErrorCode PETSc error code (PETSC_SUCCESS if successful)
 */
static PetscErrorCode create_linear_triangle_coo_matrix(Omega_h::Mesh& mesh,
                                                        Mat* A)
{
  PetscInt* coo_rows = nullptr;
  PetscInt* coo_cols = nullptr;
  PetscInt matSize = 0;
  PetscFunctionBeginUser;
  PetscCall(createSeqAIJMat(PETSC_COMM_WORLD, mesh.nverts(), mesh.nverts(), 0,
                            nullptr, A));
  PetscCall(
    build_linear_triangle_coo_rows_cols(mesh, &coo_rows, &coo_cols, &matSize));
  PetscCall(MatSetPreallocationCOO(*A, matSize, coo_rows, coo_cols));
  PetscCall(PetscFree2(coo_rows, coo_cols));
  PetscFunctionReturn(PETSC_SUCCESS);
}
PetscErrorCode calculateMassMatrix(Omega_h::Mesh& mesh, Mat* mass_out)
{
  PetscFunctionBeginUser;

  MeshField::OmegahMeshField<DefaultExecutionSpace, 2,
                             MeshField::KokkosController>
    omf(mesh);

  const auto ShapeOrder = 1;
  auto coordField = omf.getCoordField();
  const auto [shp, map] =
    MeshField::Omegah::getTriangleElement<ShapeOrder>(mesh);
  MeshField::FieldElement coordFe(mesh.nelems(), coordField, shp, map);

  auto elmMassMatrix = buildMassMatrix(mesh, coordFe);
  const PetscInt expected_nnz = 9 * mesh.nelems();
  PetscCheck(static_cast<PetscInt>(elmMassMatrix.extent(0)) == expected_nnz,
             PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ,
             "Element mass data size (%d) does not match expected linear COO "
             "entries (%d)",
             static_cast<PetscInt>(elmMassMatrix.extent(0)), expected_nnz);

  Mat mass;
  PetscCall(create_linear_triangle_coo_matrix(mesh, &mass));
  PetscCall(MatZeroEntries(mass));
  PetscCall(
    MatSetValuesCOO(mass, elmMassMatrix.data(),
                    INSERT_VALUES)); // FIXME fails here on gpu, calls into host
  // implementation... AFAIK, petsc checks
  // the type of the input array of values to
  // decide which backend to use...
  //
  if (mesh.nelems() < 10) {
    PetscCall(MatView(mass, PETSC_VIEWER_STDOUT_WORLD));
  }

  *mass_out = mass;
  PetscFunctionReturn(PETSC_SUCCESS);
}
} // namespace pcms
