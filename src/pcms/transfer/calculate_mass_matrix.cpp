#include "pcms/transfer/calculate_mass_matrix.hpp"

namespace pcms {
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
static PetscErrorCode CreateMatrix(Omega_h::Mesh &mesh, Mat*A) {
  const auto numNodesPerTri = 3; // FIXME query the mesh
  const auto matSize = numNodesPerTri * numNodesPerTri * mesh.nelems();
  auto elmVerts = Omega_h::HostRead(mesh.ask_elem_verts());
  PetscInt *oor, *ooc, cnt = 0;
  PetscFunctionBeginUser;
  PetscCall(MatCreate(PETSC_COMM_WORLD, A));
  PetscCall(
    MatSetSizes(*A, mesh.nverts(), mesh.nverts(), PETSC_DECIDE, PETSC_DECIDE));
  PetscCall(MatSetFromOptions(*A));
  /* determine for each entry in each element stiffness matrix the global row
   * and column */
  /* since the element is triangular with piecewise linear basis functions there
   * are three degrees of freedom per element, one for each vertex */
  PetscCall(PetscMalloc2(matSize, &oor, matSize, &ooc));
  for (PetscInt e = 0; e < mesh.nelems(); e++) {
    for (PetscInt vi = 0; vi < numNodesPerTri; vi++) {
      for (PetscInt vj = 0; vj < numNodesPerTri; vj++) {
        oor[cnt] = elmVerts[numNodesPerTri * e + vi];
        ooc[cnt++] = elmVerts[numNodesPerTri * e + vj];
      }
    }
  }
  PetscCall(MatSetPreallocationCOO(*A, matSize, oor, ooc));
  PetscCall(PetscFree2(oor, ooc));
  PetscFunctionReturn(PETSC_SUCCESS);
}
PetscErrorCode calculateMassMatrix(Omega_h::Mesh &mesh, Mat*mass_out) {
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

  auto host_elmMassMatrix = Kokkos::create_mirror_view(elmMassMatrix);

  Mat mass;
  PetscCall(CreateMatrix(mesh, &mass));
  PetscBool is_kokkos;
  PetscCall(
    PetscObjectBaseTypeCompare((PetscObject)mass, MATSEQAIJKOKKOS, &is_kokkos));
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
}
