#include "pcms/transfer/coo_assembly_utils.hpp"

namespace pcms
{
namespace
{
constexpr PetscInt num_nodes_per_tri = 3;
}

PetscErrorCode build_linear_triangle_coo_rows(Omega_h::Mesh& mesh,
                                              PetscInt** rows_out,
                                              PetscInt* nnz_out)
{
  PetscFunctionBeginUser;
  PetscCheck(rows_out != nullptr, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
             "rows_out must not be null");
  PetscCheck(nnz_out != nullptr, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
             "nnz_out must not be null");
  PetscCheck(mesh.dim() == 2, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG,
             "Expected 2D mesh, got dim=%d", mesh.dim());
  PetscCheck(mesh.family() == OMEGA_H_SIMPLEX, PETSC_COMM_SELF,
             PETSC_ERR_ARG_WRONG,
             "Expected simplex mesh family for linear-triangle COO assembly");

  const PetscInt nnz = mesh.nelems() * num_nodes_per_tri;

  auto elm_verts = Omega_h::HostRead(mesh.ask_elem_verts());
  PetscInt* rows = nullptr;
  PetscCall(PetscMalloc1(nnz, &rows));

  PetscInt idx = 0;
  for (PetscInt e = 0; e < mesh.nelems(); ++e) {
    for (PetscInt vi = 0; vi < num_nodes_per_tri; ++vi) {
      rows[idx++] = elm_verts[num_nodes_per_tri * e + vi];
    }
  }
  PetscCheck(idx == nnz, PETSC_COMM_SELF, PETSC_ERR_PLIB,
             "COO row fill count mismatch: idx=%d nnz=%d", idx, nnz);

  *rows_out = rows;
  *nnz_out = nnz;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode build_linear_triangle_coo_rows_cols(Omega_h::Mesh& mesh,
                                                   PetscInt** rows_out,
                                                   PetscInt** cols_out,
                                                   PetscInt* nnz_out)
{
  PetscFunctionBeginUser;
  PetscCheck(rows_out != nullptr, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
             "rows_out must not be null");
  PetscCheck(cols_out != nullptr, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
             "cols_out must not be null");
  PetscCheck(nnz_out != nullptr, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
             "nnz_out must not be null");
  PetscCheck(mesh.dim() == 2, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG,
             "Expected 2D mesh, got dim=%d", mesh.dim());
  PetscCheck(mesh.family() == OMEGA_H_SIMPLEX, PETSC_COMM_SELF,
             PETSC_ERR_ARG_WRONG,
             "Expected simplex mesh family for linear-triangle COO assembly");

  const PetscInt nnz = mesh.nelems() * num_nodes_per_tri * num_nodes_per_tri;

  auto elm_verts = Omega_h::HostRead(mesh.ask_elem_verts());
  PetscInt* rows = nullptr;
  PetscInt* cols = nullptr;
  PetscCall(PetscMalloc2(nnz, &rows, nnz, &cols));

  /* determine for each entry in each element stiffness matrix the global row
   * and column */
  /* since the element is triangular with piecewise linear basis functions
   * there are three degrees of freedom per element, one for each vertex */
  PetscInt idx = 0;
  for (PetscInt e = 0; e < mesh.nelems(); ++e) {
    for (PetscInt vi = 0; vi < num_nodes_per_tri; ++vi) {
      for (PetscInt vj = 0; vj < num_nodes_per_tri; ++vj) {
        rows[idx] = elm_verts[num_nodes_per_tri * e + vi];
        cols[idx] = elm_verts[num_nodes_per_tri * e + vj];
        ++idx;
      }
    }
  }
  PetscCheck(idx == nnz, PETSC_COMM_SELF, PETSC_ERR_PLIB,
             "COO row/col fill count mismatch: idx=%d nnz=%d", idx, nnz);

  *rows_out = rows;
  *cols_out = cols;
  *nnz_out = nnz;
  PetscFunctionReturn(PETSC_SUCCESS);
}

} // namespace pcms
