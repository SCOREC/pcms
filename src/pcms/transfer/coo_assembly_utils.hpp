#ifndef PCMS_TRANSFER_COO_ASSEMBLY_UTILS_HPP
#define PCMS_TRANSFER_COO_ASSEMBLY_UTILS_HPP

#include <Omega_h_mesh.hpp>
#include <petscmat.h>

namespace pcms
{

// Builds COO row indices for linear triangles (3 vertices per element).
// Order is deterministic by (element, local_vertex).
// Allocates `rows_out` with PetscMalloc1; caller owns and must PetscFree.
PetscErrorCode build_linear_triangle_coo_rows(Omega_h::Mesh& mesh,
                                              PetscInt** rows_out,
                                              PetscInt* nnz_out);

// Builds COO row/column index pairs for linear-triangle element matrices.
// Order is deterministic by (element, local_row_vertex, local_col_vertex).
// Allocates `rows_out` and `cols_out` with PetscMalloc2; caller owns and must
// PetscFree2.
PetscErrorCode build_linear_triangle_coo_rows_cols(Omega_h::Mesh& mesh,
                                                   PetscInt** rows_out,
                                                   PetscInt** cols_out,
                                                   PetscInt* nnz_out);

} // namespace pcms

#endif // PCMS_TRANSFER_COO_ASSEMBLY_UTILS_HPP
