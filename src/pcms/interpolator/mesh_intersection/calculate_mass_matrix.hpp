/**
 * @file calculateMassMatrix.hpp
 * @brief Functions for calculating mass matrices on finite element meshes
 * @author [Cameron Smith]
 * @date April, 2025
 *
 * This file contains functions for creating and computing mass matrices
 * for finite element calculations using Omega_h mesh structures and PETSc.
 */

#ifndef PCMS_INTERPOLATOR_CALCULATE_MASS_MATRIX_HPP
#define PCMS_INTERPOLATOR_CALCULATE_MASS_MATRIX_HPP

#include <Omega_h_adapt.hpp>
#include <Omega_h_array_ops.hpp>
#include <Omega_h_atomics.hpp> //Omega_h::atomic_fetch_add
#include <Omega_h_build.hpp>
#include <Omega_h_class.hpp>
#include <Omega_h_compare.hpp>
#include <Omega_h_dbg.hpp>
#include <Omega_h_file.hpp> //Omega_h::binary
#include <Omega_h_for.hpp>
#include <Omega_h_recover.hpp> //project_by_fit
#include <Omega_h_shape.hpp>
#include <Omega_h_timer.hpp>
#include <iomanip> //precision
#include <iostream>
#include <petscvec_kokkos.hpp>
#include <sstream> //ostringstream

#include <pcms/interpolator/mesh_intersection/mass_matrix_integrator.hpp>
#include <MeshField.hpp>

#include <petscmat.h>

// detect floating point exceptions
#include <fenv.h>

using ExecutionSpace = Kokkos::DefaultExecutionSpace;
using MemorySpace = Kokkos::DefaultExecutionSpace::memory_space;

/**
 * @brief Sets a constant value for a field at all mesh vertices
 *
 * @tparam ShapeField Type of the shape field
 * @param mesh The Omega_h mesh
 * @param field The field to set values for
 * @param val The constant value to set at all vertices
 */
template <typename ShapeField>
void setFieldAtVertices(Omega_h::Mesh& mesh, ShapeField field,
                        const MeshField::Real val)
{
  const auto MeshDim = mesh.dim();
  auto setFieldAtVertices = KOKKOS_LAMBDA(const int& vtx)
  {
    field(0, 0, vtx, MeshField::Vertex) = val;
  };
  MeshField::parallel_for(ExecutionSpace(), {0}, {mesh.nverts()},
                          setFieldAtVertices, "setFieldAtVertices");
}

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
static PetscErrorCode CreateMatrix(Omega_h::Mesh& mesh, Mat* A)
{
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

/**
 * @brief Calculates the mass matrix for a given mesh
 *
 * This function constructs a mass matrix based on the provided mesh using
 * a finite element approach. It creates coordinate field elements, builds
 * the mass matrix using the massMatrixIntegrator, and sets up the PETSc matrix
 * with appropriate values.
 *
 * @param mesh The Omega_h mesh to calculate the mass matrix for
 * @param[out] mass_out Pointer to the resulting mass matrix
 * @return PetscErrorCode PETSc error code (PETSC_SUCCESS if successful)
 */

namespace pcms
{
inline PetscErrorCode calculateMassMatrix(Omega_h::Mesh& mesh, Mat* mass_out)
{
  PetscFunctionBeginUser;

  MeshField::OmegahMeshField<ExecutionSpace, MeshField::KokkosController> omf(
    mesh);

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
} // namespace pcms
#endif
