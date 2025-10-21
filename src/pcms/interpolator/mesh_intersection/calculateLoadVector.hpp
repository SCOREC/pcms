/**
 * @file calculateLoadVector.hpp
 * @brief Functions for calculating load vector on finite element meshes
 *
 * This file contains functions for creating and computing load vectors
 * for finite element calculations using Omega_h mesh structures and PETSc.
 */

#ifndef COMPUTING_AT_SCALE_CALCULATE_LOAD_VECTOR_HPP
#define COMPUTING_AT_SCALE_CALCULATE_LOAD_VECTOR_HPP
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

#include "loadVectorIntegrator.hpp"
#include <MeshField.hpp>

#include <petscmat.h>

// detect floating point exceptions
#include <fenv.h>

using ExecutionSpace = Kokkos::DefaultExecutionSpace;
using MemorySpace = Kokkos::DefaultExecutionSpace::memory_space;

/**
 * @brief Calculates the load
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
inline PetscErrorCode calculateLoadVector(
  Omega_h::Mesh& target_mesh, Omega_h::Mesh& source_mesh,
  const IntersectionResults& intersection, const Omega_h::Reals& source_values,
  Vec* loadVec_out)
{

  PetscFunctionBeginUser;
  const int numNodesPerTri = 3;

  const int nnz = target_mesh.nelems() * numNodesPerTri;

  // Allocate COO indices and values
  PetscInt* coo_i;
  PetscScalar* coo_vals;
  PetscCall(PetscMalloc2(nnz, &coo_i, nnz, &coo_vals));

  // Fill COO global indices and values
  auto elmVerts = Omega_h::HostRead(target_mesh.ask_elem_verts());
  auto elmLoadVector =
    buildLoadVector(target_mesh, source_mesh, intersection, source_values);

  auto hostElmLoadVector = Kokkos::create_mirror_view(elmLoadVector);
  Kokkos::deep_copy(hostElmLoadVector, elmLoadVector);

  /*
  std::cerr << " DEBUG: printing the load vector \n";
  for (int i = 0; i < hostElmLoadVector.size(); ++i) {
    std::cout << hostElmLoadVector[i] << " ";
  }
  std::cout << "\n";

  target_mesh.add_tag(2, "elmLoadVector", 3,
               Omega_h::read(Omega_h::Write<MeshField::Real>(elmLoadVector)));

  Omega_h::vtk::write_parallel("loadVector.vtk", &target_mesh, 2);
  */
  PetscInt idx = 0;
  for (PetscInt e = 0; e < target_mesh.nelems(); ++e) {
    for (PetscInt vi = 0; vi < numNodesPerTri; ++vi) {
      coo_i[idx] = elmVerts[numNodesPerTri * e + vi];
      coo_vals[idx] = hostElmLoadVector(numNodesPerTri * e + vi);
      ++idx;
    }
  }

  // create vector with preallocated COO structure
  Vec vec;
  PetscCall(VecCreate(PETSC_COMM_WORLD, &vec));
  PetscCall(VecSetSizes(vec, target_mesh.nverts(), PETSC_DECIDE));
  PetscCall(VecSetFromOptions(vec));
  PetscCall(VecSetPreallocationCOO(vec, nnz, coo_i));
  PetscCall(VecSetValuesCOO(vec, coo_vals, ADD_VALUES));
  PetscCall(PetscFree2(coo_i, coo_vals));

  if (target_mesh.nelems() < 10) {
    PetscCall(VecView(vec, PETSC_VIEWER_STDOUT_WORLD));
  }

  *loadVec_out = vec;
  PetscFunctionReturn(PETSC_SUCCESS);
}

#endif
