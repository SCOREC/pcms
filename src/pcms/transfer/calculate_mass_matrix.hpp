/**
 * @file calculateMassMatrix.hpp
 * @brief Functions for calculating mass matrices on finite element meshes
 * @author [Cameron Smith]
 * @date April, 2025
 *
 * This file contains functions for creating and computing mass matrices
 * for finite element calculations using Omega_h mesh structures and PETSc.
 */

#ifndef PCMS_TRANSFER_CALCULATE_MASS_MATRIX_HPP
#define PCMS_TRANSFER_CALCULATE_MASS_MATRIX_HPP

#include <Omega_h_mesh.hpp>
#include <petscmat.h>

namespace pcms
{
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

PetscErrorCode calculateMassMatrix(Omega_h::Mesh& mesh, Mat* mass_out);
} // namespace pcms
#endif // PCMS_TRANSFER_CALCULATE_MASS_MATRIX_HPP
