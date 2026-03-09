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

#include <pcms/transfer/mass_matrix_integrator.hpp>
#include <pcms/utility/memory_spaces.h>
#include <MeshField.hpp>

#include <petscmat.h>

// detect floating point exceptions
#include <fenv.h>

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
#endif
