#ifndef PCMS_INTERPOLATOR_ARRAY_OPS_HPP
#define PCMS_INTERPOLATOR_ARRAY_OPS_HPP

#include <pcms/interpolator/pcms_interpolator_aliases.hpp>
#include <Omega_h_array_ops.hpp>
#include <cmath>
#include <Omega_h_fail.hpp>
#include <Kokkos_MathematicalFunctions.hpp>

namespace pcms
{
namespace detail
{
/**
 * @brief Fills the kokkos scratch view
 *
 * @param value The value to be populated
 * @param team The team member
 * @param matrix The scratch matrix
 *
 */
KOKKOS_INLINE_FUNCTION
void fill(double value, member_type team, ScratchMatView matrix)
{

  int row = matrix.extent(0);
  int col = matrix.extent(1);
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, row), [=](int j) {
    for (int k = 0; k < col; ++k) {
      matrix(j, k) = value;
    }
  });
}

/**
 * @brief Fills the kokkos scratch view
 *
 * @param value The value to be populated
 * @param team The team member
 * @param vector The scratch vector
 *
 */
KOKKOS_INLINE_FUNCTION
void fill(double value, member_type team, ScratchVecView vector)
{

  int size = vector.extent(0);
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, size),
                       [=](int j) { vector(j) = value; });
}

KOKKOS_INLINE_FUNCTION
void find_sq_root_each(member_type team, ScratchVecView& array)
{
  int size = array.size();

  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, size), [=](int i) {
    OMEGA_H_CHECK_PRINTF(
      array(i) >= 0,
      "[Error:] Square root of a negative number is invalid!\n"
      "value is %12.6f\n",
      array(i));
    array(i) = Kokkos::sqrt(array(i));
  });
}

/**
 *
 *@brief Finds the inverse of each element in an array
 *
 *@param team The team member
 *@param array The scratch vector view
 *
 */
KOKKOS_INLINE_FUNCTION
void find_inverse_each(member_type team, ScratchVecView& array)
{
  int size = array.size();
  double eps = 1e-8;
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, size), [=](int i) {
    if (array(i) > eps) {
      array(i) = 1 / array(i);
    } else {
      array(i) = 0.0;
    }
  });
}

/**
 *@brief Evaluates the shrinkage factor in SVD (S/(S^2 + lambda)
 *
 *@param team The team member
 *@param lambda The regularization parameter
 *@param sigma The scratch vector view of singular values
 *
 */
KOKKOS_INLINE_FUNCTION
void calculate_shrinkage_factor(member_type team, double lambda,
                                ScratchVecView& sigma)
{
  int size = sigma.size();
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, size), [=](int i) {
    sigma(i) /= (sigma(i) * sigma(i) + lambda);
  });
}

/**
 *
 *@brief Finds the transpose of the scratch matrix
 *
 *@param team The team member
 *@param array The scratch matrix view
 *
 */
KOKKOS_INLINE_FUNCTION
ScratchMatView find_transpose(member_type team, const ScratchMatView& matrix)
{
  int row = matrix.extent(0);
  int column = matrix.extent(1);

  ScratchMatView transMatrix(team.team_scratch(1), column, row);
  fill(0.0, team, transMatrix);
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, row), [=](int i) {
    for (int j = 0; j < column; ++j) {
      transMatrix(i, j) = matrix(j, i);
    }
  });
  return transMatrix;
}

/**
 *
 *@brief Fills the diagonal entries of a matrix
 *
 *@param team The team member
 *@param diagonal_entries The scratch vector view of the entries thats to be
 *filled in a matrix diagonal
 *@param array The scratch matrix whose digonal entries to be filled
 *
 */
KOKKOS_INLINE_FUNCTION
void fill_diagonal(member_type team, const ScratchVecView& diagonal_entries,
                   ScratchMatView& matrix)
{
  int size = diagonal_entries.size();
  OMEGA_H_CHECK(size == matrix.extent(0) || size == matrix.extent(1));
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, size),
                       [=](int i) { matrix(i, i) = diagonal_entries(i); });
}

/**
 * @brief Scales each row of the input matrix by corresponding diagonal values.
 *
 * This function performs **row scaling** of the matrix, meaning each row `i`
 * of the matrix is multiplied by the scalar value `diagonal_entries(i)`.
 *
 * Mathematically, if the matrix is `A` and the diagonal scaling matrix is `D_A
 * = diag(diagonal_entries)`, the operation performed is:
 *
 *     A_scaled = D_A * A
 *
 * where `D_A` is a diagonal matrix that scales rows of `A`.
 *
 * @param team Kokkos team member used for team-level parallelism
 * @param diagonal_entries A vector of diagonal values used for scaling (should
 * match or be smaller than the number of rows)
 * @param matrix The 2D matrix whose rows will be scaled
 */
KOKKOS_INLINE_FUNCTION
void eval_row_scaling(member_type team, ScratchVecView diagonal_entries,
                      ScratchMatView matrix)
{

  int row = matrix.extent(0);
  int column = matrix.extent(1);
  int vector_size = diagonal_entries.size();

  OMEGA_H_CHECK_PRINTF(
    vector_size <= row,
    "[ERROR]: for row scaling the size of diagonal entries vector should be "
    "equal or less than the row of the matrix which is to be scaled\n"
    "size of vector = %d, row of a matrix = %d\n",
    vector_size, row);
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, row), [=](int i) {
    for (int j = 0; j < column; ++j) {
      matrix(i, j) *=
        diagonal_entries(i); // scales each row of a matrix with a corresponding
    } // value in a digonal entries
  });
}

/**
 * @brief Scales the right-hand-side (RHS) vector using the given diagonal
 * weights.
 *
 * This function performs element-wise scaling of the RHS vector. Each entry
 * `rhs_values(i)` is multiplied by the corresponding entry
 * `diagonal_entries(i)`.
 *
 * Mathematically:
 *     rhs_scaled(i) = diagonal_entries(i) * rhs_values(i)
 *
 * @param team Kokkos team member used for team-level parallelism
 * @param diagonal_entries A vector of scaling weights (same size as rhs_values)
 * @param rhs_values The RHS vector to be scaled
 */
KOKKOS_INLINE_FUNCTION
void eval_rhs_scaling(member_type team, ScratchVecView diagonal_entries,
                      ScratchVecView rhs_values)
{

  int weight_size = diagonal_entries.size();
  int rhs_size = rhs_values.size();

  OMEGA_H_CHECK_PRINTF(
    weight_size == rhs_size,
    "[ERROR]: for row scaling the size of diagonal entries vector should be "
    "equal or less than the row of the matrix which is to be scaled\n"
    "size of weight vector = %d, size of a rhs vector = %d\n",
    weight_size, rhs_size);
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, weight_size), [=](int i) {
    rhs_values(i) *= diagonal_entries(
      i); // scales each element  of a rhs with a corresponding weights
  });
}

/**
 * @brief Scales the input matrix row-wise and writes the result into an
 * adjusted matrix.
 *
 * This function first scales the rows of `matrixToScale` using the
 * corresponding entries from `diagonal_entries` (i.e., performs D * A, where D
 * is a diagonal matrix), and then copies the scaled values into
 * `adjustedMatrix`.
 *
 * Preconditions:
 * - The number of columns in `adjustedMatrix` must match the number of rows in
 * `matrixToScale` (colB == rowA)
 * - `diagonal_entries` must have size ≤ rowA (validated inside
 * `eval_row_scaling`)
 *
 * Typical use case: used in evaluating SV in U^TSV for SVD solver
 *
 * @param team Kokkos team member used for team-level parallelism
 * @param diagonal_entries A vector of scaling values for each row
 * @param matrixToScale The matrix to be row-scaled (modified in place)
 * @param adjustedMatrix The output matrix to store the scaled version
 */
KOKKOS_INLINE_FUNCTION
void scale_and_adjust(member_type team, ScratchVecView& diagonal_entries,
                      ScratchMatView& matrixToScale,
                      ScratchMatView& adjustedMatrix)
{
  size_t rowA = matrixToScale.extent(0);

  size_t rowB = adjustedMatrix.extent(0);
  size_t colB = adjustedMatrix.extent(1);
  OMEGA_H_CHECK(colB == rowA);
  eval_row_scaling(team, diagonal_entries, matrixToScale);
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, rowB), [=](int i) {
    for (int j = 0; j < colB; ++j) {
      adjustedMatrix(i, j) = matrixToScale(i, j);
    }
  });
}

} // namespace detail
} // namespace pcms

#endif
