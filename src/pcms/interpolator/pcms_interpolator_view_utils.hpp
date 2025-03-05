#ifndef PCMS_INTERPOLATOR_ARRAY_OPS_HPP
#define PCMS_INTERPOLATOR_ARRAY_OPS_HPP

#include "pcms_interpolator_aliases.hpp"
#include <Omega_h_array_ops.hpp>
#include <cmath>
#include <Omega_h_fail.hpp>

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
    array(i) = sqrt(array(i));
  });
}

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

KOKKOS_INLINE_FUNCTION
ScratchMatView find_transpose(member_type team, const ScratchMatView& matrix)
{
  int row = matrix.extent(0);
  int column = matrix.extent(1);

  ScratchMatView transMatrix(team.team_scratch(0), column, row);
  fill(0.0, team, transMatrix);
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, row), [=](int i) {
    for (int j = 0; j < column; ++j) {
      transMatrix(j, i) = matrix(j, i);
    }
  });
  return transMatrix;
}

KOKKOS_INLINE_FUNCTION
void fill_diagonal(member_type team, const ScratchVecView& diagonal_entries,
                   ScratchMatView& matrix)
{
  int size = diagonal_entries.size();
  OMEGA_H_CHECK(size == matrix.extent(0) || size == matrix.extent(1));
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, size),
                       [=](int i) { matrix(i, i) = diagonal_entries(i); });
}

KOKKOS_INLINE_FUNCTION
void eval_row_scaling(member_type team, const ScratchVecView& diagonal_entries,
                      ScratchMatView& matrix)
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
                             // value in a digonal entries
    });

}

KOKKOS_INLINE_FUNCTION
void scale_and_adjust(member_type team, const ScratchVecView& diagonal_entries, ScratchMatView& matrixToScale, ScratchMatView& adjustedMatrix){
    size_t rowA = matrixToScale.extent(0);
    size_t columnA = matrixToScale.extent(1);

    size_t rowB = adjustedMatrix.extent(0);
    size_t columnB = adjustedMatrix.extent(0);

    OMEGA_CHECK(columnB == rowA)
    eval_row_scaling(team, diagonal_entries, matrixtoScale);
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, rowB), [=](int i) {
      for (int j = 0; j < columnB; == j) {
        adjustedMatrix(i, j) = matrixtoScale(i, j);
      }
    });
}
} // end namespace pcms
} // end namespace detail
#endif
