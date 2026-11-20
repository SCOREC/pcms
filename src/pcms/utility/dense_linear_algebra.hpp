#ifndef PCMS_UTILITY_DENSE_LINEAR_ALGEBRA_HPP
#define PCMS_UTILITY_DENSE_LINEAR_ALGEBRA_HPP

// Private: dense linear algebra of the n x n matrices (n <= 3) of coordinate
// systems, maps and cells, usable in kernels.

#include "pcms/utility/types.hpp"
#include <Kokkos_Core.hpp>

namespace pcms::detail
{

KOKKOS_INLINE_FUNCTION Real Determinant(const Real (&m)[3][3], int n)
{
  if (n == 1) {
    return m[0][0];
  }
  if (n == 2) {
    return m[0][0] * m[1][1] - m[0][1] * m[1][0];
  }
  return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
         m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
         m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}

// Solves m x = b by Gaussian elimination with partial pivoting; b becomes x
// and m is overwritten.
KOKKOS_INLINE_FUNCTION void Solve(Real (&m)[3][3], int n, Real (&b)[3])
{
  for (int c = 0; c < n; ++c) {
    int pivot = c;
    for (int r = c + 1; r < n; ++r) {
      if (Kokkos::fabs(m[r][c]) > Kokkos::fabs(m[pivot][c])) {
        pivot = r;
      }
    }
    for (int k = 0; k < n; ++k) {
      const Real t = m[c][k];
      m[c][k] = m[pivot][k];
      m[pivot][k] = t;
    }
    const Real t = b[c];
    b[c] = b[pivot];
    b[pivot] = t;
    for (int r = c + 1; r < n; ++r) {
      const Real f = m[r][c] / m[c][c];
      for (int k = c; k < n; ++k) {
        m[r][k] -= f * m[c][k];
      }
      b[r] -= f * b[c];
    }
  }
  for (int c = n - 1; c >= 0; --c) {
    for (int k = c + 1; k < n; ++k) {
      b[c] -= m[c][k] * b[k];
    }
    b[c] /= m[c][c];
  }
}

// The inverse of m, column by column.
KOKKOS_INLINE_FUNCTION void Invert(const Real (&m)[3][3], int n,
                                   Real (&inverse)[3][3])
{
  for (int c = 0; c < n; ++c) {
    Real a[3][3] = {};
    Real column[3] = {};
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        a[i][j] = m[i][j];
      }
    }
    column[c] = 1;
    Solve(a, n, column);
    for (int i = 0; i < n; ++i) {
      inverse[i][c] = column[i];
    }
  }
}

} // namespace pcms::detail

#endif // PCMS_UTILITY_DENSE_LINEAR_ALGEBRA_HPP
