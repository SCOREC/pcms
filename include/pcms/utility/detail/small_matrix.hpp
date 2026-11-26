#ifndef PCMS_UTILITY_DETAIL_SMALL_MATRIX_HPP
#define PCMS_UTILITY_DETAIL_SMALL_MATRIX_HPP

#include "pcms/utility/assert.hpp"
#include "pcms/utility/types.hpp"
#include <Kokkos_Core.hpp>
#include <type_traits>

namespace pcms::detail
{

/// @brief A Dim x Dim matrix stored as its columns, m[j][i] the entry in row
///        i of column j, as Omega_h::Matrix.
template <int Dim>
using SquareMatrix = Kokkos::Array<Kokkos::Array<Real, Dim>, Dim>;

/// @brief f(std::integral_constant<int, Dim>{}) with Dim = dim.
/// @throws pcms_error unless dim is 1, 2 or 3
template <typename F>
decltype(auto) VisitDimension(int dim, F&& f)
{
  switch (dim) {
    case 1: return f(std::integral_constant<int, 1>{});
    case 2: return f(std::integral_constant<int, 2>{});
    case 3: return f(std::integral_constant<int, 3>{});
    default: throw pcms_error("VisitDimension: the dimension must be 1 to 3");
  }
}

[[nodiscard]] KOKKOS_INLINE_FUNCTION Kokkos::Array<Real, 3> Cross(
  const Kokkos::Array<Real, 3>& u, const Kokkos::Array<Real, 3>& v)
{
  return {u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
          u[0] * v[1] - u[1] * v[0]};
}

template <int Dim>
[[nodiscard]] KOKKOS_INLINE_FUNCTION Real
Determinant(const SquareMatrix<Dim>& m)
{
  static_assert(Dim >= 1 && Dim <= 3, "Determinant: 1 to 3 dimensions");
  if constexpr (Dim == 1) {
    return m[0][0];
  } else if constexpr (Dim == 2) {
    return m[0][0] * m[1][1] - m[1][0] * m[0][1];
  } else {
    const auto c = Cross(m[1], m[2]);
    return m[0][0] * c[0] + m[0][1] * c[1] + m[0][2] * c[2];
  }
}

/// @brief The inverse, in closed form, as Omega_h::invert.
template <int Dim>
[[nodiscard]] KOKKOS_INLINE_FUNCTION SquareMatrix<Dim> Invert(
  const SquareMatrix<Dim>& m)
{
  static_assert(Dim >= 1 && Dim <= 3, "Invert: 1 to 3 dimensions");
  const Real det = Determinant<Dim>(m);
  if constexpr (Dim == 1) {
    return {{{1.0 / det}}};
  } else if constexpr (Dim == 2) {
    return {{{m[1][1] / det, -m[0][1] / det}, {-m[1][0] / det, m[0][0] / det}}};
  } else {
    // row k of the inverse is the cross product of the other two columns
    const SquareMatrix<3> rows{Cross(m[1], m[2]), Cross(m[2], m[0]),
                               Cross(m[0], m[1])};
    SquareMatrix<3> inverse;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        inverse[j][i] = rows[i][j] / det;
      }
    }
    return inverse;
  }
}

/// @brief The solution x of m x = b, by Gaussian elimination with partial
///        pivoting.
template <int Dim>
[[nodiscard]] KOKKOS_INLINE_FUNCTION Kokkos::Array<Real, Dim> Solve(
  SquareMatrix<Dim> m, Kokkos::Array<Real, Dim> b)
{
  for (int c = 0; c < Dim; ++c) {
    int pivot = c;
    for (int r = c + 1; r < Dim; ++r) {
      if (Kokkos::fabs(m[c][r]) > Kokkos::fabs(m[c][pivot])) {
        pivot = r;
      }
    }
    for (int k = c; k < Dim; ++k) {
      const Real t = m[k][c];
      m[k][c] = m[k][pivot];
      m[k][pivot] = t;
    }
    const Real t = b[c];
    b[c] = b[pivot];
    b[pivot] = t;
    for (int r = c + 1; r < Dim; ++r) {
      const Real f = m[c][r] / m[c][c];
      for (int k = c; k < Dim; ++k) {
        m[k][r] -= f * m[k][c];
      }
      b[r] -= f * b[c];
    }
  }
  for (int c = Dim - 1; c >= 0; --c) {
    for (int k = c + 1; k < Dim; ++k) {
      b[c] -= m[k][c] * b[k];
    }
    b[c] /= m[c][c];
  }
  return b;
}

/// @brief m v.
template <int Dim>
[[nodiscard]] KOKKOS_INLINE_FUNCTION Kokkos::Array<Real, Dim> Multiply(
  const SquareMatrix<Dim>& m, const Kokkos::Array<Real, Dim>& v)
{
  Kokkos::Array<Real, Dim> out;
  for (int i = 0; i < Dim; ++i) {
    out[i] = m[0][i] * v[0];
  }
  for (int j = 1; j < Dim; ++j) {
    for (int i = 0; i < Dim; ++i) {
      out[i] = out[i] + m[j][i] * v[j];
    }
  }
  return out;
}

} // namespace pcms::detail

#endif // PCMS_UTILITY_DETAIL_SMALL_MATRIX_HPP
