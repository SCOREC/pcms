#ifndef MLS_RBF_OPTIONS_HPP
#define MLS_RBF_OPTIONS_HPP

#include "mdspan/mdspan.hpp"
#include "pcms/arrays.h"
#include "pcms/assert.h"
#include "pcms/memory_spaces.h"
#include "pcms/types.h"
#include <Kokkos_Core.hpp>
#include <cmath>
#include <iostream>

namespace pcms
{

namespace detail
{

template <typename T, typename MemorySpace>
T get_element_from_span(Rank1View<T, MemorySpace> span, LO index)
{
  T element;
  Kokkos::View<T*, MemorySpace> data_view(span.data_handle() + index, 1);
  auto element_host = Kokkos::create_mirror_view(data_view);
  Kokkos::deep_copy(element_host, data_view);
  element = element_host(0);
  return element;
}

/** \brief Check if the grid is uniformly spaced and ascending
 *
 * \param x input grid points
 * \param ztol tolerance for uniformity check
 * \return true if the grid is uniformly spaced and ascending, false otherwise
 *
 * \details This function checks if the input grid points are uniformly spaced
 * and in ascending order. It calculates the average spacing and verifies that
 * each interval between consecutive points is within a specified tolerance of
 * this average spacing. The tolerance is defined as a fraction (ztol) of the
 * average spacing.
 */
template <typename T, typename MemorySpace>
bool isUniformAscending(Rank1View<T, MemorySpace> x, T ztol)
{
  using execution_space = typename MemorySpace::execution_space;
  LO inx = static_cast<LO>(x.extent(0));
  T dxavg = (get_element_from_span(x, inx - 1) - get_element_from_span(x, 0)) /
            (inx - 1);
  T zeps = std::abs(ztol * dxavg);
  bool is_uniform_ascending = true;

  Kokkos::parallel_reduce(
    "check_uniform_ascending", Kokkos::RangePolicy<execution_space>(1, inx),
    KOKKOS_LAMBDA(const LO ix, bool& local_flag) {
      T zdiffx = x(ix) - x(ix - 1);
      T zdiff = zdiffx - dxavg;
      if (zdiffx < 0.0 || std::abs(zdiff) > zeps) {
        local_flag = false;
      }
    },
    Kokkos::LAnd<bool>(is_uniform_ascending));
  return is_uniform_ascending;
}

template <typename T, typename MemorySpace>
void grid_valid(Rank1View<T, MemorySpace> x)
{
  PCMS_ALWAYS_ASSERT(x.extent(0) >= 2);
  bool is_uniform_ascending = isUniformAscending<T, MemorySpace>(x, 1.0e-3);
  PCMS_ALWAYS_ASSERT(is_uniform_ascending);
}

/** \brief Clamp xget to the grid range [x[0], x[nx-1]] with tolerance
 *
 * \param xget input coordinate
 * \param zxget output clamped coordinate
 * \param x input grid points
 * \param nx number of grid points
 * \param ier error flag (0 = no error, 1 = out of range)
 *
 * \details If xget is outside the range [x[0], x[nx-1]], then check if it is
 * within a small tolerance of the range. If it is within the tolerance, then
 * clamp it to the nearest boundary value. If it is outside the tolerance, then
 * set the error flag ier to 1.
 */
template <typename T, typename MemorySpace>
KOKKOS_INLINE_FUNCTION void clamp_to_grid_range(T xget, T& zxget,
                                                Rank1View<T, MemorySpace> x,
                                                LO nx, LO& ier, T tol)
{
  if (xget < x[0] || xget > x[nx - 1]) {
    T zxtol = tol * std::max(std::abs(x[0]), std::abs(x[nx - 1]));

    if (xget < x[0] - zxtol || xget > x[nx - 1] + zxtol) {
      ier = 1; // Error code for out of range
      printf("lookup:  xget=%.6f out of range %.6f to %.6f\n", xget, x[0],
             x[nx - 1]);

      return;
    } else {
      printf("lookup:  xget=%.6f beyond range %.6f to %.6f (fixup applied)\n",
             xget, x[0], x[nx - 1]);
      zxget = (xget < x[0]) ? x[0] : x[nx - 1];
    }
  }
}

/** \brief Detect if the boundary conditions are homogeneous
 *
 * \param bcmin boundary condition data
 * \param ibcmin boundary condition flag
 * \param iflg2 flag indicating if the boundary conditions are
 * \param nx number of grid points in x direction
 *
 * \details Check if given derivatives at the boundary are all zero to identify
 * homogeneous boundary conditions. A homogeneous boundary condition  set A has
 * the property that if s1(x,y) satisfies A, then (for any real number c) c*S1
 * also satisfies A, and, if both s1(x,y) and s2(x,y) satisfy A then so also
 * does s1+s2. A non-homogeneous boundary condition need extra correction steps
 * during spline setup.
 */
template <typename T, typename MemorySpace>
void detect_inhomogeneous_bc(Rank1View<T, MemorySpace> bcmin, LO ibcmin,
                             LO& iflg2, LO nx)
{
  using execution_space = typename MemorySpace::execution_space;
  if (ibcmin == 1 || ibcmin == 2) {
    Kokkos::parallel_reduce(
      "check_nonzero", Kokkos::RangePolicy<execution_space>(0, nx),
      KOKKOS_LAMBDA(const LO ix, LO& local_flag) {
        if (bcmin(ix) != 0.0) {
          local_flag = 1;
        }
      },
      Kokkos::Max<LO>(iflg2));
  }
}

/** \brief Solve for 1d spline coefficients.
 *
 * \param k_bc1 boundary condition flag at x(0)
 * \param k_bcn boundary condition flag at x(-1)
 * \param n number of grid points
 * \param x input grid points
 * \param f output coefficients of the cubic spline
 * \param wk working space for the spline solver
 *
 * \details Reference:
 *   [1] G. E. Forsythe, M. A. Malcolm, and C. B. Moler, Computer Methods for
 * Mathematical Computations. Prentice-Hall, 1977, p. 76. [2] G. Engeln-Müllges
 * and F. Uhlig, Numerical Algorithms with Fortran. Berlin: Springer-Verlag,
 * 1996, p. 251. [3] W. A. Houlberg and D. McCune, March 2000.
 */
template <typename T, typename MemorySpace>
KOKKOS_FUNCTION void compute_spline_coefficients(LO k_bc1, LO k_bcn, LO n,
                                                 Rank1View<T, MemorySpace> x,
                                                 Rank2View<T, MemorySpace> f,
                                                 Rank1View<T, MemorySpace> wk)
{
  LO i_bc1 = k_bc1;
  LO i_bcn = k_bcn;

  // Clip to allowed ranges
  if (i_bc1 < -1 || i_bc1 > 7)
    i_bc1 = 0; // outside [-1,7] -> not-a-knot
  if (i_bcn < 0 || i_bcn > 7)
    i_bcn = 0; // outside [0,7]  -> not-a-knot

  // Periodic BC handling
  if (i_bc1 == -1)
    i_bcn = -1;

  LO i3knots = 0;
  LO i3perio = 0;

  if (n == 3) {
    i_bc1 = std::min(6, i_bc1);
    i_bcn = std::min(6, i_bcn);

    if (i_bc1 == 0)
      i3knots += 1;
    if (i_bcn == 0)
      i3knots += 1;
    if (i_bc1 == -1)
      i3perio = 1;
  }

  LO iord1;
  LO iord2;

  if (n == 2) {
    if (i_bc1 == -1) {
      i_bc1 = 5;
      i_bcn = 5;
    }

    if (i_bc1 == 0 || i_bc1 > 5)
      i_bc1 = 5;
    if (i_bcn == 0 || i_bcn > 5)
      i_bcn = 5;

    // LHS match
    if (i_bc1 == 1 || i_bc1 == 3 || i_bc1 == 5)
      iord1 = 1; // first derivative match
    else
      iord1 = 2; // second derivative match
    // RHS match
    if (i_bcn == 1 || i_bcn == 3 || i_bcn == 5)
      iord2 = 1;
    else
      iord2 = 2;
  }

  T imin = 0;
  T imax = n - 1;
  T a1 = 0.0;
  T b1 = 0.0;
  T an = 0.0;
  T bn = 0.0;
  T f0, fh, h;

  const T dx = x(1) - x(0);

  if (i_bc1 == 1) {
    a1 = f(1, 0);
  } else if (i_bc1 == 2) {
    b1 = f(2, 0);
  } else if (i_bc1 == 5) {
    a1 = (f(0, 1) - f(0, 0)) / dx;
  } else if (i_bc1 == 6) {
    b1 = 2.0 * ((f(0, 2) - f(0, 1)) / dx - (f(0, 1) - f(0, 0)) / dx) / (2 * dx);
  }

  if (i_bcn == 1) {
    an = f(1, n - 1);
  } else if (i_bcn == 2) {
    bn = f(2, n - 1);
  } else if (i_bcn == 5) {
    an = (f(0, n - 1) - f(0, n - 2)) / dx;
  } else if (i_bcn == 6) {
    bn = 2.0 *
         ((f(0, n - 1) - f(0, n - 2)) / dx - (f(0, n - 2) - f(0, n - 3)) / dx) /
         (2 * dx);
  }
  f(1, n - 1) = 0.0;
  f(2, n - 1) = 0.0;
  f(3, n - 1) = 0.0;
  if (n == 2) {
    if (i_bc1 == 5 && i_bcn == 5) {
      // Coefficients for n = 2
      f(1, 0) = (f(0, 1) - f(0, 0)) / dx;
      f(2, 0) = 0.0;
      f(3, 0) = 0.0;
      f(1, 1) = f(1, 0);
      f(2, 1) = 0.0;
      f(3, 1) = 0.0;
    } else if (iord1 == 1 && iord2 == 1) {
      T a1 = f(1, 0), an = f(1, 1);
      f(1, 0) = a1;
      f(1, 1) = an;
      h = dx;
      f0 = f(0, 0);
      fh = f(0, 1);

      f(2, 0) = (3 * (fh - f0) / (h * h) - (2 * a1 + an) / h) * 2; // 2*c2
      f(3, 0) =
        (-2 * (fh - f0) / (h * h * h) + (a1 + an) / (h * h)) * 6; // 6*c1

      f(3, 1) = f(3, 0);
      f(2, 1) = f(3, 0) * h + f(2, 0);

    } else if (iord1 == 1 && iord2 == 2) {
      T a1 = f(1, 0), bn = f(2, 1);
      f(1, 0) = a1;
      f(2, 1) = bn;
      h = dx;
      f0 = f(0, 0);
      fh = f(0, 1);

      f(2, 0) = (-bn / 4 + 3 * (fh - f0) / (2 * h * h) - 3 * a1 / (2 * h)) * 2;
      f(3, 0) =
        (bn / (4 * h) - (fh - f0) / (2 * h * h * h) + a1 / (2 * h * h)) * 6;

      f(3, 1) = f(3, 0);
      f(1, 1) = f(3, 0) * h * h / 2 + f(2, 0) * h + a1;

    } else if (iord1 == 2 && iord2 == 1) {
      T b1 = f(2, 0), an = f(1, 1);
      f(2, 0) = b1;
      f(1, 1) = an;
      h = dx;
      f0 = f(0, 0);
      fh = f(0, 1);

      f(1, 0) = 3 * (fh - f0) / (2 * h) - b1 * h / 4 - an / 2; // c3
      f(3, 0) =
        (an / (2 * h * h) - (fh - f0) / (2 * h * h * h) - b1 / (4 * h)) * 6;

      f(3, 1) = f(3, 0);
      f(2, 1) = f(3, 0) * h + f(2, 0);

    } else if (iord1 == 2 && iord2 == 2) {
      T b1 = f(2, 0), bn = f(2, 1);
      f(2, 0) = b1;
      f(2, 1) = bn;
      h = dx;
      f0 = f(0, 0);
      fh = f(0, 1);

      f(1, 0) = (fh - f0) / h - b1 * h / 3 - bn * h / 6;
      f(3, 0) = (bn - b1) / h;

      f(3, 1) = f(3, 0);
      f(1, 1) = f(3, 0) * h * h / 2 + b1 * h + f(1, 0);
    }
  } // n == 2
  else if (i3perio == 1) {
    h = 2 * dx;

    T dels = (f(0, 2) - f(0, 1)) / dx - (f(0, 1) - f(0, 0)) / dx;

    f(1, 0) = (f(0, 1) - f(0, 0)) / dx + (dx * dels) / h;
    f(2, 0) = -6.0 * dels / h;
    f(3, 0) = 12.0 * dels / (dx * h);

    f(1, 1) = (f(0, 2) - f(0, 1)) / dx - (dx * dels) / h;
    f(2, 1) = 6.0 * dels / h;
    f(3, 1) = -12.0 * dels / (dx * h);

    f(1, 2) = f(1, 0);
    f(2, 2) = f(2, 0);
    f(3, 2) = f(3, 0);
  } // i3perio == 1
  else if (i3knots == 2) {
    // Special case: nx = 3, not-a-knot on both sides
    h = 2 * dx;

    T f1 = f(0, 0) - f(0, 1);
    T f2 = f(0, 2) - f(0, 1);

    // Solve quadratic through 3 points centered at x[1]
    T aa = (f2 * dx + f1 * dx) / (dx * dx * h);
    T bb = (f2 * dx * dx - f1 * dx * dx) / (dx * dx * h);

    // Third derivative = 0 (quadratic)
    f(3, 0) = 0.0;
    f(3, 1) = 0.0;
    f(3, 2) = 0.0;

    // Second derivative is constant
    f(2, 0) = 2 * aa;
    f(2, 1) = 2 * aa;
    f(2, 2) = 2 * aa;

    // First derivatives
    f(1, 0) = bb - 2 * aa * dx;
    f(1, 1) = bb;
    f(1, 2) = bb + 2 * aa * dx;
  } // i3knots == 2
  else if (i3knots == 1) {
    if (i_bc1 == 1 || i_bc1 == 3 || i_bc1 == 5) {
      // f' LHS condition; not-a-knot RHS
      T h3 = 2 * dx;

      T f2 = f(0, 1) - f(0, 0);
      T f3 = f(0, 2) - f(0, 0);

      T aa = a1 / (dx * h3) + f3 / (h3 * h3 * (h3 - dx)) -
             f2 / (dx * dx * (h3 - dx));
      T bb = -a1 * (h3 * h3 - dx * dx) / (dx * h3 * (h3 - dx)) +
             f2 * h3 / (dx * dx * (h3 - dx)) - f3 * dx / (h3 * h3 * (h3 - dx));

      f(1, 0) = a1;
      f(2, 0) = 2 * bb;
      f(3, 0) = 6 * aa;

      f(1, 1) = 3 * aa * dx * dx + 2 * bb * dx + a1;
      f(2, 1) = 6 * aa * dx + 2 * bb;
      f(3, 1) = 6 * aa;

      f(1, 2) = 3 * aa * h3 * h3 + 2 * bb * h3 + a1;
      f(2, 2) = 6 * aa * h3 + 2 * bb;
      f(3, 2) = 6 * aa;

    } else if (i_bc1 == 2 || i_bc1 == 4 || i_bc1 == 6) {
      // f'' LHS condition; not-a-knot RHS
      T h3 = 2 * dx;

      T f2 = f(0, 1) - f(0, 0);
      T f3 = f(0, 2) - f(0, 0);

      T aa = -(b1 / 2.0) * (h3 - dx) / (h3 * h3 - dx * dx) -
             f2 / (dx * (h3 * h3 - dx * dx)) + f3 / (h3 * (h3 * h3 - dx * dx));
      T bb = -(b1 / 2.0) * dx * h3 * (h3 - dx) / (h3 * h3 - dx * dx) +
             f2 * h3 * h3 / (dx * (h3 * h3 - dx * dx)) -
             f3 * dx * dx / (h3 * (h3 * h3 - dx * dx));

      f(1, 0) = bb;
      f(2, 0) = b1;
      f(3, 0) = 6 * aa;

      f(1, 1) = 3 * aa * dx * dx + b1 * dx + bb;
      f(2, 1) = 6 * aa * dx + b1;
      f(3, 1) = 6 * aa;

      f(1, 2) = 3 * aa * h3 * h3 + b1 * h3 + bb;
      f(2, 2) = 6 * aa * h3 + b1;
      f(3, 2) = 6 * aa;

    } else if (i_bcn == 1 || i_bcn == 3 || i_bcn == 5) {
      // f' RHS condition; not-a-knot LHS
      T h2 = -dx;
      T h3 = -2 * dx;

      T f2 = f(0, 1) - f(0, 2);
      T f3 = f(0, 0) - f(0, 2);

      T aa = an / (h2 * h3) + f3 / (h3 * h3 * (h3 - h2)) -
             f2 / (h2 * h2 * (h3 - h2));
      T bb = -an * (h3 * h3 - h2 * h2) / (h2 * h3 * (h3 - h2)) +
             f2 * h3 / (h2 * h2 * (h3 - h2)) - f3 * h2 / (h3 * h3 * (h3 - h2));

      f(1, 2) = an;
      f(2, 2) = 2 * bb;
      f(3, 2) = 6 * aa;

      f(1, 1) = 3 * aa * h2 * h2 + 2 * bb * h2 + an;
      f(2, 1) = 6 * aa * h2 + 2 * bb;
      f(3, 1) = 6 * aa;

      f(1, 0) = 3 * aa * h3 * h3 + 2 * bb * h3 + an;
      f(2, 0) = 6 * aa * h3 + 2 * bb;
      f(3, 0) = 6 * aa;

    } else if (i_bcn == 2 || i_bcn == 4 || i_bcn == 6) {
      // f'' RHS condition; not-a-knot LHS
      T h2 = -dx;
      T h3 = -2 * dx;

      T f2 = f(0, 1) - f(0, 2);
      T f3 = f(0, 0) - f(0, 2);

      T aa = -(bn / 2.0) * (h3 - h2) / (h3 * h3 - h2 * h2) -
             f2 / (h2 * (h3 * h3 - h2 * h2)) + f3 / (h3 * (h3 * h3 - h2 * h2));
      T bb = -(bn / 2.0) * h2 * h3 * (h3 - h2) / (h3 * h3 - h2 * h2) +
             f2 * h3 * h3 / (h2 * (h3 * h3 - h2 * h2)) -
             f3 * h2 * h2 / (h3 * (h3 * h3 - h2 * h2));

      f(1, 2) = bb;
      f(2, 2) = bn;
      f(3, 2) = 6 * aa;

      f(1, 1) = 3 * aa * h2 * h2 + bn * h2 + bb;
      f(2, 1) = 6 * aa * h2 + bn;
      f(3, 1) = 6 * aa;

      f(1, 0) = 3 * aa * h3 * h3 + bn * h3 + bb;
      f(2, 0) = 6 * aa * h3 + bn;
      f(3, 0) = 6 * aa;
    }
  } // i3knots == 1
  else if (n > 2) {
    f(3, 0) = dx;
    f(2, 1) = (f(0, 1) - f(0, 0)) / f(3, 0);

    for (LO i = 1; i < n - 1; ++i) {
      f(3, i) = dx;
      f(1, i) = 2.0 * (f(3, i - 1) + f(3, i));
      f(2, i + 1) = (f(0, i + 1) - f(0, i)) / f(3, i);
      f(2, i) = f(2, i + 1) - f(2, i);
    }

    T elem21 = f(3, 0);
    T elemnn1 = f(3, n - 2);

    // Left boundary conditions
    if (i_bc1 == -1) {
      f(1, 0) = 2.0 * (f(3, 0) + f(3, n - 2));
      f(2, 0) = (f(0, 1) - f(0, 0)) / f(3, 0) -
                (f(0, n - 1) - f(0, n - 2)) / f(3, n - 2);
      wk[0] = f(3, n - 2);
      for (LO i = 1; i <= n - 4; ++i)
        wk[i] = 0.0;
      wk[n - 3] = f(3, n - 3);
      wk[n - 2] = f(3, n - 2);
    } else if (i_bc1 == 1 || i_bc1 == 3 || i_bc1 == 5) {
      f(1, 0) = 2.0 * f(3, 0);
      f(2, 0) = (f(0, 1) - f(0, 0)) / f(3, 0) - a1;
    } else if (i_bc1 == 2 || i_bc1 == 4 || i_bc1 == 6) {
      f(1, 0) = 2.0 * f(3, 0);
      f(2, 0) = f(3, 0) * b1 / 3.0;
      f(3, 0) = 0.0;
    } else if (i_bc1 == 7) {
      f(1, 0) = -f(3, 0);
      f(2, 0) = f(2, 2) / (2 * dx) - f(2, 1) / (2 * dx);
      f(2, 0) *= f(3, 0) * f(3, 0) / (3 * dx);
    } else {
      imin = 1;
      f(1, 1) = f(3, 0) + 2.0 * f(3, 1);
      f(2, 1) = f(2, 1) * f(3, 1) / (f(3, 0) + f(3, 1));
    }

    // Right boundary conditions
    if (i_bcn == 1 || i_bcn == 3 || i_bcn == 5) {
      f(1, n - 1) = 2.0 * f(3, n - 2);
      f(2, n - 1) = -(f(0, n - 1) - f(0, n - 2)) / f(3, n - 2) + an;
    } else if (i_bcn == 2 || i_bcn == 4 || i_bcn == 6) {
      f(1, n - 1) = 2.0 * f(3, n - 2);
      f(2, n - 1) = f(3, n - 2) * bn / 3.0;
      elemnn1 = 0.0;
    } else if (i_bcn == 7) {
      f(1, n - 1) = -f(3, n - 2);
      f(2, n - 1) = f(2, n - 2) / (2 * dx) - f(2, n - 3) / (2 * dx);
      f(2, n - 1) = -f(2, n - 1) * f(3, n - 2) * f(3, n - 2) / (3 * dx);
    } else if (i_bc1 != -1) {
      imax = n - 2;
      f(1, n - 2) = 2.0 * f(3, n - 3) + f(3, n - 2);
      f(2, n - 2) = f(2, n - 2) * f(3, n - 3) / (f(3, n - 2) + f(3, n - 3));
    }
    if (i_bc1 == -1) {
      // Periodic Boundary Conditions

      // Forward elimination
      for (LO i = 1; i <= n - 3; ++i) {
        T t = f(3, i - 1) / f(1, i - 1);
        f(1, i) -= t * f(3, i - 1);
        f(2, i) -= t * f(2, i - 1);
        wk[i] -= t * wk[i - 1];

        T q = wk[n - 2] / f(1, i - 1);
        wk[n - 2] = -q * f(3, i - 1);
        f(1, n - 2) -= q * wk[i - 1];
        f(2, n - 2) -= q * f(2, i - 1);
      }

      wk[n - 2] += f(3, n - 3);

      // Complete forward elimination
      T t = wk[n - 2] / f(1, n - 3);
      f(1, n - 2) -= t * wk[n - 3];
      f(2, n - 2) -= t * f(2, n - 3);

      // Back substitution
      f(2, n - 2) /= f(1, n - 2);
      f(2, n - 3) = (f(2, n - 3) - wk[n - 3] * f(2, n - 2)) / f(1, n - 3);
      for (LO ib = 2; ib <= n - 2; ++ib) {
        LO i = n - ib - 2;
        f(2, i) =
          (f(2, i) - f(3, i) * f(2, i + 1) - wk[i] * f(2, n - 2)) / f(1, i);
      }

      f(2, n - 1) = f(2, 0);

    } else {
      // Non-periodic Boundary Conditions

      for (LO i = imin + 1; i <= imax; ++i) {
        T t;
        if ((i == n - 2) && (imax == n - 2)) {
          t = (f(3, i - 1) - f(3, i)) / f(1, i - 1);
        } else if (i == 1) {
          t = elem21 / f(1, i - 1);
        } else if (i == n - 1) {
          t = elemnn1 / f(1, i - 1);
        } else {
          t = f(3, i - 1) / f(1, i - 1);
        }

        if ((i == imin + 1) && (imin == 1)) {
          f(1, i) -= t * (f(3, i - 1) - f(3, i - 2));
        } else {
          f(1, i) -= t * f(3, i - 1);
        }

        f(2, i) -= t * f(2, i - 1);
      }

      f(2, imax) /= f(1, imax);
      for (LO ib = 0; ib < imax - imin; ++ib) {
        LO i = imax - 1 - ib;
        if ((i == 1) && (imin == 1)) {
          f(2, i) = (f(2, i) - (f(3, i) - f(3, i - 1)) * f(2, i + 1)) / f(1, i);
        } else {
          f(2, i) = (f(2, i) - f(3, i) * f(2, i + 1)) / f(1, i);
        }
      }

      f(3, 0) = dx;
      f(3, n - 2) = dx;

      if (i_bc1 <= 0 || i_bc1 > 7) {
        f(2, 0) = (f(2, 1) * (f(3, 0) + f(3, 1)) - f(2, 2) * f(3, 0)) / f(3, 1);
      }

      if (i_bcn <= 0 || i_bcn > 7) {
        f(2, n - 1) =
          f(2, n - 2) + (f(2, n - 2) - f(2, n - 3)) * f(3, n - 2) / f(3, n - 3);
      }
    }

    // Polynomial coefficient computation
    for (LO i = 0; i < n - 1; ++i) {
      f(1, i) = (f(0, i + 1) - f(0, i)) / f(3, i) -
                f(3, i) * (f(2, i + 1) + 2.0 * f(2, i));
      f(3, i) = (f(2, i + 1) - f(2, i)) / f(3, i);
      f(2, i) *= 6.0;
      f(3, i) *= 6.0;
    }

    if (i_bc1 == -1) {
      f(1, n - 1) = f(1, 0);
      f(2, n - 1) = f(2, 0);
      f(3, n - 1) = f(3, 0);
    } else {
      f(1, n - 1) = f(1, n - 2) + dx * (f(2, n - 2) + 0.5 * dx * f(3, n - 2));
      f(2, n - 1) = f(2, n - 2) + dx * f(3, n - 2);
      f(3, n - 1) = f(3, n - 2);

      if (i_bcn == 1 || i_bcn == 3 || i_bcn == 5) {
        f(1, n - 1) = an;
      } else if (i_bcn == 2 || i_bcn == 4 || i_bcn == 6) {
        f(2, n - 1) = bn;
      }
    }
  }
} // compute_spline_coefficients

template <typename T, typename MemorySpace>
KOKKOS_FUNCTION void solve_cubic_spline_explicit(
  Rank1View<T, MemorySpace> x, LO nx, Rank2View<T, MemorySpace> fspl,
  LO ibcxmin, T bcxmin, LO ibcxmax, T bcxmax, Rank1View<T, MemorySpace> wk)
{

  const T half = 1.0 / 2.0;
  const T sixth = 1.0 / 6.0;

  if (ibcxmin == 1) {
    fspl(1, 0) = bcxmin;
  } else if (ibcxmin == 2) {
    fspl(2, 0) = bcxmin;
  }

  if (ibcxmax == 1) {
    fspl(1, nx - 1) = bcxmax;
  } else if (ibcxmax == 2) {
    fspl(2, nx - 1) = bcxmax;
  }
  compute_spline_coefficients<T, MemorySpace>(ibcxmin, ibcxmax, nx, x, fspl,
                                              wk);
  for (LO i = 0; i < nx; ++i) {
    fspl(2, i) = half * fspl(2, i);
    fspl(3, i) = sixth * fspl(3, i);
  }
}

/** \brief Evaluate the point given location in the grid
 *
 * \param selector selector array indicating which derivatives to compute
 * \param fval output values of the spline at given points
 * \param i index of the grid point
 * \param dx displacement from the grid point
 * \param fspl spline coefficients
 */
template <typename T, typename MemorySpace>
KOKKOS_INLINE_FUNCTION void cubic_evalfn_explicit(
  Kokkos::View<LO*, MemorySpace> selector, Rank1View<T, MemorySpace> fval, LO i,
  T dx, Rank2View<T, MemorySpace> fspl)
{

  LO iaval = 0;

  if (selector[0] == 3) {
    // Third derivative only
    iaval++;
    fval(iaval - 1) = 6.0 * fspl(3, i);
  } else {
    if (selector[0] > 0) {
      // Evaluate f
      iaval++;
      fval(iaval - 1) =
        fspl(0, i) + dx * (fspl(1, i) + dx * (fspl(2, i) + dx * fspl(3, i)));
    }

    if (selector[1] > 0) {
      // Evaluate df/dx
      iaval++;
      fval(iaval - 1) =
        fspl(1, i) + dx * (2.0 * fspl(2, i) + dx * 3.0 * fspl(3, i));
    }

    if (selector[2] > 0) {
      // Evaluate d2f/dx2
      iaval++;
      fval(iaval - 1) = 2.0 * fspl(2, i) + dx * 6.0 * fspl(3, i);
    }
  }
} // end evalfn

template <typename T, typename MemorySpace>
KOKKOS_INLINE_FUNCTION LO lookup_in_grid(T zxget, Rank1View<T, MemorySpace> x,
                                         LO nx)
{
  LO nxm = nx - 1;
  LO ii = static_cast<LO>(nxm * (zxget - x[0]) / (x[nx - 1] - x[0]));
  LO i = std::min(nxm - 1, ii);
  if (zxget < x[i]) {
    i = std::max(0, i - 1);
  } else if (zxget > x[i + 1]) {
    i = std::min(nxm - 1, i + 1);
  }
  return i;
}

/** \brief Lookup the grid point and compute the displacement
 *
 * \param xget input point where the spline is evaluated
 * \param x input grid points
 * \param nx number of grid points
 * \param i index of the grid point
 * \param dx displacement from the grid point
 * \param ier error flag, 0 if successful, non-zero otherwise
 * \param lookup_tol tolerance for out-of-bound checking
 */
template <typename T, typename MemorySpace>
KOKKOS_INLINE_FUNCTION void cubic_lookup_explicit(T xget,
                                                  Rank1View<T, MemorySpace> x,
                                                  LO nx, LO& i, T& dx, LO& ier,
                                                  T lookup_tol)
{
  T zxget = xget;

  // Range check
  clamp_to_grid_range(xget, zxget, x, nx, ier, lookup_tol);

  i = lookup_in_grid<T, MemorySpace>(zxget, x, nx);

  dx = zxget - x[i];
}

/** \brief Evaluate the cubic spline at a single point with selector
 *
 * \param xget input point where the spline is evaluated
 * \param selector selector array indicating which derivatives to compute
 * \param fval output values of the spline at given point
 * \param x input grid points
 * \param nx number of grid points
 * \param fspl spline coefficients
 * \param ier error flag, 0 if successful, non-zero otherwise
 * \param lookup_tol tolerance for out-of-bound checking
 *
 * \details A static function that can be used to evaluate the cubic spline at a
 * single point. Note that this function does nothing when tje input point is
 * out of the grid range (ier is set to non-zero).
 */
template <typename T, typename MemorySpace>
KOKKOS_INLINE_FUNCTION void cubic_eval_explicit(
  T xget, Kokkos::View<LO*, MemorySpace> iselect,
  Rank1View<T, MemorySpace> fval, Rank1View<T, MemorySpace> x, LO nx,
  Rank2View<T, MemorySpace> fspl, LO& ier, T lookup_tol)
{
  LO ia = 0;
  T dxa = 0.0;
  cubic_lookup_explicit(xget, x, nx, ia, dxa, ier, lookup_tol);
  if (ier != 0) {
    printf("error in evaluation, ier = %d\n", ier);
    return;
  }

  cubic_evalfn_explicit(iselect, fval, ia, dxa, fspl);
}

// TODO: take member type as parameter for further optimization
template <typename T, typename MemorySpace>
KOKKOS_FUNCTION void solve_cubic_spline_compact(
  Rank1View<T, MemorySpace> x, LO nx, Rank2View<T, MemorySpace> fspl,
  Rank2View<T, MemorySpace> fspl4, LO ibcxmin, T bcxmin, LO ibcxmax, T bcxmax,
  Rank1View<T, MemorySpace> wk)
{
  // Copy f data to fspl4 and zero out second derivative output
  for (LO i = 0; i < nx; ++i) {
    fspl4(0, i) = fspl(0, i);
    fspl(1, i) = 0.0;
  }

  // Call traditional spline generator
  solve_cubic_spline_explicit(x, nx, fspl4, ibcxmin, bcxmin, ibcxmax, bcxmax,
                              wk);

  for (LO i = 0; i < nx - 1; ++i) {
    fspl(1, i) = 2.0 * fspl4(2, i);
  }

  fspl(1, nx - 1) =
    2.0 * fspl4(2, nx - 2) + (x[nx - 1] - x[nx - 2]) * 6.0 * fspl4(3, nx - 2);
}

/** \brief Evaluate the point given location in the grid
 *
 * \param selector selector array indicating which derivatives to compute
 * \param fval output values of the spline at given points
 * \param i index of the grid point
 * \param xparam displacement from the grid point
 * \param hx step size in the x direction
 * \param hxi inverse step size in the x direction
 * \param fs2 compact spline coefficients
 */
template <typename T, typename MemorySpace>
KOKKOS_INLINE_FUNCTION void cubic_evalfn_compact(
  Kokkos::View<LO*, MemorySpace> selector, Rank1View<T, MemorySpace> fval, LO i,
  T xparam, T hx, T hxi, Rank2View<T, MemorySpace> fs2)
{
  const T sixth = 1.0 / 6.0;
  LO iadr = 0;

  if (selector[0] <= 2) {
    T xp = xparam;
    T xpi = 1.0 - xp;
    if (selector[0] == 1) {
      // Function value f(x)
      ++iadr;
      // TODO: too much temp varibale created?
      T xp2 = xp * xp;
      T xpi2 = xpi * xpi;
      T cx = xp * (xp2 - 1.0);
      T cxi = xpi * (xpi2 - 1.0);
      T hx2 = hx * hx;

      T sum = xpi * fs2(0, i) + xp * fs2(0, i + 1);
      sum += sixth * hx2 * (cxi * fs2(1, i) + cx * fs2(1, i + 1));

      fval(iadr - 1) = sum; // zero-based indexing
    }

    if (selector[1] == 1) {
      // First derivative df/dx
      ++iadr;
      T xp2 = xp * xp;
      T xpi2 = xpi * xpi;

      T cxd = 3.0 * xp2 - 1.0;
      T cxdi = -3.0 * xpi2 + 1.0;

      T sum = hxi * (fs2(0, i + 1) - fs2(0, i));
      sum += sixth * hx * (cxdi * fs2(1, i) + cxd * fs2(1, i + 1));

      fval(iadr - 1) = sum;
    }

    if (selector[2] == 1) {
      // Second derivative d2f/dx2
      ++iadr;

      T sum = xpi * fs2(1, i) + xp * fs2(1, i + 1);
      fval(iadr - 1) = sum;
    }

  } else {
    // Third derivative d3f/dx3
    iadr = 1;
    fval(iadr - 1) = hxi * (fs2(1, i + 1) - fs2(1, i));
  }
}

/** \brief Lookup the grid point and compute the displacement
 *
 * \param xget input point where the spline is evaluated
 * \param x input grid points
 * \param nx number of grid points
 * \param i index of the grid point
 * \param xparam displacement from the grid point
 * \param hx step size in the x direction
 * \param hxi inverse step size in the x direction
 * \param ier error flag, 0 if successful, non-zero otherwise
 * \param lookup_tol tolerance for looking up the grid point
 */
template <typename T, typename MemorySpace>
KOKKOS_INLINE_FUNCTION void cubic_lookup_compact(T xget,
                                                 Rank1View<T, MemorySpace> x,
                                                 LO nx, LO& i, T& xparam, T& hx,
                                                 T& hxi, LO& ier, T lookup_tol)
{
  T zxget = xget;

  clamp_to_grid_range(xget, zxget, x, nx, ier, lookup_tol);

  i = lookup_in_grid<T, MemorySpace>(zxget, x, nx);

  hx = x[i + 1] - x[i];
  hxi = 1.0 / hx;
  xparam = (zxget - x[i]) * hxi;
}

/** \brief Evaluate the compact cubic spline at a single point with selector
 *
 * \param xget input point where the spline is evaluated
 * \param ict selector array indicating which derivatives to compute
 * \param fval output values of the spline at given point
 * \param x input grid points
 * \param nx number of grid points
 * \param fs2 compact spline coefficients
 * \param ier error flag, 0 if successful, non-zero otherwise
 * \param lookup_tol tolerance for out-of-bound checking
 */
template <typename T, typename MemorySpace>
KOKKOS_INLINE_FUNCTION void cubic_eval_compact(
  T xget, Kokkos::View<LO*, MemorySpace> ict, Rank1View<T, MemorySpace> fval,
  Rank1View<T, MemorySpace> x, LO nx, Rank2View<T, MemorySpace> fspl, LO& ier,
  T lookup_tol)
{

  // Initialize output zone info
  LO i = 0;
  T xparam = 0.0;
  T hx = 0.0;
  T hxi = 0.0;

  // Find the interval containing xget
  cubic_lookup_compact(xget, x, nx, i, xparam, hx, hxi, ier, lookup_tol);
  if (ier != 0)
    return;

  // Evaluate spline at the points
  cubic_evalfn_compact(ict, fval, i, xparam, hx, hxi, fspl);
}

/** \brief Lookup the grid point and compute the displacement
 *
 * \param xget input point in x direction where the spline is evaluated
 * \param yget input point in y direction where the spline is evaluated
 * \param x input grid points in x direction
 * \param nx number of grid points in x direction
 * \param y input grid points in y direction
 * \param ny number of grid points in y direction
 * \param i index of the grid point in x direction
 * \param j index of the grid point in y direction
 * \param dx displacement from the grid point in x direction
 * \param dy displacement from the grid point in y direction
 * \param ier error flag, 0 if successful, non-zero otherwise
 * \param lookup_tol tolerance for looking up the grid point
 */
template <typename T, typename MemorySpace>
KOKKOS_INLINE_FUNCTION void bicubic_lookup_explicit(
  T xget, T yget, Rank1View<T, MemorySpace> x, LO nx,
  Rank1View<T, MemorySpace> y, LO ny, LO& i, LO& j, T& dx, T& dy, LO& ier,
  T lookup_tol)
{
  T zxget = xget;
  T zyget = yget;

  clamp_to_grid_range(xget, zxget, x, nx, ier, lookup_tol);
  clamp_to_grid_range(yget, zyget, y, ny, ier, lookup_tol);

  if (ier != 0)
    return;

  i = lookup_in_grid<T, MemorySpace>(zxget, x, nx);
  j = lookup_in_grid<T, MemorySpace>(zyget, y, ny);

  dx = zxget - x[i];
  dy = zyget - y[j];
}

/** \brief Evaluate the bicubic spline at given grid cell
 *
 * \param ict selector array indicating which derivatives to compute
 * \param fval output values of the spline at given grid cell
 * \param i index of the grid cell in x direction
 * \param j index of the grid cell in y direction
 * \param dx displacement from the grid cell in x direction
 * \param dy displacement from the grid cell in y direction
 * \param fspl bicubic spline coefficients
 */
template <typename T, typename MemorySpace>
KOKKOS_INLINE_FUNCTION void bicubic_evalfn_explicit(
  Kokkos::View<LO*, MemorySpace> ict, Rank1View<T, MemorySpace> fval, LO i,
  LO j, T dx, T dy, Rank4View<T, MemorySpace> fspl)
{
  LO iaval = 0; // Index for fval
  if (ict[0] <= 2) {
    if ((ict[0] > 0) || (ict[0] == -1)) {
      // Evaluate f
      iaval++;
      fval(iaval - 1) =
        fspl(0, 0, i, j) +
        dy *
          (fspl(0, 1, i, j) + dy * (fspl(0, 2, i, j) + dy * fspl(0, 3, i, j))) +
        dx *
          (fspl(1, 0, i, j) +
           dy * (fspl(1, 1, i, j) +
                 dy * (fspl(1, 2, i, j) + dy * fspl(1, 3, i, j))) +
           dx *
             (fspl(2, 0, i, j) +
              dy * (fspl(2, 1, i, j) +
                    dy * (fspl(2, 2, i, j) + dy * fspl(2, 3, i, j))) +
              dx * (fspl(3, 0, i, j) +
                    dy * (fspl(3, 1, i, j) +
                          dy * (fspl(3, 2, i, j) + dy * fspl(3, 3, i, j))))));
    }

    if ((ict[1] > 0) && (ict[0] != -1)) {
      // Evaluate df/dx
      iaval++;
      fval(iaval - 1) =
        fspl(1, 0, i, j) +
        dy *
          (fspl(1, 1, i, j) + dy * (fspl(1, 2, i, j) + dy * fspl(1, 3, i, j))) +
        2.0 * dx *
          (fspl(2, 0, i, j) +
           dy * (fspl(2, 1, i, j) +
                 dy * (fspl(2, 2, i, j) + dy * fspl(2, 3, i, j))) +
           1.5 * dx *
             (fspl(3, 0, i, j) +
              dy * (fspl(3, 1, i, j) +
                    dy * (fspl(3, 2, i, j) + dy * fspl(3, 3, i, j)))));
    }

    if ((ict[2] > 0) && (ict[0] != -1)) {
      // Evaluate df/dy
      iaval++;
      fval(iaval - 1) =
        fspl(0, 1, i, j) +
        dy * (2.0 * fspl(0, 2, i, j) + dy * 3.0 * fspl(0, 3, i, j)) +
        dx *
          (fspl(1, 1, i, j) +
           dy * (2.0 * fspl(1, 2, i, j) + dy * 3.0 * fspl(1, 3, i, j)) +
           dx * (fspl(2, 1, i, j) +
                 dy * (2.0 * fspl(2, 2, i, j) + dy * 3.0 * fspl(2, 3, i, j)) +
                 dx * (fspl(3, 1, i, j) + dy * (2.0 * fspl(3, 2, i, j) +
                                                dy * 3.0 * fspl(3, 3, i, j)))));
    }

    if ((ict[3] > 0) || (ict[0] == -1)) {
      // Evaluate d2f/dx2
      iaval++;
      fval(iaval - 1) =
        2.0 * (fspl(2, 0, i, j) +
               dy * (fspl(2, 1, i, j) +
                     dy * (fspl(2, 2, i, j) + dy * fspl(2, 3, i, j)))) +
        6.0 * dx *
          (fspl(3, 0, i, j) +
           dy * (fspl(3, 1, i, j) +
                 dy * (fspl(3, 2, i, j) + dy * fspl(3, 3, i, j))));
    }

    if ((ict[4] > 0) || (ict[0] == -1)) {
      // Evaluate d2f/dy2
      iaval++;
      fval(iaval - 1) =
        2.0 * fspl(0, 2, i, j) + 6.0 * dy * fspl(0, 3, i, j) +
        dx *
          (2.0 * fspl(1, 2, i, j) + 6.0 * dy * fspl(1, 3, i, j) +
           dx * (2.0 * fspl(2, 2, i, j) + 6.0 * dy * fspl(2, 3, i, j) +
                 dx * (2.0 * fspl(3, 2, i, j) + 6.0 * dy * fspl(3, 3, i, j))));
    }

    if ((ict[5] > 0) && (ict[0] != -1)) {
      // Evaluate d2f/dxdy
      iaval++;
      fval(iaval - 1) =
        fspl(1, 1, i, j) +
        dy * (2.0 * fspl(1, 2, i, j) + dy * 3.0 * fspl(1, 3, i, j)) +
        2.0 * dx *
          (fspl(2, 1, i, j) +
           dy * (2.0 * fspl(2, 2, i, j) + dy * 3.0 * fspl(2, 3, i, j)) +
           1.5 * dx *
             (fspl(3, 1, i, j) +
              dy * (2.0 * fspl(3, 2, i, j) + dy * 3.0 * fspl(3, 3, i, j))));
    }

    if (ict[0] == -1) {
      // Evaluate d4f/dx2dy2
      iaval++;
      fval(iaval - 1) =
        4.0 * fspl(2, 2, i, j) + 12.0 * dy * fspl(2, 3, i, j) +
        dx * (12.0 * fspl(3, 2, i, j) + 36.0 * dy * fspl(3, 3, i, j));
    }
  } else if (ict[0] == 3) {
    if (ict[1] == 1) {
      // d³f/dx³ (not continuous)
      iaval++;
      fval(iaval - 1) =
        6.0 * (fspl(3, 0, i, j) +
               dy * (fspl(3, 1, i, j) +
                     dy * (fspl(3, 2, i, j) + dy * fspl(3, 3, i, j))));
    }

    if (ict[2] == 1) {
      // d³f/dx²dy
      iaval++;
      fval(iaval - 1) =
        2.0 * (fspl(2, 1, i, j) +
               dy * (2.0 * fspl(2, 2, i, j) + dy * 3.0 * fspl(2, 3, i, j))) +
        6.0 * dx *
          (fspl(3, 1, i, j) +
           dy * (2.0 * fspl(3, 2, i, j) + dy * 3.0 * fspl(3, 3, i, j)));
    }

    if (ict[3] == 1) {
      // d³f/dxdy²
      iaval++;
      fval(iaval - 1) =
        2.0 * fspl(1, 2, i, j) + 6.0 * dy * fspl(1, 3, i, j) +
        2.0 * dx *
          (2.0 * fspl(2, 2, i, j) + 6.0 * dy * fspl(2, 3, i, j) +
           1.5 * dx * (2.0 * fspl(3, 2, i, j) + 6.0 * dy * fspl(3, 3, i, j)));
    }

    if (ict[4] == 1) {
      // d³f/dy³ (not continuous)
      iaval++;
      fval(iaval - 1) =
        6.0 * (fspl(0, 3, i, j) +
               dx * (fspl(1, 3, i, j) +
                     dx * (fspl(2, 3, i, j) + dx * fspl(3, 3, i, j))));
    }

  } else if (ict[0] == 4) {
    if (ict[1] == 1) {
      // d⁴f/dx³dy
      iaval++;
      fval(iaval - 1) =
        6.0 * (fspl(3, 1, i, j) +
               dy * 2.0 * (fspl(3, 2, i, j) + dy * 1.5 * fspl(3, 3, i, j)));
    }

    if (ict[2] == 1) {
      // d⁴f/dx²dy²
      iaval++;
      fval(iaval - 1) =
        4.0 * fspl(2, 2, i, j) + 12.0 * dy * fspl(2, 3, i, j) +
        dx * (12.0 * fspl(3, 2, i, j) + 36.0 * dy * fspl(3, 3, i, j));
    }

    if (ict[3] == 1) {
      // d⁴f/dxdy³ (not continuous)
      iaval++;
      fval(iaval - 1) =
        6.0 * (fspl(1, 3, i, j) +
               2.0 * dx * (fspl(2, 3, i, j) + 1.5 * dx * fspl(3, 3, i, j)));
    }

  } else if (ict[0] == 5) {
    if (ict[1] == 1) {
      // d⁵f/dx³dy² (not continuous)
      iaval++;
      fval(iaval - 1) = 12.0 * (fspl(3, 2, i, j) + dy * 3.0 * fspl(3, 3, i, j));
    }

    if (ict[2] == 1) {
      // d⁵f/dx²dy³ (not continuous)
      iaval++;
      fval(iaval - 1) = 12.0 * (fspl(2, 3, i, j) + dx * 3.0 * fspl(3, 3, i, j));
    }

  } else if (ict[0] == 6) {
    // d⁶f/dx³dy³ (not continuous)
    iaval++;
    fval(iaval - 1) = 36.0 * fspl(3, 3, i, j);
  }
} // end evalfn

/** \brief Evaluate the bicubic spline at a single point with selector
 *
 * \param xget input point in x direction where the spline is evaluated
 * \param yget input point in y direction where the spline is evaluated
 * \param iselect selector array indicating which derivatives to compute
 * \param fval output values of the spline at given point
 * \param x input grid points in x direction
 * \param nx number of grid points in x direction
 * \param y input grid points in y direction
 * \param ny number of grid points in y direction
 * \param fspl bicubic spline coefficients
 * \param ier error flag, 0 if successful, non-zero otherwise
 * \param lookup_tol tolerance for looking up the grid point
 */
template <typename T, typename MemorySpace>
KOKKOS_INLINE_FUNCTION void bicubic_eval_explicit(
  T xget, T yget, Kokkos::View<LO*, MemorySpace> iselect,
  Rank1View<T, MemorySpace> fval, Rank1View<T, MemorySpace> x, LO nx,
  Rank1View<T, MemorySpace> y, LO ny, Rank4View<T, MemorySpace> fspl, LO& ier,
  T lookup_tol)
{
  LO i = 0;
  LO j = 0;
  T dx = 0.0;
  T dy = 0.0;

  // Range finding
  bicubic_lookup_explicit(xget, yget, x, nx, y, ny, i, j, dx, dy, ier,
                          lookup_tol);
  if (ier != 0)
    return;

  // Evaluate spline function
  bicubic_evalfn_explicit(iselect, fval, i, j, dx, dy, fspl);
}

template <typename T, typename MemorySpace>
void solve_bicubic_spline_explicit(Rank1View<T, MemorySpace> x, LO inx,
                                   Rank1View<T, MemorySpace> y, LO iny,
                                   Rank4View<T, MemorySpace> fspl, LO ibcxmin,
                                   Rank1View<T, MemorySpace> bcxmin, LO ibcxmax,
                                   Rank1View<T, MemorySpace> bcxmax, LO ibcymin,
                                   Rank1View<T, MemorySpace> bcymin, LO ibcymax,
                                   Rank1View<T, MemorySpace> bcymax,
                                   Rank1View<T, MemorySpace> wk)
{
  using execution_space = typename MemorySpace::execution_space;
  using member_type = typename Kokkos::TeamPolicy<execution_space>::member_type;

  LO iflg2 = 0;
  if (ibcymin != -1) {
    detect_inhomogeneous_bc(bcymin, ibcymin, iflg2, inx);
    detect_inhomogeneous_bc(bcymax, ibcymax, iflg2, inx);
  }

  const T xo2 = 0.5;
  const T xo6 = 1.0 / 6.0;

  auto fspl_l_x =
    Rank2View<T, MemorySpace>(wk.data_handle() + 4 * inx * iny, 4 * inx, iny);
  auto wk_l =
    Rank1View<T, MemorySpace>(wk.data_handle() + 2 * 4 * inx * iny, inx * iny);
  auto fspl_l_y = Rank2View<T, MemorySpace>(wk.data_handle(), 4 * inx, iny);

  Kokkos::parallel_for(
    Kokkos::TeamPolicy<execution_space>(iny, Kokkos::AUTO),
    KOKKOS_LAMBDA(const member_type& team) {
      const LO iy = team.league_rank();

      auto fspl_x_view = Rank2View<T, MemorySpace>(
        fspl_l_x.data_handle() + 4 * inx * iy, 4, inx);
      auto wk_x_view =
        Rank1View<T, MemorySpace>(wk_l.data_handle() + iy * inx, inx);

      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, inx), [=](LO ix) {
        fspl_x_view(0, ix) = fspl(0, 0, ix, iy);
      });

      if (team.team_rank() == 0) {
        if (ibcxmin == 1)
          fspl_x_view(1, 0) = bcxmin[iy];
        else if (ibcxmin == 2)
          fspl_x_view(2, 0) = bcxmin[iy];

        if (ibcxmax == 1)
          fspl_x_view(1, inx - 1) = bcxmax[iy];
        else if (ibcxmax == 2)
          fspl_x_view(2, inx - 1) = bcxmax[iy];

        compute_spline_coefficients<T, MemorySpace>(ibcxmin, ibcxmax, inx, x,
                                                    fspl_x_view, wk_x_view);
      }

      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, inx), [=](LO ix) {
        fspl(1, 0, ix, iy) = fspl_x_view(1, ix);
        fspl(2, 0, ix, iy) = fspl_x_view(2, ix) * xo2;
        fspl(3, 0, ix, iy) = fspl_x_view(3, ix) * xo6;
      });
    });

  Kokkos::parallel_for(
    Kokkos::TeamPolicy<execution_space>(inx, Kokkos::AUTO),
    KOKKOS_LAMBDA(const member_type& team) {
      const LO ix = team.league_rank();

      auto fspl_y_view = Rank2View<T, MemorySpace>(
        fspl_l_y.data_handle() + 4 * ix * iny, 4, iny);
      auto wk_y_view =
        Rank1View<T, MemorySpace>(wk_l.data_handle() + ix * iny, iny);

      for (LO ic = 0; ic < 4; ++ic) {
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, iny), [=](LO iy) {
          fspl_y_view(0, iy) = fspl(ic, 0, ix, iy);
        });

        // Set linear BCs initially
        if (team.team_rank() == 0) {
          fspl_y_view(1, 0) = 0.0;
          fspl_y_view(2, 0) = 0.0;
          fspl_y_view(1, iny - 1) = 0.0;
          fspl_y_view(2, iny - 1) = 0.0;

          LO ibcymina = ibcymin;
          LO ibcymaxa = ibcymax;
          if (iflg2 == 1) {
            if (ibcymin == 1 || ibcymin == 2)
              ibcymina = 0;
            if (ibcymax == 1 || ibcymax == 2)
              ibcymaxa = 0;
          }

          compute_spline_coefficients<T, MemorySpace>(
            ibcymina, ibcymaxa, iny, y, fspl_y_view, wk_y_view);
        }

        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, iny), [=](LO iy) {
          fspl(ic, 1, ix, iy) = fspl_y_view(1, iy);
          fspl(ic, 2, ix, iy) = fspl_y_view(2, iy) * xo2;
          fspl(ic, 3, ix, iy) = fspl_y_view(3, iy) * xo6;
        });
      }
    });

  if (iflg2 == 1) {
    LO iasc = 0;       // Workspace base for correction splines
    LO iinc = 4 * iny; // Spacing between correction splines

    T zhxn =
      get_element_from_span(x, inx - 1) - get_element_from_span(x, inx - 2);
    T zhy =
      get_element_from_span(y, iny - 1) - get_element_from_span(y, iny - 2);

    LO jx = inx - 2;
    LO jy = iny - 2;

    LO iselect1_arr[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    LO iselect2_arr[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    if (ibcymin == 1)
      iselect1_arr[2] = 1;
    if (ibcymin == 2)
      iselect1_arr[4] = 1;
    if (ibcymax == 1)
      iselect2_arr[2] = 1;
    if (ibcymax == 2)
      iselect2_arr[4] = 1;
    Kokkos::View<LO*, MemorySpace> iselect1("iselect1", 10);
    Kokkos::parallel_for(
      Kokkos::RangePolicy<execution_space>(0, iselect1.size()),
      KOKKOS_LAMBDA(const LO i) { iselect1(i) = iselect1_arr[i]; });
    Kokkos::View<LO*, MemorySpace> iselect2("iselect2", 10);
    Kokkos::parallel_for(
      Kokkos::RangePolicy<execution_space>(0, iselect2.size()),
      KOKKOS_LAMBDA(const LO i) { iselect2(i) = iselect2_arr[i]; });

    auto zcur_shared =
      Rank1View<T, MemorySpace>(fspl_l_x.data_handle(), inx * 3);

    Kokkos::parallel_for(
      Kokkos::RangePolicy<execution_space>(0, inx), KOKKOS_LAMBDA(const LO ix) {
        LO ier = 0;
        auto zcur =
          Rank1View<T, MemorySpace>(zcur_shared.data_handle() + ix * 3, 3);
        T zdiff1 = 0.0, zdiff2 = 0.0;
        auto wk_y = Rank1View<T, MemorySpace>(wk_l.data_handle(), iny);

        if (ibcymin == 1) {
          zcur[0] =
            (ix < inx - 1)
              ? fspl(0, 1, ix, 0)
              : fspl(0, 1, jx, 0) + zhxn * (fspl(1, 1, jx, 0) +
                                            zhxn * (fspl(2, 1, jx, 0) +
                                                    zhxn * fspl(3, 1, jx, 0)));
          zdiff1 = bcymin[ix] - zcur[0];
        } else if (ibcymin == 2) {
          zcur[0] = (ix < inx - 1)
                      ? 2.0 * fspl(0, 2, ix, 0)
                      : 2.0 * (fspl(0, 2, jx, 0) +
                               zhxn * (fspl(1, 2, jx, 0) +
                                       zhxn * (fspl(2, 2, jx, 0) +
                                               zhxn * fspl(3, 2, jx, 0))));
          zdiff1 = bcymin[ix] - zcur[0];
        }

        // auto zcur = Rank1View<T, MemorySpace>(zcur_vec.data(), 3);
        if (ibcymax == 1) {
          if (ix < inx - 1) {
            zcur(0) =
              fspl(0, 1, ix, jy) +
              zhy * (2.0 * fspl(0, 2, ix, jy) + zhy * 3.0 * fspl(0, 3, ix, jy));
          } else {
            bicubic_eval_explicit(x[inx - 1], y[iny - 1], iselect2, zcur, x,
                                  inx, y, iny, fspl, ier, 4.0e-7);
            if (ier != 0)
              return;
          }
          zdiff2 = bcymax[ix] - zcur(0);
        } else if (ibcymax == 2) {
          if (ix < inx - 1) {
            zcur(0) = 2.0 * fspl(0, 2, ix, jy) + 6.0 * zhy * fspl(0, 3, ix, jy);
          } else {
            bicubic_eval_explicit(x[inx - 1], y[iny - 1], iselect2, zcur, x,
                                  inx, y, iny, fspl, ier, 4.0e-7);
            if (ier != 0)
              return;
          }
          zdiff2 = bcymax[ix] - zcur(0);
        }

        LO iadr = iasc + ix * iinc;
        for (LO iy = 0; iy < iny; ++iy)
          fspl_l_y(ix * 4, iy) = 0.0;

        fspl_l_y(1 + ix * 4, 0) = 0.0;
        fspl_l_y(2 + ix * 4, 0) = 0.0;
        fspl_l_y(1 + ix * 4, iny - 1) = 0.0;
        fspl_l_y(2 + ix * 4, iny - 1) = 0.0;

        if (ibcymin == 1)
          fspl_l_y(1 + ix * 4, 0) = zdiff1;
        else if (ibcymin == 2)
          fspl_l_y(2 + ix * 4, 0) = zdiff1;

        if (ibcymax == 1)
          fspl_l_y(1 + ix * 4, iny - 1) = zdiff2;
        else if (ibcymax == 2)
          fspl_l_y(2 + ix * 4, iny - 1) = zdiff2;

        auto fspl_s_y =
          Rank2View<T, MemorySpace>(fspl_l_y.data_handle() + iadr, 4, iny);
        compute_spline_coefficients(ibcymin, ibcymax, iny, y, fspl_s_y, wk_y);
      });

    Kokkos::parallel_for(
      Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {inx, iny - 1}),
      KOKKOS_LAMBDA(const LO ix, const LO iy) {
        // Multiply local coefficients
        fspl_l_y(2 + ix * 4, iy) *= xo2;
        fspl_l_y(3 + ix * 4, iy) *= xo6;

        // Conditional accumulation LOo global fspl
        if (ix < inx - 1) {
          fspl(0, 1, ix, iy) += fspl_l_y(1 + ix * 4, iy);
          fspl(0, 2, ix, iy) += fspl_l_y(2 + ix * 4, iy);
          fspl(0, 3, ix, iy) += fspl_l_y(3 + ix * 4, iy);
        }
      });
    Kokkos::parallel_for(
      Kokkos::TeamPolicy<execution_space>(iny - 1, Kokkos::AUTO),
      KOKKOS_LAMBDA(const member_type& team) {
        const LO iy = team.league_rank();

        auto fspl_x_view = Rank2View<T, MemorySpace>(
          fspl_l_x.data_handle() + 4 * inx * iy, 4, inx);
        auto wk_x_view =
          Rank1View<T, MemorySpace>(wk_l.data_handle() + iy * inx, inx);

        for (LO ic = 1; ic < 4; ++ic) {
          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, inx), [=](LO ix) {
            fspl_x_view(0, ix) = fspl_l_y(ic + ix * 4, iy);
          });

          if (team.team_rank() == 0) {
            fspl_x_view(1, 0) = 0.0;
            fspl_x_view(2, 0) = 0.0;
            fspl_x_view(1, inx - 1) = 0.0;
            fspl_x_view(2, inx - 1) = 0.0;

            compute_spline_coefficients(ibcxmin, ibcxmax, inx, x, fspl_x_view,
                                        wk_x_view);
          }

          Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team, inx - 1), [=](LO ix) {
              fspl(1, ic, ix, iy) += fspl_x_view(1, ix);
              fspl(2, ic, ix, iy) += fspl_x_view(2, ix) * xo2;
              fspl(3, ic, ix, iy) += fspl_x_view(3, ix) * xo6;
            });
        }
      });
  }
}

/** \brief Lookup the grid point and compute the displacement
 *
 * \param xget input point in x direction where the spline is evaluated
 * \param yget input point in y direction where the spline is evaluated
 * \param x input grid points in x direction
 * \param nx number of grid points in x direction
 * \param y input grid points in y direction
 * \param ny number of grid points in y direction
 * \param i index of the grid point in x direction
 * \param j index of the grid point in y direction
 * \param xparam displacement from the grid point in x direction
 * \param hx step size in the x direction
 * \param hxi inverse step size in the x direction
 * \param yparam displacement from the grid point in y direction
 * \param hy step size in the y direction
 * \param hyi inverse step size in the y direction
 * \param ier error flag, 0 if successful, non-zero otherwise
 * \param lookup_tol tolerance for the lookup
 */
template <typename T, typename MemorySpace>
KOKKOS_INLINE_FUNCTION void bicubic_lookup_compact(
  T xget, T yget, Rank1View<T, MemorySpace> x, LO nx,
  Rank1View<T, MemorySpace> y, LO ny, LO& i, LO& j, T& xparam, T& yparam, T& hx,
  T& hxi, T& hy, T& hyi, LO& ier, T lookup_tol)
{
  T zxget = xget;
  T zyget = yget;

  clamp_to_grid_range(xget, zxget, x, nx, ier, lookup_tol);
  clamp_to_grid_range(yget, zyget, y, ny, ier, lookup_tol);

  i = lookup_in_grid<T, MemorySpace>(zxget, x, nx);
  j = lookup_in_grid<T, MemorySpace>(zyget, y, ny);

  hx = x[i + 1] - x[i];
  hy = y[j + 1] - y[j];

  hxi = 1.0 / hx;
  hyi = 1.0 / hy;

  xparam = (zxget - x[i]) * hxi;
  yparam = (zyget - y[j]) * hyi;
}

/** \brief Evaluate the bicubic spline at given grid cell
 *
 * /param ict selector array indicating which derivatives to compute
 * \param fval output values of the spline at given grid cell
 * \param i index of the grid cell in x direction
 * \param j index of the grid cell in y direction
 * \param xparam displacement from the grid cell in x direction
 * \param yparam displacement from the grid cell in y direction
 * \param hx step size in the x direction
 * \param hxi inverse step size in the x direction
 * \param hy step size in the y direction
 * \param hyi inverse step size in the y direction
 * \param f bicubic spline coefficients
 */
template <typename T, typename MemorySpace>
KOKKOS_INLINE_FUNCTION void bicubic_evalfn_compact(
  Kokkos::View<LO*, MemorySpace> ict, LO ivec, LO ivecd,
  Rank1View<T, MemorySpace> fval, LO i, LO j, T xparam, T yparam, T hx, T hxi,
  T hy, T hyi, Rank3View<T, MemorySpace> f)
{
  constexpr T sixth = 1.0 / 6.0;
  const T z36th = sixth * sixth;
  LO iadr = 0;

  // ict[0] = 1 or 2  →   f, df/dx, df/dy, d²f/dx², d²f/dy², d²f/dxdy
  if (ict[0] <= 2) {
    T xp = xparam, xpi = 1.0 - xp;
    T yp = yparam, ypi = 1.0 - yp;
    // f (function value)
    if (ict[0] == 1) {
      iadr++;
      T xp2 = xp * xp, xpi2 = xpi * xpi;
      T cx = xp * (xp2 - 1.0);
      T cxi = xpi * (xpi2 - 1.0);
      T hx2 = hx * hx;
      T yp2 = yp * yp, ypi2 = ypi * ypi;
      T cy = yp * (yp2 - 1.0);
      T cyi = ypi * (ypi2 - 1.0);
      T hy2 = hy * hy;

      T sum = xpi * (ypi * f(0, i, j) + yp * f(0, i, j + 1)) +
              xp * (ypi * f(0, i + 1, j) + yp * f(0, i + 1, j + 1));
      sum += sixth * hx2 *
             (cxi * (ypi * f(1, i, j) + yp * f(1, i, j + 1)) +
              cx * (ypi * f(1, i + 1, j) + yp * f(1, i + 1, j + 1)));
      sum += sixth * hy2 *
             (xpi * (cyi * f(2, i, j) + cy * f(2, i, j + 1)) +
              xp * (cyi * f(2, i + 1, j) + cy * f(2, i + 1, j + 1)));
      sum += z36th * hx2 * hy2 *
             (cxi * (cyi * f(3, i, j) + cy * f(3, i, j + 1)) +
              cx * (cyi * f(3, i + 1, j) + cy * f(3, i + 1, j + 1)));
      fval(iadr - 1) = sum;
    }

    // df/dx
    if (ict[1] == 1) {
      iadr++;
      T xp2 = xp * xp, xpi2 = xpi * xpi;
      T cxd = 3.0 * xp2 - 1.0;
      T cxdi = -3.0 * xpi2 + 1.0;
      T yp2 = yp * yp, ypi2 = ypi * ypi;
      T cy = yp * (yp2 - 1.0);
      T cyi = ypi * (ypi2 - 1.0);
      T hy2 = hy * hy;

      T sum = hxi * (-(ypi * f(0, i, j) + yp * f(0, i, j + 1)) +
                     (ypi * f(0, i + 1, j) + yp * f(0, i + 1, j + 1)));
      sum += sixth * hx *
             (cxdi * (ypi * f(1, i, j) + yp * f(1, i, j + 1)) +
              cxd * (ypi * f(1, i + 1, j) + yp * f(1, i + 1, j + 1)));
      sum += sixth * hxi * hy2 *
             (-(cyi * f(2, i, j) + cy * f(2, i, j + 1)) +
              (cyi * f(2, i + 1, j) + cy * f(2, i + 1, j + 1)));
      sum += z36th * hx * hy2 *
             (cxdi * (cyi * f(3, i, j) + cy * f(3, i, j + 1)) +
              cxd * (cyi * f(3, i + 1, j) + cy * f(3, i + 1, j + 1)));
      fval(iadr - 1) = sum;
    }

    //  df/dy
    if (ict[2] == 1) {
      iadr++;
      T xp2 = xp * xp, xpi2 = xpi * xpi;
      T cx = xp * (xp2 - 1.0);
      T cxi = xpi * (xpi2 - 1.0);
      T hx2 = hx * hx;
      T yp2 = yp * yp, ypi2 = ypi * ypi;
      T cyd = 3.0 * yp2 - 1.0;
      T cydi = -3.0 * ypi2 + 1.0;

      T sum = hyi * (xpi * (-f(0, i, j) + f(0, i, j + 1)) +
                     xp * (-f(0, i + 1, j) + f(0, i + 1, j + 1)));
      sum += sixth * hx2 * hyi *
             (cxi * (-f(1, i, j) + f(1, i, j + 1)) +
              cx * (-f(1, i + 1, j) + f(1, i + 1, j + 1)));
      sum += sixth * hy *
             (xpi * (cydi * f(2, i, j) + cyd * f(2, i, j + 1)) +
              xp * (cydi * f(2, i + 1, j) + cyd * f(2, i + 1, j + 1)));
      sum += z36th * hx2 * hy *
             (cxi * (cydi * f(3, i, j) + cyd * f(3, i, j + 1)) +
              cx * (cydi * f(3, i + 1, j) + cyd * f(3, i + 1, j + 1)));
      fval(iadr - 1) = sum;
    }

    // d2f/dx2
    if (ict[3] == 1) {
      iadr++;
      T yp2 = yp * yp, ypi2 = ypi * ypi;
      T cy = yp * (yp2 - 1.0);
      T cyi = ypi * (ypi2 - 1.0);
      T hy2 = hy * hy;

      T sum = xpi * (ypi * f(1, i, j) + yp * f(1, i, j + 1)) +
              xp * (ypi * f(1, i + 1, j) + yp * f(1, i + 1, j + 1));
      sum += sixth * hy2 *
             (xpi * (cyi * f(3, i, j) + cy * f(3, i, j + 1)) +
              xp * (cyi * f(3, i + 1, j) + cy * f(3, i + 1, j + 1)));
      fval(iadr - 1) = sum;
    }

    // d2f/dy2
    if (ict[4] == 1) {
      iadr++;
      T xp2 = xp * xp, xpi2 = xpi * xpi;
      T cx = xp * (xp2 - 1.0);
      T cxi = xpi * (xpi2 - 1.0);
      T hx2 = hx * hx;

      T sum = xpi * (ypi * f(2, i, j) + yp * f(2, i, j + 1)) +
              xp * (ypi * f(2, i + 1, j) + yp * f(2, i + 1, j + 1));
      sum += sixth * hx2 *
             (cxi * (ypi * f(3, i, j) + yp * f(3, i, j + 1)) +
              cx * (ypi * f(3, i + 1, j) + yp * f(3, i + 1, j + 1)));
      fval(iadr - 1) = sum;
    }

    // d2f/dxdy
    if (ict[5] == 1) {
      iadr++;
      T xp2 = xp * xp, xpi2 = xpi * xpi;
      T cxd = 3.0 * xp2 - 1.0;
      T cxdi = -3.0 * xpi2 + 1.0;
      T yp2 = yp * yp, ypi2 = ypi * ypi;
      T cyd = 3.0 * yp2 - 1.0;
      T cydi = -3.0 * ypi2 + 1.0;

      T sum =
        hxi * hyi *
        (f(0, i, j) - f(0, i, j + 1) - f(0, i + 1, j) + f(0, i + 1, j + 1));
      sum += sixth * hx * hyi *
             (cxdi * (-f(1, i, j) + f(1, i, j + 1)) +
              cxd * (-f(1, i + 1, j) + f(1, i + 1, j + 1)));
      sum += sixth * hxi * hy *
             (-(cydi * f(2, i, j) + cyd * f(2, i, j + 1)) +
              (cydi * f(2, i + 1, j) + cyd * f(2, i + 1, j + 1)));
      sum += z36th * hx * hy *
             (cxdi * (cydi * f(3, i, j) + cyd * f(3, i, j + 1)) +
              cxd * (cydi * f(3, i + 1, j) + cyd * f(3, i + 1, j + 1)));
      fval(iadr - 1) = sum;
    }
  }

  // ict[0] = 3  → 3rd‑order derivative combinations
  else if (ict[0] == 3) {
    // d3f/dx3
    if (ict[1] == 1) {
      iadr++;
      T yp = yparam, ypi = 1.0 - yp;
      T yp2 = yp * yp, ypi2 = ypi * ypi;
      T cy = yp * (yp2 - 1.0);
      T cyi = ypi * (ypi2 - 1.0);
      T hy2 = hy * hy;
      T sum = hxi * (-(ypi * f(1, i, j) + yp * f(1, i, j + 1)) +
                     (ypi * f(1, i + 1, j) + yp * f(1, i + 1, j + 1)));
      sum += sixth * hy2 * hxi *
             (-(cyi * f(3, i, j) + cy * f(3, i, j + 1)) +
              (cyi * f(3, i + 1, j) + cy * f(3, i + 1, j + 1)));
      fval(iadr - 1) = sum;
    }

    // d3f/dx2dy
    if (ict[2] == 1) {
      iadr++;
      T xp = xparam, xpi = 1.0 - xp;
      T yp = yparam, ypi = 1.0 - yp;
      T yp2 = yp * yp, ypi2 = ypi * ypi;
      T cyd = 3.0 * yp2 - 1.0;
      T cydi = -3.0 * ypi2 + 1.0;
      T sum = hyi * (xpi * (-f(1, i, j) + f(1, i, j + 1)) +
                     xp * (-f(1, i + 1, j) + f(1, i + 1, j + 1)));
      sum += sixth * hy *
             (xpi * (cydi * f(3, i, j) + cyd * f(3, i, j + 1)) +
              xp * (cydi * f(3, i + 1, j) + cyd * f(3, i + 1, j + 1)));
      fval(iadr - 1) = sum;
    }

    // d3f/dxdy2
    if (ict[3] == 1) {
      iadr++;
      T xp = xparam, xpi = 1.0 - xp;
      T xp2 = xp * xp, xpi2 = xpi * xpi;
      T cxd = 3.0 * xp2 - 1.0;
      T cxdi = -3.0 * xpi2 + 1.0;
      T yp = yparam, ypi = 1.0 - yp;
      T sum = hxi * (-(ypi * f(2, i, j) + yp * f(2, i, j + 1)) +
                     (ypi * f(2, i + 1, j) + yp * f(2, i + 1, j + 1)));
      sum += sixth * hx *
             (cxdi * (ypi * f(3, i, j) + yp * f(3, i, j + 1)) +
              cxd * (ypi * f(3, i + 1, j) + yp * f(3, i + 1, j + 1)));
      fval(iadr - 1) = sum;
    }

    // d3f/dy3
    if (ict[4] == 1) {
      iadr++;
      T xp = xparam, xpi = 1.0 - xp;
      T xp2 = xp * xp, xpi2 = xpi * xpi;
      T cx = xp * (xp2 - 1.0);
      T cxi = xpi * (xpi2 - 1.0);
      T hx2 = hx * hx;
      T sum = hyi * (xpi * (-f(2, i, j) + f(2, i, j + 1)) +
                     xp * (-f(2, i + 1, j) + f(2, i + 1, j + 1)));
      sum += sixth * hx2 * hyi *
             (cxi * (-f(3, i, j) + f(3, i, j + 1)) +
              cx * (-f(3, i + 1, j) + f(3, i + 1, j + 1)));
      fval(iadr - 1) = sum;
    }
  }

  // ict[0] = 4  → 4th‑order derivative combinations
  else if (ict[0] == 4) {
    // d4f/dx3dy
    if (ict[1] == 1) {
      iadr++;
      T yp = yparam, ypi = 1.0 - yp;
      T yp2 = yp * yp, ypi2 = ypi * ypi;
      T cyd = 3.0 * yp2 - 1.0;
      T cydi = -3.0 * ypi2 + 1.0;
      T sum =
        hxi * hyi *
        (f(1, i, j) - f(1, i, j + 1) - f(1, i + 1, j) + f(1, i + 1, j + 1));
      sum += sixth * hy * hxi *
             (-(cydi * f(3, i, j) + cyd * f(3, i, j + 1)) +
              (cydi * f(3, i + 1, j) + cyd * f(3, i + 1, j + 1)));
      fval(iadr - 1) = sum;
    }

    // d4f/dx2dy2
    if (ict[2] == 1) {
      iadr++;
      T xp = xparam, xpi = 1.0 - xp;
      T yp = yparam, ypi = 1.0 - yp;
      T sum = xpi * (ypi * f(3, i, j) + yp * f(3, i, j + 1)) +
              xp * (ypi * f(3, i + 1, j) + yp * f(3, i + 1, j + 1));
      fval(iadr - 1) = sum;
    }

    // d4f/dxdy3
    if (ict[3] == 1) {
      iadr++;
      T xp = xparam, xpi = 1.0 - xp;
      T xp2 = xp * xp, xpi2 = xpi * xpi;
      T cxd = 3.0 * xp2 - 1.0;
      T cxdi = -3.0 * xpi2 + 1.0;
      T sum =
        hxi * hyi *
        (f(2, i, j) - f(2, i, j + 1) - f(2, i + 1, j) + f(2, i + 1, j + 1));
      sum += sixth * hx * hyi *
             (cxdi * (-f(3, i, j) + f(3, i, j + 1)) +
              cxd * (-f(3, i + 1, j) + f(3, i + 1, j + 1)));
      fval(iadr - 1) = sum;
    }
  }

  // ict[0] = 5  → 5th‑order derivative combinations
  else if (ict[0] == 5) {
    // d5f/dx3dy2
    if (ict[1] == 1) {
      iadr++;
      T yp = yparam, ypi = 1.0 - yp;
      T sum = hxi * (-(ypi * f(3, i, j) + yp * f(3, i, j + 1)) +
                     (ypi * f(3, i + 1, j) + yp * f(3, i + 1, j + 1)));
      fval(iadr - 1) = sum;
    }

    // d5f/dx2dy3
    if (ict[2] == 1) {
      iadr++;
      T xp = xparam, xpi = 1.0 - xp;
      T sum = hyi * (xpi * (-f(3, i, j) + f(3, i, j + 1)) +
                     xp * (-f(3, i + 1, j) + f(3, i + 1, j + 1)));
      fval(iadr - 1) = sum;
    }
  }

  // ict[0] = 6  → 6th‑order derivative (dx³dy³)
  else if (ict[0] == 6) {
    iadr++;
    T sum = hxi * hyi *
            (f(3, i, j) - f(3, i, j + 1) - f(3, i + 1, j) + f(3, i + 1, j + 1));
    fval(iadr - 1) = sum;
  }
}

/** \brief Evaluate the bicubic spline at a single point with selector
 *
 * \param xget input point in x direction where the spline is evaluated
 * \param yget input point in y direction where the spline is evaluated
 * \param ict selector array indicating which derivatives to compute
 * \param fval output values of the spline at given point
 * \param x input grid points in x direction
 * \param nx number of grid points in x direction
 * \param y input grid points in y direction
 * \param ny number of grid points in y direction
 * \param f bicubic spline coefficients
 * \param ier error flag, 0 if successful, non-zero otherwise
 * \param lookup_tol tolerance for the lookup
 */
template <typename T, typename MemorySpace>
KOKKOS_INLINE_FUNCTION void bicubic_eval_compact(
  T xget, T yget, Kokkos::View<LO*, MemorySpace> ict,
  Rank1View<T, MemorySpace> fval, Rank1View<T, MemorySpace> x, LO nx,
  Rank1View<T, MemorySpace> y, LO ny, Rank3View<T, MemorySpace> f, LO& ier,
  T lookup_tol)
{
  // Local variables
  LO i = 0, j = 0;
  T xparam = 0.0, yparam = 0.0;
  T hx = 0.0, hy = 0.0;
  T hxi = 0.0, hyi = 0.0;

  // Call herm2xy to locate cell and compute params
  bicubic_lookup_compact(xget, yget, x, nx, y, ny, i, j, xparam, yparam, hx,
                         hxi, hy, hyi, ier, lookup_tol);

  if (ier != 0)
    return;

  // Call evalfn with scalar-vector LOerface
  bicubic_evalfn_compact(ict, 1, 1, fval, i, j, xparam, yparam, hx, hxi, hy,
                         hyi, f);
}

template <typename T, typename MemorySpace>
void solve_bicubic_spline_compact(Rank1View<T, MemorySpace> x, LO nx,
                                  Rank1View<T, MemorySpace> y, LO ny,
                                  Rank3View<T, MemorySpace> f, LO ibcxmin,
                                  Rank1View<T, MemorySpace> bcxmin, LO ibcxmax,
                                  Rank1View<T, MemorySpace> bcxmax, LO ibcymin,
                                  Rank1View<T, MemorySpace> bcymin, LO ibcymax,
                                  Rank1View<T, MemorySpace> bcymax,
                                  Rank1View<T, MemorySpace> wk)
{
  using execution_space = typename MemorySpace::execution_space;
  using member_type = typename Kokkos::TeamPolicy<execution_space>::member_type;
  LO iflg2 = 0;

  // Check if inhomogeneous y-boundary conditions exist
  if (ibcymin != -1) {
    detect_inhomogeneous_bc(bcymin, ibcymin, iflg2, nx);
    detect_inhomogeneous_bc(bcymax, ibcymax, iflg2, nx);
  }

  auto fspl_l_x = Rank2View<T, MemorySpace>(wk.data_handle(), 2 * ny, nx);
  auto wk_l =
    Rank1View<T, MemorySpace>(wk.data_handle() + 2 * ny * nx, nx * ny);
  auto fwk4_l_x =
    Rank2View<T, MemorySpace>(wk.data_handle() + 3 * ny * nx, 4 * ny, nx);

  Kokkos::parallel_for(
    Kokkos::TeamPolicy<execution_space>(ny, Kokkos::AUTO),
    KOKKOS_LAMBDA(const member_type& team) {
      const LO iy = team.league_rank();
      T zbcmin = 0.0, zbcmax = 0.0;

      auto fwk_x_view =
        Rank2View<T, MemorySpace>(fspl_l_x.data_handle() + 2 * nx * iy, 2, nx);
      auto wk_x_view =
        Rank1View<T, MemorySpace>(wk_l.data_handle() + nx * iy, nx);
      auto fwk4_x_view =
        Rank2View<T, MemorySpace>(fwk4_l_x.data_handle() + 4 * nx * iy, 4, nx);

      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nx),
                           [=](LO ix) { fwk_x_view(0, ix) = f(0, ix, iy); });

      if (team.team_rank() == 0) {
        if (ibcxmin == 1 || ibcxmin == 2)
          zbcmin = bcxmin[iy];
        if (ibcxmax == 1 || ibcxmax == 2)
          zbcmax = bcxmax[iy];

        solve_cubic_spline_compact(x, nx, fwk_x_view, fwk4_x_view, ibcxmin,
                                   zbcmin, ibcxmax, zbcmax, wk_x_view);
      }

      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nx),
                           [=](LO ix) { f(1, ix, iy) = fwk_x_view(1, ix); });
    });

  Kokkos::parallel_for(
    Kokkos::TeamPolicy<execution_space>(nx, Kokkos::AUTO),
    KOKKOS_LAMBDA(const member_type& team) {
      const LO ix = team.league_rank();
      LO ibcmin = ibcymin, ibcmax = ibcymax;

      auto fwk_y_view =
        Rank2View<T, MemorySpace>(fspl_l_x.data_handle() + 2 * ny * ix, 2, ny);
      auto wk_y_view =
        Rank1View<T, MemorySpace>(wk_l.data_handle() + ny * ix, ny);
      auto fwk4_y_view =
        Rank2View<T, MemorySpace>(fwk4_l_x.data_handle() + 4 * ny * ix, 4, ny);

      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, ny),
                           [=](LO iy) { fwk_y_view(0, iy) = f(0, ix, iy); });

      if (team.team_rank() == 0) {
        if (iflg2 == 1) {
          if (ibcymin == 1 || ibcymin == 2)
            ibcmin = 0;
          if (ibcymax == 1 || ibcymax == 2)
            ibcmax = 0;
        }

        solve_cubic_spline_compact(y, ny, fwk_y_view, fwk4_y_view, ibcmin, 0.0,
                                   ibcmax, 0.0, wk_y_view);
      }
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, ny),
                           [=](LO iy) { f(2, ix, iy) = fwk_y_view(1, iy); });
    });

  Kokkos::parallel_for(
    Kokkos::TeamPolicy<execution_space>(nx, Kokkos::AUTO),
    KOKKOS_LAMBDA(const member_type& team) {
      const LO ix = team.league_rank();
      LO ibcmin = ibcymin, ibcmax = ibcymax;

      auto fwk_y_view =
        Rank2View<T, MemorySpace>(fspl_l_x.data_handle() + 2 * ny * ix, 2, ny);
      auto wk_y_view =
        Rank1View<T, MemorySpace>(wk_l.data_handle() + ny * ix, ny);
      auto fwk4_y_view =
        Rank2View<T, MemorySpace>(fwk4_l_x.data_handle() + 4 * ny * ix, 4, ny);

      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, ny),
                           [=](LO iy) { fwk_y_view(0, iy) = f(1, ix, iy); });

      if (team.team_rank() == 0) {
        if (iflg2 == 1) {
          if (ibcymin == 1 || ibcymin == 2)
            ibcmin = 0;
          if (ibcymax == 1 || ibcymax == 2)
            ibcmax = 0;
        }

        solve_cubic_spline_compact(y, ny, fwk_y_view, fwk4_y_view, ibcmin, 0.0,
                                   ibcmax, 0.0, wk_y_view);
      }
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, ny),
                           [=](LO iy) { f(3, ix, iy) = fwk_y_view(1, iy); });
    });

  // Correct for inhomogeneous boundary conditions if needed
  if (iflg2 == 1) {
    T zbcmin = 0.0;
    T zbcmax = 0.0;

    auto fcorr =
      Rank3View<T, MemorySpace>(wk.data_handle() + 7 * ny * nx, 2, nx, ny);

    Kokkos::parallel_for(
      Kokkos::TeamPolicy<execution_space>(nx, Kokkos::AUTO),
      KOKKOS_LAMBDA(const member_type& team) {
        const LO ix = team.league_rank();
        T zdiff1 = 0.0, zdiff2 = 0.0;

        auto fwk_y_view = Rank2View<T, MemorySpace>(
          fspl_l_x.data_handle() + 2 * ny * ix, 2, ny);
        auto wk_y_view =
          Rank1View<T, MemorySpace>(wk_l.data_handle() + ny * ix, ny);
        auto fwk4_y_view = Rank2View<T, MemorySpace>(
          fwk4_l_x.data_handle() + 4 * ny * ix, 4, ny);

        if (ibcymin == 1) {
          zdiff1 = bcymin[ix] -
                   ((f(0, ix, 1) - f(0, ix, 0)) / (y[1] - y[0]) +
                    (y[1] - y[0]) * (-2.0 * f(2, ix, 0) - f(2, ix, 1)) / 6.0);
        } else if (ibcymin == 2) {
          zdiff1 = bcymin[ix] - f(2, ix, 0);
        } else {
          zdiff1 = 0.0;
        }

        if (ibcymax == 1) {
          zdiff2 =
            bcymax[ix] -
            ((f(0, ix, ny - 1) - f(0, ix, ny - 2)) / (y[ny - 1] - y[ny - 2]) +
             (y[ny - 1] - y[ny - 2]) *
               (2.0 * f(2, ix, ny - 1) + f(2, ix, ny - 2)) / 6.0);
        } else if (ibcymax == 2) {
          zdiff2 = bcymax[ix] - f(2, ix, ny - 1);
        } else {
          zdiff2 = 0.0;
        }

        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, ny),
                             [=](LO iy) { fwk_y_view(0, iy) = 0.0; });
        if (team.team_rank() == 0) {
          solve_cubic_spline_compact(y, ny, fwk_y_view, fwk4_y_view, ibcymin,
                                     zdiff1, ibcymax, zdiff2, wk_y_view);
        }
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, ny), [=](LO iy) {
          fcorr(0, ix, iy) = fwk_y_view(1, iy);
        });
      });

    Kokkos::parallel_for(
      Kokkos::TeamPolicy<execution_space>(ny, Kokkos::AUTO),
      KOKKOS_LAMBDA(const member_type& team) {
        const LO iy = team.league_rank();
        auto fwk_x_view = Rank2View<T, MemorySpace>(
          fspl_l_x.data_handle() + 2 * nx * iy, 2, nx);
        auto wk_x_view =
          Rank1View<T, MemorySpace>(wk_l.data_handle() + nx * iy, nx);
        auto fwk4_x_view = Rank2View<T, MemorySpace>(
          fwk4_l_x.data_handle() + 4 * nx * iy, 4, nx);

        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nx), [=](LO ix) {
          fwk_x_view(0, ix) = fcorr(0, ix, iy);
        });

        if (team.team_rank() == 0) {
          solve_cubic_spline_compact(x, nx, fwk_x_view, fwk4_x_view, ibcxmin,
                                     zbcmin, ibcxmax, zbcmax, wk_x_view);
        }

        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nx), [=](LO ix) {
          fcorr(1, ix, iy) = fwk_x_view(1, ix);
        });
      });

    Kokkos::parallel_for(
      Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {nx, ny}),
      KOKKOS_LAMBDA(const LO ix, const LO iy) {
        f(2, ix, iy) += fcorr(0, ix, iy);
        f(3, ix, iy) += fcorr(1, ix, iy);
      });
  }
}

} // namespace detail

/** \brief Boundary condition options for spline interpolation
 *
 * This enum defines the various boundary conditions that can be applied
 * during spline interpolation. The boundary condition flags can take the
 * following values: PERIODIC | ibcxmin=-1 | ibcxmax=-1  --  periodic boundary
 * condition, where s'(x(0))=s'(x(-1)) and s''(x(0))=s''(x(-1)). NOT_A_KNOT |
 * ibcxmin=0 | ibcxmax=0  --  not-a-knot boundary condition.
 * FIRST_DERIVATIVE_MATCH | ibcxmin=1 | ibcxmax=1  --  first derivative match at
 * x(0) or x(-1). SECOND_DERIVATIVE_MATCH | ibcxmin=2 | ibcxmax=2  --  second
 * derivative match at x(0) or x(-1). FIRST_DERIVATIVE_ZERO | ibcxmin=3 |
 * ibcxmax=3  --  first derivative equal to zero at x(0) or x(-1).
 * SECOND_DERIVATIVE_ZERO | ibcxmin=4 | ibcxmax=4  --  second derivative equal
 * to zero at x(0) or x(-1). FIRST_DIV_DIFF_MATCH | ibcxmin=5 | ibcxmax=5  --
 * first derivative matches the first divided difference. SECOND_DIV_DIFF_MATCH
 * | ibcxmin=6 | ibcxmax=6  --  second derivative matches the second divided
 * difference. THIRD_DIV_DIFF_MATCH | ibcxmin=7 | ibcxmax=7  --  third
 * derivative matches the third divided difference.
 */
enum class BoundaryCondition
{
  PERIODIC = -1,  ///< Periodic boundary condition: s'(x(0))=s'(x(-1)) and
                  ///< s''(x(0))=s''(x(-1))
  NOT_A_KNOT = 0, ///< Not-a-knot boundary condition
  FIRST_DERIVATIVE_MATCH = 1,  ///< First derivative match at boundary
  SECOND_DERIVATIVE_MATCH = 2, ///< Second derivative match at boundary
  FIRST_DERIVATIVE_ZERO = 3,   ///< First derivative equal to zero at boundary
  SECOND_DERIVATIVE_ZERO = 4,  ///< Second derivative equal to zero at boundary
  FIRST_DIV_DIFF_MATCH =
    5, ///< First derivative matches first divided difference
  SECOND_DIV_DIFF_MATCH =
    6, ///< Second derivative matches second divided difference
  THIRD_DIV_DIFF_MATCH =
    7 ///< Third derivative matches third divided difference
};

/** \brief CubicSplineInterpolator class
 *
 * \details Derived from from Princeton Splines library
 *   Original developers: PPPL Computational Plasma Physics Group
 *   Original release: February 5, 1999
 */
template <typename T, typename MemorySpace>
class CubicSplineInterpolator
{
public:
  virtual void evaluate(Rank1View<T, MemorySpace> xvec,
                        Rank2View<T, MemorySpace> fval) = 0;
  virtual void evaluate(Kokkos::View<LO*, MemorySpace> selector,
                        Rank1View<T, MemorySpace> xvec,
                        Rank2View<T, MemorySpace> fval, T lookup_tol) = 0;
};

/** \brief ExplicitCubicSplineInterpolator class
 *
 * \details This class provides methods for 1d explicit cubic spline
 * interpolation and serves as a base class for higher dimensions and compact
 * cubic spline interpolator implementations.
 *
 * The 1d explicit cubic spline are evaluated using the following formula:
 *     s(x)=C(1,i)+(x-x(i))*C(2,i)+(x-x(i))**2*C(3,i)+(x-x(i))**3*C(4,i)
 *         for x(i).le.x.le.x(i+1)
 * where C(1,i), C(2,i), C(3,i), and C(4,i) are the coefficients of the cubic
 * spline.
 *
 * Performance considerations: The memory cost of coeeficients are:
 * - explicit -- 4**dimension*(N1*N2*...*Nn)
 * - compact  -- 2**dimension*(N1*N2*...*Nn)
 * Explicit spline representation offers better performance in 1D cases.
 * (Compact representations actually call the explicit spline setup functions.)
 * Compact spline representation provides performance advantages in higher
 * dimensions due to memory efficiency gains that outweigh the modest increase
 * in computational cost.
 */
template <typename T, typename MemorySpace>
class ExplicitCubicSplineInterpolator
  : public CubicSplineInterpolator<T, MemorySpace>
{
public:
  using execution_space = typename MemorySpace::execution_space;
  using member_type = typename Kokkos::TeamPolicy<execution_space>::member_type;

  /** \brief Construct the interpolator object and set up the coefficients
   *
   * \param x input grid points
   * \param values values of given grid points
   * \param ibcxmin bounday condition flag at x(0)
   * \param bcxmin bounday condition data at x(0)
   * \param ibcxmax bounday condition flag at x(-1)
   * \param bcxmax bounday condition data at x(-1)
   */
  ExplicitCubicSplineInterpolator(Rank1View<T, MemorySpace> x,
                                  Rank1View<T, MemorySpace> values,
                                  BoundaryCondition ibcxmin, T bcxmin,
                                  BoundaryCondition ibcxmax, T bcxmax);

  /** \brief Construct the interpolator object and set up the coefficients with
   * default boundary conditions
   *
   * \param x input grid points
   * \param values values of given grid points
   *
   * \details The boundary condition flags are set to 0, which corresponds to
   * not-a-knot boundary condition.
   */
  ExplicitCubicSplineInterpolator(Rank1View<T, MemorySpace> x,
                                  Rank1View<T, MemorySpace> values);

  /** \brief Evaluate the function value of cubic spline at given points
   *
   * \param xvec input points where the spline is evaluated
   * \param fval output values of the spline at given points
   *
   * \details The size of fval should match the size of xvec.
   */
  void evaluate(Rank1View<T, MemorySpace> xvec,
                Rank2View<T, MemorySpace> fval) override;

  /** \brief Evaluate the cubic spline at given points with selector
   *
   * \param selector selector array indicating which derivatives to compute
   * \param xvec input points where the spline is evaluated
   * \param fval output values of the spline at given points
   * \param lookup_tol tolerance for the lookup
   *
   * \details The size of fval should match the size of xvec and the size of
   * selector. The selector array can be used to specify which derivatives to
   * compute. For example, if selector = {1, 1, 0}, then the function value and
   * first derivative will be computed.
   */
  void evaluate(Kokkos::View<LO*, MemorySpace> selector,
                Rank1View<T, MemorySpace> xvec, Rank2View<T, MemorySpace> fval,
                T lookup_tol = 4.0e-7) override;

private:
  Rank1View<T, MemorySpace> x_;
  Kokkos::View<T*, MemorySpace> coefficients_;
  Rank2View<T, MemorySpace> get_coefficients_view();
};

/**
 * CompactCubicSplineInterpolator class
 * This class provides methods for 1d compact cubic spline interpolation.
 * It inherits from ExplicitCubicSplineInterpolator and implements the
 * compact representation of cubic spline coefficients.
 *
 * The 1d compact cubic spline are evaluated using the following formula:
 * If x(i).le.x.le.x(i+1), let
 *   h = (x(i+1)-x(i))
 *   p = (x-x(i))/h
 *   q = 1-p
 *   s(x)= q*D(1,i) + p*D(1,i+1) +
 *       (h**2/6)*{[q**3-q]*D(2,i) + [p**3-p]*D(2,i+1)}
 * where D(1,i) and D(2,i) are the coefficients of the compact cubic spline.
 */
template <typename T, typename MemorySpace>
class CompactCubicSplineInterpolator
  : public CubicSplineInterpolator<T, MemorySpace>
{
public:
  using execution_space = typename MemorySpace::execution_space;
  using member_type = typename Kokkos::TeamPolicy<execution_space>::member_type;

  /** \brief Construct the interpolator object and set up the coefficients
   *
   * \param x input grid points
   * \param values values of given grid points
   * \param ibcxmin bounday condition flag at x(0)
   * \param bcxmin bounday condition data at x(0)
   * \param ibcxmax bounday condition flag at x(-1)
   * \param bcxmax bounday condition data at x(-1)
   */
  CompactCubicSplineInterpolator(Rank1View<T, MemorySpace> x,
                                 Rank1View<T, MemorySpace> values,
                                 BoundaryCondition ibcxmin, T bcxmin,
                                 BoundaryCondition ibcxmax, T bcxmax);

  /** \brief Construct the interpolator object and set up the coefficients with
   * default boundary conditions
   *
   * \param x input grid points
   * \param values values of given grid points
   */
  CompactCubicSplineInterpolator(Rank1View<T, MemorySpace> x,
                                 Rank1View<T, MemorySpace> values);

  /** \brief Evaluate the function value of compact cubic spline at given points
   *
   * \param xvec input points where the spline is evaluated
   * \param fval output values of the spline at given points
   */
  void evaluate(Rank1View<T, MemorySpace> xvec,
                Rank2View<T, MemorySpace> fval) override;

  /** \brief Evaluate the compact cubic spline at given points with selector
   *
   * \param selector selector array indicating which derivatives to compute
   * \param xvec input points where the spline is evaluated
   * \param fval output values of the spline at given points
   * \param lookup_tol tolerance for the lookup
   *
   * \details The selector follows the same rules as in
   * ExplicitCubicSplineInterpolator.
   */
  void evaluate(Kokkos::View<LO*, MemorySpace> selector,
                Rank1View<T, MemorySpace> xvec, Rank2View<T, MemorySpace> fval,
                T lookup_tol = 4.0e-7) override;

private:
  Rank1View<T, MemorySpace> x_;
  Kokkos::View<T*, MemorySpace> coefficients_;
  Rank2View<T, MemorySpace> get_coefficients_view();
};

template <typename T, typename MemorySpace>
class BiCubicSplineInterpolator
{
public:
  virtual void evaluate(Rank1View<T, MemorySpace> xvec,
                        Rank1View<T, MemorySpace> yvec,
                        Rank2View<T, MemorySpace> fval) = 0;
  virtual void evaluate(Kokkos::View<LO*, MemorySpace> selector,
                        Rank1View<T, MemorySpace> xvec,
                        Rank1View<T, MemorySpace> yvec,
                        Rank2View<T, MemorySpace> fval, T lookup_tol) = 0;
};

/** \brief BiCubicSplineInterpolator class
 * This class provides methods for 2D explicit bicubic spline interpolation.
 * It inherits from ExplicitCubicSplineInterpolator and implements the bicubic
 * representation of cubic spline coefficients.
 *
 * The 2D explicit bicubic spline are evaluated using the following formula:
 *     s(x,y)=  C11 + dx*C21 + dx**2*C31 + dx**3*C41
 *         +dy*(C12 + dx*C22 + dx**2*C32 + dx**3*C42)
 *      +dy**2*(C13 + dx*C23 + dx**2*C33 + dx**3*C43)
 *      +dy**3*(C14 + dx*C24 + dx**2*C34 + dx**3*C44)
 * where Cij are the coefficients of the bicubic spline.
 */
template <typename T, typename MemorySpace>
class ExplicitBiCubicSplineInterpolator
  : public BiCubicSplineInterpolator<T, MemorySpace>
{
public:
  using execution_space = typename MemorySpace::execution_space;
  using member_type = typename Kokkos::TeamPolicy<execution_space>::member_type;

  /** \brief Construct the interpolator object and set up the coefficients
   *
   * \param x input grid points in x direction
   * \param y input grid points in y direction
   * \param values values of given grid points
   * \param ibcxmin boundary condition flag at x(0)
   * \param bcxmin boundary condition data at x(0)
   * \param ibcxmax boundary condition flag at x(-1)
   * \param bcxmax boundary condition data at x(-1)
   * \param ibcymin boundary condition flag at y(0)
   * \param bcymin boundary condition data at y(0)
   * \param ibcymax boundary condition flag at y(-1)
   * \param bcymax boundary condition data at y(-1)
   *
   * \details The boundary condition flags take the same values as in
   * CubicSplineInterpolator at each direction.
   */
  ExplicitBiCubicSplineInterpolator(
    Rank1View<T, MemorySpace> x, Rank1View<T, MemorySpace> y,
    Rank1View<T, MemorySpace> values, BoundaryCondition ibcxmin,
    Rank1View<T, MemorySpace> bcxmin, BoundaryCondition ibcxmax,
    Rank1View<T, MemorySpace> bcxmax, BoundaryCondition ibcymin,
    Rank1View<T, MemorySpace> bcymin, BoundaryCondition ibcymax,
    Rank1View<T, MemorySpace> bcymax);

  /** \brief Construct the interpolator object and set up the coefficients with
   * default boundary conditions
   *
   * \param x input grid points in x direction
   * \param y input grid points in y direction
   * \param values values of given grid points
   */
  ExplicitBiCubicSplineInterpolator(Rank1View<T, MemorySpace> x,
                                    Rank1View<T, MemorySpace> y,
                                    Rank1View<T, MemorySpace> values);

  /** \brief Evaluate the bicubic spline at given points
   *
   * \param xvec input points in x direction where the spline is evaluated
   * \param yvec input points in y direction where the spline is evaluated
   * \param fval output values of the spline at given points
   */
  void evaluate(Rank1View<T, MemorySpace> xvec, Rank1View<T, MemorySpace> yvec,
                Rank2View<T, MemorySpace> fval) override;

  /** \brief Evaluate the bicubic spline at given points with selector
   *
   * \param iselect selector array indicating which derivatives to compute
   * \param xvec input points in x direction where the spline is evaluated
   * \param yvec input points in y direction where the spline is evaluated
   * \param fval output values of the spline at given points
   * \param lookup_tol tolerance for the lookup
   *
   * \details The selector follows the the following rules:
   *   iselect(1)=1 -- evaluate f
   *   iselect(2)=1 -- evaluate df/dx
   *   iselect(3)=1 -- evaluate df/dy
   *   iselect(4)=1 -- evaluate d2f/dx2
   *   iselect(5)=1 -- evaluate d2f/dy2
   *   iselect(6)=1 -- evaluate d2f/dxdy
   */
  void evaluate(Kokkos::View<LO*, MemorySpace> iselect,
                Rank1View<T, MemorySpace> xvec, Rank1View<T, MemorySpace> yvec,
                Rank2View<T, MemorySpace> fval, T lookup_tol = 4.0e-7) override;

private:
  Rank1View<T, MemorySpace> x_;
  Rank1View<T, MemorySpace> y_;
  Kokkos::View<T*, MemorySpace> coefficients_;
  Rank4View<T, MemorySpace> get_coefficients_view();
};

/** \brief CompactBiCubicSplineInterpolator class
 * This class provides methods for 2D compact bicubic spline interpolation.
 * It inherits from ExplicitBiCubicSplineInterpolator and implements the compact
 * representation of bicubic spline coefficients.
 *
 * The 2D compact bicubic spline are evaluated using the following formula:
 * for x(i).le.x.le.x(i+1) and y(j).le.y.le.y(j+1), let
 *   hx = x(i+1)-x(i)    hy = y(j+1)-y(j)
 *   px = (x-x(i))/hx    py = (y-y(j))/hy
 *   qx = 1-px           qy = 1-py
 *   rx = (px**3-px)     ry = (py**3-py)
 *   sx = (qx**3-qx)     sy = (qy**3-qy)
 *
 *  s(x,y) = qx*(qy*D(1,i,j)+py*D(1,i,j+1)) +
 *        px*(qy*D(1,i+1,j)+py*D(1,i+1,j+1)) +
 *        (hx**2/6)*{sx*(qy*D(2,i,j)+py*D(2,i,j+1)) +
 *                   rx*(qy*D(2,i+1,j)+py*D(2,i+1,j+1))} +
 *        (hy**2/6)*{qx*(sy*D(3,i,j)+ry*D(3,i,j+1)) +
 *                   px*(sy*D(3,i+1,j)+ry*D(3,i+1,j+1))} +
 *        (hx**2*hy**2/36)*{sx*(sy*D(4,i,j)+ry*D(4,i,j+1)) +
 *                   rx*(sy*D(4,i+1,j)+ry*D(4,i+1,j+1))}
 */
template <typename T, typename MemorySpace>
class CompactBiCubicSplineInterpolator
  : public BiCubicSplineInterpolator<T, MemorySpace>
{
public:
  using execution_space = typename MemorySpace::execution_space;
  using member_type = typename Kokkos::TeamPolicy<execution_space>::member_type;

  /** \brief Construct the interpolator object and set up the coefficients
   *
   * \param x input grid points in x direction
   * \param y input grid points in y direction
   * \param values values of given grid points
   * \param ibcxmin boundary condition flag at x(0)
   * \param bcxmin boundary condition data at x(0)
   * \param ibcxmax boundary condition flag at x(-1)
   * \param bcxmax boundary condition data at x(-1)
   * \param ibcymin boundary condition flag at y(0)
   * \param bcymin boundary condition data at y(0)
   * \param ibcymax boundary condition flag at y(-1)
   * \param bcymax boundary condition data at y(-1)
   *
   * \details The boundary condition flags take the same values as in
   * CubicSplineInterpolator at each direction.
   */
  CompactBiCubicSplineInterpolator(
    Rank1View<T, MemorySpace> x, Rank1View<T, MemorySpace> y,
    Rank1View<T, MemorySpace> values, BoundaryCondition ibcxmin,
    Rank1View<T, MemorySpace> bcxmin, BoundaryCondition ibcxmax,
    Rank1View<T, MemorySpace> bcxmax, BoundaryCondition ibcymin,
    Rank1View<T, MemorySpace> bcymin, BoundaryCondition ibcymax,
    Rank1View<T, MemorySpace> bcymax);

  /** \brief Construct the interpolator object and set up the coefficients with
   * default boundary conditions \param x input grid points in x direction
   * \param y input grid points in y direction
   * \param values values of given grid points
   */
  CompactBiCubicSplineInterpolator(Rank1View<T, MemorySpace> x,
                                   Rank1View<T, MemorySpace> y,
                                   Rank1View<T, MemorySpace> values);

  /** \brief Evaluate the bicubic spline at given points
   *
   * \param xvec input points in x direction where the spline is evaluated
   * \param yvec input points in y direction where the spline is evaluated
   * \param fval output values of the spline at given points
   */
  void evaluate(Rank1View<T, MemorySpace> xvec, Rank1View<T, MemorySpace> yvec,
                Rank2View<T, MemorySpace> fval) override;

  /** \brief Evaluate the bicubic spline at given points with selector
   *
   * \param iselect selector array indicating which derivatives to compute
   * \param xvec input points in x direction where the spline is evaluated
   * \param yvec input points in y direction where the spline is evaluated
   * \param fval output values of the spline at given points
   * \param lookup_tol tolerance for locating the interval in the grid
   *
   * \details The selector follows the same rules as in
   * ExplicitBiCubicSplineInterpolator.
   */
  void evaluate(Kokkos::View<LO*, MemorySpace> iselect,
                Rank1View<T, MemorySpace> xvec, Rank1View<T, MemorySpace> yvec,
                Rank2View<T, MemorySpace> fval, T lookup_tol = 4.0e-7) override;

private:
  Rank1View<T, MemorySpace> x_;
  Rank1View<T, MemorySpace> y_;
  Kokkos::View<T*, MemorySpace> coefficients_;
  Rank3View<T, MemorySpace> get_coefficients_view();
};

template <typename T, typename MemorySpace>
Rank2View<T, MemorySpace>
ExplicitCubicSplineInterpolator<T, MemorySpace>::get_coefficients_view()
{
  return Rank2View<T, MemorySpace>(coefficients_.data(), 4, x_.extent(0));
}

template <typename T, typename MemorySpace>
Rank2View<T, MemorySpace>
CompactCubicSplineInterpolator<T, MemorySpace>::get_coefficients_view()
{
  return Rank2View<T, MemorySpace>(coefficients_.data(), 2, x_.extent(0));
}

template <typename T, typename MemorySpace>
Rank4View<T, MemorySpace>
ExplicitBiCubicSplineInterpolator<T, MemorySpace>::get_coefficients_view()
{
  return Rank4View<T, MemorySpace>(coefficients_.data(), 4, 4, x_.extent(0),
                                   y_.extent(0));
}

template <typename T, typename MemorySpace>
Rank3View<T, MemorySpace>
CompactBiCubicSplineInterpolator<T, MemorySpace>::get_coefficients_view()
{
  return Rank3View<T, MemorySpace>(coefficients_.data(), 4, x_.extent(0),
                                   y_.extent(0));
}

template <typename T, typename MemorySpace>
struct InitCubicCoeffFunctor
{

  Rank2View<T, MemorySpace> fspl;
  Rank1View<T, MemorySpace> values;

  // Constructor to initialize captured variables
  InitCubicCoeffFunctor(Rank2View<T, MemorySpace> fspl_,
                        Rank1View<T, MemorySpace> values_)
    : fspl(fspl_), values(values_)
  {
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const LO i) const { fspl(0, i) = values(i); }
};

template <typename T, typename MemorySpace>
struct SolveCompactCubicSplineFunctor
{
  using execution_space = typename MemorySpace::execution_space;
  using member_type = typename Kokkos::TeamPolicy<execution_space>::member_type;

  Rank1View<T, MemorySpace> x;
  LO nx;
  Rank2View<T, MemorySpace> fspl;
  Rank2View<T, MemorySpace> fspl4;
  LO ibcxmin;
  T bcxmin;
  LO ibcxmax;
  T bcxmax;
  Rank1View<T, MemorySpace> wk;

  SolveCompactCubicSplineFunctor(Rank1View<T, MemorySpace> x_, LO nx_,
                                 Rank2View<T, MemorySpace> fspl_,
                                 Rank2View<T, MemorySpace> fspl4_, LO ibcxmin_,
                                 T bcxmin_, LO ibcxmax_, T bcxmax_,
                                 Rank1View<T, MemorySpace> wk_)
    : x(x_),
      nx(nx_),
      fspl(fspl_),
      fspl4(fspl4_),
      ibcxmin(ibcxmin_),
      bcxmin(bcxmin_),
      ibcxmax(ibcxmax_),
      bcxmax(bcxmax_),
      wk(wk_)
  {
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const member_type& team_member) const
  {
    Kokkos::single(Kokkos::PerTeam(team_member), [=]() {
      detail::solve_cubic_spline_compact(x, nx, fspl, fspl4, ibcxmin, bcxmin,
                                         ibcxmax, bcxmax, wk);
    });
  }
};

template <typename T, typename MemorySpace>
struct SolveExplicitCubicSplineFunctor
{
  using execution_space = typename MemorySpace::execution_space;
  using member_type = typename Kokkos::TeamPolicy<execution_space>::member_type;

  Rank1View<T, MemorySpace> x;
  LO nx;
  Rank2View<T, MemorySpace> fspl;
  LO ibcxmin;
  T bcxmin;
  LO ibcxmax;
  T bcxmax;
  Rank1View<T, MemorySpace> wk;

  SolveExplicitCubicSplineFunctor(Rank1View<T, MemorySpace> x_, LO nx_,
                                  Rank2View<T, MemorySpace> fspl_, LO ibcxmin_,
                                  T bcxmin_, LO ibcxmax_, T bcxmax_,
                                  Rank1View<T, MemorySpace> wk_)
    : x(x_),
      nx(nx_),
      fspl(fspl_),
      ibcxmin(ibcxmin_),
      bcxmin(bcxmin_),
      ibcxmax(ibcxmax_),
      bcxmax(bcxmax_),
      wk(wk_)
  {
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const member_type& team_member) const
  {
    Kokkos::single(Kokkos::PerTeam(team_member), [=]() {
      detail::solve_cubic_spline_explicit(x, nx, fspl, ibcxmin, bcxmin, ibcxmax,
                                          bcxmax, wk);
    });
  }
};

template <typename T, typename MemorySpace>
ExplicitCubicSplineInterpolator<T, MemorySpace>::
  ExplicitCubicSplineInterpolator(Rank1View<T, MemorySpace> x,
                                  Rank1View<T, MemorySpace> values,
                                  BoundaryCondition ibcxmin, T bcxmin,
                                  BoundaryCondition ibcxmax, T bcxmax)
{
  detail::grid_valid<T, MemorySpace>(x);
  x_ = x;
  LO nx = x.extent(0);
  Kokkos::View<T*, MemorySpace> fspl_view("coefficients", 4 * nx);
  coefficients_ = fspl_view;
  auto fspl = get_coefficients_view();
  InitCubicCoeffFunctor<T, MemorySpace> init_functor(fspl, values);
  Kokkos::parallel_for("initialize_coefficients", nx, init_functor);

  Kokkos::View<T*, MemorySpace> wk_view("working space", 10);
  auto wk = Rank1View<T, MemorySpace>(wk_view.data(), 10);
  // TODO: pass member type to solve_spline
  SolveExplicitCubicSplineFunctor<T, MemorySpace> functor(
    x, nx, fspl, static_cast<LO>(ibcxmin), bcxmin, static_cast<LO>(ibcxmax),
    bcxmax, wk);
  Kokkos::parallel_for("solve_cubic_spline_explicit",
                       Kokkos::TeamPolicy<execution_space>(1, Kokkos::AUTO),
                       functor);
}

template <typename T, typename MemorySpace>
ExplicitCubicSplineInterpolator<T, MemorySpace>::
  ExplicitCubicSplineInterpolator(Rank1View<T, MemorySpace> x,
                                  Rank1View<T, MemorySpace> values)
  : ExplicitCubicSplineInterpolator(x, values, BoundaryCondition::NOT_A_KNOT, 0,
                                    BoundaryCondition::NOT_A_KNOT, 0)
{
  std::cout << "Not a knot boundary condition is applied by default"
            << std::endl;
}

template <typename T, typename MemorySpace>
CompactCubicSplineInterpolator<T, MemorySpace>::CompactCubicSplineInterpolator(
  Rank1View<T, MemorySpace> x, Rank1View<T, MemorySpace> values,
  BoundaryCondition ibcxmin, T bcxmin, BoundaryCondition ibcxmax, T bcxmax)
{
  detail::grid_valid<T, MemorySpace>(x);
  x_ = x;
  LO nx = x.extent(0);
  Kokkos::View<T*, MemorySpace> fspl_view("coefficients", 2 * nx);
  coefficients_ = fspl_view;
  auto fspl = get_coefficients_view();
  InitCubicCoeffFunctor<T, MemorySpace> init_functor(fspl, values);
  Kokkos::parallel_for("initialize_coefficients", nx, init_functor);

  Kokkos::View<T*, MemorySpace> fspl4_view("explicit coefficients", 4 * nx);
  auto fspl4 = Rank2View<T, MemorySpace>(fspl4_view.data(), 4, nx);
  Kokkos::View<T*, MemorySpace> wk_view("working space", 10);
  auto wk = Rank1View<T, MemorySpace>(wk_view.data(), 10);

  SolveCompactCubicSplineFunctor<T, MemorySpace> functor(
    x, nx, fspl, fspl4, static_cast<LO>(ibcxmin), bcxmin,
    static_cast<LO>(ibcxmax), bcxmax, wk);
  Kokkos::parallel_for("solve_cubic_spline_compact",
                       Kokkos::TeamPolicy<execution_space>(1, Kokkos::AUTO),
                       functor);
}

template <typename T, typename MemorySpace>
CompactCubicSplineInterpolator<T, MemorySpace>::CompactCubicSplineInterpolator(
  Rank1View<T, MemorySpace> x, Rank1View<T, MemorySpace> fx)
  : CompactCubicSplineInterpolator(x, fx, BoundaryCondition::NOT_A_KNOT, 0,
                                   BoundaryCondition::NOT_A_KNOT, 0)
{
  std::cout << "Not a knot boundary condition is applied by default"
            << std::endl;
}

template <typename T, typename MemorySpace>
void ExplicitCubicSplineInterpolator<T, MemorySpace>::evaluate(
  const Kokkos::View<LO*, MemorySpace> selector, Rank1View<T, MemorySpace> xvec,
  Rank2View<T, MemorySpace> fval, T lookup_tol)
{
  LO ivd = fval.extent(1);
  auto fspl = get_coefficients_view();
  Kokkos::parallel_for(
    Kokkos::RangePolicy<execution_space>(0, xvec.extent(0)),
    KOKKOS_CLASS_LAMBDA(const LO i) {
      LO ier = 0;
      auto fval_view =
        Rank1View<T, MemorySpace>(fval.data_handle() + i * ivd, ivd);
      detail::cubic_eval_explicit(xvec(i), selector, fval_view, x_,
                                  x_.extent(0), fspl, ier, lookup_tol);
    });
}

template <typename T, typename MemorySpace>
void ExplicitCubicSplineInterpolator<T, MemorySpace>::evaluate(
  Rank1View<T, MemorySpace> xvec, Rank2View<T, MemorySpace> fval)
{
  std::array<LO, 3> ict_arr = {1, 0, 0};
  evaluate(Kokkos::View<LO*, MemorySpace>(ict_arr.data(), 3), xvec, fval);
};

template <typename T, typename MemorySpace>
void CompactCubicSplineInterpolator<T, MemorySpace>::evaluate(
  Kokkos::View<LO*, MemorySpace> selector, Rank1View<T, MemorySpace> xvec,
  Rank2View<T, MemorySpace> fval, T lookup_tol)
{
  LO ivd = fval.extent(1);
  auto fspl = get_coefficients_view();
  Kokkos::parallel_for(
    Kokkos::RangePolicy<execution_space>(0, xvec.extent(0)),
    KOKKOS_CLASS_LAMBDA(const LO i) {
      LO ier = 0;
      auto fval_view =
        Rank1View<T, MemorySpace>(fval.data_handle() + i * ivd, ivd);
      detail::cubic_eval_compact(xvec(i), selector, fval_view, x_, x_.extent(0),
                                 fspl, ier, lookup_tol);
    });
}

// TODO: evluation functions are exactly the same, so put them in base class and
// override?
//       need to figure out if this is possible with the Kokkos
template <typename T, typename MemorySpace>
void CompactCubicSplineInterpolator<T, MemorySpace>::evaluate(
  Rank1View<T, MemorySpace> xvec, Rank2View<T, MemorySpace> fval)
{
  std::array<LO, 3> ict_arr = {1, 0, 0};
  evaluate(Kokkos::View<LO*, MemorySpace>(ict_arr.data(), 3), xvec, fval);
};

template <typename T, typename MemorySpace>
void ExplicitBiCubicSplineInterpolator<T, MemorySpace>::evaluate(
  Kokkos::View<LO*, MemorySpace> iselect, Rank1View<T, MemorySpace> xvec,
  Rank1View<T, MemorySpace> yvec, Rank2View<T, MemorySpace> fval, T lookup_tol)
{
  PCMS_ALWAYS_ASSERT(xvec.extent(0) == yvec.extent(0) &&
                     xvec.extent(0) == fval.extent(0));

  LO ivd = fval.extent(1);
  auto fspl = get_coefficients_view();
  Kokkos::parallel_for(
    Kokkos::RangePolicy<execution_space>(0, xvec.extent(0)),
    KOKKOS_CLASS_LAMBDA(const LO i) {
      LO ier = 0;
      auto fval_view =
        Rank1View<T, MemorySpace>(fval.data_handle() + i * ivd, ivd);
      detail::bicubic_eval_explicit(xvec(i), yvec(i), iselect, fval_view, x_,
                                    x_.extent(0), y_, y_.extent(0), fspl, ier,
                                    lookup_tol);
    });
}

template <typename T, typename MemorySpace>
void ExplicitBiCubicSplineInterpolator<T, MemorySpace>::evaluate(
  Rank1View<T, MemorySpace> xvec, Rank1View<T, MemorySpace> yvec,
  Rank2View<T, MemorySpace> fval)
{
  LO ict_arr[6] = {1, 0, 0, 0, 0, 0};
  Kokkos::View<LO*, MemorySpace> ict("select", 6);
  Kokkos::parallel_for(
    Kokkos::RangePolicy<execution_space>(0, ict.size()),
    KOKKOS_LAMBDA(const LO i) { ict(i) = ict_arr[i]; });
  evaluate(ict, xvec, yvec, fval);
}

template <typename T, typename MemorySpace>
void CompactBiCubicSplineInterpolator<T, MemorySpace>::evaluate(
  Kokkos::View<LO*, MemorySpace> iselect, Rank1View<T, MemorySpace> xvec,
  Rank1View<T, MemorySpace> yvec, Rank2View<T, MemorySpace> fval, T lookup_tol)
{

  PCMS_ALWAYS_ASSERT(xvec.extent(0) == yvec.extent(0) &&
                     xvec.extent(0) == fval.extent(0));
  LO ivd = fval.extent(1);
  auto fspl = get_coefficients_view();
  Kokkos::parallel_for(
    Kokkos::RangePolicy<execution_space>(0, xvec.extent(0)),
    KOKKOS_CLASS_LAMBDA(const LO i) {
      LO ier = 0;
      auto fval_view =
        Rank1View<T, MemorySpace>(fval.data_handle() + i * ivd, ivd);
      detail::bicubic_eval_compact(xvec(i), yvec(i), iselect, fval_view, x_,
                                   x_.extent(0), y_, y_.extent(0), fspl, ier,
                                   lookup_tol);
    });
}

template <typename T, typename MemorySpace>
void CompactBiCubicSplineInterpolator<T, MemorySpace>::evaluate(
  Rank1View<T, MemorySpace> xvec, Rank1View<T, MemorySpace> yvec,
  Rank2View<T, MemorySpace> fval)
{
  LO ict_arr[6] = {1, 0, 0, 0, 0, 0};
  Kokkos::View<LO*, MemorySpace> ict("select", 6);
  Kokkos::parallel_for(
    Kokkos::RangePolicy<execution_space>(0, ict.size()),
    KOKKOS_LAMBDA(const LO i) { ict(i) = ict_arr[i]; });
  evaluate(ict, xvec, yvec, fval);
}

template <typename T, typename MemorySpace>
struct InitExplicitBiCubicCoeffFunctor
{
  Rank4View<T, MemorySpace> fspl;
  Rank1View<T, MemorySpace> values;
  LO ny;

  InitExplicitBiCubicCoeffFunctor(Rank4View<T, MemorySpace> fspl_,
                                  Rank1View<T, MemorySpace> values_, LO ny_)
    : fspl(fspl_), values(values_), ny(ny_)
  {
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const LO iy, const LO ix) const
  {
    fspl(0, 0, ix, iy) = values(ix * ny + iy);
  }
};

template <typename T, typename MemorySpace>
struct InitCompactBiCubicCoeffFunctor
{
  Rank3View<T, MemorySpace> fspl;
  Rank1View<T, MemorySpace> values;
  LO ny;

  InitCompactBiCubicCoeffFunctor(Rank3View<T, MemorySpace> fspl_,
                                 Rank1View<T, MemorySpace> values_, LO ny_)
    : fspl(fspl_), values(values_), ny(ny_)
  {
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const LO iy, const LO ix) const
  {
    fspl(0, ix, iy) = values(ix * ny + iy);
  }
};

template <typename T, typename MemorySpace>
ExplicitBiCubicSplineInterpolator<T, MemorySpace>::
  ExplicitBiCubicSplineInterpolator(
    Rank1View<T, MemorySpace> x, Rank1View<T, MemorySpace> y,
    Rank1View<T, MemorySpace> values, BoundaryCondition ibcxmin,
    Rank1View<T, MemorySpace> bcxmin, BoundaryCondition ibcxmax,
    Rank1View<T, MemorySpace> bcxmax, BoundaryCondition ibcymin,
    Rank1View<T, MemorySpace> bcymin, BoundaryCondition ibcymax,
    Rank1View<T, MemorySpace> bcymax)
{
  detail::grid_valid<T, MemorySpace>(x);
  detail::grid_valid<T, MemorySpace>(y);

  x_ = x;
  y_ = y;
  LO ny = y.extent(0);
  LO nx = x.extent(0);
  Kokkos::View<T*, MemorySpace> fspl_view("fspl view", 4 * 4 * nx * ny);
  coefficients_ = fspl_view;
  auto fspl = get_coefficients_view();
  InitExplicitBiCubicCoeffFunctor<T, MemorySpace> init_functor(fspl, values,
                                                               ny);
  Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ny, nx}),
                       init_functor);

  Kokkos::View<T*, MemorySpace> wk_view("wk", 9 * nx * ny);
  auto wk = Rank1View<T, MemorySpace>(wk_view.data(), 9 * nx * ny);

  detail::solve_bicubic_spline_explicit(
    x, nx, y, ny, fspl, static_cast<LO>(ibcxmin), bcxmin,
    static_cast<LO>(ibcxmax), bcxmax, static_cast<LO>(ibcymin), bcymin,
    static_cast<LO>(ibcymax), bcymax, wk);
}

template <typename T, typename MemorySpace>
ExplicitBiCubicSplineInterpolator<T, MemorySpace>::
  ExplicitBiCubicSplineInterpolator(Rank1View<T, MemorySpace> x,
                                    Rank1View<T, MemorySpace> y,
                                    Rank1View<T, MemorySpace> values)
  : ExplicitBiCubicSplineInterpolator(
      x, y, values, BoundaryCondition::NOT_A_KNOT, {},
      BoundaryCondition::NOT_A_KNOT, {}, BoundaryCondition::NOT_A_KNOT, {},
      BoundaryCondition::NOT_A_KNOT, {})
{
  std::cout << "Not a knot boundary condition is applied by default"
            << std::endl;
}

template <typename T, typename MemorySpace>
CompactBiCubicSplineInterpolator<T, MemorySpace>::
  CompactBiCubicSplineInterpolator(
    Rank1View<T, MemorySpace> x, Rank1View<T, MemorySpace> y,
    Rank1View<T, MemorySpace> values, BoundaryCondition ibcxmin,
    Rank1View<T, MemorySpace> bcxmin, BoundaryCondition ibcxmax,
    Rank1View<T, MemorySpace> bcxmax, BoundaryCondition ibcymin,
    Rank1View<T, MemorySpace> bcymin, BoundaryCondition ibcymax,
    Rank1View<T, MemorySpace> bcymax)
{
  detail::grid_valid<T, MemorySpace>(x);
  detail::grid_valid<T, MemorySpace>(y);

  x_ = x;
  y_ = y;
  LO ny = y.extent(0);
  LO nx = x.extent(0);
  Kokkos::View<T*, MemorySpace> fspl_view("fspl view", 4 * nx * ny);
  coefficients_ = fspl_view;
  auto fspl = get_coefficients_view();

  InitCompactBiCubicCoeffFunctor<T, MemorySpace> init_functor(fspl, values, ny);
  Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ny, nx}),
                       init_functor);

  Kokkos::View<T*, MemorySpace> wk_view("wk", 9 * nx * ny);
  auto wk = Rank1View<T, MemorySpace>(wk_view.data(), 9 * nx * ny);

  detail::solve_bicubic_spline_compact(
    x, nx, y, ny, fspl, static_cast<LO>(ibcxmin), bcxmin,
    static_cast<LO>(ibcxmax), bcxmax, static_cast<LO>(ibcymin), bcymin,
    static_cast<LO>(ibcymax), bcymax, wk);
}

template <typename T, typename MemorySpace>
CompactBiCubicSplineInterpolator<T, MemorySpace>::
  CompactBiCubicSplineInterpolator(Rank1View<T, MemorySpace> x,
                                   Rank1View<T, MemorySpace> y,
                                   Rank1View<T, MemorySpace> values)
  : CompactBiCubicSplineInterpolator(
      x, y, values, BoundaryCondition::NOT_A_KNOT, {},
      BoundaryCondition::NOT_A_KNOT, {}, BoundaryCondition::NOT_A_KNOT, {},
      BoundaryCondition::NOT_A_KNOT, {})
{
  std::cout << "Not a knot boundary condition is applied by default"
            << std::endl;
}

} // namespace pcms
#endif
