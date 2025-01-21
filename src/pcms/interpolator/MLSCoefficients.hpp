#ifndef PCMS_INTERPOLATOR_MLS_COEFFICIENTS_HPP
#define PCMS_INTERPOLATOR_MLS_COEFFICIENTS_HPP

#include <cmath>
#include <Omega_h_fail.hpp>

#include <KokkosBatched_Util.hpp>
#include <KokkosBatched_LU_Decl.hpp>
#include <KokkosBatched_LU_Team_Impl.hpp>
#include <KokkosBatched_SolveLU_Decl.hpp>
#include <KokkosBatched_Gemm_Decl.hpp>
#include <KokkosBlas.hpp>
#include "points.hpp"
#include <type_traits>

static constexpr int MAX_DIM = 3;

/**
 * basisSliceLengths, basissSize and BasisPoly are needed to evaluate the
 * polynomial basis for any degree and dimension For instance, polynomial basis
 * vector for dim = 2 and degree = 3 at the point A(x,y) looks like {1, x, y,
 * xx, xy, yy, xxx, xxy, xyy,yyy}. The slices can be written as [1] degree 0 [x]
 * & [y]                degree 1 [xx] &  [xy, yy]         degree 2 [xxx] & [xxy,
 * xyy, yyy]  degree 3
 *
 * the recursive pattern becomes:
 * [1]                      degree 0
 * x*[1] & y*[1]            degree 1
 * x*[x] & y*[x, y]         degree 2
 * xx*[x] & y*[xx,xy,yy]    degree 3
 *
 * lengths of the slices:
 * Degree \ Dim | x | y |
 * =========================
 *        1     | 1 | 1 |
 *        2     | 1 | 2 |
 *        3     | 1 | 3 |
 */

/**
 * @brief computes the slice lengths of the polynomial basis
 * @param array: takes kokkos view array as an input and computes the lengths of
 * slices and fills the array
 *
 * @notes it takes a host array as an input, could have done for device
 * but there was a race conditions for a device
 */
KOKKOS_INLINE_FUNCTION
void basisSliceLengths(Kokkos::View<int**, Kokkos::HostSpace>& array)
{
  int degree = array.extent(0);
  int dim = array.extent(1);

  for (int j = 0; j < dim; ++j) {
    array(0, j) = 1;
  }

  for (int i = 0; i < degree; ++i) {
    array(i, 0) = 1;
  }

  for (int i = 1; i < degree; ++i) {
    for (int j = 1; j < dim; ++j) {
      array(i, j) = array(i, j - 1) + array(i - 1, j);
    }
  }
}

/**
 * @brief computes the size of the polynomial basis vector
 *
 * @param array: it is the array of the slices length
 * @return sum: sum of each element of slice lengths array gives the basis
 * vector size
 *
 * @note it takes the host array
 */

KOKKOS_INLINE_FUNCTION
int basisSize(const Kokkos::View<int**, Kokkos::HostSpace>& array)
{
  int sum = 1;
  int degree = array.extent(0);
  int dim = array.extent(1);

  for (int i = 0; i < degree; ++i) {
    for (int j = 0; j < dim; ++j) {
      sum += array(i, j);
    }
  }

  return sum;
}

/**
 * @brief evaluates the polynomial basis
 *   for example, if it dim = 2 and degree = 2
 *   the basis_vector looks like
 *   basis_vector = {1,x,y,xx,xy,yy}
 *
 *   @param basis_vector: basis_vector 1D array is initialised and passed as
 * input
 *   @param slice_length: slice length is the 2D array that consists of the
 * length of slices
 *   @param p: instance of an class type Coord
 *
 *   @return basis_vector filled with polynomial basis
 */
KOKKOS_INLINE_FUNCTION
void BasisPoly(ScratchVecView& basis_vector, const MatViewType& slice_length,
               Coord& p)
{
  basis_vector(0) = 1;
  int dim = slice_length.extent(1);
  int degree = slice_length.extent(0);

  int prev_col = 0;
  int curr_col = 1;

  double point[MAX_DIM];
  point[0] = p.x;
  point[1] = p.y;

  if (dim == 3) {
    point[2] = p.z;
  }

  for (int i = 0; i < degree; ++i) {
    int offset = curr_col;
    for (int j = 0; j < dim; ++j) {
      for (int k = 0; k < slice_length(i, j); ++k) {
        basis_vector(offset + k) = basis_vector(prev_col + k) * point[j];
      }

      offset += slice_length(i, j);
    }

    prev_col = curr_col;
    curr_col = offset;
  }
}

// create vandermondeMatrix
KOKKOS_INLINE_FUNCTION
void CreateVandermondeMatrix(ScratchMatView vandermonde_matrix,
                             const ScratchMatView& local_source_points, int j,
                             const MatViewType& slice_length)
{
  int N = local_source_points.extent(0);
  int dim = local_source_points.extent(1);

  Coord source_point;
  source_point.x = local_source_points(j, 0);
  source_point.y = local_source_points(j, 1);
  if (dim == 3) {
    source_point.z = local_source_points(j, 2);
  }
  ScratchVecView basis_vector_supports =
    Kokkos::subview(vandermonde_matrix, j, Kokkos::ALL());
  BasisPoly(basis_vector_supports, slice_length, source_point);
}

// radial basis function vector (phi(s_ti))
template <typename Func,
          std::enable_if_t<std::is_invocable_r_v<double, Func, double, double>,
                           bool> = true>
KOKKOS_INLINE_FUNCTION void PhiVector(ScratchVecView Phi,
                                      const Coord target_point,
                                      const ScratchMatView local_source_points,
                                      int j, double cuttoff_dis_sq,
                                      Func rbf_func)
{
  int N = local_source_points.extent(0);
  double dx = target_point.x - local_source_points(j, 0);
  double dy = target_point.y - local_source_points(j, 1);
  double ds_sq = dx * dx + dy * dy;
  Phi(j) = rbf_func(ds_sq, cuttoff_dis_sq);
  OMEGA_H_CHECK_PRINTF(!std::isnan(Phi(j)),
                       "ERROR: Phi(j) in PhiVector is NaN for j = %d "
                       "ds_sq=%.16f, cuttoff_dis_sq=%.16f",
                       j, ds_sq, cuttoff_dis_sq);
}

// takes a matrix A(m,n) and vector b(m) and calulates
// A^Tb --> A^T(n,m) b(m)
KOKKOS_INLINE_FUNCTION
void ScaleColumnTransMatrix(ScratchMatView result_matrix,
                            const ScratchMatView& matrix,
                            const ScratchVecView& vector, member_type team,
                            int j)
{

  int M = matrix.extent(0);
  int N = matrix.extent(1);

  ScratchVecView matrix_row = Kokkos::subview(matrix, j, Kokkos::ALL());
  for (int k = 0; k < N; k++) {
    OMEGA_H_CHECK_PRINTF(!std::isnan(matrix_row(k)),
                         "ERROR: given matrix is NaN for k = %d\n", k);
    OMEGA_H_CHECK_PRINTF(!std::isnan(vector(j)),
                         "ERROR: given vector is NaN for j = %d\n", j);
    result_matrix(k, j) = matrix_row(k) * vector(j);
    OMEGA_H_CHECK_PRINTF(!std::isnan(result_matrix(k, j)),
                         "ERROR: result_matrix is NaN for k = %d, j = %d\n", k,
                         j);
  }
}

// dot product
KOKKOS_INLINE_FUNCTION
void dot_product(member_type team, const ScratchVecView& result_sub,
                 const ScratchVecView& SupportValues_sub, double& target_value)
{
  int N = result_sub.extent(0);
  for (int j = 0; j < N; ++j) {
    target_value += result_sub(j) * SupportValues_sub[j];
  }
  OMEGA_H_CHECK_PRINTF(!std::isnan(target_value),
                       "ERROR: NaN found in dot_product: N: %d\n", N);
}

// TODO: Implement QR decomposition to solve the linear system
// convert normal equation P^T Q P x = P^T Q b to Ax = b';
// A = P^T Q P & b' = P^T Q b

/**
 * @struct ResultConvertNormal
 * @brief Represents the results from P^T Q(scaled matrix), P^T Q P (square
 * matrix) and P^T Q b (transformed rhs)
 *
 * scaled_matrix is the matrix obtained after scaling the column of the given
 * matrix square_matrix is the matrix obtained after P^T Q P operations
 * transformed_rhs is obtained after P^T b operation
 *
 * This class is used to store the scaled matrix, square matrix and transformed
 * rhs
 */

struct ResultConvertNormal
{
  ScratchMatView scaled_matrix;
  ScratchMatView square_matrix;
  ScratchVecView transformed_rhs;
};

KOKKOS_INLINE_FUNCTION
ResultConvertNormal ConvertNormalEq(const ScratchMatView matrix,
                                    const ScratchVecView weight_vector,
                                    const ScratchVecView rhs, member_type team)
{

  int m = matrix.extent(0);
  int n = matrix.extent(1);

  ResultConvertNormal result;

  result.scaled_matrix = ScratchMatView(team.team_scratch(0), n, m);

  result.square_matrix = ScratchMatView(team.team_scratch(0), n, n);

  result.transformed_rhs = ScratchVecView(team.team_scratch(0), n);

  // performing P^T Q
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, n), [=](int j) {
    for (int k = 0; k < m; ++k) {
      result.scaled_matrix(j, k) = 0;
    }
  });

  team.team_barrier();

  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, m), [=](int j) {
    ScaleColumnTransMatrix(result.scaled_matrix, matrix, weight_vector, team,
                           j);
  });

  team.team_barrier();

  // performing (P^T Q P)
  KokkosBatched::TeamGemm<
    member_type, KokkosBatched::Trans::NoTranspose,
    KokkosBatched::Trans::NoTranspose,
    KokkosBatched::Algo::Gemm::Unblocked>::invoke(team, 1.0,
                                                  result.scaled_matrix, matrix,
                                                  0.0, result.square_matrix);

  team.team_barrier();

  KokkosBlas::Experimental::
    Gemv<KokkosBlas::Mode::Team, KokkosBlas::Algo::Gemv::Unblocked>::invoke(
      team, 'N', 1.0, result.scaled_matrix, rhs, 0.0, result.transformed_rhs);

  team.team_barrier();

  return result;
}

/**
 * @brief Solves the matrix equation Ax = b
 *
 * This function uses LU decomposition to solve the matrix equation
 *
 * @param square_matrix The given matrix A (must be square)
 * @param rhs The known vector b
 * @return rhs The unknown vector x after solving matrix equation is overwritten
 * in rhs
 *
 */
KOKKOS_INLINE_FUNCTION
void SolveMatrix(const ScratchMatView square_matrix, const ScratchVecView rhs,
                 member_type team)
{

  // Perform LU decomposition
  KokkosBatched::TeamLU<
    member_type, KokkosBatched::Algo::LU::Unblocked>::invoke(team,
                                                             square_matrix);

  team.team_barrier();

  //   Solve the equation with forward and backward solves
  KokkosBatched::TeamSolveLU<
    member_type, KokkosBatched::Trans::NoTranspose,
    KokkosBatched::Algo::SolveLU::Unblocked>::invoke(team, square_matrix, rhs);
}

#endif
