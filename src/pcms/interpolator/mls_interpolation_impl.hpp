#ifndef PCMS_INTERPOLATOR_MLS_INTERPOLATION_IMP_HPP
#define PCMS_INTERPOLATOR_MLS_INTERPOLATION_IMP_HPP

#include <cmath>
#include <Omega_h_fail.hpp>
#include <type_traits>

#include <KokkosBatched_Util.hpp>
#include <KokkosBatched_LU_Decl.hpp>
#include <KokkosBatched_LU_Team_Impl.hpp>
#include <KokkosBatched_SolveLU_Decl.hpp>
#include <KokkosBlas.hpp>
#include <KokkosBlas1_team_dot.hpp>
#include <Kokkos_StdAlgorithms.hpp>
#include "pcms_interpolator_aliases.hpp"
#include "adj_search.hpp"
#include "pcms/assert.h"
#include "pcms/profile.h"
#include "pcms_interpolator_view_utils.hpp"
#include "pcms_interpolator_logger.hpp"

#include <KokkosBlas2_gemv.hpp> //KokkosBlas::gemv
#include "KokkosBatched_SVD_Decl.hpp"
#include "KokkosBatched_SVD_Serial_Impl.hpp"
#include <KokkosBlas2_serial_gemv_impl.hpp>
#include <KokkosBatched_Gemm_Decl.hpp>

static constexpr int MAX_DIM = 6;

/**
 * calculate_basis_slice_lengths, calculate_basis_vector_size and
 * eval_basis_vector are needed to evaluate the polynomial basis for any degree
 * and dimension For instance, polynomial basis vector for dim = 2 and degree =
 * 3 at the point (x,y) looks like {1, x, y, xx, xy, yy, xxx, xxy, xyy,yyy}. The
 * slices can be written as [1]                      degree 0 [x] & [y] degree 1
 * [xx] &  [xy, yy]         degree 2
 * [xxx] & [xxy,xyy, yyy]   degree 3
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

namespace pcms
{

namespace detail
{
/**
 * @brief Computes the slice lengths of the polynomial basis
 *
 * This function takes kokkos view array as an input and computes the lengths of
 * slices and fills the array
 *
 * @param[in,out] array The 2D array
 *
 * @note It takes a host array as an input, could have done for device
 * but there was a race conditions for a device
 */
void calculate_basis_slice_lengths(IntHostMatView& array);

/**
 * @brief Computes the size of the polynomial basis vector
 *
 * @param array The 2D array of the slices length
 * @return The basis vector size
 *
 * @note It takes the host array
 */
int calculate_basis_vector_size(const IntHostMatView& array);

/**
 * @brief Calculates the scratch size required
 *
 * This function uses the size of the vectors and matrices and
 * calculates the scratch size required in each team and extracts the
 * maximum shared size
 *
 * @param None
 *
 * @return The scratch size
 */
int calculate_scratch_shared_size(const SupportResults& support,
                                  const int ntargets, int basis_size, int dim);
/**
 *
 * @brief Performs minmax normalization
 *
 * This function takes coordinates and dimension
 * and outputs the normalized coordinates
 *
 * @param coordinates The Coordinates
 * @param dim The Dimension
 *
 * @return normalized coordinates The normlaised coordinates
 */
Reals min_max_normalization(Reals& coordinates, int dim);

/**
 * @brief Evaluates the polynomial basis
 *
 *   for example, if it dim = 2 and degree = 2, the basis_vector looks like
 *   basis_vector = {1,x,y,xx,xy,yy}
 *
 *   @param[in] slice_length A 2D array of slices length
 *   @param[in] p A reference to the coordinate struct
 *   @param[in,out] basis_vector The polynomial basis vector
 *
 *
 */
KOKKOS_INLINE_FUNCTION
void eval_basis_vector(const IntDeviceMatView& slice_length, const double* p,
                       ScratchVecView& basis_vector)
{
  basis_vector(0) = 1;
  int dim = slice_length.extent(1);
  int degree = slice_length.extent(0);

  int prev_col = 0;
  int curr_col = 1;

  double point[MAX_DIM];

  for (int i = 0; i < dim; ++i) {
    point[i] = p[i];
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

/**
 *  @Brief Normalizes the support coordinates about the given target point
 *
 *  This function takes the pivot point Coord object
 *  and support coordinates as input and gives normalized coordinates
 *
 *  @param[in] team The team member
 *  @param[in, out] pivot The Coord object that stores the coordinate of the
 *target pivot
 *  @param[in,out] support_coordinates The orginal coordinates when in, when out
 *  		normalized coordinates
 **
 */
KOKKOS_INLINE_FUNCTION
void normalize_supports(const member_type& team, double* target_point,
                        ScratchMatView& support_coordinates)
{
  int nsupports = support_coordinates.extent(0);
  int dim = support_coordinates.extent(1);

  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nsupports), [=](int i) {
    Real pivot_point[MAX_DIM];

    for (int j = 0; j < dim; ++j) {
      support_coordinates(i, j) -= target_point[j];
    }
  });

  for (int j = 0; j < dim; ++j) {
    target_point[j] = 0.0;
  }
}

/**
 * @brief Creates a vandermonde matrix
 *
 * @param[in] local_source_points The coordinates of the local supports
 * @param[in] j The jth row
 * @param[in] slice_length The slice lengths of the polynomial basis
 * @param[in,out] vandermonde_matrix The Vandermonde Matrix
 *
 */
KOKKOS_INLINE_FUNCTION
void create_vandermonde_matrix(const ScratchMatView& local_source_points, int j,
                               const IntDeviceMatView& slice_length,
                               ScratchMatView vandermonde_matrix)
{
  int N = local_source_points.extent(0);
  int dim = local_source_points.extent(1);

  double source_point[MAX_DIM] = {};
  for (int i = 0; i < dim; ++i) {
    source_point[i] = local_source_points(j, i);
  }

  ScratchVecView basis_vector_supports =
    Kokkos::subview(vandermonde_matrix, j, Kokkos::ALL());
  eval_basis_vector(slice_length, source_point, basis_vector_supports);
}

/**
 *@brief Computes basis function vector
 *
 *This function takes a basis function the user wants and calculates the value
 *for each local source supports
 *
 * @tparam Func The type of the radial basis function
 * @param[in] target_point The coordinate of the target point where the
 *interpolation is carried out
 * @param[in] source_points The coordinates of the source support points
 * @param[in] j The jth row
 * @param[in] cutoff_dis_sq The cutoff radius squared
 * @param[in] rbf_func The radial basis function of choice
 * @param[in, out] phi A radial basis function value
 *
 */
template <typename Func,
          std::enable_if_t<std::is_invocable_r_v<double, Func, double, double>,
                           bool> = true>
KOKKOS_INLINE_FUNCTION void compute_phi_vector(
  const double* target_point, const ScratchMatView& local_source_points, int j,
  double cuttoff_dis_sq, Func rbf_func, ScratchVecView phi)
{
  int N = local_source_points.extent(0);
  int dim = local_source_points.extent(1);

  double ds_sq = 0;

  for (int i = 0; i < dim; ++i) {
    double temp = target_point[i] - local_source_points(j, i);
    ds_sq += temp * temp;
  }
  phi(j) = rbf_func(ds_sq, cuttoff_dis_sq);
  OMEGA_H_CHECK_PRINTF(!std::isnan(phi(j)),
                       "ERROR: Phi(j) in compute_phi_vector is NaN for j = %d "
                       "ds_sq=%.16f, cuttoff_dis_sq=%.16f",
                       j, ds_sq, cuttoff_dis_sq);
}

/**
 * @brief Scales the column of the transpose of the given matrix
 *
 * This function takes a matrix P(m,n) and vector phi(m) and scales the
 * each column of the transpose of P by each corresponding element in vector
 * phi
 *
 * @param[in] matrix The matrix
 * @param[in] vector The vector
 * @param[in] team The team member
 * @param[in, out] result_matrix
 *
 */
KOKKOS_INLINE_FUNCTION
void scale_column_trans_matrix(const ScratchMatView& matrix,
                               const ScratchVecView& vector, member_type team,
                               int j, ScratchMatView result_matrix)
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

KOKKOS_INLINE_FUNCTION
void add_regularization(const member_type& team, ScratchMatView& square_matrix,
                        Real lambda_factor)
{
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, square_matrix.extent(0)),
                       [=](int i) { square_matrix(i, i) += lambda_factor; });
}

/**
 * @struct ResultConvertNormal
 * @brief Stores the results of matrix and vector transformations.
 *
 * This struct represents the results of the following operations:
 * - `P^T Q` (scaled matrix): The matrix obtained after scaling the columns of
 * the given matrix.
 * - `P^T Q P` (square matrix): The matrix obtained after performing the `P^T Q
 * P` operation.
 * - `P^T b` (transformed RHS): The vector obtained after applying the `P^T b`
 * operation.
 *
 * The struct is used to store:
 * - The scaled matrix (`scaled_matrix`).
 * - The square matrix (`square_matrix`).
 * - The transformed right-hand side vector (`transformed_rhs`).
 */
struct ResultConvertNormal
{
  ScratchMatView scaled_matrix;
  ScratchMatView square_matrix;
  ScratchVecView transformed_rhs;
};

/**
 * @brief Converts a normal equation into a simplified form.
 *
 * This function takes a matrix, a basis/weight vector, and a right-hand side
 * (RHS) vector and computes the scaled matrix, square matrix, and transformed
 * RHS.
 *
 * For example, the normal equation `P^T Q P x = P^T Q b` is converted into the
 * form `Ax = b'`, where:
 * - `A` (square matrix) = `P^T Q P`
 * - `b'` (transformed RHS) = `P^T Q b`
 *
 * @param matrix The input matrix.
 * @param weight_vector The weight or basis vector.
 * @param rhs The right-hand side vector of the equation.
 * @param team The team member executing the task (optional, if applicable).
 *
 * @return The result of the normal equation conversion as a struct containing:
 * - The scaled matrix.
 * - The square matrix.
 * - The transformed RHS.
 *
 * @todo Implement QR decomposition to solve the linear system.
 */
KOKKOS_INLINE_FUNCTION
ResultConvertNormal convert_normal_equation(const ScratchMatView& matrix,
                                            const ScratchVecView& weight_vector,
                                            const ScratchVecView& rhs,
                                            member_type team,
                                            double lambda_factor)
{

  int m = matrix.extent(0);
  int n = matrix.extent(1);

  ScratchMatView scaled_matrix(team.team_scratch(1), n, m);

  ScratchMatView square_matrix(team.team_scratch(1), n, n);

  ScratchVecView transformed_rhs(team.team_scratch(1), n);

  // performing P^T Q
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, n), [=](int j) {
    for (int k = 0; k < m; ++k) {
      scaled_matrix(j, k) = 0;
    }
  });

  team.team_barrier();

  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, m), [=](int j) {
    scale_column_trans_matrix(matrix, weight_vector, team, j, scaled_matrix);
  });

  team.team_barrier();

  // performing (P^T Q P)
  KokkosBatched::TeamGemm<
    member_type, KokkosBatched::Trans::NoTranspose,
    KokkosBatched::Trans::NoTranspose,
    KokkosBatched::Algo::Gemm::Unblocked>::invoke(team, 1.0, scaled_matrix,
                                                  matrix, 0.0, square_matrix);

  team.team_barrier();

  add_regularization(team, square_matrix, lambda_factor);

  team.team_barrier();

  KokkosBlas::Experimental::
    Gemv<KokkosBlas::Mode::Team, KokkosBlas::Algo::Gemv::Unblocked>::invoke(
      team, 'N', 1.0, scaled_matrix, rhs, 0.0, transformed_rhs);

  team.team_barrier();

  return ResultConvertNormal{scaled_matrix, square_matrix, transformed_rhs};
}

/**
 * @brief Solves the matrix equation Ax = b' using LU decomposition.
 *
 * This function solves the given matrix equation using LU decomposition.
 * The solution vector `x'` overwrites the input `rhs` vector.
 *
 * @param square_matrix The input matrix A (must be square).
 * @param rhs The right-hand side vector b. After execution, it is overwritten
 * with the solution vector x.
 * @param team The team member executing the task.
 *
 * @return None. The solution is directly written to the `rhs` vector.
 */
KOKKOS_INLINE_FUNCTION
void solve_matrix(const ScratchMatView& square_matrix,
                  const ScratchVecView& rhs, member_type team)
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

KOKKOS_INLINE_FUNCTION
void solve_matrix_svd(member_type team, const ScratchVecView& weight,
                      ScratchVecView rhs_values, ScratchMatView matrix,
                      ScratchVecView solution_vector, double lambda)
{

  int row = matrix.extent(0);

  int column = matrix.extent(1);

  int weight_size = weight.size();

  OMEGA_H_CHECK_PRINTF(
    weight_size == row,
    "the size of the weight vector should be equal to the row of the matrix\n"
    "weight vector size = %d, row of matrix = %d\n",
    weight_size, row);

  eval_row_scaling(team, weight, matrix);

  team.team_barrier();

  eval_rhs_scaling(team, weight, rhs_values);
  team.team_barrier();

  // initilize U (orthogonal matrices)  array
  ScratchMatView U(team.team_scratch(1), row, row);
  fill(-5.0, team, U);

  // initialize Vt
  ScratchMatView Vt(team.team_scratch(1), column, column);
  fill(-5.0, team, Vt);

  // initialize sigma
  ScratchVecView sigma(team.team_scratch(1), column);
  fill(-5.0, team, sigma);

  // initialize work
  ScratchVecView work(team.team_scratch(1), row);
  fill(-5.0, team, work);

  // initialise sigma_mat_inv  S^-1
  ScratchMatView temp_matrix(team.team_scratch(1), column, row);
  fill(0, team, temp_matrix);

  // initialise V S^-1 U^T
  ScratchMatView vsigmaInvUtMul(team.team_scratch(1), column, row);
  fill(0, team, vsigmaInvUtMul);

  if (team.team_rank() == 0) {
    KokkosBatched::SerialSVD::invoke(KokkosBatched::SVD_USV_Tag(), matrix, U,
                                     sigma, Vt, work, 1e-6);
  }
  team.team_barrier();

  calculate_shrinkage_factor(team, lambda, sigma);

  auto Ut = find_transpose(team, U);

  scale_and_adjust(team, sigma, Ut, temp_matrix); // S^-1 U^T

  KokkosBatched::TeamGemm<
    member_type, KokkosBatched::Trans::Transpose,
    KokkosBatched::Trans::NoTranspose,
    KokkosBatched::Algo::Gemm::Unblocked>::invoke(team, 1.0, Vt, temp_matrix,
                                                  0.0, vsigmaInvUtMul);

  KokkosBlas::TeamGemv<
    member_type, KokkosBlas::Trans::NoTranspose,
    KokkosBlas::Algo::Gemv::Unblocked>::invoke(team, 1.0, vsigmaInvUtMul,
                                               rhs_values, 0.0,
                                               solution_vector);
}

/**
 * @brief Maps the data from source mesh to target mesh
 *
 * @param source_values Scalar array view of source field values
 * @param source_coordinates Scalar array view of the coordinates of control
 * points of source field
 * @param target_coordinates Scalar array view of the coordinates of control
 * points of target field
 * @param support The object that encapsulates support info
 * @param dim The dimension of the simulations
 * @param degree The degree of the interpolation order
 * @param rbf_func The radial basis function choice
 * @param[in/out] Scalar array view of the interpolated field in target mesh
 *
 *
 */
template <typename Func>
void mls_interpolation(RealConstDefaultScalarArrayView source_values,
                       RealConstDefaultScalarArrayView source_coordinates,
                       RealConstDefaultScalarArrayView target_coordinates,
                       const SupportResults& support, const LO& dim,
                       const LO& degree, Func rbf_func,
                       RealDefaultScalarArrayView approx_target_values,
                       double lambda)
{
  PCMS_FUNCTION_TIMER;
  static_assert(std::is_invocable_r_v<double, Func, double, double>,
                "rbf_func, takes radius and cutoff, returns weight");
  static_assert(!std::is_pointer_v<Func>,
                "function pointer will fail in GPU execution context");

  int nsources = source_coordinates.size() / dim;

  OMEGA_H_CHECK_PRINTF(
    source_values.size() == nsources,
    "[ERROR] The size of the source values and source coordinates is not "
    "same. "
    "The current sizes are :\n"
    "source_values_size = %d, source_coordinates_size = %d\n",
    source_values.size(), nsources);

  const auto ntargets = target_coordinates.size() / dim;

  OMEGA_H_CHECK_PRINTF(approx_target_values.size() == ntargets,
                       "[ERROR] The size of the approx target values and the "
                       "number of targets is "
                       "not same. The current numbers are :\n"
                       "approx_target_values = %d, ntargets = %d\n",
                       approx_target_values.size(), ntargets);

  IntHostMatView host_slice_length(
    "stores slice length of  polynomial basis in host", degree, dim);

  Kokkos::deep_copy(host_slice_length, 0);

  calculate_basis_slice_lengths(host_slice_length);

  auto basis_size = calculate_basis_vector_size(host_slice_length);

  IntDeviceMatView slice_length(
    "stores slice length of polynomial basis in device", degree, dim);

  auto slice_length_hd = Kokkos::create_mirror_view(slice_length);
  Kokkos::deep_copy(slice_length_hd, host_slice_length);
  Kokkos::deep_copy(slice_length, slice_length_hd);

  int shared_size =
    calculate_scratch_shared_size(support, ntargets, basis_size, dim);

  team_policy tp(ntargets, Kokkos::AUTO);

  int scratch_size = tp.scratch_size_max(1);
  printf("Scratch Size = %d\n", scratch_size);
  printf("Shared Size = %d\n", shared_size);
  PCMS_ALWAYS_ASSERT(scratch_size > shared_size);

  // calculates the interpolated values
  Kokkos::parallel_for(
    "MLS coefficients", tp.set_scratch_size(1, Kokkos::PerTeam(scratch_size)),
    KOKKOS_LAMBDA(const member_type& team) {
      int league_rank = team.league_rank();
      int start_ptr = support.supports_ptr[league_rank];
      int end_ptr = support.supports_ptr[league_rank + 1];
      int nsupports = end_ptr - start_ptr;

      //  local_source_point stores the coordinates of source supports of a
      //  given target
      ScratchMatView local_source_points(team.team_scratch(1), nsupports, dim);

      // rbf function values of source supports Phi(n,n)
      ScratchVecView phi_vector(team.team_scratch(1), nsupports);

      //  vondermonde matrix P from the vectors of basis vector of supports
      ScratchMatView vandermonde_matrix(team.team_scratch(1), nsupports,
                                        basis_size);
      // stores known vector (b)
      ScratchVecView support_values(team.team_scratch(1), nsupports);

      // basis of target
      ScratchVecView target_basis_vector(team.team_scratch(1), basis_size);

      // solution coefficients (solution vector)
      ScratchVecView solution_coefficients(team.team_scratch(1), basis_size);

      // Initialize the scratch matrices and  vectors
      fill(0.0, team, local_source_points);
      fill(0.0, team, vandermonde_matrix);
      fill(0.0, team, phi_vector);
      fill(0.0, team, support_values);
      fill(0.0, team, target_basis_vector);
      fill(0.0, team, solution_coefficients);

      Logger logger(10);
      // storing the coords of local supports
      int count = -1;
      for (int j = start_ptr; j < end_ptr; ++j) {
        count++;
        auto index = support.supports_idx[j];

        for (int i = 0; i < dim; ++i) {
          local_source_points(count, i) = source_coordinates[index * dim + i];
        }
      }

      logger.logMatrix(team, LogLevel::DEBUG, local_source_points,
                       "Support Coordinates");
      double target_point[MAX_DIM] = {};

      for (int i = 0; i < dim; ++i) {
        target_point[i] = target_coordinates[league_rank * dim + i];
      }

      logger.logArray(team, LogLevel::DEBUG, target_point, dim,
                      "Target points");

      /** phi(nsupports) is the array of rbf functions evaluated at the
       * source supports In the actual implementation, Phi(nsupports,
       * nsupports) is the diagonal matrix & each diagonal element is the phi
       * evaluated at each source points
       *
       * step 1: evaluate the phi vector with the original dimension
       */

      Kokkos::parallel_for(
        Kokkos::TeamThreadRange(team, nsupports), [=](int j) {
          compute_phi_vector(target_point, local_source_points, j,
                             support.radii2[league_rank], rbf_func, phi_vector);
        });

      team.team_barrier();

      /** support_values(nsupports) (or known rhs vector b) is the vector of
       * the quantity that we want interpolate
       *
       *
       * step 4: find local supports function values
       */
      Kokkos::parallel_for(
        Kokkos::TeamThreadRange(team, nsupports), [=](const int j) {
          support_values(j) =
            source_values[support.supports_idx[start_ptr + j]];
          OMEGA_H_CHECK_PRINTF(!std::isnan(support_values(j)),
                               "ERROR: NaN found: at support %d\n", j);
        });

      team.team_barrier();

      logger.log(team, LogLevel::DEBUG, "The search  starts");
      logger.logVector(team, LogLevel::DEBUG, support_values, "Support values");

      /**
       *
       * the local_source_points is of the type ScratchMatView with
       * coordinates information; row is number of local supports
       * & column is dim
       *
       * step 2: normalize local source supports and target point
       */

      normalize_supports(team, target_point, local_source_points);
      team.team_barrier();
      /**
       *
       * evaluates the basis vector of a given target point Coord target_point;
       * this can evaluate monomial basis vector for any degree of polynomial
       * step 3: call basis vector evaluation function (eval_basis_vector);
       */
      eval_basis_vector(slice_length, target_point, target_basis_vector);

      /** vandermonde_matrix(nsupports, basis_size) vandermonde Matrix is
       * created with the basis vector of source supports stacking on top of
       * each other
       *
       * step 4: create vandermonde matrix
       */
      Kokkos::parallel_for(
        Kokkos::TeamThreadRange(team, nsupports), [=](int j) {
          create_vandermonde_matrix(local_source_points, j, slice_length,
                                    vandermonde_matrix);
        });
>>>>>>> mls_interpolation_svd_solver

      team.team_barrier();

      logger.logMatrix(team, LogLevel::DEBUG, vandermonde_matrix,
                       "vandermonde matrix");

      OMEGA_H_CHECK_PRINTF(
        support.radii2[league_rank] > 0,
        "ERROR: radius2 has to be positive but found to be %.16f\n",
        support.radii2[league_rank]);

      solve_matrix_svd(team, phi_vector, support_values, vandermonde_matrix,
                       solution_coefficients, lambda);
      team.team_barrier();

      double target_value = KokkosBlas::Experimental::dot(
        team, solution_coefficients, target_basis_vector);
      logger.logScalar(team, LogLevel::DEBUG, target_value,
                       "interpolated value");
      if (team.team_rank() == 0) {
        OMEGA_H_CHECK_PRINTF(!std::isnan(target_value), "Nan at %d\n",
                             league_rank);
        approx_target_values[league_rank] = target_value;
      }
    });
}

/**
 * @brief \overload
 * Maps the data from source mesh to target mesh
 *
 * @param source_values Read array of Source field values
 * @param source_coordinates Read array of the coordinates of control points
 * of source field
 * @param target_coordinates Read array of the coordinates of control points
 * of target field
 * @param support The object that enpasulates support info
 * @param dim The dimension of the simulations
 * @param degree The degree of the interpolation order
 * @param rbf_func The radial basis function choice
 * @return  Write array of the interpolated field in target mesh
 *
 *
 */
template <typename Func>
Write<Real> mls_interpolation(const Reals source_values,
                              const Reals source_coordinates,
                              const Reals target_coordinates,
                              const SupportResults& support, const LO& dim,
                              const LO& degree, Func rbf_func, double lambda)
{
  const auto nsources = source_coordinates.size() / dim;

  const auto ntargets = target_coordinates.size() / dim;

  RealConstDefaultScalarArrayView source_values_array_view(
    source_values.data(), source_values.size());

  RealConstDefaultScalarArrayView source_coordinates_array_view(
    source_coordinates.data(), source_coordinates.size());

  RealConstDefaultScalarArrayView target_coordinates_array_view(
    target_coordinates.data(), target_coordinates.size());

  RealDefaultScalarArrayView radii2_array_view(support.radii2.data(),
                                               support.radii2.size());

  Write<Real> interpolated_values(ntargets, 0, "approximated target values");

  RealDefaultScalarArrayView interpolated_values_array_view(
    interpolated_values.data(), interpolated_values.size());

  mls_interpolation(source_values_array_view, source_coordinates_array_view,
                    target_coordinates_array_view, support, dim, degree,
                    rbf_func, interpolated_values_array_view, lambda);

  return interpolated_values;
}

} // namespace detail
} // namespace pcms
#endif
