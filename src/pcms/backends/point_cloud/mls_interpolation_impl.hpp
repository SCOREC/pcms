#ifndef PCMS_FIELD_EVALUATOR_MLS_INTERPOLATION_IMPL_HPP
#define PCMS_FIELD_EVALUATOR_MLS_INTERPOLATION_IMPL_HPP

#include <cmath>
#include <type_traits>

#include <KokkosBatched_Util.hpp>
#include <KokkosBatched_LU_Decl.hpp>
#include <KokkosBatched_LU_Team_Impl.hpp>
#include <KokkosBatched_SolveLU_Decl.hpp>
#include <KokkosBlas.hpp>
#include <KokkosBlas1_team_dot.hpp>
#include <Kokkos_StdAlgorithms.hpp>
#include <pcms/backends/point_cloud/mls_interpolation.hpp>
#include <pcms/backends/point_cloud/interpolator_aliases.hpp>
#include <pcms/utility/check.hpp>
#include <pcms/utility/assert.hpp>
#include <pcms/utility/detail/profile.hpp>
#include <pcms/backends/point_cloud/interpolator_view_utils.hpp>
#include <pcms/backends/point_cloud/interpolator_logger.hpp>

#include <KokkosBlas2_gemv.hpp>
#include <KokkosBatched_SVD_Decl.hpp>
#include <KokkosBatched_SVD_Serial_Impl.hpp>
#include <KokkosBlas2_serial_gemv_impl.hpp>
#include <KokkosBatched_Gemm_Decl.hpp>

static constexpr int MAX_DIM = 6;

/**
 * calculate_basis_slice_lengths, calculate_basis_vector_size and
 * eval_basis_vector are needed to evaluate the polynomial basis for any degree
 * and dimension For instance, polynomial basis vector for dim = 2 and degree =
 * 3 at the point (x,y) looks like {1, x, y, xx, xy, yy, xxx, xxy, xyy,yyy}. The
 * slices can be written as
 * [1]                      degree 0
 * [x] & [y]         degree 1
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

namespace pcms::field
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
 *      normalized coordinates
 **
 */
KOKKOS_INLINE_FUNCTION
void normalize_supports(const member_type& team, double* target_point,
                        ScratchMatView& support_coordinates)
{
  int nsupports = support_coordinates.extent(0);
  int dim = support_coordinates.extent(1);

  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nsupports), [=](int i) {
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
  int dim = local_source_points.extent(1);

  double ds_sq = 0;
  for (int i = 0; i < dim; ++i) {
    double temp = target_point[i] - local_source_points(j, i);
    ds_sq += temp * temp;
  }
  phi(j) = rbf_func(ds_sq, cuttoff_dis_sq);
  PCMS_CHECK_PRINTF(!std::isnan(phi(j)),
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
                               const ScratchVecView& vector,
                               member_type /*unused*/, int j,
                               ScratchMatView result_matrix)
{
  int N = matrix.extent(1);

  ScratchVecView matrix_row = Kokkos::subview(matrix, j, Kokkos::ALL());
  for (int k = 0; k < N; k++) {
    PCMS_CHECK_PRINTF(!std::isnan(matrix_row(k)),
                      "ERROR: given matrix is NaN for k = %d\n", k);
    PCMS_CHECK_PRINTF(!std::isnan(vector(j)),
                      "ERROR: given vector is NaN for j = %d\n", j);
    result_matrix(k, j) = matrix_row(k) * vector(j);
    PCMS_CHECK_PRINTF(!std::isnan(result_matrix(k, j)),
                      "ERROR: result_matrix is NaN for k = %d, j = %d\n", k, j);
  }
}

/**
 * @brief Solves a regularized weighted linear system WAx = Wb using SVD
 * decomposition within a Kokkos parallel team.
 *
 * This function solves a (possibly overdetermined) linear system using Singular
 * Value Decomposition (SVD), The regularization term (lambda) helps stabilize
 * the solution especially in ill-conditioned cases.
 *
 *
 * @param team The Kokkos team member executing this function.
 *
 * @param weight A scratch-space vector of RBF weights (diagonal elements of
 * weight matrix W). size is m.
 *
 * @param rhs_values A scratch-space vector representing the right-hand side
 * vector b (m). It may be modified during the computation.
 *
 * @param matrix A scratch-space matrix representing the system matrix A (m x
 * n).
 *
 * @param solution_vector The output vector x (n). It will contain the computed
 * solution after the function completes.
 *
 * @param lambda The regularization coefficient. Set to 0.0 for pure SVD; use a
 * positive value for Tikhonov regularization.
 *
 * @param tol  The  tolerance threshold used to discard small singular values in
 * the decomposition
 *
 *
 * @note This function does not return a value; the solution is written directly
 * to `solution_vector`.
 */
KOKKOS_INLINE_FUNCTION
void pseudo_inverse_svd(member_type team, const ScratchVecView& weight,
                        ScratchMatView matrix, ScratchMatView pseudo_inverse,
                        double lambda, double tol)
{
  int row = matrix.extent(0);
  int column = matrix.extent(1);
  int weight_size = weight.size();

  PCMS_CHECK_PRINTF(
    weight_size == row,
    "the size of the weight vector should be equal to the row of the matrix\n"
    "weight vector size = %d, row of matrix = %d\n",
    weight_size, row);

  eval_row_scaling(team, weight, matrix);
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

  fill(0, team, pseudo_inverse);

  if (team.team_rank() == 0) {
    KokkosBatched::SerialSVD::invoke(KokkosBatched::SVD_USV_Tag(), matrix, U,
                                     sigma, Vt, work, tol);
  }
  team.team_barrier();

  calculate_shrinkage_factor(team, lambda, sigma);
  auto Ut = find_transpose(team, U);

  scale_and_adjust(team, sigma, Ut, temp_matrix); // S^-1 U^T

  KokkosBatched::TeamGemm<
    member_type, KokkosBatched::Trans::Transpose,
    KokkosBatched::Trans::NoTranspose,
    KokkosBatched::Algo::Gemm::Unblocked>::invoke(team, 1.0, Vt, temp_matrix,
                                                  0.0, pseudo_inverse);
}

KOKKOS_INLINE_FUNCTION
void solve_matrix_svd(member_type team, const ScratchVecView& weight,
                      ScratchVecView rhs_values, ScratchMatView matrix,
                      ScratchVecView solution_vector, double lambda, double tol)
{
  ScratchMatView vsigmaInvUtMul(team.team_scratch(1), matrix.extent(1),
                                matrix.extent(0));
  pseudo_inverse_svd(team, weight, matrix, vsigmaInvUtMul, lambda, tol);
  eval_rhs_scaling(team, weight, rhs_values);
  team.team_barrier();

  KokkosBlas::TeamGemv<
    member_type, KokkosBlas::Trans::NoTranspose,
    KokkosBlas::Algo::Gemv::Unblocked>::invoke(team, 1.0, vsigmaInvUtMul,
                                               rhs_values, 0.0,
                                               solution_vector);
}

/**
 * @brief The MLS coefficient rows: for each target point, the weights c_j with
 * which the MLS value at the target is sum_j c_j source_value[idx_j] over its
 * supports (the CSR rows of the supports): c_j = phi_j sum_k t_k M_kj, where
 * M = V S^-1 U^T is the (regularized) pseudo-inverse of the weighted
 * Vandermonde matrix and t is the polynomial basis at the target.
 *
 * @param source_coordinates Scalar array view of the coordinates of the source
 * points
 * @param target_coordinates Scalar array view of the coordinates of the target
 * points
 * @param support The object that encapsulates support info
 * @param dim The dimension of the simulations
 * @param degree The degree of the interpolation order
 * @param rbf_func The radial basis function choice
 * @param[out] coefficients one entry per support entry (support.supports_idx);
 * zero for targets without supports
 */
template <typename Func>
void mls_coefficients(RealConstDefaultRank1View source_coordinates,
                      RealConstDefaultRank1View target_coordinates,
                      const SupportResults& support, const LO& dim,
                      const LO& degree, Func rbf_func,
                      RealDefaultRank1View coefficients, double lambda,
                      double tol)
{
  PCMS_FUNCTION_TIMER;
  static_assert(std::is_invocable_r_v<double, Func, double, double>,
                "rbf_func, takes radius and cutoff, returns weight");
  static_assert(!std::is_pointer_v<Func>,
                "function pointer will fail in GPU execution context");

  const auto ntargets = target_coordinates.size() / dim;

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
  PCMS_ALWAYS_ASSERT(scratch_size > shared_size);

  Kokkos::parallel_for(
    "MLS coefficients", tp.set_scratch_size(1, Kokkos::PerTeam(shared_size)),
    KOKKOS_LAMBDA(const member_type& team) {
      int league_rank = team.league_rank();
      int start_ptr = support.supports_ptr[league_rank];
      int end_ptr = support.supports_ptr[league_rank + 1];
      int nsupports = end_ptr - start_ptr;

      // A missing query has no MLS system to solve and no coefficients. All
      // threads in this team take the same branch, before allocating scratch
      // or entering the SVD. The point evaluator applies its missing-point
      // policy afterwards.
      if (nsupports == 0) {
        return;
      }

      //  local_source_point stores the coordinates of source supports of a
      //  given target
      ScratchMatView local_source_points(team.team_scratch(1), nsupports, dim);

      // rbf function values of source supports Phi(n,n)
      ScratchVecView phi_vector(team.team_scratch(1), nsupports);

      //  vondermonde matrix P from the vectors of basis vector of supports
      ScratchMatView vandermonde_matrix(team.team_scratch(1), nsupports,
                                        basis_size);

      // basis of target
      ScratchVecView target_basis_vector(team.team_scratch(1), basis_size);

      // the pseudo-inverse M (basis_size x nsupports)
      ScratchMatView pseudo_inverse(team.team_scratch(1), basis_size,
                                    nsupports);

      // Initialize the scratch matrices and  vectors
      fill(0.0, team, local_source_points);
      fill(0.0, team, vandermonde_matrix);
      fill(0.0, team, phi_vector);
      fill(0.0, team, target_basis_vector);

      /**
       *
       * the local_source_points is of the type ScratchMatView with
       * coordinates information; row is number of local supports
       * & column is dim
       */

      int count = -1;
      for (int j = start_ptr; j < end_ptr; ++j) {
        count++;
        auto index = support.supports_idx[j];
        for (int i = 0; i < dim; ++i) {
          local_source_points(count, i) = source_coordinates[index * dim + i];
        }
      }

      double target_point[MAX_DIM] = {};
      for (int i = 0; i < dim; ++i) {
        target_point[i] = target_coordinates[league_rank * dim + i];
      }

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

      /**
       * step 2: normalize local source supports and target point
       */

      normalize_supports(team, target_point, local_source_points);
      team.team_barrier();
      /**
       *
       * evaluates the basis vector of a given target point Coord target_point;
       * this can evaluate monomial basis vector for any degree of polynomial
       *
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
      team.team_barrier();

      PCMS_CHECK_PRINTF(
        support.radii2[league_rank] > 0,
        "ERROR: radius2 has to be positive but found to be %.16f\n",
        support.radii2[league_rank]);

      /**
       * step 5: the pseudo-inverse of the weighted system by SVD
       */
      pseudo_inverse_svd(team, phi_vector, vandermonde_matrix, pseudo_inverse,
                         lambda, tol);
      team.team_barrier();

      /**
       * step 6: c_j = phi_j sum_k t_k M_kj
       */
      Kokkos::parallel_for(
        Kokkos::TeamThreadRange(team, nsupports), [=](int j) {
          double c = 0;
          for (int k = 0; k < basis_size; ++k) {
            c += target_basis_vector(k) * pseudo_inverse(k, j);
          }
          c *= phi_vector(j);
          PCMS_CHECK_PRINTF(!std::isnan(c), "Nan at %d\n", league_rank);
          coefficients[start_ptr + j] = c;
        });
    });
}

/**
 * @brief \overload
 * The coefficient rows of the targets, as a device array aligned with the
 * entries of `support` (see the array form).
 */
template <typename Func>
Kokkos::View<Real*, DeviceMemorySpace> mls_coefficients(
  Kokkos::View<const Real*, DeviceMemorySpace> source_coordinates,
  Kokkos::View<const Real*, DeviceMemorySpace> target_coordinates,
  const SupportResults& support, LO dim, LO degree, Func rbf_func,
  double lambda, double tol)
{
  Kokkos::View<Real*, DeviceMemorySpace> coefficients(
    "mls_coefficients", support.supports_idx.extent(0));
  mls_coefficients(
    RealConstDefaultRank1View(source_coordinates.data(),
                              source_coordinates.size()),
    RealConstDefaultRank1View(target_coordinates.data(),
                              target_coordinates.size()),
    support, dim, degree, rbf_func,
    RealDefaultRank1View(coefficients.data(), coefficients.size()), lambda,
    tol);
  return coefficients;
}

/**
 * @brief The MLS values at the targets from the coefficient rows:
 * values[i] = sum_j c_j source_values[idx_j] over the supports of target i.
 */
inline Kokkos::View<Real*, DeviceMemorySpace> apply_mls_coefficients(
  Kokkos::View<const Real*, DeviceMemorySpace> coefficients,
  Kokkos::View<const Real*, DeviceMemorySpace> source_values,
  const SupportResults& support, LO ntargets)
{
  Kokkos::View<Real*, DeviceMemorySpace> values("approximated target values",
                                                ntargets);
  Kokkos::parallel_for(
    "MLS apply coefficients", ntargets, KOKKOS_LAMBDA(const LO i) {
      Real sum = 0;
      for (LO j = support.supports_ptr[i]; j < support.supports_ptr[i + 1];
           ++j) {
        sum += coefficients[j] * source_values[support.supports_idx[j]];
      }
      values[i] = sum;
    });
  return values;
}

/**
 * @brief Maps the data from source mesh to target mesh: the coefficient rows
 * applied to the source values (a per-call solve; an evaluator that
 * evaluates repeatedly keeps the rows instead).
 */
template <typename Func>
Kokkos::View<Real*, DeviceMemorySpace> mls_interpolation(
  Kokkos::View<const Real*, DeviceMemorySpace> source_values,
  Kokkos::View<const Real*, DeviceMemorySpace> source_coordinates,
  Kokkos::View<const Real*, DeviceMemorySpace> target_coordinates,
  const SupportResults& support, LO dim, LO degree, Func rbf_func,
  double lambda, double tol)
{
  PCMS_CHECK_PRINTF(
    static_cast<LO>(source_values.size()) ==
      static_cast<LO>(source_coordinates.size()) / dim,
    "[ERROR] The size of the source values and source coordinates is not "
    "same. The current sizes are :\nsource_values_size = %d, "
    "source_coordinates_size = %d\n",
    static_cast<int>(source_values.size()),
    static_cast<int>(source_coordinates.size() / dim));
  const auto coefficients =
    mls_coefficients(source_coordinates, target_coordinates, support, dim,
                     degree, rbf_func, lambda, tol);
  return apply_mls_coefficients(coefficients, source_values, support,
                                static_cast<LO>(target_coordinates.size()) /
                                  dim);
}

} // namespace detail
} // namespace pcms::field

#endif // PCMS_FIELD_EVALUATOR_MLS_INTERPOLATION_IMPL_HPP
