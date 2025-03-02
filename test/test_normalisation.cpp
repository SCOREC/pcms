#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <pcms/interpolator/mls_interpolation.hpp>
#include <pcms/interpolator/pcms_interpolator_aliases.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <vector>
#include <iostream>
#include <Omega_h_array_ops.hpp>
using namespace std;
using namespace Omega_h;
using namespace pcms;

Reals fillFromMatrix(const Matrix<2, 10>& matrix)
{
  int num_points = 10;
  Write<Real> array(20, 0, "array to store coordinates");
  parallel_for(
    "fills the array view from matrix", num_points, KOKKOS_LAMBDA(const int i) {
      int index = i * 2;

      array[index] = matrix[i][0];
      array[index + 1] = matrix[i][1];
    });
  return read(array);
}

TEST_CASE("test_normalization_routine")
{

  SECTION("test_minmax_normalization")
  {

    printf("-------------the test for min_max starts-----------\n");

    Real tolerance = 1E-5;
    int num_points = 10;

    Omega_h::Matrix<2, 10> matrix{{0.45, 3.18},   {0.12, -1.2},  {-2.04, -1.45},
                                  {1.52, -0.98},  {-0.22, 1.45}, {0.96, 3.62},
                                  {-0.26, -1.62}, {0.82, 4.22},  {5.6, 3.62},
                                  {1.2, 1.2}};

    auto coordinates = fillFromMatrix(matrix);

    auto results = pcms::detail::min_max_normalization(coordinates, 2);

    printf("results after normalization");
    HostRead<Real> host_results(results);

    double expected_results[10][2] = {
      {0.325916, 0.821918}, {0.282723, 0.071918}, {0.000000, 0.029110},
      {0.465969, 0.109589}, {0.238220, 0.525685}, {0.392670, 0.897260},
      {0.232984, 0.000000}, {0.374346, 1.000000}, {1.000000, 0.897260},
      {0.424084, 0.482877}};

    printf("<scaled coordinates>\n");

    for (int i = 0; i < 20; ++i) {
      printf("%12.6f\n", host_results[i]);
    }

    for (int i = 0; i < 10; ++i) {

      int index = i * 2;

      CAPTURE(i, host_results[index], expected_results[i][0], tolerance);
      CHECK_THAT(host_results[index],
                 Catch::Matchers::WithinAbs(expected_results[i][0], tolerance));

      CAPTURE(i, host_results[index + 1], expected_results[i][1], tolerance);
      CHECK_THAT(host_results[index + 1],
                 Catch::Matchers::WithinAbs(expected_results[i][1], tolerance));

      printf("-------------the test for min_max starts-----------\n");
    }

  } // end SECTION

  SECTION("test_normalization_relative_distance")
  {
    Real tolerance = 1E-5;
    int nvertices_target = 3;

    Omega_h::Matrix<3, 3> target_coordinates{
      {0.5, 2.0, 1.4}, {0.1, 1.2, 0.3}, {0.9, 2.8, 2.2}};

    Omega_h::Matrix<3, 6> supports_target_1{
      {0.45, 3.18, 0.38}, {0.12, 1.2, 1.0},  {2.04, 1.45, 3.2},
      {1.52, 0.98, 3.2},  {0.22, 1.45, 3.2}, {0.96, 3.62, 2.2}};

    Omega_h::Matrix<3, 8> supports_target_2{
      {0.88, 0.98, -0.8}, {1.2, 0.23, -1.21}, {-1.23, 2.85, 1.45},
      {0.63, 0.98, 0.2},  {-0.5, 1.45, 0.33}, {1.4, 2.61, -0.98},
      {0.88, 1.72, 0.38}, {1.2, 0.50, 2.1}};

    Omega_h::Matrix<3, 4> supports_target_3{{3.15, 3.12, 2.91},
                                            {2.24, 1.12, 1.12},
                                            {0.12, 4.29, 1.98},
                                            {1.24, 1.54, 4.2}};

    LOs nsupports({6, 8, 4}, "number of supports");
    int dim = 3;

    auto total_coordinates = Omega_h::get_sum(nsupports) * dim;

    REQUIRE(total_coordinates == 54);
    Write<Real> normalized_support_coordinates(
      total_coordinates, 0, "coordinates after normalization");

    Write<Real> normalized_target_coordinates(
      total_coordinates, 0, "coordinates after normalization");

    Write<Real> all_support_coordinates(total_coordinates, 0,
                                        "support coordinates all");
    Write<Real> all_target_coordinates(total_coordinates, 0,
                                       "target coordinates all");

    team_policy tp(nvertices_target, Kokkos::AUTO);
    Kokkos::parallel_for(
      "inside team", tp.set_scratch_size(0, Kokkos::PerTeam(200)),
      KOKKOS_LAMBDA(const member_type& team) {
        int league_rank = team.league_rank();

        int num_supports = nsupports[league_rank];
        ScratchMatView local_supports(team.team_scratch(0), num_supports, dim);
        detail::fill(0, team, local_supports);

        Coord target_point;
        target_point.x = target_coordinates[league_rank][0];
        target_point.y = target_coordinates[league_rank][1];
        target_point.z = target_coordinates[league_rank][2];

        Kokkos::parallel_for(
          Kokkos::TeamThreadRange(team, num_supports), [=](int i) {
            for (int j = 0; j < dim; ++j) {
              int index = league_rank * num_supports * dim + i * dim + j;

              all_target_coordinates[index] =
                target_coordinates[league_rank][j];
              if (league_rank == 0) {
                local_supports(i, j) = supports_target_1[i][j];
                all_support_coordinates[index] = supports_target_1[i][j];
              } else if (league_rank == 1) {
                local_supports(i, j) = supports_target_2[i][j];
                all_support_coordinates[index] = supports_target_2[i][j];
              } else {
                local_supports(i, j) = supports_target_3[i][j];
                all_support_coordinates[index] = supports_target_3[i][j];
              }
            }
          });

        team.team_barrier();
        detail::normalize_supports(team, target_point, local_supports);

        Kokkos::parallel_for(
          Kokkos::TeamThreadRange(team, num_supports), [=](int i) {
            for (int j = 0; j < dim; ++j) {
              int index = league_rank * num_supports * dim + i * dim + j;
              normalized_support_coordinates[index] = local_supports(i, j);
              if (j == 0) {
                normalized_target_coordinates[index] = target_point.x;
              } else if (j == 1) {
                normalized_target_coordinates[index] = target_point.y;
              } else {
                normalized_target_coordinates[index] = target_point.z;
              }
            }
          });
      });

    auto host_all_support_coordinates =
      HostRead<Real>(read(all_support_coordinates));

    auto host_all_target_coordinates =
      HostRead<Real>(read(all_target_coordinates));

    auto host_normalized_support_coordinates =
      HostRead<Real>(read(normalized_support_coordinates));

    auto host_normalized_target_coordinates =
      HostRead<Real>(read(normalized_target_coordinates));

    for (int i = 0; i < total_coordinates; ++i) {
      auto result_support =
        host_all_support_coordinates[i] - host_all_target_coordinates[i];

      auto result_target =
        host_all_target_coordinates[i] - host_all_target_coordinates[i];

      CHECK_THAT(host_normalized_support_coordinates[i],
                 Catch::Matchers::WithinAbs(result_support, tolerance));

      CHECK_THAT(host_normalized_target_coordinates[i],
                 Catch::Matchers::WithinAbs(result_target, tolerance));
    }
  } // end SECTION

} // end TESTCASE
