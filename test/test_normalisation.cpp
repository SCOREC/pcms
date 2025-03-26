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

        double target_point[MAX_DIM] = {};

        for (int i = 0; i < dim; ++i) {
          target_point[i] = target_coordinates[league_rank][i];
        }

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
              normalized_target_coordinates[index] = target_point[j];
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
