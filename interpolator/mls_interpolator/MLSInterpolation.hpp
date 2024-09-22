#ifndef MLS_INTERPOLATION_HPP
#define MLS_INTERPOLATION_HPP

#include "MLSCoefficients.hpp"
#include "adj_search_dega2.hpp"
#include "adj_search.hpp"
#include "points.hpp"

using namespace Omega_h;
using namespace pcms;

Write<Real> mls_interpolation(const Reals source_values,
                              const Reals source_coordinates,
                              const Reals target_coordinates,
                              const SupportResults& support, const LO& dim,
                              Write<Real> radii2) {
  const auto nvertices_source = source_coordinates.size() / dim;
  const auto nvertices_target = target_coordinates.size() / dim;

  // Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace> range_policy(1,
  // nvertices_target);

  Kokkos::View<size_t*> shmem_each_team(
      "stores the size required for each team", nvertices_target);
  // Write<int> shmem_each_team(nvertices_target, 0, "stores the size required
  // for each team member");

  Kokkos::parallel_for(
      "calculate the size required for scratch for each team", nvertices_target,
      KOKKOS_LAMBDA(const int i) {
        int start_ptr = support.supports_ptr[i];
        int end_ptr = support.supports_ptr[i + 1];
        int nsupports = end_ptr - start_ptr;

        size_t total_shared_size = 0;

        total_shared_size += ScratchMatView::shmem_size(6, 6) * 4;
        total_shared_size += ScratchMatView::shmem_size(6, nsupports) * 2;
        total_shared_size += ScratchMatView::shmem_size(nsupports, 6);
        total_shared_size += ScratchVecView::shmem_size(6);
        total_shared_size += ScratchVecView::shmem_size(nsupports) * 3;
        total_shared_size += ScratchMatView::shmem_size(nsupports, 2);
        shmem_each_team(i) = total_shared_size;
      });

  size_t shared_size;
  Kokkos::parallel_reduce(
      "find_max", nvertices_target,
      KOKKOS_LAMBDA(const int i, size_t& max_val_temp) {
        if (shmem_each_team(i) > max_val_temp) {
          max_val_temp = shmem_each_team(i);
        }
      },
      Kokkos::Max<size_t>(shared_size));
  printf("shared size = %d \n", shared_size);

  Write<Real> approx_target_values(nvertices_target, 0,
                                   "approximated target values");

  team_policy tp(nvertices_target, Kokkos::AUTO);

  int scratch_size = tp.scratch_size_max(0);
  printf("scratch size is %d \n", scratch_size);

  Kokkos::parallel_for(
      "MLS coefficients", tp.set_scratch_size(0, Kokkos::PerTeam(shared_size)),
      KOKKOS_LAMBDA(const member_type& team) {
        int i = team.league_rank();
        int start_ptr = support.supports_ptr[i];
        int end_ptr = support.supports_ptr[i + 1];

        int nsupports = end_ptr - start_ptr;

        ScratchMatView local_source_points(team.team_scratch(0), nsupports, 2);
        int count = -1;
        for (int j = start_ptr; j < end_ptr; ++j) {
          count++;
          auto index = support.supports_idx[j];
          local_source_points(count, 0) = source_coordinates[index * dim];
          local_source_points(count, 1) = source_coordinates[index * dim + 1];
        }

        ScratchMatView lower(team.team_scratch(0), 6, 6);

        ScratchMatView forward_matrix(team.team_scratch(0), 6, 6);

        ScratchMatView moment_matrix(team.team_scratch(0), 6, 6);

        ScratchMatView inv_mat(team.team_scratch(0), 6, 6);

        ScratchMatView V(team.team_scratch(0), nsupports, 6);

        ScratchMatView Ptphi(team.team_scratch(0), 6, nsupports);

        ScratchMatView resultant_matrix(team.team_scratch(0), 6, nsupports);

        ScratchVecView targetMonomialVec(team.team_scratch(0), 6);

        ScratchVecView SupportValues(team.team_scratch(0), nsupports);

        ScratchVecView result(team.team_scratch(0), nsupports);

        ScratchVecView Phi(team.team_scratch(0), nsupports);

        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 6), [=](int j) {
          for (int k = 0; k < 6; ++k) {
            lower(j, k) = 0;
            forward_matrix(j, k) = 0;
            moment_matrix(j, k) = 0;
            inv_mat(j, k) = 0;
          }

          targetMonomialVec(j) = 0;
          for (int k = 0; k < nsupports; ++k) {
            resultant_matrix(j, k) = 0;

            Ptphi(j, k) = 0;
          }
        });

        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nsupports),
                             [=](int j) {
                               for (int k = 0; k < 6; ++k) {
                                 V(j, k) = 0;
                               }

                               SupportValues(j) = 0;
                               result(j) = 0;
                               Phi(j) = 0;
                             });

        Coord target_point;

        target_point.x = target_coordinates[i * dim];

        target_point.y = target_coordinates[i * dim + 1];

        BasisPoly(targetMonomialVec, target_point);

        Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team, nsupports),
            [=](int j) { VandermondeMatrix(V, local_source_points, j); });

        team.team_barrier();

        Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team, nsupports), [=](int j) {
              PhiVector(Phi, target_point, local_source_points, j, radii2[i]);
            });

        team.team_barrier();

        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nsupports),
                             [=](int j) { PTphiMatrix(Ptphi, V, Phi, j); });

        team.team_barrier();

        MatMatMul(team, moment_matrix, Ptphi, V);
        team.team_barrier();

        inverse_matrix(team, moment_matrix, lower, forward_matrix, inv_mat);

        team.team_barrier();

        MatMatMul(team, resultant_matrix, inv_mat, Ptphi);
        team.team_barrier();

        MatVecMul(team, targetMonomialVec, resultant_matrix, result);
        team.team_barrier();

        Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team, nsupports), [=](const int i) {
              SupportValues(i) =
                  source_values[support.supports_idx[start_ptr + i]];
            });

        double tgt_value = 0;
        dot_product(team, result, SupportValues, tgt_value);
        if (team.team_rank() == 0) {
          approx_target_values[i] = tgt_value;
        }
      });

  return approx_target_values;
}

#endif
