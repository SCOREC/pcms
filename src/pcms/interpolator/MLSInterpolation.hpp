#ifndef MLS_INTERPOLATION_HPP
#define MLS_INTERPOLATION_HPP

#include "MLSCoefficients.hpp"
#include "adj_search.hpp"
#include <cassert>
#include <cmath>
using namespace Omega_h;

Write<Real> mls_interpolation(const Reals source_values,
                              const Reals source_coordinates,
                              const Reals target_coordinates,
                              const SupportResults& support, const LO& dim,
                              const LO& degree, Write<Real> radii2)
{
  const auto nvertices_source = source_coordinates.size() / dim;
  const auto nvertices_target = target_coordinates.size() / dim;

  // Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace> range_policy(1,
  // nvertices_target);

  // static_assert(degree > 0," the degree of polynomial basis should be atleast
  // 1");

  Kokkos::View<size_t*> shmem_each_team(
    "stores the size required for each team", nvertices_target);

  Kokkos::View<int**, Kokkos::HostSpace> host_slice_length(
    "stores slice length of  polynomial basis in host", degree, dim);
  Kokkos::deep_copy(host_slice_length, 0);
  basisSliceLengths(host_slice_length);

  auto basis_size = basisSize(host_slice_length);
  printf("basis_size is %d\n", basis_size);
  MatViewType slice_length("stores slice length of polynomial basis in device",
                           degree, dim);
  auto slice_length_hd = Kokkos::create_mirror_view(slice_length);
  Kokkos::deep_copy(slice_length_hd, host_slice_length);
  Kokkos::deep_copy(slice_length, slice_length_hd);

  Kokkos::parallel_for(
    "calculate the size required for scratch for each team", nvertices_target,
    KOKKOS_LAMBDA(const int i) {
      int start_ptr = support.supports_ptr[i];
      int end_ptr = support.supports_ptr[i + 1];
      int nsupports = end_ptr - start_ptr;

      size_t total_shared_size = 0;

      total_shared_size +=
        ScratchMatView::shmem_size(basis_size, basis_size) * 4;
      total_shared_size +=
        ScratchMatView::shmem_size(basis_size, nsupports) * 2;
      total_shared_size += ScratchMatView::shmem_size(nsupports, basis_size);
      total_shared_size += ScratchVecView::shmem_size(basis_size);
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

  Write<Real> approx_target_values(nvertices_target, 0,
                                   "approximated target values");

  team_policy tp(nvertices_target, Kokkos::AUTO);

  int scratch_size = tp.scratch_size_max(0);

  assert(scratch_size > shared_size &&
         "The required scratch size exceeds the max available scratch size");

  Kokkos::parallel_for(
    "MLS coefficients", tp.set_scratch_size(0, Kokkos::PerTeam(shared_size)),
    KOKKOS_LAMBDA(const member_type& team) {
      int i = team.league_rank();
      int start_ptr = support.supports_ptr[i];
      int end_ptr = support.supports_ptr[i + 1];
      int nsupports = end_ptr - start_ptr;

      ScratchMatView local_source_points(team.team_scratch(0), nsupports, dim);
      int count = -1;
      for (int j = start_ptr; j < end_ptr; ++j) {
        count++;
        auto index = support.supports_idx[j];
        local_source_points(count, 0) = source_coordinates[index * dim];
        local_source_points(count, 1) = source_coordinates[index * dim + 1];
        if (dim == 3) {
          local_source_points(count, 2) = source_coordinates[index * dim + 2];
        }
      }

      ScratchMatView lower(team.team_scratch(0), basis_size, basis_size);

      ScratchMatView forward_matrix(team.team_scratch(0), basis_size,
                                    basis_size);

      ScratchMatView moment_matrix(team.team_scratch(0), basis_size,
                                   basis_size);

      ScratchMatView inv_mat(team.team_scratch(0), basis_size, basis_size);

      ScratchMatView V(team.team_scratch(0), nsupports, basis_size);

      ScratchMatView Ptphi(team.team_scratch(0), basis_size, nsupports);

      ScratchMatView resultant_matrix(team.team_scratch(0), basis_size,
                                      nsupports);

      ScratchVecView targetMonomialVec(team.team_scratch(0), basis_size);

      ScratchVecView SupportValues(team.team_scratch(0), nsupports);

      ScratchVecView result(team.team_scratch(0), nsupports);

      ScratchVecView Phi(team.team_scratch(0), nsupports);

      //	Kokkos::deep_copy(lower, 0.0);
      //	Kokkos::deep_copy(forward_matrix, 0.0);
      //	Kokkos::deep_copy(moment_matrix, 0.0);
      //	Kokkos::deep_copy(inv_mat, 0.0);
      //	Kokkos::deep_copy(V, 0.0);
      //	Kokkos::deep_copy(Ptphi, 0.0);
      //	Kokkos::deep_copy(resultant_matrix, 0.0);
      //	Kokkos::deep_copy(targetMonomialVec, 0.0);
      //	Kokkos::deep_copy(SupportValues, 0.0);
      //	Kokkos::deep_copy(result, 0.0);
      //	Kokkos::deep_copy(Phi, 0.0);
      //
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, basis_size),
                           [=](int j) {
                             for (int k = 0; k < basis_size; ++k) {
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
                             for (int k = 0; k < basis_size; ++k) {
                               V(j, k) = 0;
                             }

                             SupportValues(j) = 0;
                             result(j) = 0;
                             Phi(j) = 0;
                           });

      Coord target_point;
      target_point.x = target_coordinates[i * dim];
      target_point.y = target_coordinates[i * dim + 1];

      if (dim == 3) {
        target_point.z = target_coordinates[i * dim + 2];
      }
      BasisPoly(targetMonomialVec, slice_length, target_point);

      Kokkos::parallel_for(
        Kokkos::TeamThreadRange(team, nsupports), [=](int j) {
          VandermondeMatrix(V, local_source_points, j, slice_length);
        });

      team.team_barrier();

      Kokkos::parallel_for(
        Kokkos::TeamThreadRange(team, nsupports), [=](int j) {
          OMEGA_H_CHECK_PRINTF(
            radii2[i] > 0,
            "ERROR: radius2 has to be positive but found to be %.16f\n",
            radii2[i]);
          PhiVector(Phi, target_point, local_source_points, j, radii2[i]);
        });

      // sum phi
      double sum_phi = 0;
      Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(team, nsupports),
        [=](const int j, double& lsum) { lsum += Phi(j); }, sum_phi);
      OMEGA_H_CHECK_PRINTF(!std::isnan(sum_phi),
                           "ERROR: sum_phi is NaN for i=%d\n", i);
      OMEGA_H_CHECK_PRINTF(sum_phi != 0, "ERROR: sum_phi is zero for i=%d\n",
                           i);

      // normalize phi with sum_phi
      Kokkos::parallel_for(
        Kokkos::TeamThreadRange(team, nsupports), [=](int j) {
          OMEGA_H_CHECK_PRINTF(
            !std::isnan(Phi(j)),
            "ERROR: Phi(j) is NaN before normalization for j = %d\n", j);
          Phi(j) = Phi(j) / sum_phi;
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
          SupportValues(i) = source_values[support.supports_idx[start_ptr + i]];
          OMEGA_H_CHECK_PRINTF(!std::isnan(SupportValues(i)),
                               "ERROR: NaN found: at support %d\n", i);
        });

      double tgt_value = 0;
      dot_product(team, result, SupportValues, tgt_value);
      if (team.team_rank() == 0) {
        OMEGA_H_CHECK_PRINTF(!std::isnan(tgt_value), "Nan at %d\n", i);
        approx_target_values[i] = tgt_value;
      }
    });

  return approx_target_values;
}

#endif
