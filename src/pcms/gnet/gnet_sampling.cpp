#include <pcms/gnet/gnet_sampling.hpp>
#include <pcms/field/evaluator/mls_interpolation.hpp>
#include <Kokkos_MathematicalFunctions.hpp>
#include <Omega_h_array_ops.hpp>
#include <Omega_h_for.hpp>
#include <algorithm>
#include <cstdio>
#include <iostream>

namespace pcms::gnet
{

namespace
{

constexpr int dim = MAX_GNET_DIM;

Omega_h::Reals volume_to_device(
  const Kokkos::View<double*, Kokkos::HostSpace>& host_volume3d)
{
  using DeviceSpace = Kokkos::DefaultExecutionSpace::memory_space;
  Kokkos::View<double*, DeviceSpace> device_volume3d("particle density",
                                                     host_volume3d.size());

  // copy from host to device view
  Kokkos::deep_copy(device_volume3d, host_volume3d);

  Omega_h::Write<Omega_h::Real> volume3d(device_volume3d.size(), 0,
                                         "stores volume3d");

  Omega_h::parallel_for(
    device_volume3d.size(),
    OMEGA_H_LAMBDA(int i) { volume3d[i] = device_volume3d(i); });

  return Omega_h::read(volume3d);
}

} // namespace

Omega_h::Reals generate_random_points_in_grid(
  const pcms::UniformGrid<MAX_GNET_DIM>& grid, Omega_h::LO num_points,
  Kokkos::Random_XorShift64_Pool<>& rand_pool)
{
  Omega_h::Write<Omega_h::Real> sampled_points(num_points * dim);

  Kokkos::parallel_for(
    "sample_points_in_grid", num_points, KOKKOS_LAMBDA(const int i) {
      auto rand_gen = rand_pool.get_state();

      for (int d = 0; d < dim; ++d) {
        Omega_h::Real t = rand_gen.drand(); // [0, 1)
        if (d == 0)
          t = Kokkos::sqrt(t);
        sampled_points[i * dim + d] = grid.bot_left[d] + t * grid.edge_length[d];
      }

      rand_pool.free_state(rand_gen);
    });

  return Omega_h::read(sampled_points);
}

Omega_h::Real get_min_volume(
  const Kokkos::View<double*, Kokkos::HostSpace>& host_volume3d)
{
  return Omega_h::get_min(volume_to_device(host_volume3d));
}

Omega_h::Real get_max_volume(
  const Kokkos::View<double*, Kokkos::HostSpace>& host_volume3d)
{
  return Omega_h::get_max(volume_to_device(host_volume3d));
}

GtcMarkers sample_gtc_markers(const Omega_h::Reals& source_values,
                              const Omega_h::Reals& source_coordinates,
                              const pcms::UniformGrid<MAX_GNET_DIM>& grid5d,
                              const GnetDims& dims,
                              const SamplingParams& params)
{
  const int total_samples = params.total_samples;
  const int batch_size = params.batch_size;
  const double threshold_value = params.threshold;
  const auto min_num_supports = compute_min_num_supports(params.degree, dim);

  // sampled points array
  Omega_h::Write<Omega_h::Real> gtc_marker_coordinates(
    total_samples * dim, 0, "stores the 5 coordinate values");
  Omega_h::Write<Omega_h::Real> gtc_f5d(total_samples, 0,
                                        "stores the 5 coordinate values");

  // initialize random pool once
  Kokkos::Random_XorShift64_Pool<> rand_pool(params.seed);

  Kokkos::View<int> accepted_count("accepted count");

  int total_accepted = 0;
  int round = 0;

  while (total_accepted < total_samples) {
    ++round;
    Kokkos::deep_copy(accepted_count, 0);

    // generate candidate points
    auto sampled_tgt_pts =
      generate_random_points_in_grid(grid5d, batch_size, rand_pool);

    // neighbor search
    auto support = searchNeighbors(sampled_tgt_pts, dims.mrho_max,
                                   dims.mtor_max, dims.mlmd_max, dims.meng_max,
                                   dims.zeta_max, dims.emax, min_num_supports);
    auto host_support_ptr = Omega_h::HostRead(support.supports_ptr);
    if (total_accepted == 0) {
      printf("Total supports : %d\n", host_support_ptr[1]);
    }
    // mapping from gnet to gtc on candidate points
    auto tgt_values = pcms::mls_interpolation(
      source_values, source_coordinates, sampled_tgt_pts, support, dim,
      params.degree, pcms::RadialBasisFunction::NO_OP, params.lambda);

    // buffer for accepted indices
    Omega_h::Write<Omega_h::LO> accepted_indices(batch_size,
                                                 "accepted indices");

    // evaluate threshold
    Omega_h::parallel_for(
      batch_size, OMEGA_H_LAMBDA(int i) {
        if (tgt_values[i] >= threshold_value) {
          int idx = Kokkos::atomic_fetch_add(&accepted_count(), 1);

          // store index of accepted sample
          accepted_indices[idx] = i;
        }
      });

    // copy only the first n valid samples until it reaches total_samples
    int accepted_this_batch = 0;
    Kokkos::deep_copy(accepted_this_batch, accepted_count);

    int remaining = total_samples - total_accepted;
    int num_copy_to_buffer = std::min(accepted_this_batch, remaining);

    Omega_h::parallel_for(
      num_copy_to_buffer, OMEGA_H_LAMBDA(int i) {
        int src_idx = accepted_indices[i];
        for (int j = 0; j < dim; ++j) {
          gtc_marker_coordinates[(total_accepted + i) * dim + j] =
            sampled_tgt_pts[src_idx * dim + j];
        }

        gtc_f5d[total_accepted + i] = tgt_values[src_idx];
      });

    total_accepted += num_copy_to_buffer;

    std::cerr << "Round #" << round
              << ": Accepted this batch = " << num_copy_to_buffer
              << ", Total = " << total_accepted << "\n";
  }

  return GtcMarkers{gtc_marker_coordinates, gtc_f5d};
}

} // namespace pcms::gnet
