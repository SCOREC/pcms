#ifndef PCMS_GNET_GNET_SAMPLING_HPP
#define PCMS_GNET_GNET_SAMPLING_HPP

#include <pcms/gnet/gnet_grid.hpp>
#include <pcms/utility/uniform_grid.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <Omega_h_array.hpp>
#include <cstdint>

namespace pcms::gnet
{

/// Settings for rejection sampling of GTC markers
struct SamplingParams
{
  int total_samples;  // markers to produce
  int batch_size;     // candidates per round
  double threshold;   // accept candidates with interpolated value >= threshold
  int degree;         // MLS polynomial degree
  double lambda;      // MLS regularization
  std::uint64_t seed; // random seed
};

/// Accepted GTC markers on the device
struct GtcMarkers
{
  Omega_h::Write<Omega_h::Real> coordinates; // total_samples * 5
  Omega_h::Write<Omega_h::Real> values;      // total_samples
};

/**
 * @brief Uniform random points in grid. Dimension 0 (rho) is drawn as
 * sqrt(u) so points are uniform in area.
 */
Omega_h::Reals generate_random_points_in_grid(
  const pcms::UniformGrid<MAX_GNET_DIM>& grid, Omega_h::LO num_points,
  Kokkos::Random_XorShift64_Pool<>& rand_pool);

Omega_h::Real get_min_volume(
  const Kokkos::View<double*, Kokkos::HostSpace>& host_volume3d);

Omega_h::Real get_max_volume(
  const Kokkos::View<double*, Kokkos::HostSpace>& host_volume3d);

/**
 * @brief Rejection sampling of GTC markers from GNET data.
 *
 * Each round draws batch_size candidates in grid5d, interpolates the GNET
 * values to them with MLS, and keeps candidates whose value is >= threshold,
 * until total_samples markers are accepted.
 */
GtcMarkers sample_gtc_markers(const Omega_h::Reals& source_values,
                              const Omega_h::Reals& source_coordinates,
                              const pcms::UniformGrid<MAX_GNET_DIM>& grid5d,
                              const GnetDims& dims,
                              const SamplingParams& params);

} // namespace pcms::gnet

#endif // PCMS_GNET_GNET_SAMPLING_HPP
