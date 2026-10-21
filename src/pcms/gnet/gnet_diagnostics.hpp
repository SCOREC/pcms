#ifndef PCMS_GNET_GNET_DIAGNOSTICS_HPP
#define PCMS_GNET_GNET_DIAGNOSTICS_HPP

#include <pcms/gnet/gnet_grid.hpp>
#include <pcms/utility/uniform_grid.h>
#include <Kokkos_Core.hpp>
#include <Omega_h_array.hpp>

namespace pcms::gnet
{

using Histogram2D = Kokkos::View<double**, Kokkos::LayoutRight>;

/**
 * @brief Energy-pitch histogram: sums f over points binned on grid5d
 * dimensions 3 (energy) and 4 (pitch).
 * @param f5d Value at each point
 * @param sampled_points Flattened 5D coordinates of the points
 */
Histogram2D calculate_velocity_distribution(
  const Omega_h::Reals& f5d, const Omega_h::Reals& sampled_points,
  const pcms::UniformGrid<MAX_GNET_DIM>& grid5d);

/**
 * @brief Rho-theta histogram: sums f over points binned on grid5d
 * dimensions 0 (rho) and 1 (theta).
 */
Histogram2D calculate_2d_distribution(
  const Omega_h::Reals& f5d, const Omega_h::Reals& sampled_points,
  const pcms::UniformGrid<MAX_GNET_DIM>& grid5d);

/**
 * @brief Packs points into an N x 6 table: rho, theta, zeta, energy, pitch,
 * value.
 */
Kokkos::View<double**, Kokkos::LayoutRight> get2DArray(
  const Omega_h::Reals& f5d, const Omega_h::Reals& sampled_points);

} // namespace pcms::gnet

#endif // PCMS_GNET_GNET_DIAGNOSTICS_HPP
