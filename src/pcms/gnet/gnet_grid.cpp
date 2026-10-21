#include <pcms/gnet/gnet_grid.hpp>
#include <pcms/field/evaluator/mls_interpolation_impl.hpp>
#include <iostream>

namespace pcms::gnet
{

Omega_h::Write<Omega_h::Real> get_gnet_coordinates(
  const size_t& mrho_max, const size_t& mtor_max, const size_t& mlmd_max,
  const size_t& meng_max, const double& zeta_max, const double& emax,
  const Omega_h::LOs& pol_max_scan_array, const size_t& npoints,
  const size_t& dim)
{

  Omega_h::Write<Omega_h::Real> gnet_coords(npoints * dim, -2,
                                            "stores gnet coordinates");

  Omega_h::parallel_for(
    npoints, OMEGA_H_LAMBDA(const int i) {
      const size_t grid = i / (meng_max * mlmd_max);
      const size_t energy = (i / mlmd_max) % meng_max;
      const size_t pitch = i % mlmd_max;
      const size_t zeta = grid % mtor_max;
      const size_t zeta_total = grid / mtor_max;
      const size_t rho = find_rho_binary_search(pol_max_scan_array, zeta_total);
      const size_t theta = zeta_total - pol_max_scan_array[rho];

      gnet_coords[i * dim + 0] = compute_radial_coordinate(mrho_max, rho);
      gnet_coords[i * dim + 1] = compute_poloidal_coordinate(rho, theta);
      gnet_coords[i * dim + 2] =
        compute_toroidal_coordinate(mtor_max, zeta_max, zeta);
      gnet_coords[i * dim + 3] =
        compute_energy_coordinate(meng_max, emax, energy);
      gnet_coords[i * dim + 4] = compute_pitch_coordinate(mlmd_max, pitch);
    });

  return gnet_coords;
}

Omega_h::LOs get_mpol_max_array(size_t rho_max)
{
  Omega_h::Write<Omega_h::LO> mpol_max_array(rho_max, 0, "polar max array");

  Omega_h::parallel_for(
    rho_max, OMEGA_H_LAMBDA(int i) { mpol_max_array[i] = 4 * (i + 1); });

  return read(mpol_max_array);
}

Omega_h::LOs get_pol_max_scan_array(size_t mrho_max)
{

  Omega_h::Write<Omega_h::LO> pol_max_scan_array(mrho_max + 1, 0,
                                                 "pol max scan array");

  Kokkos::parallel_scan(
    mrho_max, KOKKOS_LAMBDA(int i, int& partial_sum, bool is_final) {
      partial_sum += 4 * (i + 1);

      if (is_final) {
        pol_max_scan_array[i + 1] = partial_sum;
      }
    });

  return Omega_h::read(pol_max_scan_array);
}

size_t compute_min_num_supports(const size_t degree, const size_t dim)
{

  pcms::IntHostMatView slice_lengths_array("array", degree, dim);
  Kokkos::deep_copy(slice_lengths_array, 0);
  pcms::detail::calculate_basis_slice_lengths(slice_lengths_array);

  return pcms::detail::calculate_basis_vector_size(slice_lengths_array);
}

pcms::SupportResults searchNeighbors(
  const Omega_h::Reals& target_points, const size_t mrho_max,
  const size_t mtor_max, const size_t mlmd_max, const size_t meng_max,
  const double zeta_max, const double emax, size_t min_num_supports)
{

  size_t num_target_points = target_points.size() / MAX_GNET_DIM;

  auto pol_max_scan_array = get_pol_max_scan_array(mrho_max);

  BinNeighbors binNeighbors{};

  binNeighbors.eval(min_num_supports);

  size_t total_neighbors =
    binNeighbors.num_rho_neighbor_bins * binNeighbors.num_theta_neighbor_bins *
    binNeighbors.num_zeta_neighbor_bins * binNeighbors.num_pitch_neighbor_bins *
    binNeighbors.num_energy_neighbor_bins;

  Omega_h::Write<Omega_h::LO> supports_ptr(
    num_target_points + 1, 0, "stores the ptr for support indices");

  int supports_size;
  Kokkos::parallel_scan(
    num_target_points,
    KOKKOS_LAMBDA(int i, int& partial_sum, bool is_final) {
      partial_sum += total_neighbors;
      if (is_final) {
        supports_ptr[i + 1] = partial_sum;
      }
    },
    supports_size);

  // Validate results
  Omega_h::HostRead<Omega_h::LO> host_supports_ptr(supports_ptr);

  if (host_supports_ptr[num_target_points] != supports_size) {
    std::cerr << "Error: Mismatch detected! supports_ptr[num_target_points] != "
                 "supports_size\n";
  }
  Omega_h::Write<Omega_h::LO> supports_idx(supports_size, 0,
                                           "stores the index of supports");

  Omega_h::parallel_for(
    num_target_points, OMEGA_H_LAMBDA(const Omega_h::LO id) {
      Omega_h::Real query_point[MAX_GNET_DIM];

      for (int i = 0; i < MAX_GNET_DIM; ++i) {
        query_point[i] = target_points[id * MAX_GNET_DIM + i];
      }

      GnetSearch search(query_point);

      Omega_h::LO support_each_target_ids[MAX_MARKER_NEIGHBORS] = {};

      search.analyticSearchGnet(mrho_max, mtor_max, mlmd_max, meng_max,
                                zeta_max, emax, min_num_supports,
                                pol_max_scan_array, support_each_target_ids);

      int support_counter = supports_ptr[id];

      for (int i = 0; i < total_neighbors; i++) {
        supports_idx[support_counter + i] = support_each_target_ids[i];
      }
    });

  Omega_h::Write<Omega_h::Real> ignored(num_target_points, 1);
  return pcms::SupportResults{Omega_h::read(supports_ptr),
                              Omega_h::read(supports_idx), ignored};
}

} // namespace pcms::gnet
