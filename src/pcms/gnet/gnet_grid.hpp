#ifndef PCMS_GNET_GNET_GRID_HPP
#define PCMS_GNET_GNET_GRID_HPP

#include <Omega_h_array.hpp>
#include <Omega_h_for.hpp>
#include <Kokkos_MathematicalFunctions.hpp>
#include <pcms/localization/adj_search.hpp>
#include <cmath>
#include <cstddef>

namespace pcms::gnet
{

constexpr static int MAX_GNET_DIM = 5;

constexpr static int MAX_NEIGHBORS_EACH_DIR = 200;

constexpr static int MAX_GRID_NEIGHBORS = 200;

constexpr static int MAX_MARKER_NEIGHBORS = 500;

/// GNET phase-space grid sizes, as stored in the GNET file
struct GnetDims
{
  size_t mrho_max; // radial bins
  size_t mtor_max; // toroidal bins
  size_t mlmd_max; // pitch bins
  size_t meng_max; // energy bins
  double zeta_max; // toroidal extent (rad)
  double emax;     // max energy (eV)
};

struct BinNeighbors
{
  size_t num_in_each_dir = 0;
  size_t num_rho_neighbor_bins = num_in_each_dir;
  size_t num_theta_neighbor_bins = num_in_each_dir + 1;
  size_t num_zeta_neighbor_bins = num_in_each_dir + 1;
  size_t num_energy_neighbor_bins = num_in_each_dir + 1;
  size_t num_pitch_neighbor_bins = num_in_each_dir + 1;

  OMEGA_H_INLINE
  void eval(size_t min_num_supports)
  {
    num_in_each_dir = std::pow(min_num_supports, 0.2); // 1/dim dim =5
    num_rho_neighbor_bins = num_in_each_dir + 1;
    num_theta_neighbor_bins = num_in_each_dir + 1;
    num_zeta_neighbor_bins = num_in_each_dir + 1;
    num_energy_neighbor_bins = num_in_each_dir + 1;
    num_pitch_neighbor_bins = num_in_each_dir + 1;
  }
};

/**
 * @brief Calculates the radial bin index for a given radial coordinate
 * @param mrho_max Maximum number of radial bins
 * @param query_radial_coord The radial coordinate to find the bin for
 * (normalized, [0,1])
 * @return The index of the radial bin containing the query point
 */
OMEGA_H_INLINE
size_t search_radial_bin(const size_t mrho_max, const double query_radial_coord)
{
  double radial_bin_id = mrho_max * query_radial_coord - 0.5;
  return static_cast<size_t>(radial_bin_id < 0.0 ? 0.0 : radial_bin_id);
}

OMEGA_H_INLINE
double compute_radial_coordinate(const size_t mrho_max,
                                 const size_t radial_bin_id)
{
  return double(0.5 + radial_bin_id) / mrho_max;
}

/**
 * @brief Calculates the poloidal bin index for a given query coordinate
 * @param radial_bin_number The input radial bin number used to determine the
 * maximum number of poloidal bins
 * @param query_poloidal_coord The poloidal coordinate for which the bin is
 * computed
 * @return The index of the poloidal bin corresponding to the given coordinate
 */
OMEGA_H_INLINE
size_t search_poloidal_bin(const size_t radial_bin_number,
                           const double query_poloidal_coord)
{
  size_t mpol_max = 4 * (radial_bin_number + 1);
  double dtheta = 2 * M_PI / mpol_max;
  double poloidal_bin_id = query_poloidal_coord / dtheta - 0.5;
  return static_cast<size_t>(poloidal_bin_id < 0.0 ? 0.0 : poloidal_bin_id);
}

OMEGA_H_INLINE
double compute_poloidal_coordinate(const size_t radial_bin_number,
                                   const size_t poloidal_bin_id)
{
  size_t mpol_max = 4 * (radial_bin_number + 1);
  double dtheta = 2 * M_PI / mpol_max;
  return (poloidal_bin_id + 0.5) * dtheta;
}

/**
 * @brief Calculates the toroidal bin index for a given query coordinate
 * @param mtor_max Maximum number of toroidal bins
 * @param zeta_max Maximum value of the toroidal coordinate
 * @param query_toroidal_coord The toroidal coordinate to find the bin for
 * @return The index of the toroidal bin containing the query point
 */
OMEGA_H_INLINE
size_t search_toroidal_bin(size_t mtor_max, double zeta_max,
                           double query_toroidal_coord)
{
  // width of each toroidal bin
  const double dzeta = zeta_max / static_cast<double>(mtor_max);

  double toroidal_bin_id = query_toroidal_coord / dzeta - 0.5;

  return static_cast<size_t>(toroidal_bin_id < 0.0 ? 0.0 : toroidal_bin_id);
}

OMEGA_H_INLINE
double compute_toroidal_coordinate(const size_t mtor_max, const double zeta_max,
                                   const size_t toroidal_bin_id)
{
  double dzeta = static_cast<double>(zeta_max / mtor_max);
  return static_cast<double>((toroidal_bin_id + 0.5) * dzeta);
}

/**
 * @brief Calculates the energy bin index for a given query coordinate
 * @param meng_max Maximum number of energy bins
 * @param emax Maximum value of the energy coordinate
 * @param query_energy_coord The energy coordinate to find the bin for
 * @return The index of the energy bin containing the query point
 */
OMEGA_H_INLINE
size_t search_energy_bin(const size_t meng_max, const double emax,
                         const double query_energy_coord)
{
  double energy_bin_id = meng_max * query_energy_coord / emax - 0.5;

  return static_cast<size_t>(energy_bin_id < 0.0 ? 0.0 : energy_bin_id);
}

OMEGA_H_INLINE
double compute_energy_coordinate(const size_t meng_max, const double emax,
                                 const size_t energy_bin_id)
{
  return emax * (energy_bin_id + 0.5) / meng_max;
}

/**
 * @brief Calculates the pitch bin index for a given query coordinate
 * @param mlmd_max Maximum number of pitch bins
 * @param query_pitch_coord The pitch coordinate to find the bin for
 * @return The index of the pitch bin containing the query point
 */
OMEGA_H_INLINE
size_t search_pitch_bin(const size_t mlmd_max, const double query_pitch_coord)
{
  double pitch_bin_id = 0.5 * mlmd_max * (query_pitch_coord + 1) - 0.5;

  return static_cast<size_t>(pitch_bin_id < 0.0 ? 0.0 : pitch_bin_id);
}

OMEGA_H_INLINE
double compute_pitch_coordinate(const size_t mlmd_max,
                                const size_t pitch_bin_id)
{
  return 2 * (pitch_bin_id + 0.5) / mlmd_max - 1.0;
}

OMEGA_H_INLINE
int find_rho_binary_search(const Omega_h::LOs& pol_max_scan_array,
                           const int theta_total)
{
  int left = 0;
  int right = pol_max_scan_array.size() - 1;

  while (left < right) {
    int mid = left + (right - left) / 2;

    if (pol_max_scan_array[mid] <= theta_total &&
        theta_total < pol_max_scan_array[mid + 1]) {
      return mid;
    }

    if (pol_max_scan_array[mid] > theta_total) {
      right = mid;
    } else {
      left = mid + 1;
    }
  }
  return left;
}

/**
 * @brief Evaluates the grid ID based on the bin indices
 * @param mpol_sum_at_radial_bin_num Sum of poloidal bins up to the given radial
 * bin number
 * @param poloidal_bin_number The poloidal bin number
 * @param toroidal_bin_number The toroidal bin number
 * @param mtor_max Maximum number of toroidal bins
 * @return The grid ID corresponding to the given bin indices
 */
OMEGA_H_INLINE
size_t eval_grid_id(size_t mpol_sum_at_radial_bin_num,
                    size_t poloidal_bin_number, size_t toroidal_bin_number,
                    size_t mtor_max)
{
  size_t grid_id =
    (mpol_sum_at_radial_bin_num + poloidal_bin_number) * mtor_max +
    toroidal_bin_number;
  return grid_id;
}

OMEGA_H_INLINE
size_t eval_marker_id(size_t mlmd_max, size_t meng_max, size_t grid_id,
                      size_t pitch_bin_number, size_t energy_bin_number)
{
  size_t marker_id = grid_id * mlmd_max * meng_max +
                     energy_bin_number * mlmd_max + pitch_bin_number;
  return marker_id;
}

// bin_type = 0 for linear; bin_type = 1 for periodic
OMEGA_H_INLINE
void get_bin_neighbors(size_t max_bin_number, size_t immediate_bin_id,
                       size_t num_neighbor_bins, int bin_type,
                       Omega_h::LO* bin_neighbors)
{
  bin_neighbors[0] = immediate_bin_id;

  size_t left = immediate_bin_id, right = immediate_bin_id;
  size_t count = 1;

  while (count < num_neighbor_bins) {
    if (count % 2 != 1) {
      if (bin_type == 1) {
        left = (left - 1 + max_bin_number) % max_bin_number;
        bin_neighbors[count] = left;

      } else {

        if (left > 0) {
          left -= 1;
          bin_neighbors[count] = left;

        } else {
          right += 1;
          bin_neighbors[count] = right;
        }
      }

    } else {
      if (bin_type == 1) {
        right = (right + 1) % max_bin_number;
        bin_neighbors[count] = right;

      } else {
        if (right < max_bin_number - 1) {
          right += 1;
          bin_neighbors[count] = right;

        } else {
          left -= 1;
          bin_neighbors[count] = left;
        }
      }
    }

    count++;
  }
}

/**
 * @brief Gets the radial bins neighboring a query point
 * @param num_rho_neighbor_bins Number of neighboring radial bins to find
 * @param mrho_max Maximum number of radial bins
 * @param query_radial_coord The radial coordinate of the query point
 * @return A list of neighboring radial bin indices
 *
 * @notes bin_type for radial bins = 0 // 0 for linear discretization
 *                                   // 1 for periodic discretization
 */
OMEGA_H_INLINE
void get_radial_bin_neighbors(size_t num_rho_neighbor_bins, size_t mrho_max,
                              double query_radial_coord,
                              Omega_h::LO* radial_bin_neighbors)
{
  size_t radial_bin_id = search_radial_bin(mrho_max, query_radial_coord);
  int bin_type = 0;

  get_bin_neighbors(mrho_max, radial_bin_id, num_rho_neighbor_bins, bin_type,
                    radial_bin_neighbors);
}

/**
 * @brief Gets the poloidal bins neighboring a query point
 * @param num_poloidal_neighbor_bins Number of neighboring poloidal bins to
 * find
 * @param radial_bin_number The radial bin number of the query point
 * @param query_poloidal_coord The poloidal coordinate of the query point
 * @return A list of neighboring poloidal bin indices
 *
 * @notes bin_type for poloidal bins = 1
 *
 */
OMEGA_H_INLINE
void get_poloidal_bin_neighbors(size_t num_poloidal_neighbor_bins,
                                size_t radial_bin_number,
                                double query_poloidal_coord,
                                Omega_h::LO* poloidal_bin_neighbors)
{
  size_t poloidal_bin_id =
    search_poloidal_bin(radial_bin_number, query_poloidal_coord);

  size_t mpol_max = 4 * (radial_bin_number + 1);

  int bin_type = 1;

  get_bin_neighbors(mpol_max, poloidal_bin_id, num_poloidal_neighbor_bins,
                    bin_type, poloidal_bin_neighbors);
}

/**
 * @brief Gets the toroidal bins neighboring a query point
 * @param num_toroidal_neighbor_bins Number of neighboring toroidal bins to
 * find
 * @param mtor_max Maximum number of toroidal bins
 * @param zeta_max Maximum value of the toroidal coordinate
 * @param query_toroidal_coord The toroidal coordinate of the query point
 * @return A list of neighboring toroidal bin indices
 *
 *
 * @notes bin_type for toroidal bins = 1
 */
OMEGA_H_INLINE
void get_toroidal_bin_neighbors(size_t num_toroidal_neighbor_bins,
                                size_t mtor_max, double zeta_max,
                                double query_toroidal_coord,
                                Omega_h::LO* toroidal_bin_neighbors)
{
  size_t toroidal_bin_id =
    search_toroidal_bin(mtor_max, zeta_max, query_toroidal_coord);

  int bin_type = 1;

  get_bin_neighbors(mtor_max, toroidal_bin_id, num_toroidal_neighbor_bins,
                    bin_type, toroidal_bin_neighbors);
}

/**
 * @brief Gets the pitch bins neighboring a query point
 * @param num_pitch_neighbor_bins Number of neighboring pitch bins to find
 * @param mlmd_max Maximum number of pitch bins
 * @param query_pitch_coord The pitch coordinate of the query point
 * @return A list of neighboring pitch bin indices
 *
 *
 * @notes bin_name_id for pitch bins = 3
 */
OMEGA_H_INLINE
void get_pitch_bin_neighbors(size_t num_pitch_neighbor_bins, size_t mlmd_max,
                             double query_pitch_coord,
                             Omega_h::LO* pitch_bin_neighbors)
{
  size_t pitch_bin_id = search_pitch_bin(mlmd_max, query_pitch_coord);

  int bin_type = 0;

  get_bin_neighbors(mlmd_max, pitch_bin_id, num_pitch_neighbor_bins, bin_type,
                    pitch_bin_neighbors);
}

/**
 * @brief Gets the energy bins neighboring a query point
 * @param num_energy_neighbor_bins Number of neighboring energy bins to find
 * @param meng_max Maximum number of energy bins
 * @param query_energy_coord The energy coordinate of the query point
 * @return A list of neighboring energy bin indices
 *
 *  @notes bin_type for energy bins = 0
 */
OMEGA_H_INLINE
void get_energy_bin_neighbors(size_t num_energy_neighbor_bins, size_t meng_max,
                              double emax, double query_energy_coord,
                              Omega_h::LO* energy_bin_neighbors)
{

  size_t energy_bin_id = search_energy_bin(meng_max, emax, query_energy_coord);

  int bin_type = 0;

  get_bin_neighbors(meng_max, energy_bin_id, num_energy_neighbor_bins, bin_type,
                    energy_bin_neighbors);
}

OMEGA_H_INLINE
void get_grid_neighbors(
  size_t num_rho_neighbor_bins, Omega_h::LO* radial_bin_neighbors,
  size_t num_theta_neighbor_bins, double query_point_theta,
  Omega_h::LO* poloidal_bin_neighbors, size_t num_zeta_neighbor_bins,
  Omega_h::LO* toroidal_bin_neighbors, Omega_h::LOs pol_max_scan_array,
  size_t mtor_max, Omega_h::LO* grid_neighbors_ids)
{

  int count = -1;

  for (size_t i = 0; i < num_rho_neighbor_bins; ++i) {

    size_t current_radial_bin_number = radial_bin_neighbors[i];

    size_t mpol_sum_at_current_radial =
      pol_max_scan_array[current_radial_bin_number];

    get_poloidal_bin_neighbors(num_theta_neighbor_bins,
                               current_radial_bin_number, query_point_theta,
                               poloidal_bin_neighbors);

    for (size_t j = 0; j < num_theta_neighbor_bins; ++j) {

      size_t current_poloidal_bin_number = poloidal_bin_neighbors[j];

      for (size_t k = 0; k < num_zeta_neighbor_bins; ++k) {

        count += 1;

        size_t current_toroidal_bin_number = toroidal_bin_neighbors[k];

        size_t current_grid_id =
          eval_grid_id(mpol_sum_at_current_radial, current_poloidal_bin_number,
                       current_toroidal_bin_number, mtor_max);

        grid_neighbors_ids[count] = current_grid_id;
      }
    }
  }
}

OMEGA_H_INLINE
void get_marker_neighbors(size_t grid_neighbors_ids_size,
                          Omega_h::LO* grid_neighbors_ids,
                          size_t num_pitch_neighbor_bins, size_t mlmd_max,
                          Omega_h::LO* pitch_bin_neighbors,
                          size_t num_energy_neighbor_bins, size_t meng_max,
                          Omega_h::LO* energy_bin_neighbors,
                          Omega_h::LO* support_idx)
{

  int count = -1;

  for (size_t i = 0; i < grid_neighbors_ids_size; ++i) {

    size_t current_grid_id = grid_neighbors_ids[i];

    for (size_t j = 0; j < num_energy_neighbor_bins; ++j) {

      size_t current_energy_bin_number = energy_bin_neighbors[j];

      for (size_t k = 0; k < num_pitch_neighbor_bins; ++k) {

        size_t current_pitch_bin_number = pitch_bin_neighbors[k];

        count += 1;

        auto current_marker_id =
          eval_marker_id(mlmd_max, meng_max, current_grid_id,
                         current_pitch_bin_number, current_energy_bin_number);

        support_idx[count] = current_marker_id;
      }
    }
  }
}

/**
 * @brief Class for performing Gnet search operations
 */
class GnetSearch
{

private:
  Omega_h::Real query_point_[MAX_GNET_DIM];

public:
  /**
   * @brief Constructor for GnetSearch
   * @param query_point The 5D query point (rho, theta, zeta, energy, pitch)
   */
  OMEGA_H_INLINE
  GnetSearch(Omega_h::Real query_point[MAX_GNET_DIM])
  {
    for (int i = 0; i < MAX_GNET_DIM; ++i) {
      query_point_[i] = query_point[i];
    }
  }

  /**
   * @brief Destructor for GnetSearch
   */
  OMEGA_H_INLINE
  ~GnetSearch() {}

  /**
   * @brief Searches for neighboring bins for the query points
   * @return A list of neighboring bin indices
   */
  OMEGA_H_INLINE
  void analyticSearchGnet(size_t mrho_max, size_t mtor_max, size_t mlmd_max,
                          size_t meng_max, double zeta_max, double emax,
                          size_t min_num_supports,
                          Omega_h::LOs pol_max_scan_array,
                          Omega_h::LO* support_idx)
  {

    BinNeighbors binNeighbors{};

    binNeighbors.eval(min_num_supports);

    double query_point_rho = query_point_[0];
    double query_point_theta = query_point_[1];
    double query_point_zeta = query_point_[2];
    double query_point_energy = query_point_[3];
    double query_point_pitch = query_point_[4];

    Omega_h::LO radial_bin_neighbors[MAX_NEIGHBORS_EACH_DIR] = {};

    Omega_h::LO toroidal_bin_neighbors[MAX_NEIGHBORS_EACH_DIR] = {};

    Omega_h::LO poloidal_bin_neighbors[MAX_NEIGHBORS_EACH_DIR] = {};

    get_radial_bin_neighbors(binNeighbors.num_rho_neighbor_bins, mrho_max,
                             query_point_rho, radial_bin_neighbors);

    get_toroidal_bin_neighbors(binNeighbors.num_zeta_neighbor_bins, mtor_max,
                               zeta_max, query_point_zeta,
                               toroidal_bin_neighbors);

    size_t grid_neighbors_ids_size = binNeighbors.num_rho_neighbor_bins *
                                     binNeighbors.num_theta_neighbor_bins *
                                     binNeighbors.num_zeta_neighbor_bins;

    Omega_h::LO grid_neighbors_ids[MAX_GRID_NEIGHBORS] = {};

    Omega_h::LO pitch_bin_neighbors[MAX_NEIGHBORS_EACH_DIR] = {};

    Omega_h::LO energy_bin_neighbors[MAX_NEIGHBORS_EACH_DIR] = {};

    get_pitch_bin_neighbors(binNeighbors.num_pitch_neighbor_bins, mlmd_max,
                            query_point_pitch, pitch_bin_neighbors);

    get_energy_bin_neighbors(binNeighbors.num_energy_neighbor_bins, meng_max,
                             emax, query_point_energy, energy_bin_neighbors);

    get_grid_neighbors(
      binNeighbors.num_rho_neighbor_bins, radial_bin_neighbors,
      binNeighbors.num_theta_neighbor_bins, query_point_theta,
      poloidal_bin_neighbors, binNeighbors.num_zeta_neighbor_bins,
      toroidal_bin_neighbors, pol_max_scan_array, mtor_max, grid_neighbors_ids);

    get_marker_neighbors(grid_neighbors_ids_size, grid_neighbors_ids,
                         binNeighbors.num_pitch_neighbor_bins, mlmd_max,
                         pitch_bin_neighbors,
                         binNeighbors.num_energy_neighbor_bins, meng_max,
                         energy_bin_neighbors, support_idx);
  }
};

/**
 * @brief Coordinates of every GNET marker (cell center), flattened as
 * [rho, theta, zeta, energy, pitch] per marker, in marker-id order.
 */
Omega_h::Write<Omega_h::Real> get_gnet_coordinates(
  const size_t& mrho_max, const size_t& mtor_max, const size_t& mlmd_max,
  const size_t& meng_max, const double& zeta_max, const double& emax,
  const Omega_h::LOs& pol_max_scan_array, const size_t& npoints,
  const size_t& dim);

/// Number of poloidal bins for each radial bin: 4 * (rho + 1)
Omega_h::LOs get_mpol_max_array(size_t rho_max);

/// Exclusive prefix sum of the poloidal bin counts (size rho_max + 1)
Omega_h::LOs get_pol_max_scan_array(size_t mrho_max);

/// Size of the MLS polynomial basis for the given degree and dimension
size_t compute_min_num_supports(const size_t degree, const size_t dim);

/**
 * @brief Finds the GNET supports (neighbor markers) of each target point
 * @param target_points Flattened 5D target coordinates
 * @return Supports in CSR format, as expected by pcms::mls_interpolation
 */
pcms::SupportResults searchNeighbors(
  const Omega_h::Reals& target_points, const size_t mrho_max,
  const size_t mtor_max, const size_t mlmd_max, const size_t meng_max,
  const double zeta_max, const double emax, size_t min_num_supports);

} // namespace pcms::gnet

#endif // PCMS_GNET_GNET_GRID_HPP
