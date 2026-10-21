#ifndef PCMS_GNET_GNET_IO_HPP
#define PCMS_GNET_GNET_IO_HPP

#include <Kokkos_Core.hpp>
#include <Omega_h_array.hpp>
#include <cstdint>
#include <string>

namespace pcms::gnet
{

/**
 * @brief Reads a GNET 5D distribution file (HDF5).
 *
 * Expected layout:
 *   /distribution5D/dimensions/{mxenergy, mxrho, mxzeta, mxpitch}  (int)
 *   /distribution5D/range/{engmax_eV, zetamax_rad}                 (double)
 *   /distribution5D/data/f5d       (grid, energy, pitch)           (double)
 *   /distribution5D/data/volume3d  (grid)                          (double)
 *
 * All data is read in the constructor; throws std::runtime_error if the file
 * or a dataset cannot be read.
 */
class GnetHDF5Reader
{
public:
  explicit GnetHDF5Reader(const std::string& filename);
  ~GnetHDF5Reader();

  GnetHDF5Reader(const GnetHDF5Reader&) = delete;
  GnetHDF5Reader& operator=(const GnetHDF5Reader&) = delete;

  int getMaxEnergyBins() const { return meng_max; }
  int getMaxRhoBins() const { return mrho_max; }
  int getMaxToroidalBins() const { return mtor_max; }
  int getMaxPitchBins() const { return mlmd_max; }
  double getMaxEnergyEV() const { return emax; }
  double getMaxZetaRad() const { return zeta_max; }

  Kokkos::View<double***, Kokkos::HostSpace> getF5D() const { return host_f5d; }
  Kokkos::View<double*, Kokkos::HostSpace> getVolume3D() const
  {
    return host_volume;
  }

private:
  std::int64_t file_id; // hid_t

  int meng_max, mrho_max, mtor_max, mlmd_max;
  double emax, zeta_max;

  Kokkos::View<double***, Kokkos::HostSpace> host_f5d;
  Kokkos::View<double*, Kokkos::HostSpace> host_volume;
};

/**
 * @brief Writes 2D datasets to a new HDF5 file (an existing file is
 * overwritten).
 */
class HDF5Writer
{
public:
  explicit HDF5Writer(const std::string& filename);
  ~HDF5Writer();

  HDF5Writer(const HDF5Writer&) = delete;
  HDF5Writer& operator=(const HDF5Writer&) = delete;

  void writef2D(const std::string& dataset_name,
                const Kokkos::View<double**, Kokkos::HostSpace>& data);

private:
  std::int64_t file_id; // hid_t
};

using HostMatrix = Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace>;

/// Writes the GTC marker file: /gtc_n6d, N x 6 (rho, theta, zeta, energy,
/// pitch, value)
void write_gtc_output(const std::string& path, const HostMatrix& markers);

/**
 * @brief Number density f5d / volume3d for every GNET marker, in marker-id
 * order (grid, energy, pitch), on the device.
 */
Omega_h::Reals calculate_number_density(
  const size_t mrho_max, const Omega_h::LOs& mpol_max, const size_t mtor_max,
  const Omega_h::LOs& pol_max_scan_array, const size_t meng_max,
  const size_t mlmd_max, const Kokkos::View<double***, Kokkos::HostSpace>& f5d,
  const Kokkos::View<double*, Kokkos::HostSpace>& volume3d,
  const int total_num_particles);

} // namespace pcms::gnet

#endif // PCMS_GNET_GNET_IO_HPP
