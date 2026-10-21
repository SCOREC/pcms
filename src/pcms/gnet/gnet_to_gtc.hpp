#ifndef PCMS_GNET_GNET_TO_GTC_HPP
#define PCMS_GNET_GNET_TO_GTC_HPP

#include <pcms/gnet/gnet_grid.hpp>
#include <pcms/utility/uniform_grid.h>
#include <Kokkos_Core.hpp>
#include <Omega_h_array.hpp>
#include <cstdint>
#include <string>

namespace pcms::gnet
{

/// Settings for a GNET -> GTC run (defaults match the reference setup)
struct GnetToGtcOptions
{
  std::string input;          // GNET HDF5 file (required)
  std::string output;         // GTC HDF5 output file (required)
  int samples = 30000000;     // GTC markers to produce
  int batch_size = 1000000;   // candidates per sampling round
  std::uint64_t seed = 14000;
  int degree = 1;               // MLS polynomial degree
  double lambda = 1e-4;         // MLS regularization
  double threshold_scale = 1e10; // threshold = threshold_scale / max(volume3d)
  bool verbose = false;
};

/// GNET data prepared as MLS source points
struct GnetData
{
  GnetDims dims;
  int marker_total;
  Omega_h::Reals source_coordinates; // marker_total * 5
  Omega_h::Reals source_values;      // number density per marker
  Kokkos::View<double*, Kokkos::HostSpace> volume3d;
};

/**
 * @brief Sampling domain for GTC markers:
 * rho [0,1], theta [0,2pi], zeta [0,2pi/5], energy [2000,200000] eV,
 * pitch [-1,1].
 */
pcms::UniformGrid<MAX_GNET_DIM> default_gtc_grid();

/// Reads a GNET file and builds the source points and number density
GnetData load_gnet(const std::string& path, bool verbose = false);

/// Runs the full GNET -> GTC mapping and writes the GTC marker file
void run_gnet_to_gtc(const GnetToGtcOptions& options);

} // namespace pcms::gnet

#endif // PCMS_GNET_GNET_TO_GTC_HPP
