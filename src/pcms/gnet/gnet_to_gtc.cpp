#include <pcms/gnet/gnet_to_gtc.hpp>
#include <pcms/gnet/gnet_diagnostics.hpp>
#include <pcms/gnet/gnet_io.hpp>
#include <pcms/gnet/gnet_sampling.hpp>
#include <Kokkos_Timer.hpp>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace pcms::gnet
{

namespace
{

constexpr int dim = MAX_GNET_DIM;

HostMatrix to_host(const Kokkos::View<double**, Kokkos::LayoutRight>& view)
{
  HostMatrix host("host copy", view.extent(0), view.extent(1));
  Kokkos::deep_copy(host, view);
  return host;
}

void print_options(const GnetToGtcOptions& o)
{
  std::cout << "GNET -> GTC\n"
            << "  input           : " << o.input << "\n"
            << "  output          : " << o.output << "\n"
            << "  gnet histogram  : "
            << (o.gnet_histogram.empty() ? "(not written)" : o.gnet_histogram)
            << "\n"
            << "  samples         : " << o.samples << "\n"
            << "  batch size      : " << o.batch_size << "\n"
            << "  seed            : " << o.seed << "\n"
            << "  MLS degree      : " << o.degree << "\n"
            << "  MLS lambda      : " << o.lambda << "\n"
            << "  threshold scale : " << o.threshold_scale << "\n";
}

} // namespace

pcms::UniformGrid<MAX_GNET_DIM> default_gtc_grid()
{
  pcms::UniformGrid<MAX_GNET_DIM> grid5d;
  grid5d.edge_length = {1.0, 2.0 * M_PI, 2.0 * M_PI / 5.0, 198000.0, 2.0};
  grid5d.bot_left = {0, 0, 0, 2000.0, -1.0};
  grid5d.divisions = {10, 20, 60, 50, 33};
  return grid5d;
}

GnetData load_gnet(const std::string& path, bool verbose)
{
  GnetHDF5Reader reader(path);

  GnetDims dims{static_cast<size_t>(reader.getMaxRhoBins()),
                static_cast<size_t>(reader.getMaxToroidalBins()),
                static_cast<size_t>(reader.getMaxPitchBins()),
                static_cast<size_t>(reader.getMaxEnergyBins()),
                reader.getMaxZetaRad(),
                reader.getMaxEnergyEV()};

  auto mpol_max = get_mpol_max_array(dims.mrho_max);
  auto pol_max_scan_array = get_pol_max_scan_array(dims.mrho_max);
  Omega_h::HostRead<Omega_h::LO> host_pol_max_scan_array(pol_max_scan_array);

  int grid_total = host_pol_max_scan_array[dims.mrho_max] * dims.mtor_max;
  int marker_total = grid_total * dims.mlmd_max * dims.meng_max;

  auto host_f5d = reader.getF5D();
  auto host_volume3d = reader.getVolume3D();
  if (static_cast<int>(host_f5d.extent(0)) != grid_total ||
      host_f5d.extent(1) != dims.meng_max ||
      host_f5d.extent(2) != dims.mlmd_max) {
    throw std::runtime_error("f5d dimensions do not match the GNET grid in " +
                             path);
  }

  std::cout << "GNET grid: rho " << dims.mrho_max << ", zeta "
            << dims.mtor_max << " (max " << dims.zeta_max << " rad), energy "
            << dims.meng_max << " (max " << dims.emax << " eV), pitch "
            << dims.mlmd_max << "\n"
            << "  grid points (rho, theta, zeta): " << grid_total
            << ", markers: " << marker_total << "\n";

  auto source_coordinates =
    get_gnet_coordinates(dims.mrho_max, dims.mtor_max, dims.mlmd_max,
                         dims.meng_max, dims.zeta_max, dims.emax,
                         pol_max_scan_array, marker_total, dim);

  auto source_values = calculate_number_density(
    dims.mrho_max, mpol_max, dims.mtor_max, pol_max_scan_array, dims.meng_max,
    dims.mlmd_max, host_f5d, host_volume3d, marker_total);

  if (verbose) {
    Omega_h::HostRead<Omega_h::Real> host_coords(source_coordinates);
    Omega_h::HostRead<Omega_h::Real> host_values(source_values);
    std::cout << "  first source points [rho, theta, zeta, energy, pitch] "
                 "-> number density:\n";
    for (int i = 0; i < std::min(10, marker_total); ++i) {
      std::cout << "    [";
      for (int d = 0; d < dim; ++d) {
        std::cout << host_coords[i * dim + d] << (d < dim - 1 ? ", " : "");
      }
      std::cout << "] -> " << host_values[i] << "\n";
    }
  }

  return GnetData{dims, marker_total, Omega_h::read(source_coordinates),
                  source_values, host_volume3d};
}

void run_gnet_to_gtc(const GnetToGtcOptions& options)
{
  print_options(options);
  Kokkos::Timer total_timer;
  Kokkos::Timer timer;

  // 1. read GNET data and build MLS source points
  auto gnet = load_gnet(options.input, options.verbose);
  Kokkos::fence();
  std::cout << "  [load]       " << timer.seconds() << " s\n";

  auto grid5d = default_gtc_grid();

  const double min_v3d = get_min_volume(gnet.volume3d);
  const double max_v3d = get_max_volume(gnet.volume3d);
  const double threshold = options.threshold_scale / max_v3d;
  std::cout << " Volume3D min: " << min_v3d << ", max: " << max_v3d
            << ", rejection threshold: " << threshold << "\n";

  // 2. optional GNET energy-pitch histogram
  if (!options.gnet_histogram.empty()) {
    timer.reset();
    auto gnet_epitch = calculate_velocity_distribution(
      gnet.source_values, gnet.source_coordinates, grid5d);
    write_gnet_histogram(options.gnet_histogram, to_host(gnet_epitch));
    std::cout << "  [gnet hist]  " << timer.seconds() << " s -> "
              << options.gnet_histogram << "\n";
  }

  // 3. rejection sampling of GTC markers
  timer.reset();
  SamplingParams params{options.samples, options.batch_size, threshold,
                        options.degree,  options.lambda,     options.seed};
  auto markers = sample_gtc_markers(gnet.source_values,
                                    gnet.source_coordinates, grid5d, gnet.dims,
                                    params);
  Kokkos::fence();
  std::cout << "  [sampling]   " << timer.seconds() << " s\n";

  // 4. GTC histograms and output table
  timer.reset();
  Omega_h::Reals values(markers.values);
  Omega_h::Reals coordinates(markers.coordinates);
  auto table = get2DArray(values, coordinates);
  auto epitch = calculate_velocity_distribution(values, coordinates, grid5d);
  auto rtheta = calculate_2d_distribution(values, coordinates, grid5d);
  write_gtc_output(options.output, to_host(table), to_host(epitch),
                   to_host(rtheta));
  std::cout << "  [write]      " << timer.seconds() << " s -> "
            << options.output << "\n";

  std::cout << "GNET to GTC Mapping Completed Successfully! ("
            << total_timer.seconds() << " s)\n";
}

} // namespace pcms::gnet
