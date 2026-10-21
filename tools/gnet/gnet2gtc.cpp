// gnet2gtc: map a GNET 5D distribution to GTC markers.
//
// Reads a GNET HDF5 file, interpolates it with MLS onto random candidate
// points, keeps candidates above a threshold (rejection sampling) and writes
// the GTC markers to an HDF5 file.

#include <pcms/gnet/gnet_to_gtc.hpp>
#include <Kokkos_Core.hpp>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>

namespace
{

void print_usage(const char* prog)
{
  const pcms::gnet::GnetToGtcOptions d;
  std::cout
    << "Usage: " << prog << " --input <gnet.h5> --output <gtc.h5> [options]\n"
    << "\n"
    << "Required:\n"
    << "  --input <file>            GNET HDF5 file\n"
    << "  --output <file>           GTC output HDF5 file (overwritten)\n"
    << "\n"
    << "Options:\n"
    << "  --samples <n>             GTC markers to produce (default " << d.samples << ")\n"
    << "  --batch <n>               candidates per sampling round (default " << d.batch_size << ")\n"
    << "  --seed <n>                random seed (default " << d.seed << ")\n"
    << "  --degree <n>              MLS polynomial degree (default " << d.degree << ")\n"
    << "  --lambda <x>              MLS regularization (default " << d.lambda << ")\n"
    << "  --threshold-scale <x>     threshold = scale / max(volume3d) (default " << d.threshold_scale << ")\n"
    << "  --verbose                 print sample source points\n"
    << "  --help                    show this message\n"
    << "\n"
    << "Kokkos options (e.g. --kokkos-device-id=<n>) are also accepted.\n";
}

// Parses arguments into options. Returns false on --help or error.
bool parse_args(int argc, char* argv[], pcms::gnet::GnetToGtcOptions& o)
{
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    auto value = [&]() -> std::string {
      if (i + 1 >= argc) {
        throw std::runtime_error("missing value for " + arg);
      }
      return argv[++i];
    };

    try {
      if (arg == "--help" || arg == "-h") {
        return false;
      } else if (arg == "--input") {
        o.input = value();
      } else if (arg == "--output") {
        o.output = value();
      } else if (arg == "--samples") {
        o.samples = std::stoi(value());
      } else if (arg == "--batch") {
        o.batch_size = std::stoi(value());
      } else if (arg == "--seed") {
        o.seed = std::stoull(value());
      } else if (arg == "--degree") {
        o.degree = std::stoi(value());
      } else if (arg == "--lambda") {
        o.lambda = std::stod(value());
      } else if (arg == "--threshold-scale") {
        o.threshold_scale = std::stod(value());
      } else if (arg == "--verbose") {
        o.verbose = true;
      } else {
        throw std::runtime_error("unknown option " + arg);
      }
    } catch (const std::logic_error&) { // std::stoi/stod failures
      throw std::runtime_error("invalid value for " + arg);
    }
  }

  if (o.input.empty() || o.output.empty()) {
    throw std::runtime_error("--input and --output are required");
  }
  // marker coordinates are indexed with int: samples * 5 must fit
  if (o.samples <= 0 || o.samples > 400000000) {
    throw std::runtime_error("--samples must be in (0, 400000000]");
  }
  if (o.batch_size <= 0 || o.batch_size > 400000000) {
    throw std::runtime_error("--batch must be in (0, 400000000]");
  }
  if (o.degree < 1) {
    throw std::runtime_error("--degree must be >= 1");
  }
  return true;
}

} // namespace

int main(int argc, char* argv[])
{
  // handle --help before Kokkos, which prints its own help for it
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--help" || arg == "-h") {
      print_usage(argv[0]);
      return EXIT_SUCCESS;
    }
  }

  Kokkos::ScopeGuard guard(argc, argv);

  pcms::gnet::GnetToGtcOptions options;
  try {
    if (!parse_args(argc, argv, options)) {
      print_usage(argv[0]);
      return EXIT_SUCCESS;
    }
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << "\n\n";
    print_usage(argv[0]);
    return EXIT_FAILURE;
  }

  try {
    pcms::gnet::run_gnet_to_gtc(options);
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << "\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
