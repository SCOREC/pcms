// Uses PCMS through its installed public headers only: builds two Lagrange
// spaces on uniform grids, interpolates a constant field from one to the
// other, and checks the result. It also includes a C API header to check that
// the installed C headers are usable.

#include <pcms/pcms.hpp>
#include <pcms/capi/interpolator.h>

#include <Kokkos_Core.hpp>
#include <mpi.h>

#include <cmath>
#include <cstddef>
#include <vector>

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  int mismatches = 0;
  {
    Kokkos::ScopeGuard kokkos(argc, argv);

    // Two uniform grids over the unit square.
    const pcms::UniformGrid<2> coarse{
      .edge_length = {1.0, 1.0}, .bot_left = {0.0, 0.0}, .divisions = {4, 4}};
    const pcms::UniformGrid<2> fine{
      .edge_length = {1.0, 1.0}, .bot_left = {0.0, 0.0}, .divisions = {9, 7}};
    const auto cartesian = pcms::csys::Cartesian::Create(2);

    auto source_space = pcms::CreateLagrangeSpace(
      pcms::CreateUniformGridDiscretization(coarse), 1, 1, cartesian);
    auto target_space = pcms::CreateLagrangeSpace(
      pcms::CreateUniformGridDiscretization(fine), 1, 1, cartesian);

    auto source = source_space->CreateFunction<pcms::Real>("source");
    auto target = target_space->CreateFunction<pcms::Real>("target");

    // A constant field is reproduced exactly by linear interpolation.
    const pcms::Real value = 2.5;
    const auto num_source =
      source.GetLayout().GetDOFHolderCoordinates().GetValues().extent(0);
    std::vector<pcms::Real> values(num_source, value);
    source.SetDOFHolderDataUncheckedHost(
      pcms::Rank2View<const pcms::Real, pcms::HostMemorySpace>(
        values.data(), num_source, 1));

    pcms::Interpolator<pcms::Real> interpolate(*source_space, *target_space);
    interpolate.Apply(source, target);

    const auto result = target.GetDOFHolderDataHost().GetValues();
    for (std::size_t i = 0; i < result.extent(0); ++i) {
      if (std::abs(result.data_handle()[i] - value) > 1e-12) {
        ++mismatches;
      }
    }
    pcms::printInfo("interpolated %d values, %d mismatches\n",
                    static_cast<int>(result.extent(0)), mismatches);
  }
  MPI_Finalize();
  return mismatches == 0 ? 0 : 1;
}
