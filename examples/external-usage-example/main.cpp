// Uses PCMS through its installed public headers only: builds the geometry of
// an Omega_h mesh and a Lagrange space on it, evaluates a linear function at
// points with QueryPoints, and interpolates it onto a uniform grid. It also
// includes a C API header to check that the installed C headers are usable.

#include <pcms/backends/uniform_grid/geometry.hpp>
#include <pcms/pcms.hpp>
#include <pcms/backends/omega_h/omega_h.hpp>
#include <pcms/capi/interpolator.h>

#include <Kokkos_Core.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_library.hpp>
#include <mpi.h>

#include <cmath>
#include <cstdio>
#include <cstddef>
#include <vector>

namespace
{
constexpr pcms::Real LinearTol = 1e-10;

pcms::Real Linear(pcms::Real x, pcms::Real y)
{
  return 1.0 + 2.0 * x + 3.0 * y;
}
} // namespace

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  int mismatches = 0;
  {
    Kokkos::ScopeGuard kokkos(argc, argv);
    Omega_h::Library lib{};
    auto mesh =
      Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1.0, 1.0, 0.0, 4, 4, 0);
    const auto cartesian = pcms::csys::CreateCartesian(2);

    // The geometry of the mesh, and a linear Lagrange space on it.
    auto geometry = pcms::omega_h::CreateGeometry(mesh, cartesian);
    auto space = pcms::CreateLagrangeSpace(geometry, 1, 1);
    auto source = space->CreateFunction<pcms::Real>("source");

    // source(x, y) = 1 + 2x + 3y at the DOF holders: linear interpolation
    // reproduces it exactly.
    const auto dof_coords = source.GetLayout().GetDOFHolderCoordinates();
    const auto num_dofs = dof_coords.GetValues().extent(0);
    Kokkos::View<pcms::Real**, pcms::HostMemorySpace> coords_host("coords_host",
                                                                  num_dofs, 2);
    auto coords_mirror =
      Kokkos::create_mirror_view(pcms::DeviceMemorySpace(), coords_host);
    pcms::ConvertMismatchLayoutView2D(coords_mirror, dof_coords.GetValues());
    Kokkos::deep_copy(coords_host, coords_mirror);
    std::vector<pcms::Real> values(num_dofs);
    for (std::size_t i = 0; i < num_dofs; ++i) {
      values[i] = Linear(coords_host(i, 0), coords_host(i, 1));
    }
    source.SetDOFHolderDataHost(
      pcms::ValueView<const pcms::Real, pcms::HostMemorySpace>(
        source.GetData().GetValueBasis(),
        pcms::Rank2View<const pcms::Real, pcms::HostMemorySpace>(values.data(),
                                                                 num_dofs, 1)));

    // Evaluate at points that move between calls.
    const std::vector<pcms::Real> points{0.3, 0.4, 0.75, 0.1, 0.9, 0.9};
    const std::size_t num_points = points.size() / 2;
    Kokkos::View<pcms::Real**, pcms::DeviceMemorySpace> point_data(
      "points", num_points, 2);
    auto point_host = Kokkos::create_mirror_view(point_data);
    for (std::size_t i = 0; i < num_points; ++i) {
      point_host(i, 0) = points[2 * i];
      point_host(i, 1) = points[2 * i + 1];
    }
    Kokkos::deep_copy(point_data, point_host);
    pcms::QueryPoints query(geometry);
    query.Locate(pcms::CoordinateView<pcms::DeviceMemorySpace>(
      cartesian,
      pcms::Rank2View<pcms::Real, pcms::DeviceMemorySpace>(
        point_data.data(), point_data.extent(0), point_data.extent(1))));
    Kokkos::View<pcms::Real**, pcms::DeviceMemorySpace> out("out", num_points,
                                                            1);
    query.Evaluate(source,
                   pcms::ValueView<pcms::Real, pcms::DeviceMemorySpace>(
                     source.GetData().GetValueBasis(),
                     pcms::Rank2View<pcms::Real, pcms::DeviceMemorySpace>(
                       out.data(), out.extent(0), out.extent(1))));
    const auto out_host =
      Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace(), out);
    for (std::size_t i = 0; i < num_points; ++i) {
      if (std::abs(out_host(i, 0) - Linear(points[2 * i], points[2 * i + 1])) >
          LinearTol) {
        ++mismatches;
      }
    }

    // Interpolate onto a uniform grid inside the mesh.
    const pcms::UniformGrid<2> grid{
      .edge_length = {0.8, 0.6}, .bot_left = {0.1, 0.2}, .divisions = {9, 7}};
    auto target_space = pcms::CreateLagrangeSpace(
      pcms::uniform_grid::CreateGeometry(grid, cartesian), 1, 1);
    auto target = target_space->CreateFunction<pcms::Real>("target");
    const auto interpolate =
      pcms::CreateInterpolator<pcms::Real>(space, target_space);
    interpolate->Apply(source, target);

    const auto target_coords = target.GetLayout().GetDOFHolderCoordinates();
    const auto num_target = target_coords.GetValues().extent(0);
    Kokkos::View<pcms::Real**, pcms::HostMemorySpace> target_coords_host(
      "target_coords_host", num_target, 2);
    auto target_mirror =
      Kokkos::create_mirror_view(pcms::DeviceMemorySpace(), target_coords_host);
    pcms::ConvertMismatchLayoutView2D(target_mirror, target_coords.GetValues());
    Kokkos::deep_copy(target_coords_host, target_mirror);
    const auto result = target.GetDOFHolderDataHost().GetValues();
    for (std::size_t i = 0; i < num_target; ++i) {
      if (std::abs(result.data_handle()[i] -
                   Linear(target_coords_host(i, 0), target_coords_host(i, 1))) >
          LinearTol) {
        ++mismatches;
      }
    }
    std::printf("evaluated and interpolated %d values, %d mismatches\n",
                static_cast<int>(num_points + num_target), mismatches);
  }
  MPI_Finalize();
  return mismatches == 0 ? 0 : 1;
}
