// Uses PCMS core through its installed public headers only, without Omega_h
// or MPI: evaluates a linear function with QueryPoints in a linear Lagrange
// space on a uniform grid and in a polynomial-reconstruction (MLS) space on a
// point cloud.

#include <pcms/backends/point_cloud/geometry.hpp>
#include <pcms/backends/uniform_grid/geometry.hpp>
#include <pcms/field/function_spaces.hpp>
#include <pcms/field/query_points.hpp>
#include <pcms/geometry/geometry.hpp>
#include <pcms/utility/coordinate_systems/cartesian.hpp>

#include <Kokkos_Core.hpp>

#include <cmath>
#include <cstdio>
#include <cstddef>
#include <memory>
#include <vector>

namespace
{
// linear interpolation reproduces a linear function exactly
constexpr pcms::Real LagrangeTol = 1e-10;
// linear MLS on a 9x9 point grid, as the MLS reproduction tests
constexpr pcms::Real MLSTol = 5e-3;

pcms::Real Linear(pcms::Real x, pcms::Real y)
{
  return 1.0 + 2.0 * x + 3.0 * y;
}

// Sets Linear at the DOF holders of a function on `space`, evaluates it at a
// few points inside [0, 1]^2 with QueryPoints and returns the number of
// values further than `tol` from Linear.
int EvaluateLinear(const std::shared_ptr<const pcms::FunctionSpace>& space,
                   const std::shared_ptr<const pcms::CoordinateSystem>& system,
                   pcms::Real tol)
{
  auto source = space->CreateFunction<pcms::Real>("source");
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

  const std::vector<pcms::Real> points{0.3, 0.4, 0.75, 0.1, 0.6, 0.6};
  const std::size_t num_points = points.size() / 2;
  Kokkos::View<pcms::Real**, pcms::DeviceMemorySpace> point_data("points",
                                                                 num_points, 2);
  auto point_host = Kokkos::create_mirror_view(point_data);
  for (std::size_t i = 0; i < num_points; ++i) {
    point_host(i, 0) = points[2 * i];
    point_host(i, 1) = points[2 * i + 1];
  }
  Kokkos::deep_copy(point_data, point_host);
  pcms::QueryPoints query(space->GetGeometry());
  query.Locate(pcms::CoordinateView<pcms::DeviceMemorySpace>(
    system, pcms::Rank2View<pcms::Real, pcms::DeviceMemorySpace>(
              point_data.data(), point_data.extent(0), point_data.extent(1))));
  Kokkos::View<pcms::Real**, pcms::DeviceMemorySpace> out("out", num_points, 1);
  query.Evaluate(source, pcms::ValueView<pcms::Real, pcms::DeviceMemorySpace>(
                           source.GetData().GetValueBasis(),
                           pcms::Rank2View<pcms::Real, pcms::DeviceMemorySpace>(
                             out.data(), out.extent(0), out.extent(1))));
  const auto out_host =
    Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace(), out);
  int mismatches = 0;
  for (std::size_t i = 0; i < num_points; ++i) {
    if (std::abs(out_host(i, 0) - Linear(points[2 * i], points[2 * i + 1])) >
        tol) {
      ++mismatches;
    }
  }
  return mismatches;
}
} // namespace

int main(int argc, char** argv)
{
  int mismatches = 0;
  {
    Kokkos::ScopeGuard kokkos(argc, argv);
    const auto cartesian = pcms::csys::CreateCartesian(2);

    // a linear Lagrange space on a uniform grid
    const pcms::UniformGrid<2> grid{
      .edge_length = {1.0, 1.0}, .bot_left = {0.0, 0.0}, .divisions = {4, 4}};
    const int lagrange = EvaluateLinear(
      pcms::CreateLagrangeSpace(
        pcms::uniform_grid::CreateGeometry(grid, cartesian), 1, 1),
      cartesian, LagrangeTol);

    // linear MLS on a point cloud (a 9x9 grid of points)
    const int n = 9;
    std::vector<pcms::Real> sites;
    for (int j = 0; j < n; ++j) {
      for (int i = 0; i < n; ++i) {
        sites.push_back(static_cast<pcms::Real>(i) / (n - 1));
        sites.push_back(static_cast<pcms::Real>(j) / (n - 1));
      }
    }
    const int mls = EvaluateLinear(
      pcms::CreatePolynomialReconstructionSpace(
        pcms::point_cloud::CreateGeometry(
          pcms::CoordinateView<pcms::HostMemorySpace>(
            cartesian, pcms::Rank2View<pcms::Real, pcms::HostMemorySpace>(
                         sites.data(), n * n, 2)),
          pcms::NeighbourhoodOptions{.radius = 0.45, .min_supports = 8}),
        pcms::MLSOptions{.degree = 1}),
      cartesian, MLSTol);

    mismatches = lagrange + mls;
    std::printf("uniform-grid Lagrange: %d mismatches; point-cloud MLS: %d "
                "mismatches\n",
                lagrange, mls);
  }
  return mismatches == 0 ? 0 : 1;
}
