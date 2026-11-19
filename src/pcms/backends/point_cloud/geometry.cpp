#include "pcms/backends/point_cloud/point_cloud_geometry.hpp"

namespace pcms::point_cloud
{

std::shared_ptr<const Geometry> CreateGeometry(
  CoordinateView<HostMemorySpace> coords, NeighbourhoodOptions options)
{
  const auto values = coords.GetValues();
  const int dim = static_cast<int>(values.extent(1));
  Kokkos::View<const Real**, Kokkos::HostSpace> host_view(
    values.data_handle(), values.extent(0), values.extent(1));
  auto device_view = Kokkos::View<Real**>("device_view", host_view.extent(0),
                                          host_view.extent(1));
  DeepCopyMismatchLayouts(device_view, host_view);
  auto layout = std::make_shared<pcms::detail::PointCloudLayout>(
    dim, device_view, coords.GetCoordinateSystem());
  return std::make_shared<PointCloudGeometry>(std::move(layout),
                                              std::move(options));
}

std::shared_ptr<const Geometry> CreateGeometry(
  std::shared_ptr<const FieldLayout> sites, NeighbourhoodOptions options)
{
  if (sites == nullptr) {
    throw pcms_error("CreateGeometry: the sites must not be null");
  }
  return std::make_shared<PointCloudGeometry>(std::move(sites),
                                              std::move(options));
}

} // namespace pcms::point_cloud
