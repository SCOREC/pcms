#include "point_cloud.h"
#include "pcms/profile.h"
#include "pcms/assert.h"

namespace pcms
{
struct PointCloudLocalizationHint
{
  PointCloudLocalizationHint(CoordinateView<HostMemorySpace> coordinate_view)
    : coordinates_("", coordinate_view.GetCoordinates().extent(0),
                   coordinate_view.GetCoordinates().extent(1))
  {
    auto coords = coordinate_view.GetCoordinates();
    int n = coords.extent(0);
    Kokkos::parallel_for(
      n, KOKKOS_LAMBDA(int i) {
        for (int j = 0; j < coords.extent(1); ++j) {
          coordinates_(i, j) = coords(i, j);
        }
      });
  }

  Kokkos::View<Real**> coordinates_;
};

PointCloud::PointCloud(std::string name, CoordinateSystem coordinate_system,
                       const PointCloudLayout& layout)
  : name_(name),
    coordinate_system_(coordinate_system),
    layout_(layout),
    data_("", layout_.GetDOFHolderCoordinates().extent(0))
{
}

const std::string& PointCloud::GetName() const
{
  return name_;
}

CoordinateSystem PointCloud::GetCoordinateSystem() const
{
  return coordinate_system_;
}

FieldDataView<const Real, HostMemorySpace> PointCloud::GetDOFHolderData() const
{
  Rank1View<const Real, pcms::HostMemorySpace> array_view{std::data(data_),
                                                          std::size(data_)};
  FieldDataView<const Real, HostMemorySpace> data_view{array_view,
                                                       GetCoordinateSystem()};
  return data_view;
}

CoordinateView<HostMemorySpace> PointCloud::GetDOFHolderCoordinates() const
{
  return CoordinateView<HostMemorySpace>{GetCoordinateSystem(),
                                         layout_.GetDOFHolderCoordinates()};
}

void PointCloud::SetDOFHolderData(
  FieldDataView<const Real, HostMemorySpace> data)
{
  PCMS_FUNCTION_TIMER;
  if (data.GetCoordinateSystem() != coordinate_system_) {
    throw std::runtime_error("Coordinate system mismatch");
  }

  Rank1View<const Real, HostMemorySpace> values = data.GetValues();
  PCMS_ALWAYS_ASSERT(static_cast<LO>(values.size()) == data_.size());
  Kokkos::parallel_for(
    values.size(), KOKKOS_LAMBDA(int i) { data_(i) = values[i]; });
}

LocalizationHint PointCloud::GetLocalizationHint(
  CoordinateView<HostMemorySpace> coordinate_view) const
{
  auto hint = std::make_shared<PointCloudLocalizationHint>(coordinate_view);
  return LocalizationHint{hint};
}

void PointCloud::Evaluate(LocalizationHint location,
                          FieldDataView<double, HostMemorySpace> results) const
{
  throw std::runtime_error("Not implemented");
}

void PointCloud::EvaluateGradient(FieldDataView<double, HostMemorySpace> results)
{
  throw std::runtime_error("Not implemented");
}

const FieldLayout& PointCloud::GetLayout() const
{
  return layout_;
}

bool PointCloud::CanEvaluateGradient()
{
  return false;
}

int PointCloud::Serialize(
  Rank1View<double, pcms::HostMemorySpace> buffer,
  Rank1View<const pcms::LO, pcms::HostMemorySpace> permutation) const
{
  PCMS_FUNCTION_TIMER;
  if (buffer.size() > 0) {
    Kokkos::parallel_for(
      data_.size(), KOKKOS_LAMBDA(int i) { buffer[permutation[i]] = data_(i); });
  }
  return data_.size();
}

void PointCloud::Deserialize(
  Rank1View<const double, pcms::HostMemorySpace> buffer,
  Rank1View<const pcms::LO, pcms::HostMemorySpace> permutation)
{
    PCMS_FUNCTION_TIMER;
    Kokkos::parallel_for(
      data_.size(), KOKKOS_LAMBDA(int i) { data_(i) = buffer[permutation[i]]; });
}

} // namespace pcms
