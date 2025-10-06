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

PointCloud::PointCloud(const PointCloudLayout& layout)
  : layout_(layout),
    data_("", layout_.GetDOFHolderCoordinates().GetCoordinates().extent(0))
{
}

Rank1View<const Real, HostMemorySpace> PointCloud::GetDOFHolderData() const
{
  return make_const_array_view(data_);
}

void PointCloud::SetDOFHolderData(Rank1View<const Real, HostMemorySpace> data)
{
  PCMS_FUNCTION_TIMER;
  PCMS_ALWAYS_ASSERT(data.size() == data_.size());
  Kokkos::parallel_for(
    data.size(), KOKKOS_LAMBDA(int i) { data_(i) = data[i]; });
}

LocalizationHint PointCloud::GetLocalizationHint(
  CoordinateView<HostMemorySpace> coordinate_view) const
{
  auto hint = std::make_shared<PointCloudLocalizationHint>(coordinate_view);
  return LocalizationHint{hint};
}

void PointCloud::Evaluate(
  LocalizationHint /* unused */,
  FieldDataView<Real, HostMemorySpace> /* unused */) const
{
  throw std::runtime_error("Not implemented");
}

void PointCloud::EvaluateGradient(
  FieldDataView<Real, HostMemorySpace> /* unused */)
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
  Rank1View<Real, pcms::HostMemorySpace> buffer,
  Rank1View<const pcms::LO, pcms::HostMemorySpace> permutation) const
{
  PCMS_FUNCTION_TIMER;
  if (buffer.size() > 0) {
    Kokkos::parallel_for(
      data_.size(),
      KOKKOS_LAMBDA(int i) { buffer[permutation[i]] = data_(i); });
  }
  return data_.size();
}

void PointCloud::Deserialize(
  Rank1View<const Real, pcms::HostMemorySpace> buffer,
  Rank1View<const pcms::LO, pcms::HostMemorySpace> permutation)
{
  PCMS_FUNCTION_TIMER;
  Kokkos::parallel_for(
    data_.size(), KOKKOS_LAMBDA(int i) { data_(i) = buffer[permutation[i]]; });
}

} // namespace pcms
